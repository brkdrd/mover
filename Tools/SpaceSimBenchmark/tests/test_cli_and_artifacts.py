#!/usr/bin/env python3
"""Executable regression tests for the SpaceSimBenchmark CLI and artifact I/O.

Usage:
    python3 test_cli_and_artifacts.py <path-to-SpaceSimBenchmark>

These tests drive the real binary as a subprocess, so they work for a native
binary and for a WSL/Linux binary launched from a POSIX shell. They were written
against the pre-fix harness and fail (RED) on it:

  * --rate inf / nan / -120 was accepted and produced a PASS run with an
    "inf"/"nan" rate in the JSON manifest and a zero timestep;
  * integer options were parsed with strtoul/strtod without an end pointer, so
    "128x", "-1" and overflow were silently truncated;
  * --json/--csv pointing into a missing directory only printed a warning and
    still exited 0 with RESULT: PASS;
  * fwrite/fprintf/fflush/fclose failures (e.g. /dev/full) were ignored;
  * the JSON manifest had no build provenance or source fingerprints;
  * dangling destination links could become aliases after the first artifact
    write, silently overwriting CSV with JSON while still reporting PASS.

After the fixes every case below passes (GREEN).

Only the Python standard library is used. The tests never write outside a
temporary directory and never assume the process runs as root or non-root.
"""

import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

BINARY = None
HASH16 = re.compile(r"^[0-9a-f]{16}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
RUN_TIMEOUT = 120


def run_bench(*args, timeout=RUN_TIMEOUT):
    return subprocess.run([BINARY, *args], capture_output=True, text=True, timeout=timeout)


class CliValidationTests(unittest.TestCase):
    """Invalid options must be rejected before any workload runs."""

    def assert_rejected(self, *args, expect_no_artifacts=True):
        with tempfile.TemporaryDirectory() as td:
            json_path = os.path.join(td, "run.json")
            csv_path = os.path.join(td, "run.csv")
            r = run_bench(*args, "--json", json_path, "--csv", csv_path)
            self.assertEqual(r.returncode, 2,
                             "expected CLI exit 2 for %r\nstdout:\n%s\nstderr:\n%s"
                             % (args, r.stdout, r.stderr))
            self.assertNotIn("RESULT: PASS", r.stdout)
            if expect_no_artifacts:
                self.assertFalse(os.path.exists(json_path), "no JSON artifact may be produced for %r" % (args,))
                self.assertFalse(os.path.exists(csv_path), "no CSV artifact may be produced for %r" % (args,))
            return r

    def test_rate_inf_is_rejected(self):
        r = self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "inf")
        self.assertIn("rate", r.stderr.lower())

    def test_rate_nan_is_rejected(self):
        self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "nan")

    def test_rate_negative_is_rejected(self):
        self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "-120")

    def test_rate_zero_is_rejected(self):
        self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "0")

    def test_rate_with_trailing_junk_is_rejected(self):
        self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "120junk")

    def test_rate_overflowing_timestep_is_rejected(self):
        # 1/1e-300 is not representable as a float; dt would be inf.
        self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "1e-300")

    def test_rate_underflowing_timestep_is_rejected(self):
        # 1/1e300 rounds to a zero float; the run would have dt = 0.
        self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--rate", "1e300")

    def test_rate_missing_value_is_rejected(self):
        r = run_bench("--scenario", "separated", "--rate")
        self.assertNotEqual(r.returncode, 0)

    def test_bodies_out_of_range_is_rejected(self):
        for bad in ("1", "0", "4097", "-2", "128x", "1e3", ""):
            with self.subTest(bodies=bad):
                self.assert_rejected("--scenario", "separated", "--bodies", bad, "--steps", "1")

    def test_steps_out_of_range_is_rejected(self):
        for bad in ("0", "-1", "1200x", "10000001", ""):
            with self.subTest(steps=bad):
                self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", bad)

    def test_threads_out_of_range_is_rejected(self):
        for bad in ("-1", "65", "4x", ""):
            with self.subTest(threads=bad):
                self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--threads", bad)

    def test_seed_trailing_junk_is_rejected(self):
        for bad in ("0x", "xyz", "0x5eed4d0000000001zz", "-1"):
            with self.subTest(seed=bad):
                self.assert_rejected("--scenario", "separated", "--bodies", "2", "--steps", "1", "--seed", bad)

    def test_missing_scenario_value_exits_2(self):
        r = run_bench("--scenario")
        self.assertEqual(r.returncode, 2)

    def test_unknown_argument_exits_2(self):
        r = run_bench("--not-an-option")
        self.assertEqual(r.returncode, 2)

    def test_help_exits_zero(self):
        r = run_bench("--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("Exit codes", r.stdout)

    def test_minimal_valid_workload_still_runs(self):
        # The smallest workload that the reviewer used to expose rate=inf must be
        # a valid run on its own: bodies=2, steps=1, rate=120.
        with tempfile.TemporaryDirectory() as td:
            json_path = os.path.join(td, "min.json")
            csv_path = os.path.join(td, "min.csv")
            r = run_bench("--scenario", "separated", "--bodies", "2", "--steps", "1",
                          "--rate", "120", "--json", json_path, "--csv", csv_path)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("RESULT: PASS", r.stdout)
            doc = load_json_strict(json_path)
            self.assertEqual(doc["manifest"]["bodies"], 2)
            self.assertEqual(doc["manifest"]["steps"], 1)
            self.assertEqual(doc["manifest"]["rate_hz"], 120.0)
            self.assertTrue(doc["manifest"]["timestep_s"] > 0.0)
            self.assertTrue(math.isfinite(doc["manifest"]["simulated_seconds"]))
            self.assertEqual(len(read_csv(csv_path)), 1)  # one measured step


class ArtifactDestinationTests(unittest.TestCase):
    """Aliased destinations are CLI errors and must not create or truncate files."""

    def snapshot(self, td):
        state = {}
        for parent, dirs, files in os.walk(td, followlinks=False):
            for name in dirs + files:
                path = os.path.join(parent, name)
                key = os.path.relpath(path, td)
                if os.path.islink(path):
                    state[key] = ("symlink", os.readlink(path))
                elif os.path.isdir(path):
                    state[key] = ("directory",)
                else:
                    with open(path, "rb") as f:
                        state[key] = ("file", f.read())
        return state

    def assert_destination_rejected(self, td, *args):
        before = self.snapshot(td)
        assert BINARY is not None
        result = subprocess.run(
            [BINARY, "--bodies", "2", "--steps", "1", *args],
            cwd=td, capture_output=True, text=True, timeout=RUN_TIMEOUT)
        # Check mutation and workload evidence even when the exit code is wrong.
        with self.subTest(check="no filesystem mutation"):
            self.assertEqual(self.snapshot(td), before,
                             "rejection must preserve files, symlinks and absent targets")
        with self.subTest(check="no workload"):
            for marker in ("instrument:", "scenario=", "state hash:", "RESULT:"):
                self.assertNotIn(marker, result.stdout)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertTrue(result.stderr.strip(), "rejection must explain the invalid destination")
        return result

    def assert_collision_rejected(self, td, json_path, csv_path):
        result = self.assert_destination_rejected(td, "--json", json_path, "--csv", csv_path)
        self.assertIn("--json", result.stderr)
        self.assertIn("--csv", result.stderr)

    def symlink(self, target, link, *, directory=False):
        try:
            os.symlink(target, link, target_is_directory=directory)
        except (OSError, NotImplementedError) as exc:
            self.skipTest("symlink unavailable: %s" % exc)

    def test_dangling_json_alias_to_csv_target_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            self.symlink("target.out", os.path.join(td, "alias.out"))
            self.assert_destination_rejected(td, "--json", "alias.out", "--csv", "target.out")

    def test_dangling_csv_alias_to_json_target_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            self.symlink("target.out", os.path.join(td, "alias.out"))
            self.assert_destination_rejected(td, "--json", "target.out", "--csv", "alias.out")

    def test_two_dangling_aliases_to_common_target_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            self.symlink("target.out", os.path.join(td, "json.out"))
            self.symlink("target.out", os.path.join(td, "csv.out"))
            self.assert_destination_rejected(td, "--json", "json.out", "--csv", "csv.out")

    def test_single_dangling_destination_is_rejected(self):
        for flag in ("--json", "--csv"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as td:
                self.symlink("absent.out", os.path.join(td, "alias.out"))
                self.assert_destination_rejected(td, flag, "alias.out")

    def test_dangling_symlink_chain_is_rejected(self):
        for flag in ("--json", "--csv"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as td:
                self.symlink("absent.out", os.path.join(td, "middle.out"))
                self.symlink("middle.out", os.path.join(td, "alias.out"))
                other = "--csv" if flag == "--json" else "--json"
                with open(os.path.join(td, "keep.out"), "wb") as f:
                    f.write(b"preserve existing other artifact\n")
                self.assert_destination_rejected(td, flag, "alias.out", other, "keep.out")

    def test_dangling_symlink_parent_is_rejected(self):
        for flag in ("--json", "--csv"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as td:
                self.symlink("absent-directory", os.path.join(td, "alias"), directory=True)
                self.assert_destination_rejected(td, flag, "alias/new.out")

    def test_symlink_loops_are_rejected(self):
        for flag in ("--json", "--csv"):
            for destination in ("loop", "loop/new.out"):
                with self.subTest(flag=flag, destination=destination), tempfile.TemporaryDirectory() as td:
                    self.symlink("loop", os.path.join(td, "loop"), directory=True)
                    self.assert_destination_rejected(td, flag, destination)

    def test_nondirectory_ancestor_is_rejected_before_workload(self):
        for flag in ("--json", "--csv"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as td:
                with open(os.path.join(td, "file"), "wb") as f:
                    f.write(b"not a directory\n")
                self.assert_destination_rejected(td, flag, "file/new.out")

    def test_distinct_existing_files_and_resolving_symlinks_remain_valid(self):
        for links in (False, True):
            with self.subTest(links=links), tempfile.TemporaryDirectory() as td:
                jp, cp = (os.path.join(td, name) for name in ("run.json", "run.csv"))
                for path in (jp, cp):
                    with open(path, "wb") as f:
                        f.write(b"replace this artifact\n")
                json_path, csv_path = jp, cp
                if links:
                    # Both a resolving leaf chain and a resolving parent are valid.
                    self.symlink("run.json", os.path.join(td, "middle.out"))
                    self.symlink("middle.out", os.path.join(td, "alias.out"))
                    os.mkdir(os.path.join(td, "parent"))
                    self.symlink("parent", os.path.join(td, "parent-alias"), directory=True)
                    json_path = os.path.join(td, "alias.out")
                    csv_path = os.path.join(td, "parent-alias", "..", "run.csv")
                result = run_bench("--bodies", "2", "--steps", "1",
                                   "--json", json_path, "--csv", csv_path)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(load_json_strict(jp)["result"], "PASS")
                self.assertEqual(len(read_csv(cp)), 1)
                if links:
                    self.assertEqual(os.readlink(json_path), "middle.out")
                    self.assertEqual(os.readlink(os.path.join(td, "middle.out")), "run.json")
                    self.assertEqual(os.readlink(os.path.join(td, "parent-alias")), "parent")

    def test_identical_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            self.assert_collision_rejected(td, "same.out", "same.out")

    def test_dot_relative_alias_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            self.assert_collision_rejected(td, os.path.join(td, "same.out"), "./same.out")

    def test_parent_normalized_alias_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            os.mkdir(os.path.join(td, "child"))
            self.assert_collision_rejected(td, "same.out", "child/../same.out")

    def check_existing_alias(self, link):
        with tempfile.TemporaryDirectory() as td:
            original = os.path.join(td, "original.out")
            alias = os.path.join(td, "alias.out")
            with open(original, "wb") as f:
                f.write(b"must remain unchanged\n")
            try:
                link(original, alias)
            except (OSError, NotImplementedError) as exc:
                self.skipTest("filesystem alias unavailable: %s" % exc)
            self.assert_collision_rejected(td, original, alias)

    def test_existing_symlink_alias_is_rejected_without_mutation(self):
        self.check_existing_alias(os.symlink)

    def test_existing_hardlink_alias_is_rejected_without_mutation(self):
        self.check_existing_alias(os.link)

    def test_symlink_parent_alias_is_rejected_before_creation(self):
        with tempfile.TemporaryDirectory() as td:
            parent = os.path.join(td, "parent")
            os.mkdir(parent)
            alias = os.path.join(td, "alias")
            try:
                os.symlink(parent, alias, target_is_directory=True)
            except (OSError, NotImplementedError) as exc:
                self.skipTest("directory symlink unavailable: %s" % exc)
            self.assert_collision_rejected(td, "parent/new.out", "alias/new.out")
            self.assertEqual(os.listdir(parent), [])


class ArtifactFailureTests(unittest.TestCase):
    """A requested artifact that cannot be written must fail the run."""

    def assert_artifact_failure(self, target, *, csv=False):
        with tempfile.TemporaryDirectory() as td:
            missing = os.path.join(td, "missing", "out.file")
            path = missing if target == "missing" else target
            other = os.path.join(td, "other.out")
            args = ["--scenario", "separated", "--bodies", "4", "--steps", "2"]
            if csv:
                args += ["--csv", path, "--json", other]
            else:
                args += ["--json", path, "--csv", other]
            r = run_bench(*args)
            self.assertNotEqual(r.returncode, 0,
                                "artifact write failure must be non-zero\nstdout:\n%s\nstderr:\n%s"
                                % (r.stdout, r.stderr))
            self.assertNotIn("RESULT: PASS", r.stdout,
                             "a failed artifact must never be reported as PASS")
            if path == missing:
                self.assertFalse(os.path.exists(path))

    def test_json_missing_parent_directory_fails(self):
        self.assert_artifact_failure("missing")

    def test_csv_missing_parent_directory_fails(self):
        self.assert_artifact_failure("missing", csv=True)

    def test_json_path_is_a_directory_fails(self):
        with tempfile.TemporaryDirectory() as td:
            self.assert_artifact_failure(td)

    def test_csv_path_is_a_directory_fails(self):
        with tempfile.TemporaryDirectory() as td:
            self.assert_artifact_failure(td, csv=True)

    @unittest.skipUnless(os.path.exists("/dev/full"), "/dev/full is Linux-only")
    def test_json_dev_full_fails(self):
        # fopen succeeds; the error surfaces on flush/close.
        self.assert_artifact_failure("/dev/full")

    @unittest.skipUnless(os.path.exists("/dev/full"), "/dev/full is Linux-only")
    def test_csv_dev_full_fails(self):
        self.assert_artifact_failure("/dev/full", csv=True)


@unittest.skipUnless(os.name == "posix" and shutil.which("bash"), "requires POSIX bash")
class CompileCommandPreflightTests(unittest.TestCase):
    """Reject untrustworthy provenance before invoking any benchmark workload."""

    def invoke_script(self, text):
        td = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(td, ignore_errors=True))
        commands = os.path.join(td, "compile_commands.json")
        if text is not None:
            with open(commands, "w", encoding="utf-8") as f:
                f.write(text)
        marker = os.path.join(td, "invoked")
        sentinel = os.path.join(td, "benchmark-sentinel")
        with open(sentinel, "w", encoding="utf-8") as f:
            # No simulated benchmark data: just record any attempted invocation.
            f.write("#!%s\nimport os, pathlib\n"
                    "pathlib.Path(os.environ['BENCH_TEST_MARKER']).touch()\n"
                    "raise SystemExit(93)\n" % sys.executable)
        os.chmod(sentinel, 0o755)
        out = os.path.join(td, "results")
        script = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "run_benchmark.sh"))
        env = dict(os.environ, SPACESIM_BENCH_COMPILE_COMMANDS=commands,
                   BENCH_TEST_MARKER=marker, PYTHON=sys.executable)
        result = subprocess.run(["bash", script, sentinel, out], env=env,
                                capture_output=True, text=True, timeout=RUN_TIMEOUT)
        return result, marker, out

    def assert_preflight_rejected(self, text):
        result, marker, out = self.invoke_script(text)
        self.assertFalse(os.path.exists(marker), "benchmark invoked before provenance rejection")
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("compile_commands", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.assertNotIn("=== separated", result.stdout)
        self.assertFalse(any(name.endswith((".json", ".csv", ".log"))
                             for name in os.listdir(out)))

    def test_missing_compile_commands_is_rejected_before_runs(self):
        self.assert_preflight_rejected(None)

    def test_malformed_or_nonstrict_commands_are_rejected_before_runs(self):
        for text in ("{", "[NaN]", "[Infinity]", "{}", "[null]"):
            with self.subTest(text=text):
                self.assert_preflight_rejected(text)

    def test_both_compile_targets_are_required_before_runs(self):
        benchmark = {"file": "/repo/Tools/SpaceSimBenchmark/src/main.cpp", "command": "c++ -O3 main.cpp"}
        library = {"file": "/repo/Jolt/Physics/PhysicsSystem.cpp", "command": "c++ -O3 -ffp-contract=off PhysicsSystem.cpp"}
        for entries in ([], [benchmark], [library],
                        [benchmark, dict(library, command="")],
                        [dict(benchmark, command=None), library]):
            with self.subTest(entries=entries):
                self.assert_preflight_rejected(json.dumps(entries))

    def test_valid_commands_record_distinct_target_metadata(self):
        result, marker, out = self.invoke_script(json.dumps([
            {"file": "/repo/Tools/SpaceSimBenchmark/src/main.cpp", "command": "c++ -O3 main.cpp"},
            {"file": "/repo/Jolt/Physics/PhysicsSystem.cpp",
             "command": "c++ -O3 -flto -ffp-contract=off PhysicsSystem.cpp"}]))
        self.assertTrue(os.path.exists(marker), result.stderr)
        self.assertNotEqual(result.returncode, 0)  # sentinel deliberately exits 93
        envfiles = [name for name in os.listdir(out) if name.startswith("environment_")]
        self.assertEqual(len(envfiles), 1)
        with open(os.path.join(out, envfiles[0]), encoding="utf-8") as f:
            metadata = dict(line.rstrip("\n").split("=", 1) for line in f if "=" in line)
        self.assertEqual(metadata["benchmark_compile_command"], "c++ -O3 main.cpp")
        self.assertEqual(metadata["library_compile_command"], "c++ -O3 -flto -ffp-contract=off PhysicsSystem.cpp")
        for label, expected in (("benchmark", "no"), ("library", "yes")):
            self.assertEqual(metadata[label + "_has_flto"], expected)
            self.assertEqual(metadata[label + "_has_ffp_contract_off"], expected)
            self.assertEqual(metadata[label + "_has_fast_math"], "no")


class ArtifactContentTests(unittest.TestCase):
    """Valid runs must produce strict, complete, self-describing artifacts."""

    def run_and_load(self, *args):
        td = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(td, ignore_errors=True))
        json_path = os.path.join(td, "run.json")
        csv_path = os.path.join(td, "run.csv")
        r = run_bench(*args, "--json", json_path, "--csv", csv_path)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RESULT: PASS", r.stdout)
        return load_json_strict(json_path), csv_path, json_path

    def test_time_fields_round_trip_the_actual_float_timestep(self):
        for rate in ("120", "1e12", "1e20", "1e40"):
            for steps in (1, 3):
                with self.subTest(rate=rate, steps=steps):
                    doc, _, _ = self.run_and_load(
                        "--bodies", "2", "--steps", str(steps), "--rate", rate)
                    # Mirror only the public timestep conversion, not decimal formatting.
                    dt = struct.unpack("f", struct.pack("f", 1.0 / float(rate)))[0]
                    manifest = doc["manifest"]
                    for key, expected in (("timestep_s", dt),
                                          ("simulated_seconds", dt * steps)):
                        self.assertIsInstance(manifest[key], (int, float))
                        self.assertGreater(manifest[key], 0)
                        self.assertEqual(manifest[key], expected, key)

    def test_json_is_strict_and_finite(self):
        doc, csv_path, _ = self.run_and_load("--scenario", "separated",
                                             "--bodies", "8", "--steps", "3", "--threads", "0")
        assert_finite_tree(doc)
        self.assertEqual(doc["result"], "PASS")
        self.assertTrue(HASH16.match(doc["state_hash"]))

    def test_manifest_records_build_provenance(self):
        doc, _, _ = self.run_and_load("--scenario", "separated",
                                      "--bodies", "8", "--steps", "3", "--threads", "0")
        m = doc["manifest"]
        for key in ("compiler", "build_type", "library_lto", "benchmark_lto",
                    "cross_platform_deterministic", "fast_math",
                    "compile_definitions", "benchmark_compile_flags",
                    "library_compile_flags", "cmake_cache",
                    "git_head_at_configure", "git_dirty_entries_at_configure"):
            self.assertIn(key, m)
            self.assertTrue(m[key] not in ("", None), key)
        for key in ("harness_source_sha256", "core_source_sha256", "build_sha256"):
            self.assertRegex(m[key], SHA256, key)
        self.assertGreater(m["core_source_files"], 0)
        self.assertEqual(m["fast_math"], "false")
        self.assertNotIn("-ffast-math", m["library_compile_flags"] + m["benchmark_compile_flags"])
        # The library and the harness TU are not compiled with identical flags.
        self.assertIn("-ffp-contract=off", m["library_compile_flags"])
        self.assertNotIn("-ffp-contract=off", m["benchmark_compile_flags"])
        self.assertEqual(m["benchmark_lto"], "off")

    def test_csv_has_one_row_per_step(self):
        doc, csv_path, _ = self.run_and_load("--scenario", "separated",
                                             "--bodies", "8", "--steps", "4", "--threads", "0")
        rows = read_csv(csv_path)
        self.assertEqual(len(rows), doc["manifest"]["steps"])
        self.assertTrue(all(ns >= 0 for ns in rows))

    def test_checks_cover_sampled_states_not_every_step(self):
        doc, _, _ = self.run_and_load("--scenario", "separated",
                                      "--bodies", "8", "--steps", "125", "--threads", "0")
        # Sampling runs every 60 steps and on the last step.
        self.assertEqual(doc["checks"]["samples"], 3)

    def test_repeat_hashes_match_and_threads_do_not_change_state(self):
        with tempfile.TemporaryDirectory() as td:
            hashes = {}
            for name, threads in (("a", "0"), ("b", "0"), ("c", "4")):
                jp = os.path.join(td, name + ".json")
                r = run_bench("--scenario", "separated", "--bodies", "16", "--steps", "4",
                              "--threads", threads, "--json", jp)
                self.assertEqual(r.returncode, 0, r.stderr)
                hashes[name] = load_json_strict(jp)["state_hash"]
            self.assertEqual(hashes["a"], hashes["b"], "repeats must be deterministic")
            self.assertEqual(hashes["a"], hashes["c"],
                             "separated workload must be bit-identical with 0 and 4 worker threads")
            self.assertTrue(HASH16.match(hashes["a"]))


def load_json_strict(path):
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    def reject(token):
        raise AssertionError("non-finite JSON token %r in %s" % (token, path))
    return json.loads(text, parse_constant=reject)


def assert_finite_tree(obj, trail="doc"):
    if isinstance(obj, float) and not math.isfinite(obj):
        raise AssertionError("non-finite number at %s" % trail)
    if isinstance(obj, dict):
        for k, v in obj.items():
            assert_finite_tree(v, "%s.%s" % (trail, k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            assert_finite_tree(v, "%s[%d]" % (trail, i))


def read_csv(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    if not lines or lines[0] != "step,ns":
        raise AssertionError("bad CSV header in %s" % path)
    rows = []
    for n, line in enumerate(lines[1:], 1):
        m = re.match(r"^([0-9]+),([0-9]+)$", line)
        if not m:
            raise AssertionError("bad CSV row %d in %s: %r" % (n, path, line))
        if int(m.group(1)) != n:
            raise AssertionError("CSV step index %s on row %d" % (m.group(1), n))
        rows.append(int(m.group(2)))
    return rows


def main():
    global BINARY
    if len(sys.argv) < 2:
        print("usage: test_cli_and_artifacts.py <path-to-SpaceSimBenchmark> [unittest args]", file=sys.stderr)
        return 2
    BINARY = os.path.abspath(sys.argv[1])
    if not os.path.isfile(BINARY) or not os.access(BINARY, os.X_OK):
        print("not an executable: %s" % BINARY, file=sys.stderr)
        return 2
    program = unittest.main(argv=[sys.argv[0], "-v", *sys.argv[2:]], exit=False)
    return 0 if program.result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
