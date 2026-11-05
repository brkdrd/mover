#!/usr/bin/env bash
# Runs the bounded space-dynamics throughput benchmark and records raw results.
#
#   ./run_benchmark.sh [path-to-SpaceSimBenchmark] [output-dir]
#
# The binary is built by Tools/SpaceSimBenchmark/CMakeLists.txt; its build
# directory must still contain compile_commands.json so the actual per-target
# compiler flags can be recorded (override with SPACESIM_BENCH_COMPILE_COMMANDS).
#
# Each configuration is run twice with the same seed; the end-state hash must be
# a 16-hex-digit value and must match, otherwise the script fails. The
# single-threaded and 4-worker separated-256 runs must also produce the same
# hash (same scenario and parameters). Every JSON is parsed strictly (nan/inf
# tokens are rejected), its manifest is checked against the requested
# configuration, its source fingerprint is checked against the checked-out
# harness sources, and every CSV must have exactly one row per step.
#
# Every configuration runs the same number of steps and the same simulated
# duration. Nothing is discarded: the first Update call is measured like any other.
#
# This script never overwrites a previously recorded results directory by
# itself; pass a new directory (e.g. results/verified-current) for a rerun.
set -u

BIN=${1:-}
OUT=${2:-results}
SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON=${PYTHON:-python3}

if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
	echo "usage: $0 <path-to-SpaceSimBenchmark> [output-dir]" >&2
	exit 2
fi

if ! mkdir -p "$OUT" || [ ! -d "$OUT" ]; then
	echo "cannot create output directory '$OUT'" >&2
	exit 2
fi

# Never silently overwrite a directory that already holds a run (the original
# capture in results/ in particular). A rerun goes into a new directory.
if compgen -G "$OUT/*.rep1.json" >/dev/null && [ "${SPACESIM_BENCH_ALLOW_OVERWRITE:-0}" != "1" ]; then
	echo "refusing to overwrite existing results in '$OUT' (set SPACESIM_BENCH_ALLOW_OVERWRITE=1 to override)" >&2
	exit 2
fi

# Content fingerprint of the compiled harness inputs. It must match the value
# embedded in the binary at configure time, otherwise the binary and the
# checked-out harness sources are not the same thing.
harness_fingerprint() {
	local f
	for f in src/main.cpp CMakeLists.txt; do
		printf '%s:%s\n' "$f" "$(sha256sum "$SELF_DIR/$f" | cut -d' ' -f1)"
	done | sha256sum | cut -d' ' -f1
}
SOURCE_FP=$(harness_fingerprint)

STAMP=$(date -u +%Y%m%dT%H%M%SZ)
ENVFILE="$OUT/environment_$STAMP.txt"
CC_JSON=${SPACESIM_BENCH_COMPILE_COMMANDS:-$(dirname "$BIN")/compile_commands.json}

failures=0
fail() {
	echo "FAIL: $*" >&2
	failures=$((failures + 1))
}

{
	echo "utc=$STAMP"
	echo "uname=$(uname -a)"
	echo "compiler=$(${CXX:-g++} --version | head -1)"
	if command -v lscpu >/dev/null 2>&1; then
		lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Core|Socket|L3 cache|Flags' | head -12
	fi
	if command -v git >/dev/null 2>&1; then
		echo "git_head=$(git -C "$SELF_DIR/../.." rev-parse HEAD 2>/dev/null)"
		echo "git_dirty_entries=$(git -C "$SELF_DIR/../.." status --porcelain 2>/dev/null | wc -l)"
		echo "git_dirty_scope_begin"
		git -C "$SELF_DIR/../.." status --porcelain 2>/dev/null | sed 's/^/dirty: /'
		echo "git_dirty_scope_end"
	fi
	echo "binary=$BIN"
	echo "binary_sha256=$(sha256sum "$BIN" | cut -d' ' -f1)"
	echo "harness_source_sha256=$SOURCE_FP"
} > "$ENVFILE"

# Record the flags that were actually used for the library and the harness TU.
if [ -f "$CC_JSON" ]; then
	echo "compile_commands=$CC_JSON" >> "$ENVFILE"
	echo "compile_commands_sha256=$(sha256sum "$CC_JSON" | cut -d' ' -f1)" >> "$ENVFILE"
	if ! "$PYTHON" - "$CC_JSON" >> "$ENVFILE" <<'PY'
import json
import sys


def reject(token):
    raise ValueError("non-finite JSON token %r" % token)


try:
    with open(sys.argv[1], "r", encoding="utf-8") as f:
        commands = json.load(f, parse_constant=reject)
    if not isinstance(commands, list) or any(
            not isinstance(entry, dict) or not isinstance(entry.get("file"), str)
            for entry in commands):
        raise ValueError("expected an array of compile-command objects with string file paths")
    targets = []
    for label, suffix in (("benchmark", "/SpaceSimBenchmark/src/main.cpp"),
                          ("library", "/Jolt/Physics/PhysicsSystem.cpp")):
        cmd = next((entry.get("command") for entry in commands
                    if entry["file"].replace("\\", "/").endswith(suffix)), None)
        if not isinstance(cmd, str) or not cmd.strip():
            raise ValueError("missing %s target command (%s); regenerate compile_commands.json from the benchmark build" %
                             (label, suffix))
        targets.append((label, cmd))
except (OSError, ValueError) as exc:
    print("compile_commands.json '%s': %s" % (sys.argv[1], exc), file=sys.stderr)
    sys.exit(1)

for label, cmd in targets:
    print("%s_compile_command=%s" % (label, cmd))
    print("%s_has_flto=%s" % (label, "yes" if "-flto" in cmd else "no"))
    print("%s_has_ffp_contract_off=%s" % (label, "yes" if "-ffp-contract=off" in cmd else "no"))
    print("%s_has_fast_math=%s" % (label, "yes" if ("-ffast-math" in cmd or "-Ofast" in cmd) else "no"))
PY
	then
		fail "cannot record actual compiler flags from '$CC_JSON'; fix or regenerate compile_commands.json before running"
		exit 2
	fi
	if grep -q '_has_fast_math=yes' "$ENVFILE"; then
		fail "compile_commands.json shows -ffast-math/-Ofast; this is not a measurement build"
	fi
else
	echo "compile_commands=missing" >> "$ENVFILE"
	fail "compile_commands.json not found at '$CC_JSON'; cannot record actual compiler flags"
	exit 2
fi
echo "environment recorded in $ENVFILE"

# Strict artifact validation. Exits non-zero on any inconsistency.
validate_artifacts() {
	local json=$1 csv=$2 scenario=$3 bodies=$4 steps=$5 threads=$6 hash=$7
	"$PYTHON" - "$json" "$csv" "$scenario" "$bodies" "$steps" "$threads" "$hash" "$SOURCE_FP" <<'PY'
import json
import math
import re
import sys

json_path, csv_path, scenario, bodies, steps, threads, state_hash, source_fp = sys.argv[1:9]
bodies, steps, threads = int(bodies), int(steps), int(threads)
errors = []
HASH16 = re.compile(r"^[0-9a-f]{16}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


def reject(token):
    raise ValueError("non-finite JSON token %r" % token)


try:
    with open(json_path, "r", encoding="utf-8") as f:
        doc = json.load(f, parse_constant=reject)
except Exception as exc:  # noqa: BLE001 - report the real parse failure
    print("%s: not strict JSON: %s" % (json_path, exc), file=sys.stderr)
    sys.exit(1)


def walk(obj, trail="doc"):
    if isinstance(obj, float) and not math.isfinite(obj):
        errors.append("non-finite number at %s" % trail)
    elif isinstance(obj, dict):
        for key, value in obj.items():
            walk(value, "%s.%s" % (trail, key))
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            walk(value, "%s[%d]" % (trail, index))


walk(doc)
manifest, timing, checks = doc["manifest"], doc["timing"], doc["checks"]

if doc.get("result") != "PASS":
    errors.append("result is %r, expected PASS (exit code and artifact disagree)" % doc.get("result"))
if manifest.get("scenario") != scenario:
    errors.append("scenario %r != %r" % (manifest.get("scenario"), scenario))
for key, expected in (("bodies", bodies), ("steps", steps), ("worker_threads", threads)):
    if manifest.get(key) != expected:
        errors.append("%s %r != %r" % (key, manifest.get(key), expected))
try:
    if not (math.isfinite(manifest["rate_hz"]) and manifest["rate_hz"] > 0):
        errors.append("rate_hz is not finite and positive")
    if not (math.isfinite(manifest["timestep_s"]) and manifest["timestep_s"] > 0):
        errors.append("timestep_s is not finite and positive")
    if not (math.isfinite(manifest["simulated_seconds"]) and manifest["simulated_seconds"] > 0):
        errors.append("simulated_seconds is not finite and positive")
    product = manifest["timestep_s"] * steps
    if abs(product - manifest["simulated_seconds"]) > 1e-6 * max(1.0, abs(product)):
        errors.append("simulated_seconds is not steps * timestep_s")
except KeyError as exc:
    errors.append("manifest key missing: %s" % exc)

expected_samples = steps // 60 + (1 if steps % 60 else 0)
if checks.get("samples") != expected_samples:
    errors.append("samples %r != %d sampled states" % (checks.get("samples"), expected_samples))

if not HASH16.match(str(doc.get("state_hash", ""))):
    errors.append("state_hash %r is not 16 lowercase hex digits" % doc.get("state_hash"))
elif doc["state_hash"] != state_hash:
    errors.append("state_hash %s != %s reported in the log" % (doc["state_hash"], state_hash))

for gate in doc.get("gates", []):
    if gate.get("gated") and not gate.get("pass"):
        errors.append("gated check %r failed" % gate.get("name"))

for key in ("harness_source_sha256", "core_source_sha256", "build_sha256"):
    if not SHA256.match(str(manifest.get(key, ""))):
        errors.append("%s is not a 64-hex fingerprint" % key)
if manifest.get("harness_source_sha256") != source_fp:
    errors.append("harness_source_sha256 does not match the checked-out harness sources "
                  "(binary is stale; rebuild it)")
if manifest.get("fast_math") != "false":
    errors.append("fast_math is %r" % manifest.get("fast_math"))
if manifest.get("benchmark_lto") != "off":
    errors.append("benchmark_lto is %r (expected the harness TU to be non-LTO)" % manifest.get("benchmark_lto"))
if not manifest.get("compiler"):
    errors.append("compiler is empty")

for key, value in timing.items():
    if isinstance(value, float) and not math.isfinite(value):
        errors.append("timing.%s is non-finite" % key)
if timing.get("sum_step_ns", 0) <= 0 or timing.get("throughput_steps_per_second", 0) <= 0:
    errors.append("timing is empty (sum_step_ns/throughput must be positive)")
if timing.get("run_wall_ns", 0) < timing.get("sum_step_ns", 0):
    errors.append("run_wall_ns is smaller than sum_step_ns")

try:
    with open(csv_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
except Exception as exc:  # noqa: BLE001
    print("%s: cannot read CSV: %s" % (csv_path, exc), file=sys.stderr)
    sys.exit(1)
if not lines or lines[0] != "step,ns":
    errors.append("%s: bad CSV header" % csv_path)
else:
    rows = lines[1:]
    if len(rows) != steps:
        errors.append("%s: %d CSV rows for %d steps" % (csv_path, len(rows), steps))
    for index, line in enumerate(rows, 1):
        m = re.match(r"^([0-9]+),([0-9]+)$", line)
        if not m:
            errors.append("%s: malformed row %d: %r" % (csv_path, index, line))
            break
        if int(m.group(1)) != index:
            errors.append("%s: step index %s on row %d" % (csv_path, m.group(1), index))
            break

if errors:
    for error in errors:
        print("%s: %s" % (json_path, error), file=sys.stderr)
    sys.exit(1)
PY
}

run_one() {
	local name=$1
	shift
	local scenario=$1 bodies=$2 steps=$3 threads=$4
	shift 4
	local rc1 rc2 h1 h2
	echo "=== $name ==="
	"$BIN" "$@" --json "$OUT/$name.rep1.json" --csv "$OUT/$name.rep1.csv" > "$OUT/$name.rep1.log" 2>&1
	rc1=$?
	"$BIN" "$@" --json "$OUT/$name.rep2.json" --csv "$OUT/$name.rep2.csv" > "$OUT/$name.rep2.log" 2>&1
	rc2=$?
	cat "$OUT/$name.rep1.log"
	h1=$(sed -n 's/^state hash: //p' "$OUT/$name.rep1.log" | head -1)
	h2=$(sed -n 's/^state hash: //p' "$OUT/$name.rep2.log" | head -1)

	if [ "$rc1" -ne 0 ] || [ "$rc2" -ne 0 ]; then
		echo "NOTE: $name reported a failed check or artifact error (exit $rc1/$rc2) -- recorded, not hidden" >&2
		failures=$((failures + 1))
	fi
	for h in "$h1" "$h2"; do
		if ! [[ "$h" =~ ^[0-9a-f]{16}$ ]]; then
			fail "$name: state hash '$h' is not 16 lowercase hex digits"
		fi
	done
	if [ "$h1" != "$h2" ]; then
		fail "$name: state hashes differ between repeats ($h1 vs $h2)"
	fi
	if [ "$rc1" -eq 0 ]; then
		validate_artifacts "$OUT/$name.rep1.json" "$OUT/$name.rep1.csv" "$scenario" "$bodies" "$steps" "$threads" "$h1" \
			|| fail "$name: repeat 1 artifacts are inconsistent"
	fi
	if [ "$rc2" -eq 0 ]; then
		validate_artifacts "$OUT/$name.rep2.json" "$OUT/$name.rep2.csv" "$scenario" "$bodies" "$steps" "$threads" "$h2" \
			|| fail "$name: repeat 2 artifacts are inconsistent"
	fi
	printf -v "HASH_${name//-/_}" '%s' "$h1"
}

# Primary results: single-threaded job system, identical step count and duration.
run_one separated-128 separated 128 1200 0 --scenario separated --bodies 128 --steps 1200 --threads 0
run_one separated-256 separated 256 1200 0 --scenario separated --bodies 256 --steps 1200 --threads 0
run_one interacting-128 interacting 128 1200 0 --scenario interacting --bodies 128 --steps 1200 --threads 0

# Worker-thread scale check on the same workload (not the primary figure).
run_one separated-256-4t separated 256 1200 4 --scenario separated --bodies 256 --steps 1200 --threads 4

# The declared determinism claim is that this workload is bit-identical across
# the two job systems. Verify it from the recorded hashes.
h_single=${HASH_separated_256:-}
h_four=${HASH_separated_256_4t:-}
if [ -z "$h_single" ] || [ -z "$h_four" ]; then
	fail "thread-invariance check could not read separated-256 hashes ('$h_single' vs '$h_four')"
elif [ "$h_single" != "$h_four" ]; then
	fail "separated-256 state hash differs between 0 and 4 worker threads ($h_single vs $h_four)"
else
	echo "thread invariance: separated-256 hash $h_single identical for 0 and 4 workers"
fi

# Harness self-test: reference math against analytic fixtures.
if "$BIN" --selftest > "$OUT/selftest.log" 2>&1 && grep -q '^SELFTEST: PASS$' "$OUT/selftest.log"; then
	cat "$OUT/selftest.log"
else
	cat "$OUT/selftest.log"
	fail "harness self-test failed"
fi

"$PYTHON" "$SELF_DIR/summarize.py" "$OUT" || fail "summarize.py failed"

if [ "$failures" -ne 0 ]; then
	echo "RESULT: $failures failure(s) recorded" >&2
	exit 1
fi
echo "RESULT: all runs completed, checks passed, artifacts validated, repeats deterministic"
