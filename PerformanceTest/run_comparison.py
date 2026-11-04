#!/usr/bin/env python3
"""Run headless Jolt binaries sequentially and retain raw timings and state hashes.

Example:
  python3 PerformanceTest/run_comparison.py --binary 3d=build/perf-3d/PerformanceTest \
    --binary before=build/perf-results/PerformanceTest-baseline \
    --binary after=build/perf-off/PerformanceTest --threads 1 6 12 --output build/comparison
"""

import argparse
import csv
import json
import math
import platform
from pathlib import Path
import re
import statistics
import subprocess


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[max(0, math.ceil(fraction * len(ordered)) - 1)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", action="append", required=True, help="label=executable")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 6, 12])
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--warmup", type=int, default=100, help="Exclude these initial frames from steady-state statistics only")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--hz", type=float, help="Optional integration frequency; requires binaries supporting -hz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.warmup < args.steps or args.repeat < 1 or min(args.threads) < 1:
        parser.error("Require 0 <= warmup < steps, repeat >= 1, threads >= 1")
    if args.hz is not None and (not math.isfinite(args.hz) or args.hz <= 0):
        parser.error("Require a finite positive integration frequency")
    binaries = []
    for value in args.binary:
        label, separator, filename = value.partition("=")
        if not separator or not re.fullmatch(r"[A-Za-z0-9_-]+", label):
            parser.error("Each binary must be label=path, with a simple unique label")
        if label in [entry[0] for entry in binaries]:
            parser.error("Binary labels must be unique")
        binaries.append((label, Path(filename).resolve(strict=True)))
    output = args.output.resolve()
    output.mkdir(parents=True)  # Refuse existing batches before writing any results.
    results = {"platform": platform.platform(), "steps": args.steps, "warmup": args.warmup, "runs": []}
    if args.hz is not None:
        results["integration_hz"] = args.hz
        results["simulated_seconds"] = args.steps / args.hz
    hashes = {}
    # Interleave binaries to reduce bias from changes in temperature or host load.
    for threads in args.threads:
        for repetition in range(args.repeat):
            for label, binary in binaries:
                directory = output / f"{label}-t{threads}-r{repetition}"
                directory.mkdir()  # Never overwrite a previous measurement.
                command = [str(binary), "-s=Pyramid", "-q=Discrete", f"-t={threads}", f"-i={args.steps}", "-f"]
                if args.hz is not None:
                    command.append(f"-hz={args.hz}")
                completed = subprocess.run(command, cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                (directory / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
                if completed.returncode:
                    raise RuntimeError(f"{label} failed; see {directory / 'stdout.txt'}")
                match = re.search(r"^Discrete,\s*\d+,\s*([0-9.]+),\s*(0x[0-9a-fA-F]+)", completed.stdout, re.M)
                if not match:
                    raise RuntimeError(f"Missing result from {label}")
                with (directory / f"per_frame_discrete_th{threads}.csv").open() as stream:
                    timings = [float(row[1]) for row in list(csv.reader(stream))[1:]]
                if len(timings) != args.steps or any(not math.isfinite(t) or t <= 0 for t in timings):
                    raise RuntimeError(f"Invalid frame timings from {label}")
                state_hash = match[2]
                if label in hashes and hashes[label] != state_hash:
                    raise RuntimeError(f"Non-deterministic final state for {label}: {hashes[label]} vs {state_hash}")
                hashes[label] = state_hash
                measured = timings[args.warmup:]
                row = dict(binary=label, threads=threads, repetition=repetition, state_hash=state_hash,
                           all_steps_per_second=float(match[1]), steady_steps_per_second=1000 / statistics.mean(measured),
                           median_ms=statistics.median(measured), p95_ms=percentile(measured, 0.95),
                           p99_ms=percentile(measured, 0.99), max_ms=max(measured),
                           fraction_within_120_hz=sum(t <= 1000 / 120 for t in measured) / len(measured))
                row["fraction_within_150_hz"] = sum(t <= 1000 / 150 for t in measured) / len(measured)
                row["all_p95_ms"] = percentile(timings, 0.95)
                row["all_p99_ms"] = percentile(timings, 0.99)
                row["all_fraction_within_120_hz"] = sum(t <= 1000 / 120 for t in timings) / len(timings)
                row["all_fraction_within_150_hz"] = sum(t <= 1000 / 150 for t in timings) / len(timings)
                results["runs"].append(row)
                (output / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
                print(f"{label}, {threads} threads, {row['all_steps_per_second']:.1f} steps/s, steady p95 {row['p95_ms']:.2f} ms", flush=True)


if __name__ == "__main__":
    main()
