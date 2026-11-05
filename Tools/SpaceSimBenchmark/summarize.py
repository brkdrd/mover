#!/usr/bin/env python3
"""Aggregate SpaceSimBenchmark per-run JSON files into a summary.json + markdown table.

Usage: summarize.py <results-dir>

Reads <dir>/*.rep1.json (the first repeat of each configuration) and writes
<dir>/summary.json. Percentiles are in milliseconds unless stated otherwise.
No values are adjusted or rounded for presentation; the JSON keeps full precision.
"""
import glob
import json
import os
import sys


def ms(v):
    return v / 1e6


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize.py <results-dir>", file=sys.stderr)
        return 2
    outdir = sys.argv[1]
    runs = []
    for path in sorted(glob.glob(os.path.join(outdir, "*.rep1.json"))):
        with open(path, "r", encoding="utf-8") as f:
            doc = json.load(f)
        m = doc["manifest"]
        t = doc["timing"]
        c = doc["checks"]
        runs.append(
            {
                "name": os.path.basename(path)[: -len(".rep1.json")],
                "scenario": m["scenario"],
                "bodies": m["bodies"],
                "steps": m["steps"],
                "worker_threads": m["worker_threads"],
                "job_system": m["job_system"],
                "simulated_seconds": m["simulated_seconds"],
                "samples": c["samples"],
                "compiler": m.get("compiler"),
                "build_type": m.get("build_type"),
                "harness_source_sha256": m.get("harness_source_sha256"),
                "core_source_sha256": m.get("core_source_sha256"),
                "build_sha256": m.get("build_sha256"),
                "git_head_at_configure": m.get("git_head_at_configure"),
                "timing_ns": t,
                "checks": c,
                "gates": doc["gates"],
                "state_hash": doc["state_hash"],
                "result": doc["result"],
            }
        )

    summary = {
        "source": "SpaceSimBenchmark raw per-run JSON",
        "note": "Single repeat shown per configuration; both repeats and the full per-step "
                "CSV are kept next to this file. Percentiles are nearest-rank over all steps.",
        "runs": runs,
    }
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    cols = ["run", "bodies", "thr", "p50 ms", "p95 ms", "p99 ms", "max ms",
            "steps/s", "setup ms", "wall ms", "miss>8.33", "E_total rel", "E_body rel",
            "spinL rel", "L rel", "p rel", "result"]
    print("| " + " | ".join(cols) + " |")
    print("|" + "|".join(["---"] * len(cols)) + "|")

    def rel(value, gated):
        return "{:.2e}".format(value) if gated else "n/a"

    for r in runs:
        t, c = r["timing_ns"], r["checks"]
        print(
            "| {name} | {bodies} | {thr} | {p50:.3f} | {p95:.3f} | {p99:.3f} | {mx:.3f} | "
            "{sps:.1f} | {setup:.3f} | {wall:.3f} | {miss} | {et:.2e} | {eb} | {l} | {lt:.2e} | {p:.2e} | {res} |".format(
                name=r["name"], bodies=r["bodies"], thr=r["worker_threads"],
                p50=ms(t["step_ns"]["p50"]), p95=ms(t["step_ns"]["p95"]),
                p99=ms(t["step_ns"]["p99"]), mx=ms(t["step_ns"]["max"]),
                sps=t["throughput_steps_per_second"], setup=ms(t["setup_ns"]), wall=ms(t["run_wall_ns"]),
                miss=t["deadline_misses"], et=c["total_energy_rel_drift"],
                eb=rel(c["energy_rel_drift"], c["energy_gated"]),
                l=rel(c["spin_angular_momentum_rel_drift"], c["spin_angular_momentum_gated"]),
                lt=c["total_angular_momentum_rel_drift"], p=c["total_linear_momentum_rel_drift"],
                res=r["result"],
            )
        )
    print()
    print("summary written to", os.path.join(outdir, "summary.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
