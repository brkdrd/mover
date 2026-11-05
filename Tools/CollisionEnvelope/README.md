# Collision-envelope acceptance diagnostic

This executable preserves the broader collision gate independently of the normal
regression suite. The recorded baseline failed; parent execution after repair
`81161f8b` passes all 18,432 queries with unchanged tolerances. Independent review
approved the repair and reproduced both baseline failure and repaired success.
The parent integrated native and whole-engine sanitizer suites also passed.
This bounded result is not universal collision/CCD proof.

It reproduces the independently reviewed 18,432-query sweep: 64 deterministic
rotations/distributions, six scales, six signed margins, four starting axes and
both argument orders. An analytical transformed-box corner reference verifies
its outward optimality conditions before judging the result. Material hit errors
are outside twice the caller's query tolerance; distance errors use that same
absolute budget. The narrow predicate repairs do not change this budget.

```sh
# From the repository root, GCC/WSL precise configuration:
cmake -S Tools/CollisionEnvelope -B out/collision-envelope -DCMAKE_BUILD_TYPE=Release
cmake --build out/collision-envelope --target CollisionEnvelope --parallel 2
./out/collision-envelope/CollisionEnvelope > out/collision-envelope/results.csv
```

Exit codes: **0** = no detected failures, **1** = collision acceptance failed,
**2** = invalid oracle/count or output-write failure. Preserve the exit code when
capturing logs. The executable writes all query rows to stdout and a summary to
stderr; it has no user-supplied numeric parameters.

Parent execution on the integrated `aee62ba5` core produced:

```text
queries=18432 material_misses=2154 material_false_hits=0 distance_failures=8733 oracle_failures=0 selected_scale_0.01_margin_0.001=382/512
```

The process exited **1**, as required. All 18,432 decisions and binary32 closest
squared distances matched the prior independent review's sweep exactly, after
normalizing CSV decimal encodings back to float bits. This is verification of
reproducibility. After the lower-dimensional projection repair, the same executable
source and budgets produced this parent-observed result, exit **0**:

```text
queries=18432 material_misses=0 material_false_hits=0 distance_failures=0 oracle_failures=0 selected_scale_0.01_margin_0.001=0/512
```

In the failed baseline, the selected scale-0.01, margin-0.001 population overlaps by roughly ten query
tolerances. Sample 0 with initial W axis and sphere first has analytical overlap
`1.00004253752e-5`, tolerance `9.99999997e-7`, a false intersection result and
closest squared distance `3.27894333907e-8`. Material misses also occurred at other
small scales; those failures were not acceptable contact ambiguity.

See [CollisionPredicateRepair.md](../../Docs/CollisionPredicateRepair.md) for the
repairs, and [MoverCompletion.md](../../Docs/MoverCompletion.md) for the remaining
release gates. A green bounded envelope and unit suite still do not establish
arbitrary-geometry CCD or full production workload acceptance.
