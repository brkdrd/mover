# Performance Test

The performance test application contains a couple of simple scenes to test performance of Jolt Physics. It will output the results to the TTY in CSV format.

## Commandline options

- -s=[scene]: This allows you to select a scene, [scene] can be;
    - Pyramid: A 4D pyramid of 1296 dynamic boxes, with sleeping disabled, to profile large island splitting. The upstream 3D scene contains 1240 boxes.
    - MaxBodies: The supported 4D body-count stress scene.
- -i=[iterations]: Number of physics steps before the test finishes.
- -hz=[frequency]: Physics integration frequency, a finite positive number (default 60). This changes the simulated timestep, not a wall-clock rate limit.
- -q=[quality]: This limits the motion quality types that the test will run on. By default it will test both. [quality] can be:
    - Discrete: Discrete collision detection
    - LinearCast: Linear cast continous collision detection
- -t=[num]: This sets the amount of threads the test will run on. By default it will test 1 .. number of virtual processors. Can be 'max' to run on as many thread as the CPU has.
- -no_sleep: Disable sleeping.
- -p: Outputs a profile snapshot every 100 iterations
- -r: Outputs a performance_test_[tag].jor file that contains a recording to be played back with JoltViewer
- -f: Outputs the time taken per frame to per_frame_[tag].csv
- -h: Displays a help text
- -rs: Record the simulation state in state_[tag].bin.
- -vs: Validate the recorded simulation state from state_[tag].bin. This will after every simulation step check that the state is the same as the recorded state and trigger a breakpoint if this is not the case. This is used to validate cross platform determinism.
- -repeat=[num]: Repeats all tests num times.
- -cache_distance=[meters]: Override the body-pair cache translation tolerance for experiments. The default physics setting is unchanged.
- -validate_hash=[hash]: Will validate that the hash of the simulation matches the supplied hash. Program terminates with return code 1 if it doesn't. Can be used to automatically validate determinism.

## Output

- Motion Quality: Shows the motion quality for the test.
- Thread Count: The amount of threads used for the test.
- Steps / Second: Average amount of physics steps / second over the entire duration of the test.
- Hash: A hash of all positions and rotations of the bodies at the end of the test. Can be used to verify that the test was deterministic.

## Results

For repeatable local comparisons, use `PerformanceTest/run_comparison.py`. It runs binaries sequentially, interleaves repetitions, retains per-frame timings and logs, and checks each binary's final-state hash across repetitions and thread counts:

```sh
python3 PerformanceTest/run_comparison.py \
  --binary 3d=build/perf-3d/PerformanceTest \
  --binary before=build/perf-results/PerformanceTest-baseline \
  --binary after=build/perf-off/PerformanceTest \
  --threads 1 6 12 --steps 600 --warmup 100 --repeat 3 \
  --output build/comparison
```

The binaries must use matching compiler/build settings. By default the runner retains the existing 1/60-second physics timestep. Its 120-Hz and 150-Hz budget statistics mean execution within 8.333 ms and 6.667 ms, respectively. To also change the integration timestep, pass `--hz 120` or `--hz 150`; every binary in that batch must support `-hz`. For an equal ten seconds of simulated time, use 600 updates at 60 Hz, 1200 at 120 Hz, or 1500 at 150 Hz. Scale the warmup count too (100, 200, or 250 updates). Changing the timestep changes the workload, so keep these results separate from an unchanged-timestep optimization comparison. See [the 4D optimization report](PerformanceOptimization.md) for the measured comparison, pipeline map, and limitations.

If you're interested in how Jolt scales with multiple CPUs and compares to other physics engines, take a look at [this document](https://jrouwe.nl/jolt/JoltPhysicsMulticoreScaling.pdf).
