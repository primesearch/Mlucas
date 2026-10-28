# Performance tuning: worked examples

This page expands on the [Performance tuning](../README.md#performance-tuning) section of the main README with worked `-cpu`/`-core` examples for various CPU topologies, and multi-instance launch examples.

Linux users: the [Mlucas install script](https://github.com/tdulcet/Distributed-Computing-Scripts#mlucas) automates everything on this page for you — it benchmarks instance/thread-count combinations and writes a `jobs.sh` launch script for the fastest one it finds, so most Linux users won't need to work through the manual steps below at all.

If you built with `use_hwloc` (see [Build](../README.md#build)), skip straight to [Using `-core`](#using--core) — it removes the logical-numbering reasoning below entirely. The rest of this page is for builds without hwloc, or for understanding what `-core` is doing for you under the hood.

## Using `-core`

With an hwloc-enabled build, `-core {lo:hi[:threads_per_core]}` addresses *physical* cores directly, so the same invocation works identically on Intel, AMD, ARM, or anything else hwloc understands:

* `-core 0:3` — one thread on each of physical cores 0-3 (4 threads total).
* `-core 0:0:4` — 4 threads on physical core 0 alone (e.g. a core with 4 SMT hardware threads).
* `-core 4:7` — one thread on each of physical cores 4-7 — useful for a second instance that doesn't overlap the first.

Because logical-core numbering no longer matters, the rest of this document is only needed for `-cpu`-based (non-hwloc) setups.

## Using `-cpu` without hwloc

`-cpu` passes *logical* CPU indices straight to the operating system's affinity call; Mlucas does not look at the physical topology. How those indices map to physical cores when hyperthreading/SMT is in play is decided by the OS and the firmware, not by the CPU vendor, so check it on your own machine before choosing ranges. On Linux, `lscpu -e` lists each logical CPU with its physical `CORE`, and `/sys/devices/system/cpu/cpu0/topology/thread_siblings_list` names the logical CPUs sharing core 0. Two layouts are common, on both Intel and AMD:

* **Split**: on an *n*-physical-core SMT CPU, physical core *k* maps to logical CPUs *k* and *n+k*. This is what Linux usually shows. `-cpu 0:1` runs one thread on each of physical cores 0 and 1, and `-cpu 0,n` puts two threads on physical core 0.
* **Adjacent**: physical core *k* maps to logical CPUs *2k* and *2k+1*. Windows numbers them this way, as do FreeBSD and some Linux systems. Here `-cpu 0:1` means both hardware threads of physical core 0, and one thread on each of physical cores 0 and 1 is the stride-2 range `-cpu 0:3:2`.

Other layouts exist, so treat these as examples. If you can't inspect the topology, a benchmark can tell the two common layouts apart (run from the object directory):

```sh
./Mlucas -s m -cpu 0:1 > test0.log 2>&1        # "2 physical cores, 1 thread each" if split / "1 physical core, both threads" if adjacent
mv mlucas.cfg mlucas.cfg.a
./Mlucas -s m -cpu 0:3:2 > test1.log 2>&1      # "2 physical cores, 1 thread each" under either layout
mv mlucas.cfg mlucas.cfg.b
```

Compare the `msec/iter` figures across FFT lengths in the two `.cfg` files: if they are about the same, the numbering is split; if `mlucas.cfg.b` is clearly faster, it is adjacent. Symlink the faster one to `mlucas.cfg`:

```sh
ln -s mlucas.cfg.a mlucas.cfg   # or mlucas.cfg.b, whichever won
```

### Non-hyperthreaded x86 (Intel or AMD)

No SMT ambiguity — just assign contiguous physical cores per instance. For a 12-core system running three 4-thread instances, from the top-level `Mlucas` directory (the one containing `obj`), with the self-test's `mlucas.cfg` in `obj`:

```sh
mkdir run0 run1 run2
for d in run0 run1 run2; do ln -s ../obj/mlucas.cfg "$d/"; done
(cd run0 && nohup ../obj/Mlucas -maxalloc 33 -cpu 0:3 &)
(cd run1 && nohup ../obj/Mlucas -maxalloc 33 -cpu 4:7 &)
(cd run2 && nohup ../obj/Mlucas -maxalloc 33 -cpu 8:11 &)
```

Each launch is in its own `( ... )` subshell so that its `cd` does not affect the next line. Substitute your object directory (e.g. `obj_avx2`) for `obj` if you built with an explicit SIMD mode.

### Hyperthreaded, split numbering, core count a multiple of 4

Benchmark 4-thread-per-instance vs. 8-thread-per-instance (4 physical cores + their hyperthreads) first:

```sh
./Mlucas -s m -cpu 0:3 > test0.log 2>&1            # 4 physical cores, no HT
mv mlucas.cfg mlucas.cfg.4c4t
./Mlucas -s m -cpu 0:3,12:15 > test1.log 2>&1      # same 4 cores + their HT siblings (n=12 physical cores here)
mv mlucas.cfg mlucas.cfg.4c8t
```

Then launch instances using whichever won, e.g. for a 12c/24t system where 4c/8t was faster (copy or rename the winning `mlucas.cfg.4c8t` to `obj/mlucas.cfg` first):

```sh
mkdir run0 run1 run2
for d in run0 run1 run2; do ln -s ../obj/mlucas.cfg "$d/"; done
(cd run0 && nohup ../obj/Mlucas -maxalloc 33 -cpu 0:3,12:15 &)
(cd run1 && nohup ../obj/Mlucas -maxalloc 33 -cpu 4:7,16:19 &)
(cd run2 && nohup ../obj/Mlucas -maxalloc 33 -cpu 8:11,20:23 &)
```

### Hyperthreaded, adjacent numbering, core count a multiple of 4

With adjacent-pair SMT numbering the "no-HT" benchmark uses a stride-2 range instead:

```sh
./Mlucas -s m -cpu 0:7:2 > test0.log 2>&1     # 4 physical cores, one thread each
mv mlucas.cfg mlucas.cfg.4c4t
./Mlucas -s m -cpu 0:7 > test1.log 2>&1       # same 4 cores, both SMT threads each
mv mlucas.cfg mlucas.cfg.4c8t
```

For a 12c/24t system with adjacent numbering where 4c/8t won:

```sh
mkdir run0 run1 run2
for d in run0 run1 run2; do ln -s ../obj/mlucas.cfg "$d/"; done
(cd run0 && nohup ../obj/Mlucas -maxalloc 33 -cpu 0:7 &)
(cd run1 && nohup ../obj/Mlucas -maxalloc 33 -cpu 8:15 &)
(cd run2 && nohup ../obj/Mlucas -maxalloc 33 -cpu 16:23 &)
```

(The odd-multiple-of-2 core-count cases follow the same pattern with narrower ranges — benchmark both interpretations and use whichever gives faster `sec/iter` numbers.)

### ARM

No hyperthreading to worry about. On a quad-core Arm CPU, use all 4 physical cores: `-cpu 0:3` (or `-core 0:3` with hwloc). On multi-socket ARM systems, run one Mlucas instance per socket with a `-maxalloc` reflecting the number of instances, e.g. two 4-thread jobs on an octo-core dual-socket board: `-cpu 0:3` and `-cpu 4:7`.

### Hybrid / big.LITTLE systems (Alder Lake, Apple Silicon, Odroid N2, etc.)

Check `/proc/cpuinfo` or `lscpu -ap` (Linux) to see which logical/physical core indices belong to which cluster. For example, on an Odroid N2, cores 0-1 are the 'little' dual-core A53 cluster and cores 2-5 are the 'big' quad-core A73 cluster; running one job per cluster (`-cpu 0:1` and `-cpu 2:5`, each with its own self-test-derived `mlucas.cfg`) generally maximizes total throughput, thermal limits permitting.

Apple Silicon is a special case: there is currently no reliable way to pin an Mlucas instance to just the performance or just the efficiency cluster, and splitting a single instance's threads across both has historically performed *better* than confining it to one cluster on some M-series chips (e.g. `-cpu 0:7` outperforming `-cpu 0:3` on an M1). Benchmark both a single wide instance and a pair of narrower instances to see which wins on your machine, since this varies by chip generation.

## Finding optimal parameters for a single FFT length

The standard self-test tiers already cover every supported FFT length, so this is not usually needed — but if you want tuning data for just one specific FFT length rather than a whole tier, run this from the object directory, adding the same `-core` or `-cpu` argument you will run with:

```sh
./Mlucas -fft {n}[K|M] -iters {100|1000|10000}
```

If at least half of the radix sets it tries at that length pass, this appends one line for that length to `mlucas.cfg` (creating the file if needed); otherwise it prints a `WARNING: ... radix-sets at FFT length ... passed - skipping it` message and writes nothing. See [help.txt](../help.txt) section `[2]` for the full FFT-length argument syntax.

## Wiring up AutoPrimeNet across multiple instances

Point [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) at all your run directories at once with repeated `--dir`/`-D` arguments (one per `run*` directory above) so it can keep every instance's `worktodo.txt` topped up and submit from every `results.txt`; see its own README for the exact invocation and worktype flags.
