# Performance tuning: worked examples

This page expands on the [Performance tuning](../README.md#performance-tuning) section of the main README with worked `-cpu`/`-core` examples for various CPU topologies, and a template multi-instance launch script.

If you built with `use_hwloc` (see [Build](../README.md#build)), skip straight to [Using `-core`](#using--core) — it removes the vendor-specific reasoning below entirely. The rest of this page is for builds without hwloc, or for understanding what `-core` is doing for you under the hood.

## Using `-core`

With an hwloc-enabled build, `-core {lo:hi[:threads_per_core]}` addresses *physical* cores directly, so the same invocation works identically on Intel, AMD, ARM, or anything else hwloc understands:

* `-core 0:3` — one thread on each of physical cores 0-3 (4 threads total).
* `-core 0:0:4` — 4 threads on physical core 0 alone (e.g. a core with 4 SMT hardware threads).
* `-core 4:7` — one thread on each of physical cores 4-7 — useful for a second instance that doesn't overlap the first.

Because vendor-specific logical-core numbering no longer matters, the rest of this document is only needed for `-cpu`-based (non-hwloc) setups.

## Using `-cpu` without hwloc

`-cpu` addresses *logical* core indices, and how those map to physical cores differs by CPU vendor when hyperthreading/SMT is in play:

* **Intel**: on an *n*-physical-core hyperthreaded CPU, physical core *k* maps to logical cores *k* and *n+k*. So `-cpu 0:1` runs one thread on each of physical cores 0 and 1 (no HT involved), and `-cpu 0,n` (physical core 0's two hardware threads) is what you'd use to double up on a single physical core.
* **AMD**: physical core *k*'s two SMT threads are the *adjacent* pair of logical indices: physical core *k* maps to logical cores *2k* and *2k+1*. So `-cpu 0:1` on AMD instead means "both hardware threads of physical core 0 alone" — to get one thread on each of physical cores 0 and 1 (the AMD equivalent of the Intel example above), use a stride-2 range: `-cpu 0:3:2` (logical cores 0 and 2).

Because the same-looking `-cpu` range means something different on the two vendors, benchmark both interpretations on unfamiliar hardware rather than assuming:

```
./Mlucas -s m -cpu 0:1 >& test0.log        # "2 physical cores, 1 thread each" on Intel / "1 physical core, both threads" on AMD
mv mlucas.cfg mlucas.cfg.a
./Mlucas -s m -cpu 0:3:2 >& test1.log      # the other interpretation
mv mlucas.cfg mlucas.cfg.b
```
Compare the `sec/iter`/`msec/iter` figures across FFT lengths in the two `.cfg` files and symlink the faster one to `mlucas.cfg`:
```
ln -s mlucas.cfg.a mlucas.cfg   # or mlucas.cfg.b, whichever won
```

### Non-hyperthreaded x86 (Intel or AMD)

No SMT ambiguity — just assign contiguous physical cores per instance. For a 12-core system running three 4-thread instances:
```
mkdir run0 run1 run2
cp obj/mlucas.cfg run0/ run1/ run2/
cd run0 && nohup ../obj/Mlucas -maxalloc 50 -cpu 0:3 &
cd ../run1 && nohup ../obj/Mlucas -maxalloc 50 -cpu 4:7 &
cd ../run2 && nohup ../obj/Mlucas -maxalloc 50 -cpu 8:11 &
```

### Hyperthreaded Intel, core count a multiple of 4

Benchmark 4-thread-per-instance vs. 8-thread-per-instance (4 physical cores + their hyperthreads) first:
```
./Mlucas -s m -cpu 0:3 >& test0.log            # 4 physical cores, no HT
mv mlucas.cfg mlucas.cfg.4c4t
./Mlucas -s m -cpu 0:3,12:15 >& test1.log      # same 4 cores + their HT siblings (n=12 physical cores here)
mv mlucas.cfg mlucas.cfg.4c8t
```
Then launch instances using whichever won, e.g. for a 12c/24t system where 4c/8t was faster:
```
mkdir run0 run1 run2
cp obj/mlucas.cfg run0/ run1/ run2/   # the winning mlucas.cfg.4c8t, renamed
cd run0 && nohup ../obj/Mlucas -maxalloc 50 -cpu 0:3,12:15 &
cd ../run1 && nohup ../obj/Mlucas -maxalloc 50 -cpu 4:7,16:19 &
cd ../run2 && nohup ../obj/Mlucas -maxalloc 50 -cpu 8:11,20:23 &
```

### Hyperthreaded AMD, core count a multiple of 4

AMD's adjacent-pair SMT numbering means the "no-HT" benchmark uses a stride-2 range instead:
```
./Mlucas -s m -cpu 0:7:2 >& test0.log     # 4 physical cores, one thread each
mv mlucas.cfg mlucas.cfg.4c4t
./Mlucas -s m -cpu 0:7 >& test1.log       # same 4 cores, both SMT threads each
mv mlucas.cfg mlucas.cfg.4c8t
```
For a 12c/24t AMD system where 4c/8t won:
```
mkdir run0 run1 run2
cp obj/mlucas.cfg run0/ run1/ run2/
cd run0 && nohup ../obj/Mlucas -maxalloc 50 -cpu 0:7 &
cd ../run1 && nohup ../obj/Mlucas -maxalloc 50 -cpu 8:15 &
cd ../run2 && nohup ../obj/Mlucas -maxalloc 50 -cpu 16:23 &
```

(The odd-multiple-of-2 core-count cases follow the same pattern with narrower ranges — benchmark both interpretations and use whichever gives faster `sec/iter` numbers.)

### ARM

No hyperthreading to worry about. On a quad-core Arm CPU, use all 4 physical cores: `-cpu 0:3` (or `-core 0:3` with hwloc). On multi-socket ARM systems, run one Mlucas instance per socket with a `-maxalloc` reflecting the number of instances, e.g. two 4-thread jobs on an octo-core dual-socket board: `-cpu 0:3` and `-cpu 4:7`.

### Hybrid / big.LITTLE systems (Alder Lake, Apple Silicon, Odroid N2, etc.)

Check `/proc/cpuinfo` or `lscpu -ap` (Linux) to see which logical/physical core indices belong to which cluster. For example, on an Odroid N2, cores 0-1 are the 'little' dual-core A53 cluster and cores 2-5 are the 'big' quad-core A72 cluster; running one job per cluster (`-cpu 0:1` and `-cpu 2:5`, each with its own self-test-derived `mlucas.cfg`) generally maximizes total throughput, thermal limits permitting.

Apple Silicon is a special case: there is currently no reliable way to pin an Mlucas instance to just the performance or just the efficiency cluster, and splitting a single instance's threads across both has historically performed *better* than confining it to one cluster on some M-series chips (e.g. `-cpu 0:7` outperforming `-cpu 0:3` on an M1). Benchmark both a single wide instance and a pair of narrower instances to see which wins on your machine, since this varies by chip generation.

## Finding optimal parameters for a single FFT length

If you need tuning data for one FFT length not covered by the standard self-test tiers:
```
./Mlucas -fft {n}[K|M] -iters {100|1000|10000} [-core or -cpu args]
```
This appends one line to `mlucas.cfg` for that length. See [help.txt](../help.txt) section `[2]` for the full FFT-length argument syntax.

## Wiring up AutoPrimeNet across multiple instances

Point [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) at all your run directories at once with repeated `--dir`/`-D` arguments (one per `run*` directory above) so it can keep every instance's `worktodo.txt` topped up and submit from every `results.txt`; see its own README for the exact invocation and worktype flags.
