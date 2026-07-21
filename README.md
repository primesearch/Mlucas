[![Actions Status](https://github.com/primesearch/Mlucas/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/primesearch/Mlucas/actions/workflows/ci.yml)

# Mlucas
Ernst Mayer's Mlucas and Mfactor programs for GIMPS

[Ernst Mayer passed away unexpectedly](https://www.mersenneforum.org/showthread.php?t=28890) on September 10, 2023. This repository contains his posthumously released Mlucas v21 code, which is now maintained by the Great Internet Mersenne Prime Search (GIMPS) community. AutoPrimeNet (the Python PrimeNet program) previously bundled with Mlucas is now maintained in a [separate repository](https://github.com/tdulcet/AutoPrimeNet).

Mlucas and Mfactor are 100% open source programs. Mlucas is for [primality](https://en.wikipedia.org/wiki/Primality_test) and [P-1](https://en.wikipedia.org/wiki/Pollard%27s_p_%E2%88%92_1_algorithm) testing of [Mersenne](https://en.wikipedia.org/wiki/Mersenne_prime) and [Fermat](https://en.wikipedia.org/wiki/Fermat_number) numbers, including support for the [Lucas-Lehmer](https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test), [Probable prime](https://en.wikipedia.org/wiki/Probable_prime) (PRP) and [Pépin](https://en.wikipedia.org/wiki/P%C3%A9pin%27s_test) tests. Mfactor is for trial factoring. They support x86 Intel and AMD, ARM and other CPUs.

Feature | | Mlucas | Prime95/MPrime
--- | --- | ---: | ---:
**Architectures** | x86 | ✔️ | ✔️
\- | ARM | ✔️ | 
\- | Other | ✔️ | 
**Worktypes** | LL | ✔️ | ✔️
\- | PRP | ✔️ | ✔️
\- | P-1 | ✔️ | ✔️
\- | P+1 | | ✔️
\- | ECM | | ✔️
\- | Pépin | ✔️ | ✔️
**PRP** | Proofs | | ✔️
\- | Certs | | ✔️
**Error Checking** | Jacobi | | ✔️
\- | Gerbicz | ✔️ | ✔️
**Random Shifts** | | ✔️ | ✔️
**Interface** | CLI | ✔️ | MPrime only
\- | GUI | | Prime95 only
**Multiple Workers** | | Separate runs | ✔️
**PrimeNet Support** | | Separate program | ✔️
**Max FFT Length** | | 256M<br>(**512M** with 0 shift) | 32M (AVX) -<br>64M (AVX512)
**Largest Exponent** | | 4,294,967,231<br>(**8,937,021,911** with 0 shift) | 595,700,000 (AVX) -<br>1,169,000,000 (AVX512)
**Performance** | | ~50-90% | **100%**
**Free** 🆓 | | **Yes**, GPL | No, EULA
**100% Open Source** | | ✔️ | Mostly
**Claim Full EFF Awards** | | ✔️ | 

The original [Mlucas README](https://mersenneforum.org/mayer/README.html) is available for posterity and contains a lot of information, but note that it is no longer up to date; this document supersedes it and is kept current with each release. For more information about Mlucas, please see the [Mlucas subforum](https://www.mersenneforum.org/node/91) on the Mersenne Forum. The source code for and information about historical versions of Mlucas can be found on:

* [mersenneforum.org](https://mersenneforum.org/mayer/README.html)
* [hogranch.com](https://web.archive.org/web/20170606082351/http://hogranch.com/mayer/README.html)
* [mersenne.org](https://web.archive.org/web/20100612202354/http://mersenne.org:80/freesoft/mlucas.html)
* [SourceForge](https://sourceforge.net/p/mlucas/code/)

## Contents

* [Prerequisites](#prerequisites)
* [Windows Users](#windows-users)
* [Automatic method](#automatic-method)
* [Manual method](#manual-method)
  * [Download](#download)
  * [Installing GMP and HWLOC](#installing-gmp-and-hwloc)
  * [Build](#build)
  * [Common Build Issues and Workarounds](#common-build-issues-and-workarounds)
* [Basic build self-test](#basic-build-self-test)
* [Performance tuning](#performance-tuning)
* [Getting exponents from PrimeNet](#getting-exponents-from-primenet)
* [Sending results to PrimeNet](#sending-results-to-primenet)
* [Savefile format](#savefile-format)
* [Tracking your contribution](#tracking-your-contribution)
* [Further documentation](#further-documentation)

## Prerequisites

Dependencies:
* Make
* GNU C or Clang compiler
* \*GNU Multiple Precision (GMP) library
* \*Portable Hardware Locality (hwloc) library
* \*Python 3 (only needed to run [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet))

\* Optional

Use `which <tool>` (e.g. `which gcc`, `which gdb`) to check whether a given piece is already installed; most Linux distributions and macOS releases already include Make and a C compiler. `gdb` is a nice-to-have if you ever need to debug a crash (see [Common Build Issues](#common-build-issues-and-workarounds)).

GMP is used to take the GCDs that extract factors from P-1 residues; build without it (`no_gmp`, see [Build](#build)) and P-1 factoring still runs, but won't identify factors it finds. HWLOC is recommended starting with v21: it lets Mlucas resolve physical-core/thread topology itself via the new `-core` flag, removing the need to know your CPU vendor's logical-core numbering scheme (see [Basic build self-test](#basic-build-self-test)).

## Windows Users

Native Windows builds via [MSYS2](https://www.msys2.org/) are experimental and covered in [Download](#download) below. For everyday use, Windows users are better served by the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) (WSL), following the Linux instructions.

One consideration if you use WSL: WSL2 runs inside a lightweight VM, so Mlucas cannot see or use all of the host's RAM, and the RAM-usage percentage Mlucas reports/uses for `-maxalloc` will be based on the VM's memory, not the host's. WSL1 has no such VM overhead and is available on hardware that lacks the virtualization extensions WSL2 requires, but is otherwise less actively developed by Microsoft. If in doubt, benchmark both — self-tests (see [Basic build self-test](#basic-build-self-test)) make this quick. Also note that on WSL, `-cpu`/`-core` affinity pinning is not always effective (the process can still migrate between cores), so throughput may be somewhat lower than on native Linux.

## Automatic method

Linux users can use the [Mlucas install script](https://github.com/tdulcet/Distributed-Computing-Scripts#mlucas) to automatically download, build, set up and run Mlucas, including an auto-tuning pass that benchmarks instance/thread-count combinations to find the one that maximizes total throughput on your system, and downloading, setting up and running [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) for automated PrimeNet assignment management.

For the do-it-yourselfers, the manual build and tune flow follows below.

## Manual method

### Download

#### Linux

1. Verify that the dependencies above are installed. On Debian and Ubuntu, run: `sudo apt update` and `sudo apt install build-essential libgmp-dev libhwloc-dev`.
2. If one has git installed, just run: `git clone https://github.com/primesearch/Mlucas.git`. Otherwise, download the latest archive: `wget https://github.com/primesearch/Mlucas/archive/main.tar.gz` and then decompress the files: `tar -xzvf main.tar.gz`.
3. To download AutoPrimeNet, run: `wget -nv https://raw.github.com/tdulcet/AutoPrimeNet/main/autoprimenet.py`.

#### macOS

1. Verify that the dependencies above are installed. Run: `brew install gmp hwloc`.
2. If one has git installed, just run: `git clone https://github.com/primesearch/Mlucas.git`. Otherwise, download the latest archive: `curl -fLO https://github.com/primesearch/Mlucas/archive/main.tar.gz` and then decompress the files: `tar -xzvf main.tar.gz`.
3. To download AutoPrimeNet, run: `curl -sSfLO https://raw.github.com/tdulcet/AutoPrimeNet/main/autoprimenet.py`.

#### Windows

Native Windows builds are experimental. For now, Windows users should use the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) (WSL) and follow the [Linux](#linux) instructions above instead.

1. Download and install [MSYS2](https://www.msys2.org/).
2. Verify that the dependencies above are installed. With the MINGW64 environment, run: `pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-hwloc`.
3. If one has git installed, just run: `git clone https://github.com/primesearch/Mlucas.git`. Otherwise, download the latest archive: `wget https://github.com/primesearch/Mlucas/archive/main.tar.gz` and then decompress the files: `tar -xzvf main.tar.gz`.
4. To download AutoPrimeNet, run: `wget -nv https://raw.github.com/tdulcet/AutoPrimeNet/main/autoprimenet.py`.

### Installing GMP and HWLOC

If your package manager doesn't have `libgmp-dev`/`gmp-devel`/`hwloc` packages (see the OS-specific commands above), or you need a newer version than your distro ships:

**GMP** (mandatory for P-1 factor extraction): download the current `.tar.xz` from [gmplib.org](https://gmplib.org/), then:
```
tar xJf gmp-*.tar.xz
cd gmp-<version>
./configure
make
make check
make install   # re-run with sudo if you get a 'permission denied' error
```

**HWLOC** (recommended since v21, for the `-core` affinity flag and for visualizing hardware topology, especially on manycore and hybrid big.LITTLE-style processors): most distros package it as `hwloc`/`libhwloc-dev`/`hwloc-devel`. If your package manager can't find the header, look for a `-dev`/`-devel` variant of the package (e.g. Ubuntu's runtime `hwloc` package is separate from the `libhwloc-dev` package that provides `hwloc.h`). Once installed, the `lstopo` command-line tool (occasionally packaged as `lstopo-no-graphics`) can render your machine's topology, e.g. `lstopo my_machine.svg`, which is useful when choosing `-core`/`-cpu` arguments by hand.

### Build

1. Change into the `Mlucas` directory. Run: `cd Mlucas` or `cd Mlucas-main` depending on which method one used to download it.
2. Run:
	* To build Mlucas: `bash makemake.sh [SIMD mode] [use_hwloc] [no_gmp]`.
	* To build Mfactor: `bash makemake.sh mfac [word]`, where `word` is optionally one of `1word`, `2word`, `3word`, `4word` or `nword`.

	(The bracketed arguments above are all optional and order-independent — mix and match as needed, e.g. `bash makemake.sh mfac avx2 2word`.)

By default `makemake.sh` auto-detects the highest SIMD mode your build host supports (AVX-512, AVX2, AVX, SSE2 on x86_64; ASIMD on Armv8) by probing `/proc/cpuinfo` on Linux, `sysctl` on macOS, or compiling a small feature-detection program on other platforms, and falls back to a scalar-double build if none is recognized. To cross-compile or target a different mode than the host's, pass one of `avx512`, `avx2`, `avx`, `sse2` (x86_64), `asimd` (Armv8), `k1om` (first-generation Xeon Phi/Knights Corner), or `nosimd` explicitly as the first argument, e.g. `bash makemake.sh avx2`. Add `use_hwloc` to link against libhwloc and enable the `-core` flag, or `no_gmp` to build without linking GMP (P-1 factoring will then run, but skip taking GCDs to extract factors from the resulting residue).

Each explicit SIMD mode gets its own object directory (the binary itself is still just named `Mlucas`), e.g. `bash makemake.sh avx2` builds `obj_avx2/Mlucas`; letting the script auto-detect (a plain `bash makemake.sh`) builds `obj/Mlucas`. This lets you keep binaries for multiple SIMD modes side by side without them overwriting each other.

To build with Clang or another GCC-compatible compiler instead of the default GCC, run `export CC=<compiler>` (e.g. `export CC=clang`) before invoking `makemake.sh`. To switch compilers on an already-built tree without a full rebuild, run `make clean` inside the `obj` directory, then `CC=<compiler> make -j "$(nproc)"`.

### Common Build Issues and Workarounds

* **Where to look first:** if `makemake.sh` reports errors, check the `build.log` file it writes inside the object directory (`obj`, `obj_avx2`, etc.) for the (case-insensitive) string "error".
* **Debug symbols:** builds include `-g` by default so that a crash can be diagnosed with `gdb -ex=r ./Mlucas` (gives you the file, line number and stack trace). If you want a smaller binary afterward, run `strip -g Mlucas`.
* **Benign warnings:** `build.log` will typically contain compiler warnings about type-punned pointers, signed/unsigned comparisons, and unused variables — these come from the quad-float emulation and inline-assembly code and are expected; they are not build failures.
* **`cc: command not found` (some WSL installs):** `export CC=gcc` before retrying `bash makemake.sh`.
* **ARM/Apple Silicon SIMD detection:** `makemake.sh` probes for `-march=native`/`-mcpu=native` support before relying on it in the ASIMD fallback path, so most Clang-on-Arm toolchains are now handled automatically. If you still hit a runtime segfault or an "unsupported architecture" compiler error on an unusual ARM target, pass the closest matching mode explicitly (`bash makemake.sh asimd` or `bash makemake.sh nosimd`), or add an explicit `-march=<version>` for your CPU to your compiler flags.
* **Rebuilding after a source update:** run `make clean` in the relevant `obj*` directory (or `rm -f obj*/*.o` from the top-level directory) before re-running `makemake.sh`, since Make does not always detect stale header-file dependents. Parallel builds are fast, so a full rebuild is usually the safe and quick choice.

If you get stuck, please [open an issue](https://github.com/primesearch/Mlucas/issues) with your `build.log`, compiler version and platform details.

## Basic build self-test

Once you have a linked binary, try a quick spot-check at a small FFT length:
```
./Mlucas -fft 192 -iters 100 -radset 0
```
Look for a line reading `INFO: System has [X] available processor cores.` — this is the number of *logical* (virtual) cores; if it's double your physical core count, your CPU supports hyperthreading (or SMT). This particular test case has an internally-tabulated 100-iteration reference residue, so a mismatch means something is wrong with the build (compiler miscompilation, hardware issue, etc.) and will produce a loud error message.

If that works, retry multithreaded:
```
./Mlucas -fft 192 -iters 100 -radset 0 -nthread 2
```

### Setting thread count and CPU affinity

Mlucas offers three mutually-exclusive ways to pick thread count and affinity:

* **`-core {lo:hi[:threads_per_core]}`** — *recommended if you built with `use_hwloc`.* HWLOC resolves physical cores and their hardware threads for you, so `-core` behaves identically regardless of CPU vendor. `-core 0:3` runs one thread on each of physical cores 0 through 3; `-core 0:0:4` runs 4 threads on physical core 0 alone (handy for e.g. a Knights Corner-style core with 4 hardware threads).
* **`-cpu {lo[:hi[:incr]]}[,{lo[:hi[:incr]]},...]`** — the traditional, portable option, and still required if you did not build with hwloc support. It sets affinity to *logical* core indices directly, so the mapping from physical cores to logical-core numbers differs by vendor: on Intel, physical core *n* out of *N* total maps to logical cores *n* and *n+N* when hyperthreaded, so `-cpu 0:1` runs one thread on each of physical cores 0 and 1; on AMD, physical core *n*'s two hardware threads are the adjacent pair *2n*,*2n+1*, so the same-looking `-cpu 0:1` instead means "both threads of physical core 0 alone" — the AMD equivalent of "one thread on each of physical cores 0 and 1" is the stride-2 range `-cpu 0:3:2`. See [Performance tuning](#performance-tuning) for worked examples and vendor-specific benchmarking recipes.
* **`-nthread {int}`** — deprecated; equivalent to `-cpu 0:{int-1}`. Kept for backward compatibility only.

For a complete option reference (all self-test tiers, P-1 flags, `mlucas.ini` options, etc.), see [help.txt](help.txt).

## Performance tuning

After a successful basic self-test, run one of the following self-test tiers to build a `mlucas.cfg` file of optimal per-FFT-length radix sets and timings for your hardware. Run these under unloaded or constant-load conditions for the most reliable data:

```
./Mlucas -s m [-core or -cpu flags] >& test.log
```
`-s m` ("medium") covers the FFT-length range most GIMPS assignments currently fall into and is the tier ordinary users are recommended to run (roughly an hour on a fast CPU). Other tiers — `-s t`(iny), `-s s`(mall), `-s l`(arge), `-s a`(ll), `-s h`(uge), `-s e`(gregious) — cover progressively larger FFT-length ranges; see [help.txt](help.txt) section `[1]` for exact ranges and expected runtimes.

Once `mlucas.cfg` exists, Mlucas reads it at startup to pick the best radix set for whatever exponent it's assigned. If you plan to run multiple simultaneous instances (recommended on 4+ core systems, generally 2-4 physical cores per instance), copy `mlucas.cfg` into one run directory per instance, and give each instance a disjoint `-core`/`-cpu` range and its own `worktodo.txt`/`results.txt`. A minimal example on a hwloc build with an 8-physical-core CPU, running two 4-core instances:
```
mkdir run0 run1
cp obj/mlucas.cfg run0/ && cp obj/mlucas.cfg run1/
cd run0 && nohup ../obj/Mlucas -maxalloc 50 -core 0:3 &
cd ../run1 && nohup ../obj/Mlucas -maxalloc 50 -core 4:7 &
```
`-maxalloc {+int}` caps the percentage of free system RAM a given instance's P-1 stage 2 may use (default 90%, 50% on macOS); when running several instances concurrently, lower this so their stage 2 memory needs don't collide.

Without hwloc, the correct `-cpu` ranges for hyperthreaded/SMT systems differ between Intel and AMD, and hybrid (big.LITTLE-style) and Apple Silicon systems need their own topology-specific reasoning. [docs/PERFORMANCE.md](docs/PERFORMANCE.md) has a full set of worked `-cpu` examples for those cases, plus a template for a `jobs.sh` boot-time launch script.

## Getting exponents from PrimeNet

The recommended way to fetch and manage GIMPS assignments is [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet), the actively-maintained successor to the `primenet.py` script that used to ship with Mlucas. It registers your machine, requests work of your preferred type, and writes/tops-up the `worktodo.txt` file in each run directory you point it at (`--dir`/`-D`, one per Mlucas instance), reading Mlucas's own progress/checkpoint files to know when to submit results and fetch more work. See its README for full usage (`--username`/`-u`, `--mlucas`/`-m` to select the Mlucas client, `--timeout`/`-t` for the polling interval, `-t 0` for a single-shot run).

If you'd rather manage assignments by hand via the [PrimeNet manual testing pages](https://www.mersenne.org/), create a PrimeNet account, check out exponents from the [Manual Test Assignments](https://www.mersenne.org/manual_assignment/) page, then paste the returned assignment lines directly into `worktodo.txt` (one directory per Mlucas instance; **note that Mlucas v21's workfile is named `worktodo.txt`, not the old `worktodo.ini`** — rename any pre-v21 workfile before starting a v21 build). Mlucas supports these assignment line formats:

```
Test={aid},{exponent},{TF bits},{P-1 done? 0|1}                                   # LL, first-time or double-check
DoubleCheck={aid},{exponent},{TF bits},{P-1 done? 0|1}
PRP={aid},1,2,{exponent},-1,{TF bits},{ignore}                                    # PRP, first-time
PRP={aid},1,2,{exponent},-1,{TF bits},{ignore},{base},{residue type}              # PRP double-check
Pminus1=[aid,]1,2,{exponent},-1,{B1},{B2}[,{TF bits}][,{B2 start}][,"{known factors}"]
Pfactor=[aid,]1,2,{exponent},-1,{TF bits},{LL/PRP tests saved if factor found}
```
`{aid}` (the assignment ID) may be omitted or given as `n/a` (case-insensitive). When an LL/PRP assignment still needs P-1 work, Mlucas automatically splits it into a `Pminus1=`/`Pfactor=` assignment followed by the original, running P-1 first. Note that as of April 2021 the server no longer hands out first-time LL assignments (requests are converted to LL double-checks) — PRP with Gerbicz error-checking has superseded first-time LL testing. You can also LL-test an arbitrary prime exponent not obtained from the server by putting the bare number on its own line in `worktodo.txt`, or construct your own PRP assignment (`PRP=n/a,1,2,{exponent},-1,{TF bits},0`) for the same purpose with hardware-error detection.

## Sending results to PrimeNet

If you're using AutoPrimeNet, results are submitted automatically. For manual submission, log into your PrimeNet account and paste the new lines from your run directory's `results.txt` file into the [Manual Test Results Check-in](https://www.mersenne.org/manual_result/) page. Results are written in JSON format, e.g.:
```
{"status":"C", "exponent":81253819, "worktype":"PRP-3", "res64":"CE9AB357C704369B", "residue-type":1, "fft-length":4718592, "shift-count":66554884, "error-code":"00000000", "program":{"name":"Mlucas", "version":"21.0.2"}, "timestamp":"2026-01-01 03:42:39 UTC", "aid":"0CF485CD87F1BAC48B6B05D3A2094579"}
```
If you need more time to finish an assignment before its default 180-day expiry, use the [Manual Test Time Extension](https://www.mersenne.org/manual_extension/) page.

## Savefile format

Mlucas writes savefiles ("checkpoints") at regular intervals — every 10,000 iterations by default (100,000 for >4-threaded runs; overridable via `CheckInterval` in `mlucas.ini`, see [help.txt](help.txt) section `[11]`) — so a run can safely resume after an interrupt, crash, or power loss. Residues are stored in bytewise, minimum-size, endian-independent form.

* LL, PRP-test, and P-1 stage 1 residues are kept as a redundant pair, `p{exponent}` and `q{exponent}` (the second guards against the rare case the first is corrupted). PRP savefiles are roughly twice the size of LL/P-1 ones, since they also carry the Gerbicz error-check residue.
* On completion of P-1 stage 1, `p{exponent}` is renamed `p{exponent}.s1`, both so a following LL/PRP test of the same exponent doesn't clobber it, and so it can be reused later for a deeper stage 2 run (a "stage 2 continuation" — see [help.txt](help.txt) section `[9b]`) without redoing stage 1.
* P-1 stage 2 uses a single `p{exponent}.s2` savefile, plus a `p{exponent}.s1_prod` file caching the precomputed stage 1 prime-powers product (useful for restarting large-B1 runs without recomputing it).
* LL/PRP runs also snapshot a persistent savefile every 10M iterations (`.10M`, `.20M`, ...), allowing a partial rerun of a suspect range without restarting from scratch.

**Interrupting a run safely:** Mlucas listens for `SIGINT` (Ctrl-C) and `SIGTERM` (the default `kill` signal) and, for LL, PRP, and P-1 stage 1 work, will finish the current iteration, write savefiles, and exit cleanly. P-1 stage 2's state is more complex and does *not* currently have a clean-exit path — interrupting stage 2 loses any work since the last stage 2 savefile write. `killall -STOP Mlucas` / `killall -CONT Mlucas` will suspend/resume all running instances without losing anything. As a last resort for an unresponsive instance, `kill -9 {pid}` terminates immediately (and loses any not-yet-checkpointed work).

## Tracking your contribution

You can track your overall progress (both automated and manually-submitted work) on the [PrimeNet server's producer page](https://www.mersenne.org/report_top_500/) once you're registered. Note this does not include any pre-v5-PrimeNet-server manual test results.

## Further documentation

* [help.txt](help.txt) — the full command-line option reference (self-test tiers, P-1 flags, `mlucas.ini` options, savefile internals, and more).
* [docs/PERFORMANCE.md](docs/PERFORMANCE.md) — worked `-cpu`/`-core` examples for hyperthreaded Intel/AMD, ARM, hybrid, and Apple Silicon systems, plus multi-instance launch-script templates.
* [docs/Fermat-testing.md](docs/Fermat-testing.md) — testing Fermat numbers with Mlucas.
* [docs/ANDROID.md](docs/ANDROID.md) — running Mlucas on a cluster of Android phones.
* [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — how Mlucas's FFT/DWT-based modular squaring works internally (the algorithmic Q&A formerly on the mersenneforum.org README page).

## Contributing

Pull requests welcome! New to the codebase? [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) is a contributor-oriented tour of the main subsystems (FFT-based modular squaring, test-type dispatch, error checking, self-test/cfg, threading, factoring).

---

Portions of this document adapt material from Ernst W. Mayer's original Mlucas README (Copyright (C) 2021 Ernst W. Mayer), used and modified here under the terms of the [GNU Free Documentation License](https://www.gnu.org/licenses/fdl-1.3.html), Version 1.3 or any later version, with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts. The Mlucas and Mfactor source code itself is licensed under the GPL — see [LICENSE](LICENSE).
