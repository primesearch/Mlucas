# Changelog

All notable changes to Mlucas are documented here. Version numbers refer to
the `VERSION` string reported in the program banner (`src/Mlucas.c`).

## [Unreleased]

The current released version is **21.0.2**. This release is a large
correctness/robustness/portability campaign of ~35 PRs (#64-#108) with no
version bump committed yet, so it is tracked here as Unreleased pending a
version-number decision.

### Correctness fixes

These fix cases where Mlucas silently computed or returned a wrong value.

- **Fixed a corrupted carryout in the hand-tuned x86-64 assembly for
  `mi64_mul_scalar_add_vec2`** that made every completed Mersenne PRP test on
  affected builds abort at the final residue check
  (`Assertion 'rmodb == 0ull' failed`). The asm read back a `"=m"`
  (write-only) operand instead of using the value it had just computed,
  picking up stack garbage as a spurious carry out of the low limbs. (#80)
- **Fixed an exact-integer-vs-floating-point bug in `mi64_div_by_scalar64`**
  for the single-word case: the remainder was already computed exactly, but
  the quotient relied on a floating-point divide plus an off-by-one
  "nudge" that could not correct a result more than one unit off. Replaced
  with a plain hardware integer divide. (#107, fixes #94, supersedes the
  partial fix in #23/#18)
- **Fixed `twopmmodq64()`'s `p < 64` fast path**, which used the comparison
  operator `<` instead of the left-shift `<<` (`1ull < p` instead of
  `1ull << p`), so it returned 0 or 1 instead of 2^p mod q for every p in
  [1, 63]. (#84)
- **Fixed a `&&`-for-`&` typo in the `LEADZ96`/`LEADZ128`/`LEADZ192`/`QLEADZ`
  leading-zero-count macros** that silently produced a wrong count when the
  high limb of the input was zero. Not currently reachable from any live
  code path, but fixed as a latent trap for future callers. (#102)
- **Fixed the worktodo-line parser to accept PrimeNet's abbreviated
  PRP/PRPDC assignment format** (optional `how_far_factored,tests_saved` and
  `base,residue_type` fields, plus a trailing quoted known-factors list),
  which previously aborted with `Assertion 'Expected ',' not found after
  TF_BITS field...'` and could not reach the known-factors extraction path
  at all. (#87, fixes #30)

### Crash and robustness fixes

These fix genuine crashes (segfaults, aborted assertions) on otherwise-healthy
runs.

- **Declared the AVX-512 opmask (`k0`-`k7`) and high-vector (`zmm16`-`zmm31`)
  registers as clobbered in inline asm.** No asm block in the codebase
  declared these, even though 51 blocks across 16 files write them; under
  `-O3 -flto -march=native` the compiler could keep C-level values live in
  those same registers across the asm and have them silently destroyed. This
  is believed to be the root cause of the AVX-512 self-test segfaults
  reported in #65, #52, #16, and #12; a carry-loop index-arithmetic hardening
  is included as defense-in-depth for the same corruption. (#68)
- **Fixed an out-of-bounds read in the AVX/AVX-512
  `radix16_wrapper_square`/`radix32_wrapper_square` twiddle-index scheduling**
  for the smallest FFT lengths, which read past the end of the `index[]`
  array on the final partial chunk — a segfault under ASan, or silent
  roundoff-driven self-test failure otherwise. (#69, fixes #7)
- **Fixed executed undefined behavior**: an out-of-range left shift
  (`<< 64`) in `shift_word()`'s `mod64==0` case, and null-pointer-offset
  arithmetic in the RADIX_256 DFT's init-mode calls. (#78, fixes #28, #9)
- **Fixed a function-pointer type mismatch between the threadpool's
  `run_func` signature and the `mers`/`fermat`/`pm1` worker functions**,
  undefined behavior flagged by UBSan and `-Wcast-function-type` across 39
  files. (#79, fixes #29)
- **Fixed `threadpool_init()` crashing with SIGSEGV instead of failing
  cleanly** when a `pthread_create()` call fails partway through pool setup
  (e.g. under memory pressure): the cleanup path joined uninitialized
  `pthread_t` handles for threads that were never created. (#89)
- **Relaxed the Gerbicz-check carry-out assertion for exponents `p == 63
  (mod 64)`**, where a nonzero carryout out of `mi64_mul_scalar` is a normal,
  legal outcome (the top residue limb has only one spare bit) — the old
  strict `== 0` assert crashed roughly one in three Gerbicz checks for such
  exponents. (#77, fixes #33)
- **Fixed `convert_res_bytewise_FP()` rejecting a legal carry-propagation
  boundary case** (all-ones top digit plus an incoming carry, producing
  `MSW == 0, carry == 1`), which aborted otherwise-healthy oversized-FFT
  runs; also fixed the diagnostic printing one word past the actual data.
  (#86, fixes #19)
- **Bounded the `mlucas.cfg` preferred-FFT-length scan to the caller's
  documented 2x-default limit**, fixing an abort
  (`Call to get_preferred_fft_radix returns out-of-range FFT length`) when
  the fastest cfg entry for an exponent exceeded that bound; a stale/hand-edited
  cfg now degrades gracefully to a fresh timing self-test instead of crashing.
  (#76, fixes #36)
- **Fixed P-1 stage 2 error-code accumulation**, which summed multiple
  distinct `ERR_*` codes from a single batch into a meaningless composite
  value and then aborted on the resulting out-of-range code — killing an
  otherwise-recoverable run. Adopted the bitmask scheme already prescribed
  in a code comment. (#85, fixes #49)
- **Fixed P-1 stage 1 restarts crashing when resumed under a different
  memory budget (or a changed `-b1`)**: the RAM-derived auto-sizing of the
  stage 1 bound B1 was never persisted, so a restart could pick a different
  B1 than the in-progress powering used, corrupting the restart-file
  iteration-count invariant. B1 is now persisted and recovered across
  restarts. (#72, fixes #31)
- **Fixed a restart livelock**: the per-checkpoint FFT-length-reversion
  check could be trivially true even when no roundoff error had occurred,
  causing affected runs to re-select the same larger FFT length and
  re-read the just-written savefile every checkpoint, forever. The
  reversion attempt is now gated on an actual roundoff-error count. (#74,
  fixes #21)
- **Fixed several memory leaks**: worker threadpools were never freed on
  re-init, leaking `NTHREADS` threads every time the runlength/radix set
  changed (#70, fixes #44); the per-FFT-length index arrays allocated in
  `mers_mod_square()` leaked on every FFT-length change during e.g. the
  self-test sweep (#105, fixes #96).
- **Added checked heap-allocation wrappers** (`MALLOC`/`CALLOC`/`REALLOC`)
  that abort cleanly with a file/line/size diagnostic on a genuine
  allocation failure, instead of dereferencing NULL and segfaulting;
  converted roughly 170 allocation call sites to use them. (#88, fixes #51)
- **Hardening pass over `util.c`**: replaced strict-aliasing-violating
  pointer-cast float/uint64 reinterpretation with `memcpy` (could
  miscompile at `-O2` and above), fixed undefined behavior from calling
  `isdigit`/`isspace` on a possibly-signed `char`, and bounded a previously
  unbounded path-construction buffer. No behavior change. (#90)
- **Hardening pass over `Mlucas.c`**: command-line argument copies used
  `strncpy` without guaranteeing NUL-termination for over-length tokens
  (an out-of-bounds read on subsequent string comparisons), fixed via
  `snprintf`; also initialized previously-uninitialized Suyama residue
  variables that could otherwise print/use garbage. (#91)
- **Bounded diagnostic-buffer (`cbuf`) writes that format variable-length
  input** (worktodo/ini lines, cofactor GCD strings, known-factors tokens,
  factor-candidate strings, filenames) with `snprintf` instead of
  unbounded `sprintf`, closing a buffer-overflow risk without touching the
  much larger set of fixed-message/integer-only call sites. (#92)

### Behavior changes

- **1K FFT lengths are now rejected with a clear error instead of
  crashing.** `-fft 1` (and any request that resolves to a 1K FFT) previously
  aborted deep in `BIT_REVERSE_INT` (`product of radices [2] != vector
  length [1]`) because the wrapper_square scheme requires at least 3
  radices and 1K can only factor into 2. 1K has long been documented as
  unsupported and was already dropped from `config-fermat.sh`; 2K is now
  the smallest supported FFT length, and the self-test skips the retired
  entry instead of printing a spurious build-options warning. (#104,
  fixes #101)
- **Savefile writing on interrupt (SIGINT/SIGTERM) is re-enabled.** This was
  disabled in Dec 2021 because the old signal handler was not
  async-signal-safe — it called `fprintf`/`exit()` from signal context, and
  since signals were unblocked on every thread, delivery to an FFT worker
  already holding the malloc-arena or a stdio lock would deadlock the
  process indefinitely. The new handler only sets `volatile sig_atomic_t`
  flags, FFT worker threads block the relevant signals so only the main
  thread ever handles them, and the actual savefile write happens
  synchronously on the main thread's next safe point. (#108, fixes #100)

### Build & portability

- **Replaced hardcoded compiler/tool version-cutoff assumptions in
  `makemake.sh` with capability probes** (for `-fdiagnostics-color`,
  `-flto` via a genuine two-object link probe, etc.), and fixed a new build
  break on current glibc (2.43 / GCC 16) where `CPU_ZERO`/`CPU_SET`/
  `sched_setaffinity` failed to declare because in-file `#define
  _GNU_SOURCE` no longer suffices — `_GNU_SOURCE` is now defined on the
  compile command line. (#67, fixes #60)
- **Fixed `'asm' operand has impossible constraints'` failures** in
  hand-tuned inline asm that clobbers all 14 non-reserved x86-64 GPRs while
  also taking memory operands, reproducible on ancient compilers (GCC 5.4)
  under `-flto` and on ASan builds on modern compilers. (#71, part of #60)
- **Added GCC 15 support** (fixed incompatible-pointer-type errors) **and
  musl libc support**, with Alpine Linux added to CI. (#64)
- **Fixed the standalone Mfactor build**, which failed to compile on GCC
  14+/C23 due to `#ifdef USE_FMADD` branches calling a 100-bit modular
  exponentiation family that was never implemented (#81); a follow-up fixed
  an LTO link failure from hardcoded (non-asm-local) labels in
  `twopmodq*.c` inline asm colliding across LTO partitions once #81 let
  Mfactor compile far enough to reach the link stage (#82, depends on #81).
- **Fixed Windows build portability**: probe-gated AVX-512 in the Windows
  CI job (MinGW's assembler lacks the AVX-512 extended register names),
  fixed ARM-vs-x86 build-path detection on `windows-11-arm` (MSYS2 bash
  under emulation reports `$HOSTTYPE=x86_64`), and fixed Mfactor
  multiword+FMADD build issues on Windows. (#83)
- **Fixed Mlucas/Mfactor failing to start on Windows** with `'printf' is
  not recognized...`: `set_mlucas_path()` shell-expanded the default path
  via `popen("printf ...")`, which invokes `cmd.exe` (no `printf`) on
  Windows/MinGW. Windows now uses the path verbatim with no shell
  involved, preserving the same length and trailing-slash checks as the
  POSIX path. (#106, fixes #50)
- **CI improvements**: gate the `AVX512` build mode out of old-toolchain CI
  jobs when the assembler can't emit AVX-512 extended register names
  (#73); drop the `continue-on-error` mask on the "Linux old" CI job now
  that the underlying old-toolchain build breaks are fixed by #67/#71/#73
  (#75); guard AVX-512 *execution* (not just build) on CI runners that can
  build but not run it, and switch the ASan job to a NOSIMD build after
  it failed to assemble hand-tuned SIMD asm under `-fsanitize=address`
  (#103).

---

### Notes on scope

This changelog covers the in-flight bug-fix/portability campaign tracked as
PRs #64-#108 in `primesearch/Mlucas`. As of this writing, only #64 and #66
(a routine `actions/checkout` version bump, not separately itemized above)
are merged; the remainder are open pull requests pending review/merge. PR
#93 is a draft integration branch combining #67-#92 for combined CI
validation and is not itself a separate change.
