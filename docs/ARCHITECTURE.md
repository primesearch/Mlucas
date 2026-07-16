# How Mlucas works: a contributor's map

Mlucas is Ernst Mayer's life's work: dense, highly-optimized C code refined
over more than two decades, with comments as the primary documentation. The
authoritative references remain Ernst's own detailed
[README webpage](https://mersenneforum.org/mayer/README.html), the hundreds
of posts he and others have written on
[mersenneforum.org](https://mersenneforum.org/), the `help.txt` file shipped
with the source, and the peer-reviewed research papers describing the
underlying algorithms. This page is a *supplementary* guided tour of the main
subsystems and where to find them in the code, aimed at getting a new
contributor to a working mental model in about 20 minutes. It intentionally
stays conceptual + points-at-code; for the deep math, follow the cross-links
to `docs/` and to the source comments themselves, which are frequently the
most detailed place a given algorithm is explained. Please keep this document
consistent with Ernst's README webpage where they overlap — that page has
been frozen and predates the v21 release, so this file is the place to note
anything that has changed since.

Top-level data flow for a normal production run:

```
worktodo.txt line
    -> command-line/self-test dispatch                          [src/Mlucas.c: main()]
    -> ernstMain() in Mlucas.c: parse the worktodo assignment
         ("Test=", "DoubleCheck=", "PRP=", "Pminus1=", "Fermat,Test=", ...)
    -> pick FFT length + radix set (mlucas.cfg, or run a timed self-test)
    -> set up FFT tables, read/create a savefile
    -> iterate: mod-square (+ periodic Gerbicz-check update)     [ernstMain()'s main loop]
    -> every ITERS_BETWEEN_CHECKPOINTS iterations: write savefile
    -> on completion: compute final residue, Res64/Selfridge-Hurwitz checksums
    -> generate_JSON_report() -> results.txt
```

Note the direction of the call: `main()` (the command-line/self-test driver) is defined
*later* in the file than `ernstMain()`, and calls into it — but `ernstMain()` is the one
that actually opens and parses `worktodo.txt` (see §2). It's easy to assume the opposite
from the names alone.

## 1. The core computation: FFT-based modular squaring

Both primality testing (LL, PRP, Pépin) and P-1 factoring boil down to
*repeated modular squaring* of a big residue: `x = x^2 mod N` for a few
hundred million to a few billion iterations, where `N = 2^p - 1` (Mersenne)
or `N = 2^(2^m) + 1` (Fermat). Doing this with schoolbook or even GMP-style
bignum multiplication would be far too slow at these sizes, so Mlucas
represents the big integer as a vector of floating-point "digits" and
squares it via FFT-based convolution — this is where nearly all the CPU
time goes, and where nearly all the SIMD/threading complexity lives.

### Balanced-digit / IBDWT representation

The p-bit (or 2^m-bit) integer is stored as an array `a[]` of `n` doubles,
where `n` is the FFT length. Each digit holds a variable number of bits
so that `n` doesn't have to evenly divide `p`: most digits hold
`bits_small = p/n` bits, and the rest hold `bits_small + 1` bits, arranged
so the total is exactly `p` bits (`mers_mod_square.c`, `base[]`/
`baseinv[]`). This is Crandall & Fagin's **irrational-base discrete weighted
transform (IBDWT)**: instead of a fixed power-of-2 base, each digit position
`j` has its own "weight" `w[j]` (arrays `wt0[]`/`wt1[]`, built in
`mers_mod_square.c`) that makes a plain cyclic convolution of the
weighted digits equal to multiplication mod `N`. Digits are signed
("balanced", roughly in `[-base/2, base/2)`) to keep rounding error
symmetric and minimize the maximum digit magnitude the FFT has to handle
in double precision.

### One call = one fused forward-FFT / square / inverse-FFT / carry pass

`mers_mod_square()` (`src/mers_mod_square.c`) is the Mersenne-mod entry
point; `fermat_mod_square()` (`src/fermat_mod_square.c`) is its Fermat-mod
counterpart. Both have near-identical header comments describing the same
4-step fused pipeline per call:

1. Passes 2..S-1 of a complex decimation-in-frequency (DIF) forward FFT.
2. The final forward-FFT radix pass, fused with the **pointwise squaring**
   step and the *first* radix pass of the inverse FFT, all in one memory
   pass. For Mersenne this is the "wrapper/square" step (see below); for
   Fermat it's a plain complex pointwise square.
3. Passes 2..S-1 of the complex decimation-in-time (DIT) inverse FFT, radices
   in reverse order.
4. The final inverse-FFT radix pass, fused with inverse-DWT-unweighting,
   **carry propagation** (renormalizing digits back into balanced-digit
   range, which is also where **round-off error is measured**), forward-DWT
   reweighting, and the first forward-FFT radix pass for the *next*
   iteration — so consecutive squarings never materialize a plain
   "integer" form of the residue in between.

Which radix-specific pass function actually runs is chosen at setup time
based on the global `RADIX_VEC[]`/`NRADICES` (see FFT length section below);
carry-propagation-with-ROE-check is done by generated per-radix functions
like `radix5_ditN_cy_dif1()`, `radix16_ditN_cy_dif1()`, etc. (dispatched
from within `mers_mod_square.c`).

**Round-off error (ROE):** each call returns a fractional round-off measure
in a local `fracmax`; the driver in `Mlucas.c` tracks running average/max
(`AME`/`MME`) and the globals `ROE_ITER`/`ROE_VAL` record the worst iteration
seen (`mers_mod_square.c`). A `fracmax` that hits exactly `0.4375`
is treated as a "dangerous" threshold that triggers retry logic (shorter
carry chains, or a Gerbicz-check rollback) rather than a hard abort — see
[§3 Error checking](#3-error-checking).

### Mersenne vs. Fermat: cyclic vs. negacyclic convolution

Squaring mod `2^p - 1` needs an ordinary **cyclic** convolution, computed via
a real-input FFT with a "wrapper/square" step: the last forward-FFT radix
pass and the pointwise-squaring math (`pair_square()`,
`src/pair_square.c`, implementing Crandall/Fagin's Hermitian-pair-squaring
identity `I[j] = H[j]^2 + {1+e^{2πij/N}}·{H[j]-H~[N-j]}^2/4`) are done inside
the generated `radix{16,32,64,...}_wrapper_square.c` routines that
`mers_mod_square()` calls for the final pass.

Squaring mod `2^n + 1` instead needs an **acyclic (negacyclic)**
convolution. Fermat's code gets this "for free" by premultiplying with the
`2N`-th roots of unity `rn(j) = e^{ijπ/N}` (comment in
`fermat_mod_square.c`) and exploiting `rn(j+N/2) = i·rn(j)` to fold
pairs `(j, j+N/2)` of the real input into a **right-angle transform**: a
plain length-N/2 complex FFT with no separate wrapper step. For power-of-2
FFT lengths that's the whole story (fast, good up to roughly F35); for
non-power-of-2 lengths this negacyclic weighting is combined with the same
IBDWT machinery as the Mersenne path (explicit "cyclic → acyclic"
weighting/unweighting in `fermat_mod_square.c`).

### FFT length and radix sets

`get_fft_radices()` (`src/get_fft_radices.c`) maps an FFT length
(expressed in "K doubles", e.g. `kblocks=4` → 4096-double transform) plus a
`radix_set` index to the ordered list of small-prime FFT radices whose
product is the complex transform length — e.g. `{12, 8, 16}`. This is a big
hand-tuned lookup table: leading radix ∈ roughly `{5..16}` (and constrained
further by SIMD build: AVX-512 needs the leading power-of-2 factor ≥ 8, AVX/
SSE2 need it divisible by 4), trailing radix ∈ `{16, 32}` (must match an
implemented `wrapper_square` routine), middle radices powers of 2. Multiple
radix sets are often available at the same FFT length — that's what the
self-test/timing step picks between (§4).

`get_default_fft_length(p)` (`get_fft_radices.c`) picks the smallest
supported FFT length whose `given_N_get_maxP(fftLen)`
(also in `get_fft_radices.c`) is ≥ the exponent being tested —
`given_N_get_maxP` encodes, via a closed-form bits-per-word budget (see the
PARI/GP script in a comment in `Mlucas.c`), the largest exponent that
keeps expected round-off error around the ~0.25 target for that transform
length. This is the concrete link between "bits per word" and "which FFT
length can safely test exponent p".

### Threading

Multithreaded builds run the middle FFT passes (steps 1 and 3 above) as
parallel chunks dispatched through a small custom threadpool,
`src/threadpool.c`. Public API: `threadpool_init()`,
`threadpool_add_task()` (can block until the task set completes —
this pool processes one batch of "chunks" of an FFT pass at a time, not a
long-lived async queue), `threadpool_free()` (all in `threadpool.c`). `mers_mod_square()`/
`fermat_mod_square()` each define a per-call thread-data struct
(`mers_thread_data_t`, `ferm_thread_data_t`) carrying the roots-of-unity
tables, block indices, etc. each worker needs, and dispatch chunks via
`mers_process_chunk()`/`fermat_process_chunk()`.

CPU affinity is handled in the same file: if built with `hwloc`, threads are
pinned via `hwloc_set_cpubind()` using logical-CPU
objects looked up from the runtime topology; without hwloc it falls back to
raw `sched_setaffinity()`/`CPU_SET` on Linux, `cpuset_setaffinity()` on
FreeBSD, or is left to the OS on OpenBSD (all in `threadpool.c`). hwloc is
optional at build time per the top-level README.

## 2. Test types and the top-level driver

`ernstMain()` in `src/Mlucas.c` is *the* driver. It is
not a thin dispatcher — worktodo parsing, FFT setup, the main iteration loop,
checkpoint I/O and final-result reporting all live inline in this one
function. `main()` (defined *after* `ernstMain()` in the
file) is a separate, higher-level entry point that parses command-line flags
(`-s`, `-fft`, `-cpu`/`-core`, etc. — see `help.txt`) and self-test tables
(`MersVec[]`/`MvecPRP[]`/`FermVec[]`, §4), then calls `ernstMain()` once per
assignment or self-test case. In default (no-flags) invocation, `main()`
calls `ernstMain()` with no exponent, and it's `ernstMain()` itself — at the
`read_next_assignment:` label — that opens `worktodo.txt` and
parses the first eligible line. Recognized prefixes (Prime95/PrimeNet
assignment syntax; parsing detail in a long comment in `Mlucas.c`):

| Worktodo prefix | Meaning |
|---|---|
| `Test=`, `DoubleCheck=` | Mersenne Lucas-Lehmer test |
| `PRP=` (`PRPDC=`) | Mersenne/Fermat PRP; append known factors for PRP-CF |
| `Fermat,Test=<m>` | Pépin test of F_m |
| `Pminus1=` | P-1 factoring, explicit stage bounds |
| `Pfactor=` | P-1 factoring, bounds chosen from TF depth / tests-saved |
| `Factor=` (`#if INCLUDE_TF`) | Trial factoring |

Assignment/test identification uses two small enums in `Mdata.h`:

- `MODULUS_TYPE_MERSENNE` / `MODULUS_TYPE_FERMAT` (plus a few rarer ones:
  `MERSMERS`, `GENFFTMUL`) — in `Mdata.h` — what kind of number is
  being tested.
- `TEST_TYPE_PRIMALITY` (LL for Mersenne, **Pépin** for Fermat — same
  test-type constant, distinguished by `MODULUS_TYPE`), `TEST_TYPE_PRP`
  (also covers PRP-CF, cofactor PRP — see below), `TEST_TYPE_PM1`,
  `TEST_TYPE_TF`, `TEST_TYPE_ECM` — also in `Mdata.h`.

Inside `ernstMain()`, a function pointer is set once per call —
`func_mod_square = (MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? mers_mod_square
: fermat_mod_square` — and the main iteration loop
(a `for(;;)` loop, with the squaring calls inside it) just calls
through that pointer `ITERS_BETWEEN_CHECKPOINTS` times per outer-loop pass,
without needing to know which modulus type it's dealing with.

- **LL / Pépin** (`TEST_TYPE_PRIMALITY`): straightforward repeated squaring
  starting from a small seed (3 for Pépin), `p-2` (LL) or `p` (Pépin)
  iterations, final residue compared against `-1`/`+1` mod N.
- **PRP** (`TEST_TYPE_PRP`): conceptually the Fermat probable-prime test
  `PRP_BASE^(N-1) ≡ 1`, but the Gerbicz check (§3) needs a *pure* sequence
  of squarings, so what's actually computed is the reformulated identity
  `b^(N+1) ≡ b^2 (mod N)` — see `help.txt` §7 for why, and the residue-type
  discussion in `Mlucas.c` for how the internally-computed "type 3"
  residue gets converted back to the standard "type 1" residue before being
  reported to PrimeNet. **PRP-CF** (cofactor PRP — testing `N/(known
  factors)` rather than `N` itself) reuses the exact same iteration path;
  the difference is a `known_factors` list parsed from the worktodo line
  and a final mod-division/GCD step against those factors folded into the
  final-residue handling via `Suyama_CF_PRP()` (in `Mlucas.c`) — it's a
  flag/post-processing difference, not a separate driver path. This GCD
  step (like the P-1 GCD below) requires Mlucas to have been built with
  the optional GMP library (`INCLUDE_GMP`); without GMP at build time no
  GCD is taken and no factor/cofactor-primality determination is possible.
- **P-1 factoring** (`TEST_TYPE_PM1`): stage 1 (a left-to-right binary
  powering over the product of small-prime-powers up to bound B1) is
  handled inline in `ernstMain()`'s same squaring loop, then a GCD checks
  for a factor — again requiring the optional GMP library at build time
  (`INCLUDE_GMP`; without it, Mlucas warns that no GCD will be done). If
  none is found and the assignment calls for it, **stage 2** hands off to
  a separate file, `src/pm1.c` (`pm1_stage2()`, a big-step/baby-step
  extension over primes in `(B1, B2]`) — which reuses the very same
  `func_mod_square` pointer (asserted explicitly in `pm1.c`). See
  [`docs/pm1.txt`](pm1.txt) for the underlying number theory and
  worked/known-factor examples (if `p-1` is `B`-smooth, `a^(p-1) mod N`
  reveals a nontrivial factor via GCD). Note: the Brent-Suyama extension
  to stage 2 (catching factors with one additional large prime factor
  cheaply) is *not* (yet) implemented in Mlucas — [`docs/brent-suyama.txt`](brent-suyama.txt)
  is reference/background material, not a description of shipped code.
  (Don't confuse this with `Suyama_CF_PRP()` above, the *Suyama
  cofactor-PRP* test, which is implemented and is a different thing.)
  See `help.txt` §9 for the operational `-maxalloc`/`-pm1_s2_nbuf`
  stage-2-memory flags and the `Pminus1=...,B2_start[,known_factors]`
  stage-2-continuation assignment syntax. Big-integer GCD/modexp needed
  here comes from `mi64.c` (§5).

**Checkpoints:** `write_ppm1_savefiles()` / `read_ppm1_savefiles()`
(both in `Mlucas.c`) persist the residue array, iteration count,
Res64/Selfridge-Hurwitz checksums, and — for PRP runs — the second
Gerbicz-checkproduct array and its own checksum triplet, every
`ITERS_BETWEEN_CHECKPOINTS` iterations.

**Final result:** on completion, `generate_JSON_report()` (declared in
`Mlucas.h`, called from `Mlucas.c`) builds the JSON blob written to
`results.txt` (`OFILE`), including the primality/PRP verdict, Res64
(and, for PRP tests, Res2048), and (for P-1) any factor found. Only
Res64/Res2048 go into the JSON line — the Selfridge-Hurwitz residues
(Res35m1/Res36m1, see §3) are internal error-check values, not part of
the JSON result, and are not typically consumed by other GIMPS client
software. Real examples, quoted in the source's own reference comment
in `Mlucas.c`:

```
{"status":"C", "exponent":86749043, "worktype":"LL", "res64":"9EEA7CAD97A07648",
 "fft-length":4718592, "shift-count":6030412, "error-code":"00000000", ...}

{"exponent":"102973951", "worktype":"PM1", "status":"F", "fft-length":5767168,
 "B1":1000000, "factors":["470377562071431809697977"], ...}
```

## 3. Error checking

Two independent, complementary layers protect against the transform
silently producing a wrong residue (hardware bit-flip, FFT-length-too-small
round-off, thread bug, etc.):

**Round-off error (ROE).** Every mod-square call measures how far each
carry-propagated digit landed from the nearest integer before rounding —
see §1. `Mlucas.c` accumulates average (`AME`) and max (`MME`) ROE per
checkpoint interval and prints them alongside Res64 in the progress line.
A single very-bad ROE doesn't necessarily kill the
run — it can trigger a switch to a shorter/safer carry chain
(`USE_SHORT_CY_CHAIN`) or a Gerbicz-check-driven rollback; a *pattern* of
bad ROEs at a given FFT length is a hint that the length is too small for
the exponent (`Mlucas.c` compares against
`get_default_fft_length`/`given_N_get_maxP`).

**Gerbicz check** (PRP runs, and Pépin runs when memory allows — see
`DO_GCHECK` in `Mlucas.c`), full math in
[`docs/gerbicz.txt`](gerbicz.txt). Idea in plain terms: alongside the
primary residue `u(t) = a^(2^t)`, maintain a second "checkproduct" `d`,
updated cheaply every `L` iterations as `d(t+1) = d(t) · u((t+1)·L)`. There's
a second, more expensive way to compute the same value:
`d(t+1) = u(0) · d(t)^(2^L)`. The two must agree; checking this identity
every `L^2` iterations (default: update every ~100-1000 iterations, verify
every ~10^6) catches almost any computational error in the intervening
squarings at roughly 0.1% time overhead, because it's *algebraically*
certain the two formulas agree absent an error — no separate trusted
reference residue is needed. In the code, the primary residue and
checkproduct are tracked in separate arrays (`a[]`/`b[]`/`c[]` in
`ernstMain()`, updated via extra `func_mod_square()` calls with special
`fwd_fft_only` mode bits — see the big comment in `mers_mod_square.c`
for the full state machine); on mismatch, `ierr = ERR_GERBICZ_CHECK`
(in `Mlucas.c`) triggers a rollback to the last-known-good checkpoint
rather than a fatal abort.

**Res64 / Selfridge-Hurwitz residues.** `Res64` is the low 64 bits of the
residue; `Res35m1`/`Res36m1` are the residue mod `2^35-1` and mod `2^36-1`
(cheap extra checksums, computed via `mi64_div_by_scalar64`-style code in
`Mlucas.c`). These are internal error-check values: they're compared
against known-good reference residues during self-test (§4) and printed in
progress/checkpoint output and stored in savefiles, but — as noted above —
only Res64 (and Res2048 for PRP) actually goes into the JSON result; the
Selfridge-Hurwitz residues are not typically transmitted to or used by
other GIMPS software. Mlucas does **not** implement a
separate Jacobi-symbol check (see the feature-comparison table in the
top-level README); Gerbicz is its only in-run correctness check.

## 4. Self-test & the `mlucas.cfg` system

Running with `-s <level>` (see `help.txt` for the self-test level presets,
e.g. `tiny`/`small`/`medium`/`large`/`all` — these are examples, not the
complete set Mlucas supports) drives `main()` through a battery of known
`(exponent, FFT length)` pairs from the `MersVec[]` / `MvecPRP[]` /
`FermVec[]` tables (all in `Mlucas.c`). Each table
entry stores an FFT length, an exponent (or Fermat index), and *reference*
Res64/mod-2^35-1/mod-2^36-1 residue triplets for 100/1000/10000-iteration
runs — these are known-correct values Ernst Mayer precomputed once and
checked in, so a self-test run doesn't need network access or an oracle.

For each FFT length, every available radix set (`get_fft_radices()`, §1) is
timed in turn; a radix set only "passes" if its residue matches the
reference triplet. Among the passing sets, the
fastest by wall-clock time is written to `mlucas.cfg` (or `fermat.cfg` for
Fermat runs — `CONFIGFILE`, set in `Mlucas.c`) as one line per FFT
length: iteration time, average/max ROE observed, and the actual radix
vector (not just an index, so the file stays meaningful even if the radix
tables in `get_fft_radices.c` change later) — see the format written in
`Mlucas.c`.

On a normal (non-self-test) run, `get_preferred_fft_radix(kblocks)`
(consulted in `Mlucas.c`) looks up the cfg file for the exponent's
default FFT length; if no entry exists yet, Mlucas asks the user to run a
self-test at that length first rather than guessing.

## 5. Factoring: trial factoring (Mfactor) and the `mi64` library

Before spending months on a full LL/PRP test, GIMPS first does cheap trial
factoring (TF): looking for small factors of `2^p - 1` directly, which (if
found) rule out the need for the expensive test entirely. This is
`src/factor.c` (~4800 lines).

**The math.** Any factor `q` of `2^p - 1` (p prime) must have the form
`q = 2·k·p + 1` and must satisfy `q ≡ ±1 (mod 8)` (comment in
`factor.c`, derived from `2` being a quadratic residue mod any such
factor). Combined with `p mod 60` / `k mod 60` residue-class constraints
that eliminate candidates divisible by 3 or 5, this cuts the candidate-`k`
search space to roughly 1/8 before any primality-style test is even run.
Surviving candidates `q` are tested by computing `2^p mod q` via
left-to-right binary exponentiation (the `twopmodq*` family, e.g.
`twopmodq96_q4`, `twopmodq78_3WORD_DOUBLE_q8`, and for full-multiword
candidates `mi64_twopmodq()` in `mi64.c`) — `q` divides `2^p-1` iff the
result is 1. Search ranges are expressed as bit levels (e.g. "TF from 2^64
to 2^65"); there's also GPU-accelerated sieving/candidate-testing in
`src/gpu_sieve.cu`/`src/gpu_iface.cu` (CUDA, optional build) mirroring the
CPU sieve logic.

**`mi64.c`/`mi64.h`** (~8800 lines) is the general multi-precision integer
library underlying all of this: arrays of `uint64` limbs with add/sub
(`mi64_add`), multiply, division (`mi64_div`,
`mi64_div_by_scalar64`), modular exponentiation
(`mi64_modpow_lr`, `mi64_twopmodq`/`mi64_twopmodq_qmmp`), and the
bit/scalar utilities the rest of the codebase
leans on for anything that doesn't fit in a native 64/128-bit integer.
Beyond trial factoring, this same library backs P-1 stage 2 arithmetic and
the GCD step used by cofactor-PRP (which, as noted in §2, requires GMP at
build time).

## 6. Savefiles / restart

Checkpointing is handled by the `write_ppm1_savefiles()`/`read_ppm1_savefiles()`
pair (in `Mlucas.c`), written every `ITERS_BETWEEN_CHECKPOINTS`
iterations (default depends on thread count; overridable via `CheckInterval=`
in `mlucas.ini` — `help.txt` §11). The on-disk format is documented in a long
comment immediately preceding these functions (search
`"Set of functions to Read/Write full-length residue data"` in `Mlucas.c`)
that traces back to Richard Crandall's original 1998 F24 "Pepin residue
file" definition. Fields, in order: a leading test-type byte and
modulus-type byte; an 8-byte iteration/squaring count; the residue itself
(`ceil(log2 N / 8)` bytes, endian-independent, renormalized from
balanced-digit to standard non-negative form, and *un*-shifted — see
`-shift`, `help.txt` §6, for the per-iteration residue-shift scheme and why
the shift count is stored and re-applied separately rather than baked into
the stored residue); the `Res64`/`Res35m1`/`Res36m1` checksums (re-verified
on read); the FFT length in use at
write time (lets a restart detect an earlier ROE-triggered auto-escalation
to a larger FFT length); and, for PRP runs only, a *second* residue array
holding the current Gerbicz checkproduct plus its own checksum triplet, so
a Gerbicz-check failure (§3) can roll back to exactly the last point both
arrays were known-consistent, plus two small cumulative-error counters
(ROE-warning count, Gerbicz-check-failure count). Note: the long
header comment preceding these functions (inherited from Crandall's
original 1998 definition) also describes a variable-length ASCII tail
field for "anything else needed to characterize the run" — in the
current Mlucas savefile format this tail is **not** actually written;
treat that part of the comment as stale/aspirational rather than a
description of the on-disk format. tdulcet's GIMPS status script,
[`gimps_status.py`](https://github.com/tdulcet/AutoPrimeNet/blob/main/gimps_status.py),
can parse and validate Mlucas savefiles and is a useful independent
reference/debugging tool if you need to check your understanding of the
format against a second implementation.

**Naming convention** (`help.txt` §12): for M(p), the checkpoint is kept as
a redundant pair, `p<p>` and `q<p>` (so an interrupted write to one doesn't
lose the other). At the end of a P-1 run, the primary file is renamed
`p<p>.s1` so a later LL/PRP run of the same exponent doesn't overwrite the
P-1 state, and so a deeper stage-2 continuation can reuse it later
(`help.txt` §9b). Stage 2 keeps its own `p<p>.s2` file, and a `p<p>.s1_prod`
file caches the (potentially expensive to recompute) stage-1
small-prime-powers product for large B1. On restart, `ernstMain()`
re-derives the FFT setup from the savefile's stored FFT length/exponent and
resumes the iteration loop at the saved iteration count.

## Further reading

- [`docs/gerbicz.txt`](gerbicz.txt) — Gerbicz check math in full.
- [`docs/Fermat-testing.md`](Fermat-testing.md) — Fermat-number-specific
  build/self-test/cfg notes (SIMD requirements, FFT length limits).
- [`docs/pm1.txt`](pm1.txt) — P-1 factoring worked examples.
- [`docs/brent-suyama.txt`](brent-suyama.txt) — background on the
  Brent-Suyama stage-2 extension; **not currently implemented** in Mlucas
  (see §2).
- [`help.txt`](../help.txt) — full command-line flag reference (`-s`,
  `-fft`, `-cpu`/`-core`, `-m`/`-f`, `-iters`, etc).
