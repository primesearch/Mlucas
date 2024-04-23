# Introduction to Fermat number testing with Mlucas

Mlucas is a powerful Mersenne primality testing and factoring package which supports Lucas–Lehmer 
testing, Fermat probable primality testing, and Pollard _p_–1 factorisation; justly, most applications of the 
software have concentrated on examination of the Mersenne numbers, however Mlucas is also capable of 
testing the primality of Fermat numbers, testing the compositeness of Fermat cofactors, and running _p_–1
factorisation on Fermat numbers. This document is aimed at filling in several lacunae in the Mlucas 
documentation with respect to Fermat number testing.

## Mlucas build notes

The Fermat modular squaring code is architecture-specific, ideally requiring Mlucas to be built with one of 
the SIMD modes (AVX, AVX2, SSE2, etc). This should be handled automatically by the `makemake.sh` script 
for most hardware, and with the exception of 64-bit ARM (for which workarounds may be available) the default
build should be capable of testing Fermat numbers.

For example, Apple Silicon is ARM which supports ASIMD, meaning a native Mlucas build cannot test Fermat numbers. One 
possible answer would be to build an SSE2 Mlucas on an Intel-based Mac, which should be capable of running under 
Rosetta emulation on Apple Silicon, at the penalty of a significant performance hit for larger exponents.

(When compiling on Intel with the intention of building a compatible version for Apple Silicon, the Makefile 
requires static linking of the GMP library, by the following alteration to the linker flags:
`LDLIBS = -Bstatic ${LD_ARGS[@]} -Bdynamic`
Static linking of libraries is usually not recommended, and this may not be a solution in all possible cases.)

## Fermat self-testing and populating `fermat.cfg`

Assuming you have a working build, Mlucas allows any Fermat number _F<sub>m</sub>_ = 2<sup>2<sup>_m_</sup></sup> + 1 with exponent _m_ in 
the range [13,63] to be entered as an argument for self-testing, by using the `-f` flag like so:

`./Mlucas -f <exponent> -iters <number> [-fft <fft_length>]`

However, the software only provides fast Fourier transforms (FFTs) up to 512M in length, which further 
limits the largest Fermat number that may be tested, up to _F_<sub>33</sub>.

The `-iters` flag is mandatory and should be followed by a number less than a million; the recommended 
values for `-iters` to be set to are 100, 1000, or 10000, since the self-testing code has correct residues 
for all Fermat numbers [14,33] saved, and any error in the test result will be immediately discovered.

Mlucas does not support Fermat testing for all possible FFT lengths. Generally, FFTs for testing Mersenne 
numbers are available for any length _k_·2<sup>_n_</sup> in kilobytes, where _k_ = [8,15], _n_ ≥ 0. 
In the case of Fermat numbers however, only the subset _k_ = {2, 7, 15} supports Fermat testing, along with
_k_ = 63 when _n_ ≥ 4.

When self-testing Mlucas for Mersenne numbers, a command such as `./Mlucas -s all` configures 
Mlucas for the majority of FFT lengths and saves the results in `mlucas.cfg`. A similar process must be 
carried out for the Fermat numbers, using the `./Mlucas -f` command above, which saves results in the 
file `fermat.cfg`.

A script [`config-fermat.sh`](https://github.com/primesearch/Mlucas/blob/main/config-fermat.sh) has been written to automate the generation of the `fermat.cfg` file using 
the set of possible FFT lengths. After the license, the next few lines declare some variables which may 
be customised as desired (for example, 10000 iterations will require significant runtime at the largest 
possible FFT sizes):
```bash
# Mlucas
MLUCAS=./Mlucas

# Number of iterations (use 100, 1000, or 10000 to match pre-computed values)
ITERS=100

# Minimum Fermat number (15 or greater)
MIN=15

# Maximum Fermat number (33 or less)
MAX=29
```
The script should be invoked from the build directory, typically `obj`, with the following command, which 
may include as parameters any cpu-specific options, such as `-core` or `-cpu`:
```
bash ../config-fermat.sh -cpu 0:3
```
Once self-testing has occurred, production work on Fermat numbers may be carried out using a `worktodo.txt` 
file (in previous versions of Mlucas this file was named `worktodo.ini`) with entries in one of the several 
formats below.

Only assignments for ECM factoring of Fermat numbers are distributed by Primenet, so work assignments 
cannot be obtained for Fermat numbers using the `primenet.py` script. However _p_–1 results for Fermat 
numbers up to _F_<sub>29</sub> may be submitted to [mersenne.org](https://www.mersenne.org/manual_result/) as Manual results.

## Primality testing: Pépin’s test

For the Pépin test, the `worktodo` format is simply:
```
Fermat,Test=<exponent>
```

Fermat numbers in the range [13,22] may select a non-optimal FFT length by default in production mode.
If this is the case, the FFT may be overridden using the command:
```
./Mlucas -fft FFT_LENGTH -shift 0
```
The optimal FFT lengths vary from 2K up to 512M, as shown in the table below under [testable Fermat numbers](#the-testable-fermat-numbers-fm--22m--1).

Currently, all Fermat numbers up to _F_<sub>30</sub> have received a [Pépin test](https://www.mersenneforum.org/showthread.php?t=18748), and moreover _F_<sub>31</sub> and _F_<sub>32</sub> are known 
to be composite. However, the character of their cofactors is unknown so a Pépin test would be a 
necessary prerequisite for testing them; and as for _F_<sub>33</sub>, this is yet to be 
tested for primality.

The Pépin test uses Gerbicz error checking (GEC) to assure the reliability of the computation, however 
residue shifting is not implemented for Fermat modular squaring with GEC. Thus when invoking Mlucas the 
flag `-shift 0` _must_ be added to the Mlucas command line:
```
./Mlucas -shift 0
```
If non-zero residue shifting is used, the Pépin test will not be able to progress beyond the millionth 
iteration (the default interval for GEC) as error checking will be unable to confirm whether or not the 
computation is correct.

## Cofactor compositeness testing: Suyama’s test
Suyama’s test is a Fermat probable primality test on the cofactor of a Fermat number. As a prerequisite, 
you _must_ already have run a Pépin test, which should have saved a final residue as a file named e.g. `f23` 
for a Pépin test of _F_<sub>23</sub>. Do **not** try to run a Suyama test as a single work type incorporating the preceding 
Pépin test.

The `worktodo` format for the Suyama test is:
```
PRP=N/A,1,2,<2^m>,+1,99,0,3,5,"known_factor_1[,known_factor_2,...]"
```

Note that this work format uses the numeric value of 2<sup>_m_</sup> of a Fermat number _F<sub>m</sub>_, rather than just the 
exponent _m_. At least one known prime factor must be supplied; composite factors are disallowed.

Prior to testing, Mlucas tries an “integrity check” on each known factor, which has the effect of testing 
whether the Pépin residue _R_ has been correctly calculated. This calculation can be performed at any point 
up to the final iteration. More information is available [here](https://github.com/primesearch/Mlucas/issues/4).

## Pollard _p_–1 factoring

Work entries for _p_–1 factoring have two variations, however only the `Pminus1` format is supported for 
Fermat numbers:
```
Pminus1=N/A,1,2,<2^m>,+1,B1,B2[,trial_factoring_bits][,B2_start][,"known_factors"]
```

The same note as regards the exponent for the Suyama test applies here also; however supplying known, prime factors
is optional. The `-shift 0` flag is again a necessary addition to the Mlucas command line.

The `B2_start` variable supports breaking up stage 2 of the algorithm among multiple instances of Mlucas.

## Fermat number `worktodo.txt` examples

An example of each of the work types, Pépin and Suyama tests, and Pollard _p_–1 factoring:
```
Fermat,Test=23
PRP=N/A,1,2,1073741824,+1,99,0,3,5,"640126220763137,1095981164658689"
Pminus1=N/A,1,2,8589934592,+1,10000000,10000000
```
Note the _p_–1 test of _F_<sub>33</sub> (as actually [performed](https://www.mersenneforum.org/showthread.php?t=29183) by Ernst Mayer) 
used the same _B1_ and _B2_ bounds to run only the first stage of the algorithm.

## The testable Fermat numbers _F<sub>m</sub>_ = 2<sup>2<sup>_m_</sup></sup> + 1

The following table lists the default FFT lengths selected by Mlucas for the Fermat numbers in the range 
[13,33] and the known factors, required for some of the test types above.

 _m_|       2<sup>_m_</sup>|Default FFT |Optimal FFTs     |known_factor(s)
--|---------:|-----------:|----------------:|--------------------------------------------------------------------
13|      8192|        1K\*|               2K|`"2710954639361,2663848877152141313,3603109844542291969,319546020820551643220672513"`
14|     16384|        1K\*|           2K, 4K|`"116928085873074369829035993834596371340386703423373313"`
15|     32768|        2K\*|           2K, 4K|`"1214251009,2327042503868417,168768817029516972383024127016961"`
16|     65536|        3K\*|               4K|`"825753601,188981757975021318420037633"`
17|    131072|  6K\*, 7K  |           7K, 8K|`"31065037602817,7751061099802522589358967058392886922693580423169"`
18|    262144|       13K\*|    14K, 15K, 16K|`"13631489,81274690703860512587777"`
19|    524288|       26K\*|    28K, 30K, 32K|`"70525124609,646730219521,37590055514133754286524446080499713"`
20|   1048576|       52K\*|    56K, 60K, 64K|
21|   2097152|      104K\*| 112K, 120K, 128K|`"4485296422913"`
22|   4194304|      208K\*| 224K, 240K, 256K|`"64658705994591851009055774868504577"`
23|   8388608|      448K  | 448K, 480K, 512K|`"167772161"`
24|  16777216|      896K  | 896K, 960K, 1008K, 1M|
25|  33554432|     1792K  | 1792K, 1920K, 2M|`"25991531462657,204393464266227713,2170072644496392193"`
26|  67108864|     3584K  | 3584K, 3840K, 4032K, 4M|`"76861124116481"`
27| 134217728|  7M, 7.5M  | 7M, 7.5M, 7.875M, 8M|`"151413703311361,231292694251438081"`
28| 268435456|       15M  | 14M, 15M, 15.75M, 16M|`"1766730974551267606529"`
29| 536870912|       30M  |  30M, 31.5M, 32M|`"2405286912458753"`
30|1073741824|       60M  |    60M, 63M, 64M|`"640126220763137,1095981164658689"`
31|2147483648|      120M  | 120M, 126M, 128M|`"46931635677864055013377"`
32|4294967296|      256M  | 240M, 252M, 256M|`"25409026523137"`
33|8589934592|      512M  |       504M, 512M|

\* These default FFTs selected by Mlucas have radices that cannot be used for Fermat modular squaring;
thus _F_<sub>13</sub> to _F_<sub>22</sub> cannot be tested with Mlucas without overriding the default FFT selected, _e.g._:
```
./Mlucas -fft 224K -shift 0
```
As mentioned above at [Mlucas build notes](#mlucas-build-notes), building Mlucas is highly architecture-specific, and some 
build types will not support all of the optimal FFT lengths.

Finally, 2K appears to be the smallest usable FFT for the smallest testable Fermat number _F_<sub>13</sub>, with the command line settings `-fft 2 -radset 8,8,16 -shift 0`.

The 4K FFT seems to be more reliable on a wider range of systems for the next larger Fermat numbers _F_<sub>14</sub> to _F_<sub>16</sub>.
