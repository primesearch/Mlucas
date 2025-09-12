[![Actions Status](https://github.com/primesearch/Mlucas/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/primesearch/Mlucas/actions/workflows/ci.yml)

# Mlucas
Ernst Mayer's Mlucas and Mfactor programs for GIMPS

[Ernst Mayer passed away unexpectedly](https://www.mersenneforum.org/showthread.php?t=28890) on September 10, 2023. This repository contains his posthumously released Mlucas v21 code, which is now maintained by the Great Internet Mersenne Prime Search (GIMPS) community. AutoPrimeNet (the Python PrimeNet program) previously bundled with Mlucas is now maintained in a [separate repository](https://github.com/tdulcet/AutoPrimeNet).

Mlucas and Mfactor are 100% open source programs. Mlucas is for [primality](https://en.wikipedia.org/wiki/Primality_test) and [P-1](https://en.wikipedia.org/wiki/Pollard%27s_p_%E2%88%92_1_algorithm) testing of [Mersenne](https://en.wikipedia.org/wiki/Mersenne_prime) and [Fermat](https://en.wikipedia.org/wiki/Fermat_number) numbers, including support for the [Lucas-Lehmer](https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test), [Probable prime](https://en.wikipedia.org/wiki/Probable_prime) (PRP) and [Pépin](https://en.wikipedia.org/wiki/P%C3%A9pin%27s_test) tests. Mfactor is for trial factoring. They support x86 Intel and AMD, ARM and other CPUs.

The original [Mlucas README](https://mersenneforum.org/mayer/README.html) is available for posterity and contains a lot of information, but note that it is no longer up to date. For more information about Mlucas v21, please see the [Ernst's Mlucas - the future](https://www.mersenneforum.org/showthread.php?t=28926) thread on the Mersenne Forum.

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

## Usage

### Automatic method

Linux users can use the [Mlucas install script](https://github.com/tdulcet/Distributed-Computing-Scripts#mlucas) to automatically download, build, setup and run Mlucas, including downloading, setting up and running the [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) for automated PrimeNet assignments.

### Manual method

Dependencies:
* Make
* GNU C or Clang compiler
* \*GNU Multiple Precision (GMP) library
* \*Portable Hardware Locality (hwloc) library
* \*Python 3

\* Optional

#### Download

##### Linux

1. Verify that the dependencies above are installed. On Debian and Ubuntu, run: `sudo apt update` and `sudo apt install build-essential libgmp-dev libhwloc-dev`.
2. If one has git installed, just run: `git clone https://github.com/primesearch/Mlucas.git`. Otherwise, download the latest archive: `wget https://github.com/primesearch/Mlucas/archive/main.tar.gz` and then decompress the files: `tar -xzvf main.tar.gz`.
3. To download AutoPrimeNet, run: `wget -nv https://raw.github.com/tdulcet/AutoPrimeNet/main/autoprimenet.py`.

##### macOS

1. Verify that the dependencies above are installed. Run: `brew install gmp hwloc`.
2. If one has git installed, just run: `git clone https://github.com/primesearch/Mlucas.git`. Otherwise, download the latest archive: `curl -fLO https://github.com/primesearch/Mlucas/archive/main.tar.gz` and then decompress the files: `tar -xzvf main.tar.gz`.
3. To download AutoPrimeNet, run: `curl -sSfLO https://raw.github.com/tdulcet/AutoPrimeNet/main/autoprimenet.py`.

##### Windows

Native Windows builds are experimental. For now, Windows users should use the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) (WSL) and follow the [Linux](#linux) instructions above instead.

1. Download and install [MSYS2](https://www.msys2.org/).
2. Verify that the dependencies above are installed. With the MINGW64 environment, run: `pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-hwloc`.
3. If one has git installed, just run: `git clone https://github.com/primesearch/Mlucas.git`. Otherwise, download the latest archive: `wget https://github.com/primesearch/Mlucas/archive/main.tar.gz` and then decompress the files: `tar -xzvf main.tar.gz`.
4. To download AutoPrimeNet, run: `wget -nv https://raw.github.com/tdulcet/AutoPrimeNet/main/autoprimenet.py`.

#### Build

1. Change into the `Mlucas` directory. Run: `cd Mlucas` or `cd Mlucas-main` depending on which method one used to download it.
2. Run:
	* To build Mlucas: `bash makemake.sh [use_hwloc]`.
	* To build Mfactor: `bash makemake.sh mfac [word]`, where  `word` is optionally one of `1word`, `2word`, `3word`, `4word` or `nword`.

To build with Clang or another compiler instead of GCC, run: `export CC=<compiler>`, for example: `export CC=clang`.

#### Setup and Run

1. Change into the `obj` directory. Run: `cd obj` or `cd obj_mfac` depending on if one built Mlucas or Mfactor respectively.

This README is still in progress. For now, see the original [Mlucas README](https://mersenneforum.org/mayer/README.html), which has more information about how to setup and run Mlucas. Also see [Help](#help) below. Note that with Mlucas v21, if built with the hwloc library, one would want to use the new `-core` option instead of `-cpu`.

## Help

The [help.txt](help.txt) file includes a variety of usage information not covered in the original [README](https://mersenneforum.org/mayer/README.html), concentrating largely on the Mlucas command line options. A separate documentation page covers [Fermat numbers](docs/Fermat-testing.md).

## Contributing

Pull requests welcome!
