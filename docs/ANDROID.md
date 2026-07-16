# Running Mlucas on Android phones

Old and cheap-but-capable Android phones can be repurposed as a low-cost ARM compute cluster. Credit to Matthew Ling for the original setup guide and [c3tools](https://github.com/sillygitter/c3tools) scripts this page is adapted from; see the [CellPhone Compute Cluster for GIMPS](https://www.mersenneforum.org/forumdisplay.php?f=155) forum for community discussion and troubleshooting.

Requirements: an ARM CPU (ideally Armv8/ASIMD-capable — a non-SIMD phone still gets roughly 2/3 the throughput of an ASIMD one), enough working screen to do basic setup, and no Factory Reset Protection lock. Rooting is not required.

## Setup

1. Enable Wi-Fi, enable Developer Options (usually `Settings → About Device → Build Number`, tap 7 times), then enable USB debugging and "Unknown Sources" for installing an APK from outside the Play Store.
2. Install [UserLAnd](https://github.com/CypherpunkArmory/UserLAnd) (available via F-Droid or side-loaded APK), then use it to install a Ubuntu environment configured for SSH access.
3. Inside the UserLAnd Ubuntu shell: `sudo apt install ssh` and note the device's IP with `hostname -I`.
4. From your PC, `ssh -p 2022 <username>@<ip>` using the credentials you set in UserLAnd.
5. Install build dependencies in the Ubuntu shell: `sudo apt -y install git wget python3` (see [Prerequisites](../README.md#prerequisites) for what each is for).
6. Clone and build Mlucas as in the main [Build](../README.md#build) instructions — `bash makemake.sh` will auto-detect ASIMD support, or use `bash makemake.sh nosimd` on a non-SIMD device.
7. Run the [basic self-test](../README.md#basic-build-self-test) and then [performance-tuning](../README.md#performance-tuning) self-tests to build `mlucas.cfg`. On most phones, a single Mlucas instance using all cores (e.g. `-cpu 0:3` on a quad-core phone) gives more reasonable completion times than splitting into several smaller instances, at a modest total-throughput cost — worthwhile since a GIMPS assignment can take months on phone-class hardware.
8. [c3tools](https://github.com/sillygitter/c3tools) can automate step 7 across many identically-configured phones: `git clone https://github.com/sillygitter/c3tools && ~/c3tools/setup 0:3` (passing your chosen `-cpu` range).

Once set up, consider leaving Wi-Fi off except for periodic single-shot [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) runs (`-t 0`) to submit results and fetch new work — this also limits the phone's own background chatter and ad/bloatware activity.

**Battery safety note:** phones run 24/7 under sustained CPU load are prone to battery swelling within a few months. If you notice the back of the case bulging, stop the affected device and retire or replace its battery through normal means; do not attempt to puncture or vent a swollen lithium battery — this carries a real fire risk.

## P-1 factoring on Android

Android's memory management (or its interaction with UserLAnd) tends to confuse Mlucas's automatic stage 2 available-memory detection. Override it explicitly with `-pm1_s2_nbuf {+int}` (buffer count must be a multiple of 24 or 40 — Mlucas will round down to the nearest valid multiple). As a starting point, assume roughly half of installed RAM is usable; at the current GIMPS testing wavefront each buffer is about one FFT-length-in-bytes (e.g. ~48MB at a 6M-double FFT length). After starting stage 2, watch for disk-swapping (`top`'s `kswapd` entries, or per-iteration time exceeding roughly 1.25x the stage 1 rate) as a sign to lower the buffer count. See [help.txt](../help.txt) section `[9a]` for the full `-pm1_s2_nbuf`/`-maxalloc` reference.
