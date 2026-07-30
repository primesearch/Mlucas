#!/usr/bin/env python3
"""pairwise-merge-build.py -- find open pull requests that merge without a
conflict and then do not compile.

With ~100 open PRs against one codebase the dangerous interaction is not the
merge conflict; git tells you about those.  It is the pair that merges with zero
conflicts and then fails to build.  Verified shapes:

  * one PR deletes a declaration whose only remaining use is added by the other
    (#210 drops 'stride'; #152 independently adds a new use of it in the same
    file -- 3-way apply rc=0, zero conflict lines, then
    "error: 'stride' undeclared");
  * a PR's new code lands in a header that a *second* translation unit also
    includes, and that TU does not declare what the new code needs (#179's
    LL-injection block is pulled into cy992_process_chunk()).

Three methodology traps this tool is built to avoid, each of which has already
produced a wrong answer by hand:

  1. Merge each PR against *current main*, never against the other PR.  A
     'git merge-base A B' between PRs cut from different mains resolves to an
     ancient base and invents conflicts -- that produced a bogus "36 files
     conflict" report for #208+#68 when the truth was one file.
  2. 'git apply' without --3way is not a merge test.  It fails where a real
     merge succeeds, which *masks* the interesting cases.  This tool uses
     'git merge-tree', a real 3-way merge.
  3. 'git merge-tree' prints the whole merged tree listing.  Only lines
     beginning with CONFLICT are conflicts; grepping it for paths over-reports.

Scope reduction: only PRs touching a common file can interact this way, so
file-touch sets are computed once and only intersecting pairs are enumerated.

"Affected translation unit" means every .c that *transitively includes* a
changed file, not just the changed .c files -- Mlucas has many .h files that are
effectively source (*_main_carry_loop.h, carry_gcc64.h, sse2_macro_gcc64.h), and
trap 2 above is precisely a header whose second includer breaks.

Usage
-----
    ./pairwise-merge-build.py --refresh              # fetch open PR heads
    ./pairwise-merge-build.py --solo                 # every PR vs current main
    ./pairwise-merge-build.py --pairs                # every intersecting pair
    ./pairwise-merge-build.py --pairs --only 152,210 # just these PRs
    ./pairwise-merge-build.py --plan                 # print counts, build nothing

Exit status is 1 if any tested combination fails to build.
"""

import argparse
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

REPO_UPSTREAM = "primesearch/Mlucas"
REFNS = "refs/prcheck"

# nosimd is the default first mode: it is the fastest, and declaration-level
# breakage -- which is the entire shape of this defect class -- is ISA
# independent.  Flags mirror makemake.sh: ARGS=(-DUSE_THREADS) plus the
# per-mode additions, and CFLAGS_PROBED=(-Wall -g -O3 -std=gnu99).
MODE_FLAGS = {
    "nosimd": [],
    "sse2":   ["-DUSE_SSE2", "-msse2"],
    "avx":    ["-DUSE_AVX", "-mavx"],
    "avx2":   ["-DUSE_AVX2", "-mavx2", "-mfma"],
    "avx512": ["-DUSE_AVX512", "-mavx512f", "-mavx512cd", "-mavx512dq",
               "-mavx512bw", "-mavx512vl", "-mfma"],
}
# -D_GNU_SOURCE is not optional: makemake.sh sets it via
#     CPPFLAGS ?= -D_GNU_SOURCE -I/usr/local/include -I/opt/homebrew/include
# and without it threadpool.c fails on CPU_ZERO/CPU_SET/sched_setaffinity in
# *every* tree, main included.  Leaving it out makes the tool report a
# spurious failure for every PR whose TU set reaches threadpool.c.
BASE_CFLAGS = ["-std=gnu99", "-D_GNU_SOURCE", "-DUSE_THREADS", "-Wall"]
if not os.path.exists("/usr/include/gmp.h"):
    BASE_CFLAGS.append("-DINCLUDE_GMP=0")


def git(*args, **kw):
    cwd = kw.pop("cwd", None)
    check = kw.pop("check", True)
    p = subprocess.run(["git"] + list(args), cwd=cwd, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise RuntimeError("git %s failed: %s" % (" ".join(args), p.stderr.strip()))
    return p


# ------------------------------------------------------------------ PR state --
def open_prs():
    p = subprocess.run(["gh", "pr", "list", "--repo", REPO_UPSTREAM, "--state", "open",
                        "--limit", "300", "--json",
                        "number,title,headRefName,isDraft"],
                       capture_output=True, text=True, check=True)
    return {x["number"]: x for x in json.loads(p.stdout)}


def refresh(repo, prs):
    """Fetch every open PR head into refs/prcheck/<N>.  Read-only with respect
    to the PR branches themselves; nothing here writes to a remote."""
    spec = ["+refs/pull/%d/head:%s/%d" % (n, REFNS, n) for n in prs]
    for i in range(0, len(spec), 60):
        git("fetch", "--quiet", "https://github.com/%s" % REPO_UPSTREAM, *spec[i:i + 60],
            cwd=repo)
    git("fetch", "--quiet", "https://github.com/%s" % REPO_UPSTREAM,
        "+refs/heads/main:refs/prcheck/main", cwd=repo)


def pr_ref(n):
    return "%s/%d" % (REFNS, n)


def touched_files(repo, main, n):
    mb = git("merge-base", main, pr_ref(n), cwd=repo).stdout.strip()
    out = git("diff", "--name-only", mb, pr_ref(n), cwd=repo).stdout
    return mb, set(x for x in out.split("\n") if x)


# ------------------------------------------------------------------- merging --
CONFLICT_RE = re.compile(r"^CONFLICT ", re.M)


def merge_onto(repo, base_commit, ref, label):
    """3-way merge `ref` onto `base_commit`, using their true merge-base.
    Returns (commit_or_None, [conflict lines])."""
    mb = git("merge-base", base_commit, ref, cwd=repo).stdout.strip()
    p = git("merge-tree", "--write-tree", "--merge-base=" + mb,
            base_commit, ref, cwd=repo, check=False)
    lines = p.stdout.split("\n")
    tree = lines[0].strip() if lines else ""
    # Only lines that literally begin with CONFLICT are conflicts; the rest of
    # merge-tree's output is the merged tree listing.
    conflicts = [l for l in lines if l.startswith("CONFLICT")]
    if p.returncode != 0 or conflicts or not tree:
        return None, conflicts or ["merge-tree rc=%d" % p.returncode]
    commit = git("commit-tree", tree, "-p", base_commit, "-p", ref,
                 "-m", "merge %s" % label, cwd=repo).stdout.strip()
    return commit, []


# ---------------------------------------------------- include graph / TU set --
INCLUDE_RE = re.compile(r'^\s*#\s*include\s*"([^"]+)"', re.M)


def include_graph(srcdir):
    """includers[h] = set of files that directly #include h."""
    includers = defaultdict(set)
    for name in os.listdir(srcdir):
        if not name.endswith((".c", ".h")):
            continue
        try:
            text = open(os.path.join(srcdir, name), encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        for inc in INCLUDE_RE.findall(text):
            includers[os.path.basename(inc)].add(name)
    return includers


def affected_tus(srcdir, changed):
    """Every .c that transitively includes any changed file (plus changed .c)."""
    includers = include_graph(srcdir)
    seen, stack, tus = set(), [], set()
    for path in changed:
        if not path.startswith("src/"):
            continue
        stack.append(os.path.basename(path))
    while stack:
        f = stack.pop()
        if f in seen:
            continue
        seen.add(f)
        if f.endswith(".c") and os.path.exists(os.path.join(srcdir, f)):
            tus.add(f)
        for up in includers.get(f, ()):
            if up not in seen:
                stack.append(up)
    return sorted(tus)


# -------------------------------------------------------------------- build --
def materialise(repo, commit, dest):
    os.makedirs(dest, exist_ok=True)
    p1 = subprocess.Popen(["git", "archive", "--format=tar", commit], cwd=repo,
                          stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["tar", "-x", "-C", dest], stdin=p1.stdout)
    p1.stdout.close()
    p2.wait()
    p1.wait()
    return p2.returncode == 0 and p1.returncode == 0


def compile_tus(workdir, tus, mode, cc, jobs, syntax_only):
    """Compile the given TUs; return the first failure as (tu, stderr) or None."""
    srcdir = os.path.join(workdir, "src")
    flags = [cc] + BASE_CFLAGS + MODE_FLAGS[mode]
    flags += ["-fsyntax-only"] if syntax_only else ["-O1", "-c", "-o", "/dev/null"]
    failures = []

    def one(tu):
        p = subprocess.run(flags + [tu], cwd=srcdir, capture_output=True, text=True)
        return tu, p.returncode, p.stderr

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for tu, rc, err in pool.map(one, tus):
            if rc != 0:
                failures.append((tu, err))
    return failures


def first_errors(stderr, n=4):
    out = [l for l in stderr.split("\n") if " error: " in l]
    return out[:n]


# ---------------------------------------------------------------- the driver --
class Checker:
    def __init__(self, args):
        self.a = args
        self.repo = args.repo
        self.main = git("rev-parse", args.main, cwd=self.repo).stdout.strip()
        self.tmp = args.workdir or tempfile.mkdtemp(prefix="pmb-")
        self.tmp_is_ours = not args.workdir
        os.makedirs(self.tmp, exist_ok=True)
        # Translation units that already fail to compile on *main* in this mode.
        # Subtracting them is not optional: src/factor.c does not compile under
        # -DUSE_AVX512 on 418f365 (twopmodq100_2WORD_DOUBLE_q2 has no
        # declaration), so without this every combination whose TU set reaches
        # factor.c would be reported as broken by the PRs rather than by main.
        self.main_dir = os.path.join(self.tmp, "_main")
        self.main_fail = {}
        materialise(self.repo, self.main, self.main_dir)

    def fails_on_main(self, tu):
        if tu not in self.main_fail:
            self.main_fail[tu] = bool(compile_tus(self.main_dir, [tu], self.a.mode,
                                                  self.a.cc, 1,
                                                  not self.a.full_objects))
        return self.main_fail[tu]

    def cleanup(self):
        if self.tmp_is_ours:
            shutil.rmtree(self.tmp, ignore_errors=True)

    def combo(self, refs, label):
        """Merge each ref onto main in turn, then build the affected TUs."""
        base, changed = self.main, set()
        for n in refs:
            c, conflicts = merge_onto(self.repo, base, pr_ref(n), "#%d" % n)
            if c is None:
                return {"label": label, "prs": refs, "status": "conflict",
                        "detail": conflicts[:6]}
            base = c
            changed |= self.touch[n]
        work = os.path.join(self.tmp, label.replace(" ", "").replace("+", "_"))
        shutil.rmtree(work, ignore_errors=True)
        if not materialise(self.repo, base, work):
            return {"label": label, "prs": refs, "status": "error",
                    "detail": ["could not materialise tree"]}
        tus = affected_tus(os.path.join(work, "src"), changed)
        if not tus:
            shutil.rmtree(work, ignore_errors=True)
            return {"label": label, "prs": refs, "status": "no-tu", "detail": []}
        fails = compile_tus(work, tus, self.a.mode, self.a.cc,
                            self.a.build_jobs, not self.a.full_objects)
        shutil.rmtree(work, ignore_errors=True)
        preexisting = [tu for tu, _ in fails if self.fails_on_main(tu)]
        fails = [(tu, err) for tu, err in fails if tu not in preexisting]
        if fails:
            det = []
            for tu, err in fails[:3]:
                det.append(tu + ": " + (first_errors(err)[0] if first_errors(err)
                                        else err.strip().split("\n")[-1]))
                det.extend("    " + e for e in first_errors(err)[1:3])
            return {"label": label, "prs": refs, "status": "BUILD-FAIL",
                    "ntu": len(tus), "detail": det, "premain": preexisting}
        return {"label": label, "prs": refs, "status": "ok", "ntu": len(tus),
                "detail": [], "premain": preexisting}


def main():
    # Line-buffer stdout: this tool runs for many minutes and its progress is
    # only useful if it arrives as it happens.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", default=here)
    ap.add_argument("--main", default="main")
    ap.add_argument("--refresh", action="store_true", help="fetch open PR heads first")
    ap.add_argument("--solo", action="store_true", help="test each PR alone onto main")
    ap.add_argument("--pairs", action="store_true", help="test intersecting pairs")
    ap.add_argument("--plan", action="store_true", help="print the pair counts and stop")
    ap.add_argument("--only", default="", help="restrict to these PR numbers")
    ap.add_argument("--skip", default="", help="exclude these PR numbers")
    ap.add_argument("--partners", action="store_true",
                    help="with --only: also consider every open PR that shares a "
                         "file with one of the named PRs, and test only the pairs "
                         "that involve a named PR. This is the per-PR CI shape: "
                         "~18 combinations for a typical PR rather than 916.")
    ap.add_argument("--mode", default="nosimd", choices=list(MODE_FLAGS))
    ap.add_argument("--cc", default=os.environ.get("CC", "cc"))
    ap.add_argument("--jobs", type=int, default=4, help="concurrent combinations")
    ap.add_argument("--build-jobs", type=int, default=4, help="concurrent TU compiles")
    ap.add_argument("--full-objects", action="store_true",
                    help="compile to objects instead of -fsyntax-only")
    ap.add_argument("--extra-ref", action="append", default=[], metavar="N=REV",
                    help="test a pinned revision under a synthetic PR number; "
                         "used to validate the tool against a known-broken state "
                         "of a PR that has since been fixed")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--json", default="", help="write the full result set here")
    args = ap.parse_args()

    prs = open_prs()
    pinned = {}
    for spec in args.extra_ref:
        n, _, rev = spec.partition("=")
        pinned[int(n)] = rev
        prs[int(n)] = {"number": int(n), "title": "pinned " + rev, "isDraft": False}
    only = {int(x) for x in args.only.split(",") if x.strip()}
    skip = {int(x) for x in args.skip.split(",") if x.strip()}
    nums = sorted(n for n in prs if (not only or n in only) and n not in skip)

    ck = Checker(args)
    for n, rev in pinned.items():
        sha = git("rev-parse", rev, cwd=args.repo).stdout.strip()
        git("update-ref", pr_ref(n), sha, cwd=args.repo)
    if args.refresh:
        t = time.time()
        refresh(args.repo, [n for n in nums if n not in pinned])
        print("fetched %d PR heads in %.0fs" % (len(nums), time.time() - t))

    # File-touch sets, computed once.
    ck.touch, dropped = {}, []
    for n in nums:
        try:
            _, files = touched_files(args.repo, ck.main, n)
        except RuntimeError:
            dropped.append(n)
            continue
        ck.touch[n] = files
    nums = [n for n in nums if n in ck.touch]
    if dropped:
        print("skipped %d PR(s) with no fetched head: %s"
              % (len(dropped), ",".join(map(str, dropped))))

    if args.partners and only:
        # Re-widen to the whole open set, then keep only the PRs that share a
        # file with something in --only.
        allnums = sorted(n for n in prs if n not in skip)
        for n in allnums:
            if n in ck.touch:
                continue
            try:
                _, files = touched_files(args.repo, ck.main, n)
            except RuntimeError:
                continue
            ck.touch[n] = files
        nums = sorted({n for n in ck.touch
                       if n in only or any(ck.touch[n] & ck.touch[m]
                                           for m in only if m in ck.touch)})

    pairs = [(a, b) for a, b in itertools.combinations(nums, 2)
             if ck.touch[a] & ck.touch[b]
             and (not (args.partners and only) or a in only or b in only)]
    print("PRs: %d   all pairs: %d   pairs sharing a file: %d (%.1f%%)"
          % (len(nums), len(nums) * (len(nums) - 1) // 2, len(pairs),
             100.0 * len(pairs) / max(1, len(nums) * (len(nums) - 1) // 2)))
    if args.plan:
        ck.cleanup()
        return 0

    jobs = []
    if args.solo:
        solo_nums = sorted(only) if (args.partners and only) else nums
        jobs += [((n,), "#%d" % n) for n in solo_nums if n in ck.touch]
    if args.pairs:
        jobs += [((a, b), "#%d+#%d" % (a, b)) for a, b in pairs]
    if not jobs:
        print("nothing to do: pass --solo and/or --pairs")
        return 0

    t0 = time.time()
    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for i, r in enumerate(pool.map(lambda j: ck.combo(*j), jobs), 1):
            results.append(r)
            if r["status"] in ("BUILD-FAIL", "error"):
                print("[%d/%d] %-16s %s" % (i, len(jobs), r["label"], r["status"]))
                for d in r["detail"]:
                    print("        " + d)
            elif i % 25 == 0:
                print("[%d/%d] ..." % (i, len(jobs)), file=sys.stderr)

    bad = [r for r in results if r["status"] == "BUILD-FAIL"]
    conf = [r for r in results if r["status"] == "conflict"]
    print("\n%d combination(s) in %.0fs: %d build-fail, %d conflict, %d ok, %d no-TU"
          % (len(results), time.time() - t0, len(bad), len(conf),
             sum(1 for r in results if r["status"] == "ok"),
             sum(1 for r in results if r["status"] == "no-tu")))
    if args.json:
        json.dump(results, open(args.json, "w"), indent=1)
    ck.cleanup()
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
