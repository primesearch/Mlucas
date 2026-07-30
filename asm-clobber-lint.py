#!/usr/bin/env python3
"""asm-clobber-lint.py -- find x86 inline-asm blocks that use a hard register
which is neither bound to an operand nor listed as a clobber.

Motivation: a GCC/clang extended-asm statement may only touch registers that
are (a) bound to one of its operands, or (b) named in its clobber list.  Any
other hard register the template writes is, from the compiler's point of view,
still live across the asm, so the compiler is free to park a value there.  The
asm destroys it.  The result is a wrong answer or a SIGSEGV that looks for all
the world like a compiler bug -- and sanitizers cannot see it, because there is
no C-level operation to instrument.

Five instances of exactly this defect have been found in this tree, four of
which were first (mis)diagnosed as compiler codegen bugs.

How it works
------------
Pass A ("where"): parse every src/*.c and src/*.h as raw text.  This gives the
    true file:line of each asm block, including blocks that live inside macro
    definitions in headers -- which is where most of them are, and which is
    information that preprocessed output has already thrown away.

Pass B ("is it live"): preprocess each src/*.c once per ISA mode and parse the
    result.  This settles which '#if' branch is actually compiled, with no
    guessing.  Findings that are live in no ISA mode are dead code and are not
    reported.

The two passes are joined on a signature over (registers used, clobbers
declared), not on line numbers, so the join survives macro argument
substitution.

Usage
-----
    python3 asm-clobber-lint.py                   # lint against the baseline
    python3 asm-clobber-lint.py --write-baseline  # regenerate the baseline
    python3 asm-clobber-lint.py --all             # report everything, ignore baseline
    python3 asm-clobber-lint.py --modes avx512    # restrict the ISA sweep

Exit status is 1 if any finding is not covered by the baseline file, 0 otherwise.
"""

import argparse
import glob
import hashlib
import os
import re
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

# ---------------------------------------------------------------- ISA modes --
# Flags mirror the per-mode 'ARGS+=(...)' lines in makemake.sh.
COMMON_FLAGS = ["-E", "-std=gnu99", "-D_GNU_SOURCE", "-DUSE_THREADS"]
ISA_FLAGS = {
    "nosimd": [],
    "sse2":   ["-DUSE_SSE2", "-msse2"],
    "avx":    ["-DUSE_AVX", "-mavx"],
    "avx2":   ["-DUSE_AVX2", "-mavx2", "-mfma"],
    "avx512": ["-DUSE_AVX512", "-mavx512f", "-mavx512cd", "-mavx512dq",
               "-mavx512bw", "-mavx512vl", "-mfma"],
}
ALL_MODES = list(ISA_FLAGS)

# ------------------------------------------------------- register canonicals --
_GPR_ALIASES = {}
for _wide, _parts in {
    "rax": "rax eax ax al ah",
    "rbx": "rbx ebx bx bl bh",
    "rcx": "rcx ecx cx cl ch",
    "rdx": "rdx edx dx dl dh",
    "rsi": "rsi esi si sil",
    "rdi": "rdi edi di dil",
    "rbp": "rbp ebp bp bpl",
    "rsp": "rsp esp sp spl",
}.items():
    for _p in _parts.split():
        _GPR_ALIASES[_p] = _wide
for _n in range(8, 16):
    for _sfx in ("", "d", "w", "b"):
        _GPR_ALIASES["r%d%s" % (_n, _sfx)] = "r%d" % _n

# Registers we never report: the compiler cannot be told about them anyway, and
# an asm that really did clobber them would be broken in a way this lint is not
# the right tool for.
_IGNORED = {"rsp", "rip"}

_REG_RE = re.compile(r"(?:[xyz]mm(?:3[01]|[12]?[0-9])|k[0-7]|mm[0-7]|"
                     r"r(?:1[0-5]|[89])[dwb]?|"
                     r"[re]?(?:ax|bx|cx|dx|si|di|bp|sp)|[abcd][hl]|[sd]il|[bs]pl)\b")


def canon_reg(name):
    """Canonical clobber-list spelling of a register name, or None if not a
    register we track.  %rax/%eax/%ax/%al/%ah all fold to 'rax'; xmm/ymm/zmmN
    all fold to 'xmmN' (they are one hard register)."""
    name = name.lower()
    m = re.fullmatch(r"[xyz]mm(\d+)", name)
    if m:
        return "xmm" + m.group(1)
    if re.fullmatch(r"k[0-7]", name):
        return name
    if re.fullmatch(r"mm[0-7]", name):
        return name
    return _GPR_ALIASES.get(name)


# Constraint letters that pin an operand to one specific hard register.  These
# registers are legitimately used by the template without appearing in the
# clobber list -- the compiler knows about them.  Getting this wrong is the
# difference between a lint people act on and a lint people switch off.
_CONSTRAINT_REGS = {
    "a": ["rax"], "b": ["rbx"], "c": ["rcx"], "d": ["rdx"],
    "S": ["rsi"], "D": ["rdi"], "A": ["rax", "rdx"],
}


def constraint_regs(constraint):
    """Hard registers pinned by an operand constraint string."""
    out = set()
    i = 0
    while i < len(constraint):
        ch = constraint[i]
        if ch == "Y":              # two-character x86 constraints: Yz, Yk, Ym...
            if i + 1 < len(constraint) and constraint[i + 1] == "z":
                out.add("xmm0")
            i += 2
            continue
        out.update(_CONSTRAINT_REGS.get(ch, ()))
        i += 1
    return out


# ------------------------------------------------------------------- parsing --
def strip_comments(text):
    """Remove C comments, preserving byte offsets (replace with spaces) so that
    line numbers and offsets computed later stay correct."""
    out = list(text)
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == '"' or c == "'":
            q = c
            i += 1
            while i < n and text[i] != q:
                if text[i] == "\\":
                    i += 1
                i += 1
            i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            for k in range(i, j):
                if out[k] != "\n":
                    out[k] = " "
            i = j
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            j = text.find("\n", i)
            j = n if j < 0 else j
            for k in range(i, j):
                out[k] = " "
            i = j
        else:
            i += 1
    return "".join(out)


_ASM_RE = re.compile(r"\b(?:__asm__|asm)\b\s*(?:__volatile__|volatile)?\s*"
                     r"(?:__inline__|inline)?\s*(?:goto)?\s*\(")


def _skip_string(text, i):
    """i points at the opening quote; return index just past the closing one."""
    q = text[i]
    i += 1
    n = len(text)
    while i < n and text[i] != q:
        if text[i] == "\\":
            i += 1
        i += 1
    return i + 1


def split_sections(body):
    """Split an extended-asm parenthesised body into its colon-separated
    sections.  Colons inside strings, parentheses or brackets do not count."""
    secs, cur, depth = [], [], 0
    i, n = 0, len(body)
    while i < n:
        c = body[i]
        if c in '"\'':
            j = _skip_string(body, i)
            cur.append(body[i:j])
            i = j
            continue
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == ":" and depth == 0:
            secs.append("".join(cur))
            cur = []
            i += 1
            continue
        cur.append(c)
        i += 1
    secs.append("".join(cur))
    return secs


_STR_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')


def string_literals(section):
    return [m.group(1) for m in _STR_RE.finditer(section)]


def parse_operands(section):
    """Yield the constraint string of each operand in an output/input section.
    An operand is  [name] "constraint" (expr)  -- the constraint is the string
    literal immediately preceding the parenthesised expression."""
    regs = set()
    i, n, pending = 0, len(section), None
    while i < n:
        c = section[i]
        if c in '"\'':
            j = _skip_string(section, i)
            if c == '"':
                pending = section[i + 1:j - 1]
            i = j
            continue
        if c == "(":
            depth, j = 1, i + 1
            while j < n and depth:
                if section[j] in '"\'':
                    j = _skip_string(section, j)
                    continue
                if section[j] == "(":
                    depth += 1
                elif section[j] == ")":
                    depth -= 1
                j += 1
            if pending is not None:
                regs |= constraint_regs(pending)
                pending = None
            i = j
            continue
        if c == "[":                       # symbolic operand name; not a constraint
            j = section.find("]", i)
            i = (n if j < 0 else j + 1)
            continue
        i += 1
    return regs


def template_regs(template, extended):
    """Hard registers named literally in an asm template.

    Extended asm escapes register names as '%%rax'; basic asm (no operands)
    writes '%rax'.  '%0', '%1' and '%[name]' are operand references, not hard
    registers, and are deliberately not matched."""
    text = template.replace("\\n", "\n").replace("\\t", "\t")
    out = set()
    pat = r"%%(" if extended else r"%(?!%)("
    for m in re.finditer(pat + _REG_RE.pattern + r")", text):
        r = canon_reg(m.group(1))
        if r and r not in _IGNORED:
            out.add(r)
    return out


class Finding:
    __slots__ = ("path", "line", "used", "declared", "bound", "missing", "sig", "modes")

    def __init__(self, path, line, used, declared, bound):
        self.path, self.line = path, line
        self.used, self.declared, self.bound = used, declared, bound
        self.missing = sorted(used - declared - bound, key=sortkey)
        self.sig = signature(used, declared)
        self.modes = []


def sortkey(r):
    m = re.fullmatch(r"([a-z]+)(\d+)", r)
    return (m.group(1), int(m.group(2))) if m else (r, -1)


def signature(used, declared):
    """Stable identity for an asm block: what it touches and what it declares.
    Deliberately excludes line numbers and macro arguments so that the raw-source
    parse and the preprocessed parse of the same block agree."""
    h = hashlib.sha1()
    h.update(",".join(sorted(used, key=sortkey)).encode())
    h.update(b"|")
    h.update(",".join(sorted(declared, key=sortkey)).encode())
    return h.hexdigest()[:8]


def scan_text(text, locate, want_comments_stripped=True):
    """Parse every asm statement in `text`.  `locate(offset)` maps a character
    offset to a (path, line) pair.  Returns a list of Finding."""
    if want_comments_stripped:
        text = strip_comments(text)
    findings = []
    for m in _ASM_RE.finditer(text):
        open_paren = m.end() - 1
        depth, j, n = 0, open_paren, len(text)
        while j < n:
            c = text[j]
            if c in '"\'':
                j = _skip_string(text, j)
                continue
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = text[open_paren + 1:j]
        secs = split_sections(body)
        extended = len(secs) > 1
        used = template_regs(" ".join(string_literals(secs[0])), extended)
        if not used:
            continue
        bound = set()
        for s in secs[1:3]:
            bound |= parse_operands(s)
        declared = set()
        if len(secs) > 3:
            for lit in string_literals(secs[3]):
                r = canon_reg(lit.strip())
                if r:
                    declared.add(r)
        path, line = locate(m.start())
        f = Finding(path, line, used, declared, bound)
        if f.missing:
            findings.append(f)
    return findings


# ------------------------------------------------------------ location maps --
def raw_locator(path, text):
    starts = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            starts.append(i + 1)

    def locate(off):
        lo, hi = 0, len(starts) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if starts[mid] <= off:
                lo = mid
            else:
                hi = mid - 1
        return path, lo + 1
    return locate


_LINEMARK_RE = re.compile(r'^# (\d+) "([^"]*)"')


def pp_locator(text, srcdir):
    """Build offset -> (original file, original line) from the '# N "file"'
    linemarkers gcc emits in -E output."""
    entries, off = [], 0
    cur_file, cur_line = "<unknown>", 0
    for raw in text.split("\n"):
        m = _LINEMARK_RE.match(raw)
        if m:
            cur_line = int(m.group(1))
            cur_file = os.path.relpath(m.group(2), srcdir) if os.path.isabs(m.group(2)) else m.group(2)
        else:
            entries.append((off, cur_file, cur_line))
            cur_line += 1
        off += len(raw) + 1

    def locate(o):
        lo, hi = 0, len(entries) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if entries[mid][0] <= o:
                lo = mid
            else:
                hi = mid - 1
        return entries[lo][1], entries[lo][2]
    return locate


# ----------------------------------------------------------------- the sweep --
def pass_a(srcdir):
    """Raw-source parse: precise file:line for every asm block in the tree."""
    out = []
    for path in sorted(glob.glob(os.path.join(srcdir, "*.c")) +
                       glob.glob(os.path.join(srcdir, "*.h"))):
        rel = os.path.basename(path)
        text = open(path, encoding="utf-8", errors="replace").read()
        out.extend(scan_text(text, raw_locator(rel, text)))
    return out


def pass_b(srcdir, modes, cc, jobs):
    """Preprocess each TU per ISA mode; return {signature: set(modes)} for the
    asm blocks that are actually live, plus any preprocessor failures."""
    live, errors = {}, []
    files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(srcdir, "*.c")))
    work = [(mode, f) for mode in modes for f in files]

    def one(job):
        mode, f = job
        flags = [cc] + COMMON_FLAGS + ISA_FLAGS[mode]
        try:
            p = subprocess.run(flags + [f], cwd=srcdir,
                               capture_output=True, text=True, timeout=300)
        except Exception as e:                                      # noqa: BLE001
            return mode, f, str(e), ()
        if p.returncode != 0:
            tail = p.stderr.strip().splitlines()
            return mode, f, (tail[-1] if tail else "rc=%d" % p.returncode), ()
        sigs = {fi.sig for fi in scan_text(p.stdout, pp_locator(p.stdout, srcdir),
                                           want_comments_stripped=False)}
        return mode, f, None, sigs

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for mode, f, err, sigs in pool.map(one, work):
            if err:
                errors.append((mode, f, err))
                continue
            for s in sigs:
                live.setdefault(s, set()).add(mode)
    return live, errors


# --------------------------------------------------------------- baseline io --
def key_of(path, missing, sig):
    return (path, ",".join(missing), sig)


def load_baseline(path):
    base = {}
    if not os.path.exists(path):
        return base
    for raw in open(path, encoding="utf-8"):
        raw = raw.rstrip("\n")
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) != 4:
            continue
        cnt, f, miss, sig = parts
        base[(f, miss, sig)] = int(cnt)
    return base


BASELINE_HEADER = """\
# asm-clobber-lint baseline -- known, accepted registers-used-vs-declared findings.
#
# Each row is:   count <TAB> file <TAB> missing-registers <TAB> block-signature
# The signature covers (registers used, clobbers declared); it deliberately
# omits line numbers so that unrelated edits do not invalidate the baseline.
#
# The lint fails when a finding appears that is not listed here, or when the
# count for a listed row goes up.  Counts going *down* is progress: rerun
#     python3 asm-clobber-lint.py --write-baseline
# and commit the smaller file.
#
# The bulk of this file is AVX-512 mask registers k0-k7, which are never
# declared clobbered anywhere in the tree.  That is a real latent hazard -- the
# compiler may allocate k1-k7 for its own predicates -- but fixing it is a
# large mechanical change tracked separately.  Do not add new rows here to
# silence a new defect; fix the asm.
"""


def write_baseline(path, findings):
    counts = Counter(key_of(f.path, f.missing, f.sig) for f in findings)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(BASELINE_HEADER)
        for (f, miss, sig), c in sorted(counts.items()):
            fh.write("%d\t%s\t%s\t%s\n" % (c, f, miss, sig))
    return counts


# ---------------------------------------------------------------- self-test --
# Each case is (name, source, expected missing registers).  The pairs marked
# "control" differ from the case above them by exactly the thing being tested,
# so a case that passes for the wrong reason shows up as its control passing too.
SELF_TESTS = [
    ("plain omission",
     '__asm__ volatile ("movq %[a],%%rdx \\n\\t" :: [a] "m" (x) : "cc","memory");',
     ["rdx"]),
    ("operand-pinned output (control for the above)",
     '__asm__ volatile ("movq %[a],%%rdx \\n\\t" : "=d" (o) : [a] "m" (x) : "cc");',
     []),
    ("operand-pinned input",
     '__asm__ volatile ("loop %%rcx \\n\\t" :: "c" (n) : "cc");',
     []),
    ("A constraint pins rax and rdx",
     '__asm__ volatile ("mulq %%rdx \\n\\t" : "=A" (o) : "m" (x));',
     []),
    ("operand references are not registers",
     '__asm__ volatile ("movq %0,%1 \\n\\t movq %[nm],%2" : "=r"(a),"=r"(b),"=r"(c) : [nm] "m"(d));',
     []),
    ("32-bit alias counts as the 64-bit register",
     '__asm__ volatile ("movl %[a],%%eax \\n\\t" :: [a] "m" (x) : "cc");',
     ["rax"]),
    ("32-bit alias satisfied by 64-bit clobber (control)",
     '__asm__ volatile ("movl %[a],%%eax \\n\\t" :: [a] "m" (x) : "cc","rax");',
     []),
    ("byte alias of an extended GPR",
     '__asm__ volatile ("setz %%r11b \\n\\t" :: [a] "m" (x) : "cc");',
     ["r11"]),
    ("zmm satisfied by an xmm clobber",
     '__asm__ volatile ("vmovaps (%%rax),%%zmm18 \\n\\t" :: [a] "m" (x) : "rax","xmm18");',
     []),
    ("zmm not satisfied by a lower-numbered clobber (control)",
     '__asm__ volatile ("vmovaps (%%rax),%%zmm18 \\n\\t" :: [a] "m" (x) : "rax","xmm8");',
     ["xmm18"]),
    ("mask register",
     '__asm__ volatile ("kmovw %%eax,%%k1 \\n\\t" :: [a] "m" (x) : "rax");',
     ["k1"]),
    ("basic asm uses a single %",
     '__asm__ volatile ("movq %rax,%rdx");',
     ["rax", "rdx"]),
    ("rsp is never reported",
     '__asm__ volatile ("movq %%rsp,%%rbx \\n\\t" :: [a] "m" (x) : "rbx");',
     []),
    ("colon inside an operand expression does not split sections",
     '__asm__ volatile ("movq %[a],%%r12 \\n\\t" :: [a] "m" (p ? q : r) : "cc","r12");',
     []),
    ("comment mentioning a register is ignored",
     '__asm__ volatile (/* uses %%rdx */ "movq %[a],%%r13 \\n\\t" :: [a] "m" (x) : "r13");',
     []),
    ("clobber spelled ymm still satisfies a ymm use",
     '__asm__ volatile ("vmovaps (%%rax),%%ymm7 \\n\\t" :: [a] "m" (x) : "rax","ymm7");',
     []),
]


def self_test():
    bad = 0
    for name, src, expect in SELF_TESTS:
        got = scan_text(src, raw_locator("<self-test>", src))
        miss = got[0].missing if got else []
        ok = miss == expect
        bad += not ok
        print("%-4s %-52s expected %-14s got %s"
              % ("ok" if ok else "FAIL", name, ",".join(expect) or "-",
                 ",".join(miss) or "-"))
    print("\nself-test: %d/%d cases pass." % (len(SELF_TESTS) - bad, len(SELF_TESTS)))
    return 1 if bad else 0


# ---------------------------------------------------------------------- main --
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--srcdir", default=os.path.join(here, "src"))
    ap.add_argument("--baseline", default=os.path.join(here, "asm-clobber-baseline.txt"))
    ap.add_argument("--modes", default=",".join(ALL_MODES))
    ap.add_argument("--cc", default=os.environ.get("CC", "cc"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--write-baseline", action="store_true")
    ap.add_argument("--self-test", action="store_true",
                    help="run the parser's own unit tests and exit")
    ap.add_argument("--all", action="store_true",
                    help="report every finding, ignoring the baseline")
    ap.add_argument("--include-dead", action="store_true",
                    help="also report blocks that no ISA mode compiles")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    modes = [m for m in args.modes.split(",") if m]
    for m in modes:
        if m not in ISA_FLAGS:
            sys.exit("unknown ISA mode %r (have: %s)" % (m, ", ".join(ALL_MODES)))

    raw = pass_a(args.srcdir)
    live, pperrs = pass_b(args.srcdir, modes, args.cc, args.jobs)

    for f in raw:
        f.modes = sorted(live.get(f.sig, ()), key=ALL_MODES.index)
    reported = raw if args.include_dead else [f for f in raw if f.modes]

    if args.write_baseline:
        counts = write_baseline(args.baseline, reported)
        print("wrote %s: %d rows, %d findings" % (args.baseline, len(counts), len(reported)))
        return 0

    base = load_baseline(args.baseline)
    counts = Counter(key_of(f.path, f.missing, f.sig) for f in reported)
    new = []
    for f in reported:
        k = key_of(f.path, f.missing, f.sig)
        if counts[k] > base.get(k, 0):
            new.append(f)

    if args.all:
        show = sorted(reported, key=lambda f: (f.path, f.line))
    else:
        show = sorted(new, key=lambda f: (f.path, f.line))

    for f in show:
        print("%s:%d: inline asm uses %s but neither binds nor clobbers %s  [%s]"
              % (f.path, f.line, ",".join(sorted(f.used, key=sortkey)),
                 ",".join(f.missing), ",".join(f.modes) or "dead"))

    if pperrs:
        # Not a warning.  A translation unit that will not preprocess is one this
        # lint cannot see into, and silently shrinking its own coverage is how a
        # CI check ends up reporting success while checking nothing.
        print("\n%d translation unit(s) failed to preprocess -- coverage is "
              "incomplete, so this run cannot be trusted:" % len(pperrs), file=sys.stderr)
        for mode, f, msg in pperrs[:20]:
            print("  %-7s %-32s %s" % (mode, f, msg), file=sys.stderr)
        return 1

    shrunk = sum(1 for k, c in base.items() if counts.get(k, 0) < c)
    if shrunk:
        print("\nnote: %d baseline row(s) now have fewer findings than recorded; "
              "rerun with --write-baseline to shrink the baseline." % shrunk)

    if args.all:
        print("\n%d finding(s) total, %d baseline row(s)." % (len(reported), len(base)))
        return 0
    if new:
        print("\n%d inline-asm clobber finding(s) not covered by %s."
              % (len(new), os.path.basename(args.baseline)))
        return 1
    print("asm-clobber-lint: %d finding(s), all covered by the baseline." % len(reported))
    return 0


if __name__ == "__main__":
    sys.exit(main())
