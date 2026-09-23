#!/usr/bin/env python3
"""Static check for integer-width mistakes in and around this tree's inline asm.

Two classes, both of which have produced silently wrong answers here:

  A. Same (base register, displacement) touched at more than one access width inside one GPR
     asm block. That is the shape of the uint96 {uint64 d0; uint32 d1;} store/load mismatch:
     the C side writes the 4-byte d1 slot with an 8-byte store while the asm reads it 4 bytes
     wide, or the reverse. Under -fstrict-aliasing the two are in different alias sets and the
     compiler may hoist the load above the store, corrupting one parallel lane and losing a
     factor with no diagnostic. Several sites here read a d1 slot 8 bytes wide deliberately,
     relying on the padding being zero; those are in the baseline, and the point of the check
     is that a *new* one has to be justified rather than appearing unnoticed.

  B. A named asm operand accessed at a width that does not match the C type of the identifier
     behind it - e.g. "movl %[__p]" where p is a uint64. The upper half is silently dropped.
     This is a wrong answer, not a warning: nothing in the compiler checks that the mnemonic
     suffix agrees with the operand's type.

No sanitizer sees either class: there is no C-level operation to instrument for A, and for B
the asm is opaque to the compiler by construction.

Usage:
    python3 type-width-lint.py                    # lint against the baseline
    python3 type-width-lint.py --all              # report everything, ignore the baseline
    python3 type-width-lint.py --write-baseline   # regenerate the baseline
    python3 type-width-lint.py --self-test        # prove the checks fire, and do not over-fire

Exit status is 1 if any finding is not covered by the baseline, 0 otherwise. Baseline entries
that no longer reproduce are reported but do not fail the run, so a fix leaves a nudge to
drop the stale line rather than silently keeping a dead exemption.

Caveats: x86 AT&T syntax only; access width comes from the mnemonic suffix; operands whose
expression is not a plain identifier declared in the same file cannot be type-resolved and are
skipped rather than guessed at.
"""
import argparse, collections, glob, os, re, sys, tempfile

WIDTH = {'q': 8, 'l': 4, 'w': 2, 'b': 1}
MNEM = re.compile(r'\b(mov|movs|movz|add|adc|sub|sbb|and|or|xor|cmp|test|imul|mul|inc|dec|neg|not'
                  r'|shl|shr|sar|rol|ror|shld|shrd|bt|bts|btr|xchg|lea|push|pop)(q|l|w|b)\b')
MEMOP = re.compile(r'(-?0x[0-9a-fA-F]+|-?\d+)?\(%%?(r[a-z0-9]+|e[a-z]{2})\s*(?:,[^)]*)?\)')
NAMED = re.compile(r'%\[\s*(\w+)\s*\]')
OPDECL = re.compile(r'\[\s*(\w+)\s*\]\s*"([^"]*)"\s*\(\s*([A-Za-z_]\w*)\s*\)')
TSIZE = {'uint64': 8, 'sint64': 8, 'int64': 8, 'uint32': 4, 'int32': 4, 'sint32': 4, 'int': 4,
         'uint': 4, 'uint16': 2, 'sint16': 2, 'uint8': 1, 'sint8': 1, 'double': 8, 'float': 4,
         'size_t': 8, 'intptr_t': 8}
DECL = re.compile(r'\b(' + '|'.join(TSIZE) + r')\b\s*([*\s]*)([^;{}()=]*)')


def asm_blocks(path):
    """Yield (1-based start line, block text, all lines) for each inline-asm block."""
    lines = open(path, errors='replace').read().split('\n')
    out, i = [], 0
    while i < len(lines):
        if re.search(r'\b__asm__\b|\basm\s*(volatile)?\s*\(', lines[i]):
            start, buf = i, []
            while i < len(lines) and i - start < 1500:
                buf.append(lines[i])
                if re.match(r'^\s*\)\s*;?\s*\\?\s*$', lines[i]) and len(buf) > 1:
                    break
                if re.search(r'\)\s*;\s*\\?\s*$', lines[i]) and len(buf) > 1 and '"' not in lines[i]:
                    break
                i += 1
            out.append((start + 1, '\n'.join(buf), lines))
        i += 1
    return out


def asm_text(block):
    return ' '.join(re.findall(r'"((?:[^"\\]|\\.)*)"', block))


def decl_type(lines, ident, before):
    """Resolve the declared type of `ident`, searching backwards from line `before`."""
    for n in range(min(before, len(lines)) - 1, -1, -1):
        line = lines[n]
        if '__asm__' in line:
            continue
        for m in DECL.finditer(line):
            ty, stars, rest = m.group(1), m.group(2), m.group(3)
            for tok in re.finditer(r'(\**)\s*\b' + re.escape(ident) + r'\b\s*(\[)?', rest):
                isptr = bool(stars.strip()) or bool(tok.group(1)) or bool(tok.group(2))
                return ty + ('*' if isptr else ''), (8 if isptr else TSIZE[ty]), n + 1
    return None


def scan(srcdir):
    """Return (findings, nblocks, nsimd). Each finding is (key, human-readable detail)."""
    findings, nblocks, nsimd = [], 0, 0
    files = sorted(glob.glob(os.path.join(srcdir, '*.c')) + glob.glob(os.path.join(srcdir, '*.h')))
    for path in files:
        rel = os.path.relpath(path, os.path.dirname(srcdir.rstrip('/')) or '.')
        for ln, block, lines in asm_blocks(path):
            asm = asm_text(block)
            if not asm.strip():
                continue
            nblocks += 1
            is_simd = bool(re.search(r'[xyz]mm', asm))
            nsimd += is_simd
            acc = collections.defaultdict(set)
            named = collections.defaultdict(set)
            for stmt in re.split(r'\\n|;', asm):
                m = MNEM.search(stmt)
                if not m:
                    continue
                w = WIDTH[m.group(2)]
                if not is_simd:
                    for mm in MEMOP.finditer(stmt):
                        off = mm.group(1) or '0'
                        off = int(off, 16) if 'x' in off.lower() else int(off)
                        acc[(mm.group(2), off)].add(w)
                for n in NAMED.finditer(stmt):
                    named[n.group(1)].add(w)

            # A: one slot, two widths.
            for (reg, off), widths in sorted(acc.items()):
                if len(widths) > 1:
                    key = f"A {rel} %{reg}+0x{off:x} {sorted(widths)}"
                    findings.append((key, f"{rel}: asm at line {ln}: %{reg}+0x{off:x} "
                                          f"accessed at widths {sorted(widths)} bytes"))

            # B: operand width vs the C type behind it.
            decls = {d.group(1): d.group(3) for d in OPDECL.finditer(block)}
            for name, widths in sorted(named.items()):
                ident = decls.get(name)
                if not ident:
                    continue
                t = decl_type(lines, ident, ln)
                if not t:
                    continue
                tyname, size, dln = t
                for w in sorted(widths):
                    if w != size:
                        key = f"B {rel} %[{name}] {ident} {tyname} {size} {w}"
                        findings.append((key, f"{rel}: asm at line {ln}: %[{name}] -> {ident}, "
                                              f"declared {tyname} ({size} bytes at line {dln}), "
                                              f"but accessed {w} bytes"))
    return findings, nblocks, nsimd


BASELINE_HEADER = """\
# type-width-lint baseline -- known, accepted findings.
#
# One key per line. Keys carry no line numbers, so unrelated edits do not invalidate them.
#
# The 'A' entries in twopmodq96.c are deliberate: the asm reads those uint96 d1 slots 8 bytes
# wide, which also touches the struct padding, and that is safe only because the local store is
# zeroed at allocation. Do not add to this list without the same guarantee.
#
# Regenerate with:
#     python3 type-width-lint.py --write-baseline
"""


def load_baseline(path):
    if not os.path.exists(path):
        return set()
    out = set()
    for line in open(path):
        line = line.strip()
        if line and not line.startswith('#'):
            out.add(line)
    return out


def write_baseline(path, findings):
    with open(path, 'w') as f:
        f.write(BASELINE_HEADER)
        for key in sorted({k for k, _ in findings}):
            f.write(key + '\n')


SELF_TEST_SOURCES = {
    # Should fire A: same slot at 4 and 8 bytes.
    'a_hit.c': '''
void f(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
''',
    # Should NOT fire A: same width throughout.
    'a_miss.c': '''
void f(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		"movq 0x8(%%rsi),%%rdx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
''',
    # Should fire B: movl on a uint64.
    'b_hit.c': '''
void f(void) {
	uint64 pp;
	__asm__ ("movl %[__p],%%eax	\\n\\t"
		:: [__p] "m" (pp) : "cc","rax");
}
''',
    # Should NOT fire B: movq on a uint64, movl on a uint32.
    'b_miss.c': '''
void f(void) {
	uint64 pp;
	uint32 qq;
	__asm__ ("movq %[__p],%%rax	\\n\\t"
		"movl %[__q],%%ecx	\\n\\t"
		:: [__p] "m" (pp), [__q] "m" (qq) : "cc","rax","rcx");
}
''',
    # Should NOT fire A: SIMD blocks are skipped (vector code legitimately mixes widths).
    'simd_skip.c': '''
void f(void *p) {
	__asm__ ("vmovaps 0x8(%%rsi),%%ymm0	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
''',
}


def self_test():
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, 'src')
        os.makedirs(src)
        for name, text in SELF_TEST_SOURCES.items():
            open(os.path.join(src, name), 'w').write(text)
        findings, _, _ = scan(src)
        keys = ' '.join(k for k, _ in findings)
        checks = [
            ('A fires on a mixed-width slot',        'a_hit.c' in keys and keys.count('A ') >= 1),
            ('A silent on a uniform-width slot',     'a_miss.c' not in keys),
            ('B fires on movl of a uint64',          'b_hit.c' in keys),
            ('B silent when widths match the types', 'b_miss.c' not in keys),
            ('SIMD blocks are skipped by A',         'simd_skip.c' not in keys),
        ]
        ok = True
        for what, passed in checks:
            print(f"  {'ok  ' if passed else 'FAIL'} - {what}")
            ok &= passed
        if not ok:
            print("\nself-test FAILED: the check does not behave as documented", file=sys.stderr)
            for _, detail in findings:
                print(f"    got: {detail}", file=sys.stderr)
        return 0 if ok else 1


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--srcdir', default=os.path.join(here, 'src'))
    ap.add_argument('--baseline', default=os.path.join(here, 'type-width-baseline.txt'))
    ap.add_argument('--write-baseline', action='store_true')
    ap.add_argument('--all', action='store_true', help='report everything, ignoring the baseline')
    ap.add_argument('--self-test', action='store_true',
                    help='check the checks against known-positive and known-negative cases')
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    findings, nblocks, nsimd = scan(args.srcdir)
    if args.write_baseline:
        write_baseline(args.baseline, findings)
        print(f"wrote {args.baseline} with {len({k for k, _ in findings})} entries")
        return 0

    baseline = set() if args.all else load_baseline(args.baseline)
    new = [(k, d) for k, d in findings if k not in baseline]
    seen = {k for k, _ in findings}
    stale = sorted(baseline - seen)

    print(f"asm blocks scanned: {nblocks} ({nsimd} SIMD, skipped by check A)")
    print(f"findings: {len(findings)} total, {len(new)} not covered by the baseline")

    if new:
        print("\nNot in the baseline:")
        for _, detail in sorted(set(new)):
            print(f"  {detail}")
        print("\nA: a slot written at one width and read at another is an aliasing hazard - either")
        print("   make the access widths agree, or guarantee the padding and add it to the baseline.")
        print("B: the asm operand's access width disagrees with the C type behind it, which")
        print("   silently drops the high bits. Use the mnemonic suffix that matches the type.")
    if stale:
        print(f"\nBaseline entries that no longer reproduce ({len(stale)}) - drop them:")
        for k in stale:
            print(f"  {k}")
    return 1 if new else 0


if __name__ == '__main__':
    sys.exit(main())
