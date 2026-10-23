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

Exit status is 1 if any finding is not covered by the baseline, 0 otherwise. The baseline
counts the asm blocks each key is found in, so the same finding in one more block is not
covered. Baseline entries found in fewer blocks than listed (or none) are reported but do not
fail the run, so a fix leaves a nudge to lower or drop the line rather than silently keeping
a dead exemption.

Caveats: x86 AT&T syntax only; access width comes from the mnemonic suffix; operands whose
expression is not a plain identifier declared in the same file cannot be type-resolved and are
skipped rather than guessed at.
"""
import argparse, bisect, collections, glob, os, re, sys, tempfile

WIDTH = {'q': 8, 'l': 4, 'w': 2, 'b': 1}
MNEM = re.compile(r'\b(mov|movs|movz|add|adc|sub|sbb|and|or|xor|cmp|test|imul|mul|inc|dec|neg|not'
                  r'|shl|shr|sar|rol|ror|shld|shrd|bt|bts|btr|xchg|lea|push|pop)(q|l|w|b)\b')
# Sign- and zero-extending moves carry two suffixes: the source's width, then the destination
# register's. The source is the only operand that can be in memory.
EXTMOV = re.compile(r'\b(?:movz|movs)(?:b[wlq]|w[lq])\b|\bmovslq\b')
MEMOP = re.compile(r'(-?0x[0-9a-fA-F]+|-?\d+)?\(%%?(r[a-z0-9]+|e[a-z]{2})\s*(?:,[^)]*)?\)')
NAMED = re.compile(r'%\[\s*(\w+)\s*\]')
OPDECL = re.compile(r'\[\s*(\w+)\s*\]\s*"([^"]*)"\s*\(\s*([A-Za-z_]\w*)\s*\)')
TSIZE = {'uint64': 8, 'sint64': 8, 'int64': 8, 'uint32': 4, 'int32': 4, 'sint32': 4, 'int': 4,
         'uint': 4, 'uint16': 2, 'sint16': 2, 'uint8': 1, 'sint8': 1, 'double': 8, 'float': 4,
         'size_t': 8, 'intptr_t': 8}
TYPE_RE = re.compile(r'\b(' + '|'.join(sorted(TSIZE, key=len, reverse=True)) + r')\b')
# After the type: any stars and qualifiers, the declared name, then '[' (an array) or '(' (a
# function, not a variable).
DECLARATOR = re.compile(r'((?:\s|\*|\b(?:const|volatile|restrict|__restrict|__restrict__)\b)*)'
                        r'([A-Za-z_]\w*)\s*([\[(])?')


ASM_KW = re.compile(r'\b(?:__asm__|__asm|asm)\b(?:\s*\b(?:__volatile__|volatile|goto)\b)*\s*\(')


def asm_blocks(path):
    """Return [(1-based start line, block text, all lines)] for each inline-asm statement.

    The block runs from the asm keyword to the parenthesis that closes the statement, found by
    counting parentheses outside comments and string literals - so a closing line such as
    ': "cc","memory");' ends it, and a '(' in a comment or in the asm text itself does not
    count. Comments are blanked out of the returned text."""
    raw = open(path, errors='replace').read()
    text = strip_c(raw, keep_strings=True)
    out, pos = [], 0
    while True:
        m = ASM_KW.search(text, pos)
        if not m:
            break
        i, depth = m.end(), 1
        while i < len(text) and depth:
            c = text[i]
            if c in '"\'':
                j = i + 1
                while j < len(text) and text[j] not in (c, '\n'):
                    j += 2 if text[j] == '\\' else 1
                i = j
            elif c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            i += 1
        line0 = text.count('\n', 0, m.start())
        bol = text.rfind('\n', 0, m.start()) + 1
        eol = text.find('\n', i)
        eol = len(text) if eol < 0 else eol
        out.append((line0 + 1, text[bol:eol], raw.split('\n')))
        pos = i
    return out


def asm_text(block):
    return ' '.join(re.findall(r'"((?:[^"\\]|\\.)*)"', block))


def strip_c(text, keep_strings=False):
    """Blank out comments (and, unless keep_strings, the contents of string and character
    literals) with spaces, keeping every newline so that offsets and line numbers still hold.
    Backslash-newline macro continuations become plain newlines."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith('/*', i):
            j = text.find('*/', i + 2)
            j = n if j < 0 else j + 2
            out.append(re.sub(r'[^\n]', ' ', text[i:j]))
        elif text.startswith('//', i):
            j = text.find('\n', i)
            j = n if j < 0 else j
            out.append(' ' * (j - i))
        elif c in '"\'':
            j = i + 1
            while j < n and text[j] not in (c, '\n'):
                j += 2 if text[j] == '\\' else 1
            j = min(j, n)
            body = text[i + 1:j]
            if not keep_strings:
                body = re.sub(r'[^\n]', ' ', body)
            close = c if j < n and text[j] == c else ''
            out.append(c + body + close)
            j += len(close)
        else:
            j = i + 1
            out.append(c)
        i = j
    return ''.join(out).replace('\\\n', ' \n')


def parse_decls(stripped, base, filescope, out):
    """Collect the declarations in one statement (text between ';', '{' and '}').

    Parenthesised groups are parsed recursively as their own lists, so each parameter of a
    function's parameter list carries its own type, and a cast or a call cannot lend its type
    to anything around it. At the top level, a comma-separated part that names a type starts a
    new declaration; a part that does not is a further declarator of the same type, as in
    'uint32 *p, v;' - and the stars belong to the declarator they precede, so v is a uint32.
    Appends (offset, ident, typename, size, filescope) to out."""
    flat, depth, start, groups = list(stripped), 0, 0, []
    for i, ch in enumerate(stripped):
        if ch == '(':
            if depth == 0:
                start = i
            depth += 1
        elif ch == ')' and depth:
            depth -= 1
            if depth == 0:
                groups.append((start, i))
    if depth:
        groups.append((start, len(stripped)))
    for s, e in groups:
        parse_decls(stripped[s + 1:e], base + s + 1, False, out)
        for k in range(s + 1, e):
            if flat[k] != '\n':
                flat[k] = ' '
    flat = ''.join(flat)
    ty, pos = None, 0
    for part in flat.split(','):
        off = pos
        pos += len(part) + 1
        tm = TYPE_RE.search(part)
        if tm:
            prev = part[:tm.start()].split()
            ty = None if prev and prev[-1] in ('long', 'short') else tm.group(1)
            d, doff = part[tm.end():], off + tm.end()
        else:
            d, doff = part, off
        if not ty:
            continue
        dm = DECLARATOR.match(d)
        if not dm or dm.group(3) == '(':
            continue
        isptr = '*' in dm.group(1) or dm.group(3) == '['
        out.append((base + doff + dm.start(2), dm.group(2), ty + ('*' if isptr else ''),
                    8 if isptr else TSIZE[ty], filescope))


def file_decls(text):
    """Every declaration in a file, comments and string literals stripped first so that prose
    and printf formats are not read as declarations. Returns (decls, ends): decls is a list of
    (line, ident, typename, size, filescope), filescope meaning outside any braces and any
    parameter list; ends lists the lines that start with '}' in column 0, which in this tree is
    how a function ends."""
    stripped = strip_c(text)
    raw, pos, depth = [], 0, 0
    for m in re.finditer(r'[;{}]', stripped + ';'):
        parse_decls(stripped[pos:m.start()], pos, depth == 0, raw)
        if m.group() == '{':
            depth += 1
        elif m.group() == '}':
            # A column-0 brace closes a function: resynchronise, in case an #if/#else pair
            # opened a brace twice.
            col0 = m.start() == 0 or stripped[m.start() - 1] == '\n'
            depth = 0 if col0 else max(depth - 1, 0)
        pos = m.end()
    nl = [i for i, ch in enumerate(stripped) if ch == '\n']
    decls = [(bisect.bisect_left(nl, off) + 1, ident, ty, size, fs)
             for off, ident, ty, size, fs in raw]
    ends = [n + 1 for n, line in enumerate(stripped.split('\n')) if line.startswith('}')]
    return decls, ends


def decl_type(fdecls, ident, before):
    """Resolve the declared type of `ident` as seen from line `before`: the nearest declaration
    in the enclosing function (everything after the previous column-0 '}', which takes in the
    function's own parameter list), else the nearest file-scope one. A same-named local of an
    earlier function is never used."""
    decls, ends = fdecls
    i = bisect.bisect_left(ends, before)
    start = ends[i - 1] if i else 0
    best = None
    for line, name, ty, size, fs in decls:
        if name != ident or line >= before:
            continue
        if line > start or fs:
            rank = (line > start, line)
            if best is None or rank >= best[0]:
                best = (rank, (ty, size, line))
    return best and best[1]


def split_first_operand(ops):
    """Split an AT&T operand list at its first top-level comma - not one inside the
    parentheses of a memory operand such as 0x8(%rsi,%rax,8)."""
    depth = 0
    for i, ch in enumerate(ops):
        depth += (ch == '(') - (ch == ')')
        if ch == ',' and depth == 0:
            return ops[:i], ops[i + 1:]
    return ops, ''


def scan(srcdir):
    """Return (findings, nblocks, nsimd). Each finding is (key, human-readable detail)."""
    findings, nblocks, nsimd = [], 0, 0
    files = sorted(glob.glob(os.path.join(srcdir, '*.c')) + glob.glob(os.path.join(srcdir, '*.h')))
    for path in files:
        rel = os.path.relpath(path, os.path.dirname(srcdir.rstrip('/')) or '.')
        fdecls = None
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
                x = EXTMOV.search(stmt)
                if x and (not m or x.start() <= m.start()):
                    # movzbl and kin: the source operand is accessed at the first suffix's
                    # width, the destination register at the second's.
                    mn = x.group()
                    srcop, dstop = split_first_operand(stmt[x.end():])
                    parts = [(srcop, WIDTH[mn[-2]]), (dstop, WIDTH[mn[-1]])]
                # lea computes an address and touches no memory, so it has no access width.
                elif m and m.group(1) != 'lea':
                    parts = [(stmt, WIDTH[m.group(2)])]
                else:
                    continue
                for text, w in parts:
                    if not is_simd:
                        for mm in MEMOP.finditer(text):
                            off = mm.group(1) or '0'
                            off = int(off, 16) if 'x' in off.lower() else int(off)
                            acc[(mm.group(2), off)].add(w)
                    for n in NAMED.finditer(text):
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
                if fdecls is None:
                    fdecls = file_decls('\n'.join(lines))
                t = decl_type(fdecls, ident, ln)
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
# Each line is '<count> <key>'. The key carries no line number, so unrelated edits do not
# invalidate it; the count is the number of asm blocks that produce it, so a new block with the
# same finding in the same file is still caught. The lint fails if a key is found in more
# blocks than its count.
#
# The 'A' entries in twopmodq96.c are deliberate: the asm reads those uint96 d1 slots 8 bytes
# wide, which also touches the struct padding, and that is safe only because the local store is
# zeroed at allocation. Do not add to this list without the same guarantee.
#
# Regenerate with:
#     python3 type-width-lint.py --write-baseline
"""


def load_baseline(path):
    """Return {key: count}. A line with no leading count allows one block."""
    out = collections.Counter()
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.strip()
        if line and not line.startswith('#'):
            m = re.match(r'(\d+)\s+(\S.*)$', line)
            out[m.group(2) if m else line] += int(m.group(1)) if m else 1
    return out


def write_baseline(path, findings):
    counts = collections.Counter(k for k, _ in findings)
    with open(path, 'w') as f:
        f.write(BASELINE_HEADER)
        for key in sorted(counts):
            f.write(f"{counts[key]} {key}\n")


def compare(findings, baseline):
    """Return (over, under). over lists (key, blocks found, blocks allowed, details) for each
    key found in more blocks than the baseline allows - a key absent from it allows none;
    under lists (key, found, allowed) for each key found in fewer."""
    found = collections.Counter(k for k, _ in findings)
    details = collections.defaultdict(list)
    for k, d in findings:
        details[k].append(d)
    over = [(k, found[k], baseline[k], details[k])
            for k in sorted(found) if found[k] > baseline[k]]
    under = [(k, found[k], n) for k, n in sorted(baseline.items()) if found[k] < n]
    return over, under


# (file name, what it shows, the exact finding keys it must produce, source). Each case is
# scanned as its own file, and must produce exactly the keys listed: no fewer, no more.
SELF_TEST_CASES = [
    ('a_hit.c', 'A fires on a mixed-width slot', ['A src/a_hit.c %rsi+0x8 [4, 8]'], '''
void f(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
'''),
    ('a_miss.c', 'A silent on a uniform-width slot', [], '''
void f(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		"movq 0x8(%%rsi),%%rdx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
'''),
    ('b_hit.c', 'B fires on movl of a uint64', ['B src/b_hit.c %[__p] pp uint64 8 4'], '''
void f(void) {
	uint64 pp;
	__asm__ ("movl %[__p],%%eax	\\n\\t"
		:: [__p] "m" (pp) : "cc","rax");
}
'''),
    ('b_miss.c', 'B silent when widths match the types', [], '''
void f(void) {
	uint64 pp;
	uint32 qq;
	__asm__ ("movq %[__p],%%rax	\\n\\t"
		"movl %[__q],%%ecx	\\n\\t"
		:: [__p] "m" (pp), [__q] "m" (qq) : "cc","rax","rcx");
}
'''),
    ('simd_skip.c', 'SIMD blocks are skipped by A', [], '''
void f(void *p) {
	__asm__ ("vmovaps 0x8(%%rsi),%%ymm0	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
'''),
    # Each parameter carries its own type: len is a uint32 even though the list began with a
    # uint64. Both a mismatch (movq) and a match (movl) on it, so the case fails either way if
    # len were resolved as 8 bytes.
    ('b_params.c', 'B types each parameter of a list separately',
     ['B src/b_params.c %[__len] len uint32 4 8'], '''
void f(const uint64 x[], uint32 len)
{
	__asm__ ("movq %[__len],%%rcx	\\n\\t"
		"movl %[__len],%%edx	\\n\\t"
		:: [__len] "m" (len) : "cc","rcx","rdx");
}
'''),
    # Prose in a comment is not a declaration: 'int x[] and a' must not make a an int.
    ('b_comment.c', 'B ignores declarations in comments', [], '''
void f(uint64 a)
{
	/* Given a vector int x[] and a dimension (len), ... */
	__asm__ ("movq %[__a],%%rax	\\n\\t"
		:: [__a] "m" (a) : "cc","rax");
}
'''),
    # The star belongs to its own declarator: v is a 4-byte uint32, not a pointer.
    ('b_declarator.c', 'B applies a pointer star to its own declarator only',
     ['B src/b_declarator.c %[__v] v uint32 4 8'], '''
void f(void)
{
	uint32 *ptr, v;
	__asm__ ("movq %[__v],%%rax	\\n\\t"
		"movq %[__ptr],%%rcx	\\n\\t"
		:: [__v] "m" (v), [__ptr] "m" (ptr) : "cc","rax","rcx");
}
'''),
    # Scope: g's local uint64 p is not visible in f, which sees the file-scope uint32 p declared
    # above g - so movq on it is a mismatch. f's own parameter q shadows the global uint64 q.
    ('b_scope.c', 'B resolves names in the enclosing function, else at file scope',
     ['B src/b_scope.c %[__p] p uint32 4 8', 'B src/b_scope.c %[__q] q uint32 4 8'], '''
uint32 p;
uint64 q;
void g(void)
{
	uint64 p = 0;
	(void)p;
}
void f(uint32 q)
{
	__asm__ ("movq %[__p],%%rax	\\n\\t"
		"movq %[__q],%%rcx	\\n\\t"
		:: [__p] "m" (p), [__q] "m" (q) : "cc","rax","rcx");
}
'''),
    # Two asm statements, each ending in a line that carries the clobber list: they are two
    # blocks, so the 8-byte read in one and the 4-byte read in the other are not one slot.
    ('a_two_blocks.c', 'A treats adjacent asm statements as separate blocks', [], '''
void f(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		:: [__x] "m" (p) : "cc", "memory");
	__asm__ ("movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc", "memory");
}
'''),
    # lea only computes an address: neither its memory operand nor a named one is an access.
    ('lea.c', 'lea is not a memory access, for A or for B', [], '''
void f(void)
{
	uint32 p;
	__asm__ ("leaq 0x8(%%rsi),%%rax	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		"leaq %[__p],%%rcx	\\n\\t"
		:: [__p] "m" (p) : "cc","rax","rcx","rdx");
}
'''),
    # Two-suffix extending moves read their source at the first suffix's width: movzbl reads
    # one byte, so it is a mismatch on a uint32 and, against an 8-byte read of the same slot,
    # an A finding; movslq reads four, which matches a uint32.
    ('extmov.c', 'two-suffix moves (movzbl, movslq, ...) access their source at the first width',
     ['A src/extmov.c %rsi+0x8 [1, 8]', 'B src/extmov.c %[__c] c uint32 4 1'], '''
void f(void)
{
	uint32 c, n;
	__asm__ ("movzbl %[__c],%%eax	\\n\\t"
		"movslq %[__n],%%rcx	\\n\\t"
		"movzbl 0x8(%%rsi),%%edx	\\n\\t"
		"movq 0x8(%%rsi),%%r8	\\n\\t"
		:: [__c] "m" (c), [__n] "m" (n) : "cc","rax","rcx","rdx","r8");
}
'''),
    # The same finding in two blocks of one file: one baseline entry must not cover both.
    ('baseline_count.c', 'a finding in two blocks is counted twice',
     ['A src/baseline_count.c %rsi+0x8 [4, 8]'], '''
void f(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
void g(void *p) {
	__asm__ ("movq 0x8(%%rsi),%%r12	\\n\\t"
		"movl 0x8(%%rsi),%%edx	\\n\\t"
		:: [__x] "m" (p) : "cc","memory");
}
'''),
]


def self_test():
    ok = True
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, 'src')
        os.makedirs(src)
        for name, what, expect, text in SELF_TEST_CASES:
            open(os.path.join(src, name), 'w').write(text)
        findings, _, _ = scan(src)
        for name, what, expect, text in SELF_TEST_CASES:
            got = sorted({k for k, _ in findings if k.split()[1] == 'src/' + name})
            passed = got == sorted(expect)
            print(f"  {'ok  ' if passed else 'FAIL'} - {what}")
            if not passed:
                print(f"         expected {sorted(expect)}\n         got      {got}")
            ok &= passed

        # Baseline: a count of 1 does not cover the key's second block; a count of 2 does.
        key = 'A src/baseline_count.c %rsi+0x8 [4, 8]'
        mine = [(k, d) for k, d in findings if k == key]
        bl = os.path.join(td, 'baseline.txt')
        for count, want_over in ((1, True), (2, False)):
            open(bl, 'w').write(f"# comment\n{count} {key}\n")
            over, _ = compare(mine, load_baseline(bl))
            passed = bool(over) == want_over
            print(f"  {'ok  ' if passed else 'FAIL'} - a baseline count of {count} "
                  f"{'does not cover' if want_over else 'covers'} a finding in two blocks")
            ok &= passed
    if not ok:
        print("\nself-test FAILED: the check does not behave as documented", file=sys.stderr)
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

    baseline = collections.Counter() if args.all else load_baseline(args.baseline)
    over, under = compare(findings, baseline)
    nnew = sum(found - allowed for _, found, allowed, _ in over)

    print(f"asm blocks scanned: {nblocks} ({nsimd} SIMD, skipped by check A)")
    print(f"findings: {len(findings)} total, {nnew} not covered by the baseline")

    if over:
        print("\nNot in the baseline:")
        for key, found, allowed, details in over:
            if allowed:
                print(f"  {key}: found in {found} blocks, the baseline allows {allowed}:")
            for d in details:
                print(f"  {'  ' if allowed else ''}{d}")
        print("\nA: a slot written at one width and read at another is an aliasing hazard - either")
        print("   make the access widths agree, or guarantee the padding and add it to the baseline.")
        print("B: the asm operand's access width disagrees with the C type behind it, which")
        print("   silently drops the high bits. Use the mnemonic suffix that matches the type.")
    if under:
        print(f"\nBaseline entries found in fewer blocks than listed ({len(under)}) - "
              f"lower the count, or drop the line at 0:")
        for key, found, allowed in under:
            print(f"  {key}: found in {found}, listed as {allowed}")
    return 1 if over else 0


if __name__ == '__main__':
    sys.exit(main())
