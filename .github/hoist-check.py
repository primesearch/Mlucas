import re,sys
# Per uint96 d1 slot (+0x8/0x18/0x28/0x38 off the per-thread local store), find the FIRST store whose
# written byte range covers that slot, and the FIRST 32-bit load of it, and compare their addresses.
# Stores may be 8-byte GPR movs OR 16/32/64-byte vector stores (gcc 15 vectorises the store group),
# so width is taken from the instruction. %rsp-based accesses are frame spills, not the local store.
WIDTH={'xmm':16,'ymm':32,'zmm':64}
SLOTS=(0x8,0x18,0x28,0x38)
# Source register: 64-bit (rax, r11) or 32-bit (eax, edi, r11d). The 32-bit forms matter - they
# are the type-correct store #312 introduces - so e?? names must be matched, not just r?? and r??d.
gpr_st=re.compile(r'^%(r[a-z0-9]+|[a-d]x|[sd]i|e[a-z]{2}),(-?0x[0-9a-f]+)\(%(r[a-z0-9]+)\)$')
vec_st=re.compile(r'^%([xyz]mm)\d+,(-?0x[0-9a-f]+)\(%(r[a-z0-9]+)\)$')
ld32  =re.compile(r'^(-?0x[0-9a-f]+)\(%(r[a-z0-9]+)\),%(e[a-z]{2}|r[0-9]+d)$')
def analyse(path,sym,label):
    lines=open(path).read().split('\n')
    i=[n for n,l in enumerate(lines) if l.startswith(sym)][0]
    stores=[]; loads={}; narrow={}
    for l in lines[i+1:]:
        # NB: the inline asm defines its own labels (twopmodq96_q4_pshiftjmp*), which objdump
        # prints as symbol headers. Stopping at the first one truncates the function mid-asm, so
        # only break on a symbol that is not one of those internal labels.
        m0=re.match(r'^[0-9a-f]+ <([^>]*)>:',l)
        if m0 and not m0.group(1).startswith(sym.split('<')[1].rstrip('>:')): break
        m=re.match(r'\s*([0-9a-f]+):\s+(\S+)\s+(.*)',l)
        if not m: continue
        a,mn,ops=int(m.group(1),16),m.group(2),m.group(3).split('#')[0].strip()
        x=gpr_st.match(ops)
        if mn=='mov' and x and x.group(3)!='rsp':
            if x.group(1).endswith('d') or x.group(1).startswith('e'):   # 32-bit: type-correct case
                if int(x.group(2),16) in SLOTS: narrow.setdefault(int(x.group(2),16),[]).append((a,x.group(3)))
            else:
                stores.append((a,int(x.group(2),16),8,x.group(3),mn))
        x=vec_st.match(ops)
        if x and x.group(3)!='rsp':
            stores.append((a,int(x.group(2),16),WIDTH[x.group(1)],x.group(3),mn))
        x=ld32.match(ops)
        if mn=='mov' and x and x.group(2)!='rsp':
            loads.setdefault(int(x.group(1),16),[]).append((a,x.group(2)))
    # the local-store base = the base used by stores covering the most distinct slot offsets
    from collections import defaultdict
    cov=defaultdict(set)
    for a,off,w,b,mn in stores:
        for s in SLOTS:
            if off<=s<off+w: cov[b].add(s)
    if not cov:
        # No wide store covering any slot. Either the build has #312's type-correct stores - a
        # 4-byte store into a 4-byte ->d1, which is exactly what removes the defect - or the probe
        # is blind. Those are opposite conclusions, so say which.
        if narrow:
            print(f"{label}")
            for s in SLOTS:
                n=narrow.get(s)
                print(f"   +{s:#5x}: " + (f"4-byte store @{min(n)[0]:x} - same width as the load" if n else "no store seen"))
            print("   VERDICT: type-correct stores (#312 applied) - no width mismatch to hoist against")
        else:
            print(f"{label}\n  PROBE FOUND NO STORES AND NO NARROW STORES - probe is blind, no conclusion")
        return
    base=max(cov,key=lambda b:len(cov[b]))
    print(f"{label}   (store base %{base}, covers {len(cov[base])}/4 slots)")
    bad=0; miss=0
    for s in SLOTS:
        cand=[(a,w,mn) for a,off,w,b,mn in stores if b==base and off<=s<off+w]
        l=loads.get(s)
        if not cand or not l: print(f"   +{s:#5x}: INCOMPLETE store={cand} load={l}"); miss+=1; continue
        sa,w,mn=min(cand); la,lb=min(l)
        h=la<sa; bad+=h
        print(f"   +{s:#5x}: store @{sa:x} ({mn}, {w}B)   first 32-bit load @{la:x}(%{lb})  -> {'HOISTED' if h else 'ok'}")
    print("   VERDICT:", "INCOMPLETE - no conclusion" if miss else (f"** {bad} LANE(S) HOISTED **" if bad else "no hoist"))
analyse(sys.argv[1],sys.argv[2],sys.argv[3])
