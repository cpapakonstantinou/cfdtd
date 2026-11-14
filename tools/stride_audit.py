#!/usr/bin/env python3
# usage: tools/stride_audit.py <repo root>, see docs/dev/Roadmap.md, "Row-contiguity audit"
# For every innermost for-loop, list array accesses whose contiguous, last, subscript does not
# use the loop variable while an earlier subscript does: those walk memory with a stride.
import re, sys, glob, collections
FOR = re.compile(r'\bfor\s*\(\s*(?:int\s+|size_t\s+)?(\w+)\s*=')
ACC = re.compile(r'\b([A-Za-z_]\w*)((?:\s*\[[^\[\]]*(?:\[[^\[\]]*\][^\[\]]*)*\])+)')
def subs(s):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == '[':
            if depth == 0: cur = ''
            else: cur += ch
            depth += 1
        elif ch == ']':
            depth -= 1
            if depth == 0: out.append(cur)
            else: cur += ch
        elif depth > 0: cur += ch
    return out
def body_of(src, start):
    # body of the loop whose 'for' is at start: either {...} or a single statement
    i = src.index('(', start); d = 0
    while True:
        if src[i] == '(': d += 1
        elif src[i] == ')':
            d -= 1
            if d == 0: break
        i += 1
    j = i + 1
    while src[j].isspace(): j += 1
    if src[j] == '{':
        d = 0; k = j
        while True:
            if src[k] == '{': d += 1
            elif src[k] == '}':
                d -= 1
                if d == 0: return src[j:k+1]
            k += 1
    return src[j:src.index(';', j)+1]
report = collections.defaultdict(list)
for f in sorted(glob.glob(sys.argv[1] + '/src/*.cc')):
    src = open(f).read()
    src = re.sub(r'//[^\n]*', '', src); src = re.sub(r'/\*.*?\*/', '', src, flags=re.S)
    for m in FOR.finditer(src):
        var = m.group(1)
        try: body = body_of(src, m.start())
        except Exception: continue
        if FOR.search(body): continue  # only innermost loops
        line = src.count('\n', 0, m.start()) + 1
        bad = set()
        for a in ACC.finditer(body):
            ss = subs(a.group(2))
            if len(ss) < 2: continue
            uses = [re.search(r'\b%s\b' % var, s) is not None for s in ss]
            if any(uses[:-1]) and not uses[-1]:
                bad.add(a.group(1))
        if bad:
            fn = None
            head = src[:m.start()]
            fm = list(re.finditer(r'\n(?:[\w:<>,\*&\s]+?)\b(\w+)\s*\([^;{]*\)\s*\{', head))
            fn = fm[-1].group(1) if fm else '?'
            report[f.split('/')[-1]].append((line, fn, var, sorted(bad)))
tot = 0
for f, items in report.items():
    for line, fn, var, bad in items:
        tot += 1
        print(f'{f}:{line} {fn}: inner loop {var}: {", ".join(bad[:6])}{" ..." if len(bad) > 6 else ""}')
print('strided innermost loops:', tot, file=sys.stderr)
