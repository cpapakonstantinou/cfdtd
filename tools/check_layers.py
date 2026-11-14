#!/usr/bin/env python3
"""Report dependencies that go against the solver's layering (docs/dev/Roadmap.md §2.5).

Layers, bottom to top:
  0 primitives   core, io (h5io), geometry/csg
  1 core library grid, fields, materials (with its nonlinear, hydro and thermal parts), boundary,
                 update, sources, monitors, geometry
  2 application  config, app, main.cc, tools/

For every file of a core module it lists
  - includes of application headers (app/, config/),
  - lines naming the application's run configuration and run state
    (run_config(), run_state(), RunConfig, RunState, fdtd::legacy),
  - field reads of them (run_cfg.<group>.<field>, run_state().<field>),
and the include edges between core modules that go against the module order
  grid < fields < materials < boundary < update < {sources, monitors, geometry}.

Usage: tools/check_layers.py [--files] [--strict MODULE ...]
  --files   list the offending lines per file, not only the counts
  --strict  exit 1 if one of the named modules has any upward dependency (for modules that are clean)
"""
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIMITIVE = {'core', 'io'}
CORE = ['grid', 'fields', 'materials', 'boundary', 'update', 'sources', 'monitors', 'geometry']
APP = {'app', 'config'}
# rank in the module order; modules with the same rank may include each other
RANK = {'grid': 0, 'fields': 1, 'materials': 2, 'boundary': 3, 'update': 4}
for m in CORE[5:]:
    RANK[m] = 5
INCLUDE = re.compile(r'^\s*#\s*include\s*"([a-z_]+)/[^"]+"')
CFG_READ = re.compile(r'\brun_cfg\.\w+\.\w+|\brun_config\(\)\.\w+\.\w+|\brun_state\(\)\.\w+')
APP_USE = re.compile(r'\b(run_config\s*\(|run_state\s*\(|RunConfig\b|RunState\b|legacy::|fdtd::legacy\b)')


def module_of(path):
    parts = path.split(os.sep)
    if parts[0] in ('inc', 'src') and len(parts) > 2:
        return parts[1]
    return None


def main(argv):
    show = '--files' in argv
    strict = argv[argv.index('--strict') + 1:] if '--strict' in argv else []
    upward = defaultdict(list)      # module -> [(file, line, text)]
    app_use = defaultdict(list)
    order = defaultdict(list)
    reads = defaultdict(int)       # module -> configuration / run-state field reads
    for top in ('inc', 'src'):
        for d, _, fs in os.walk(os.path.join(ROOT, top)):
            for f in fs:
                if not f.endswith(('.h', '.cc')):
                    continue
                full = os.path.join(d, f)
                rel = os.path.relpath(full, ROOT)
                mod = module_of(rel)
                if mod not in CORE and mod not in PRIMITIVE:
                    continue
                text = open(full, errors='ignore').read()
                code = re.sub(r'/\*.*?\*/', lambda m: '\n' * m.group().count('\n'), text, flags=re.S)
                for n, line in enumerate(code.split('\n'), 1):
                    line = re.sub(r'//.*', '', line)
                    reads[mod] += len(CFG_READ.findall(line))
                    m = INCLUDE.match(line)
                    if m:
                        dep = m.group(1)
                        if dep in APP or (mod in PRIMITIVE and dep in CORE):
                            upward[mod].append((rel, n, line.strip()))
                        elif mod in RANK and dep in RANK and RANK[dep] > RANK[mod]:
                            order[mod].append((rel, n, line.strip()))
                    elif APP_USE.search(line):
                        app_use[mod].append((rel, n, line.strip()))
    print(f"{'module':<10} {'app includes':>12} {'app uses':>9} {'field reads':>11} {'order':>6}  files")
    bad = False
    for mod in sorted(PRIMITIVE) + CORE:
        files = {r for r, _, _ in upward[mod] + app_use[mod] + order[mod]}
        print(f"{mod:<10} {len(upward[mod]):>12} {len(app_use[mod]):>9} {reads[mod]:>11} {len(order[mod]):>6}  {len(files)}")
        if mod in strict and files:
            bad = True
    total = sum(len(v) for v in upward.values()), sum(len(v) for v in app_use.values()), sum(reads.values()), sum(len(v) for v in order.values())
    print(f"{'total':<10} {total[0]:>12} {total[1]:>9} {total[2]:>11} {total[3]:>6}")
    if show:
        for title, table in (('application includes', upward), ('application reads', app_use), ('against the module order', order)):
            print(f"\n== {title}")
            for mod in sorted(table):
                for r, n, t in table[mod]:
                    print(f"{r}:{n}: {t[:110]}")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
