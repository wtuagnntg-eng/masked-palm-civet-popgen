#!/usr/bin/env python3
import argparse, re
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--root', required=True)
p.add_argument('--output', required=True)
a=p.parse_args()
rows=[]
for log in sorted(Path(a.root).glob('K*/rep*/admixture.log')):
    txt=log.read_text(errors='ignore')
    m=re.search(r'CV error \(K=(\d+)\):\s*([0-9.eE+-]+)', txt)
    if m:
        rows.append((int(m.group(1)), log.parent.name, float(m.group(2)), str(log.parent)))
with open(a.output,'w') as f:
    f.write('K\treplicate\tCV_error\trun_dir\n')
    for r in sorted(rows): f.write(f'{r[0]}\t{r[1]}\t{r[2]:.10g}\t{r[3]}\n')
