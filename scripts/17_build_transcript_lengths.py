#!/usr/bin/env python3
import argparse,re,collections
p=argparse.ArgumentParser(); p.add_argument('--gff',required=True); p.add_argument('--output',required=True); a=p.parse_args()
lengths=collections.Counter()
with open(a.gff) as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        x=line.rstrip().split('\t')
        if len(x)<9 or x[2] != 'CDS': continue
        attrs=x[8]
        # Exclude only CDS features explicitly annotated as partial.
        if re.search(r'(?:^|;)partial=true(?:;|$)', attrs, flags=re.I):
            continue
        m=re.search(r'(?:^|;)Parent=([^;]+)',attrs)
        if not m: continue
        L=int(x[4])-int(x[3])+1
        for tr in m.group(1).split(','):
            lengths[tr]+=L
with open(a.output,'w') as o:
    o.write('transcript_id\tcds_length_bp\n')
    for tr in sorted(lengths):
        if lengths[tr] <= 50: continue
        o.write(f'{tr}\t{lengths[tr]}\n')
