#!/usr/bin/env python3
import argparse,csv
p=argparse.ArgumentParser(); p.add_argument('--het',required=True); p.add_argument('--metadata',required=True); p.add_argument('--output',required=True); a=p.parse_args()
meta={}
with open(a.metadata) as f:
    r=csv.DictReader(f,delimiter='\t')
    for x in r: meta[x['sample_id']]=x['ecotype']
with open(a.het) as f, open(a.output,'w') as o:
    h=f.readline().split(); idx={k:i for i,k in enumerate(h)}
    o.write('sample_id\tecotype\tHo\n')
    for line in f:
        x=line.split(); sid=x[idx['IID']]; ohom=float(x[idx['O(HOM)']]); nnm=float(x[idx['N(NM)']]); ho=1-ohom/nnm if nnm else float('nan')
        o.write(f'{sid}\t{meta.get(sid,"NA")}\t{ho:.10g}\n')
