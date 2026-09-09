#!/usr/bin/env python3
import argparse,csv
p=argparse.ArgumentParser(); p.add_argument('--transcripts',required=True); p.add_argument('--mapping',required=True); p.add_argument('--output',required=True); a=p.parse_args()
mp={}
with open(a.mapping) as f:
    for r in csv.DictReader(f,delimiter='\t'):
        mp[r['transcript_id']]=r['gene_name']
missing=[]; genes=[]
with open(a.transcripts) as f:
    for line in f:
        tr=line.strip()
        if not tr: continue
        if tr in mp and mp[tr]: genes.append(mp[tr])
        else: missing.append(tr)
with open(a.output,'w') as o:
    for g in sorted(set(genes)): o.write(g+'\n')
if missing:
    print(f'WARNING: {len(missing)} transcript IDs were absent from transcript_to_gene.tsv', flush=True)
