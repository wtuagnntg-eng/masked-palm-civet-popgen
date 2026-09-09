#!/usr/bin/env python3
import argparse,csv,collections
p=argparse.ArgumentParser(); p.add_argument('--hom',required=True); p.add_argument('--metadata',required=True); p.add_argument('--fai',required=True); p.add_argument('--rename',required=True); p.add_argument('--output',required=True); a=p.parse_args()
meta={}
with open(a.metadata) as f:
    for r in csv.DictReader(f,delimiter='\t'): meta[r['sample_id']]=r['ecotype']
old=set()
with open(a.rename) as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        old.add(line.split('\t')[0])
auto_len=0
with open(a.fai) as f:
    for line in f:
        x=line.split('\t')
        if x[0] in old: auto_len += int(x[1])
D=collections.defaultdict(lambda: {'n':0,'bp':0,'c1':0,'c2':0,'c3':0,'c4':0,'bp1':0,'bp2':0,'bp3':0,'bp4':0})
with open(a.hom) as f:
    h=f.readline().split(); idx={k:i for i,k in enumerate(h)}
    for line in f:
        x=line.split(); sid=x[idx['IID']]
        kb=float(x[idx['KB']]); bp=int(round(kb*1000)); d=D[sid]; d['n']+=1; d['bp']+=bp
        if kb < 1000: c=1
        elif kb < 5000: c=2
        elif kb <= 10000: c=3
        else: c=4
        d[f'c{c}']+=1; d[f'bp{c}']+=bp
with open(a.output,'w') as o:
    o.write('sample_id\tecotype\tN_ROH\tROH_Mb\tFROH\tN_0.1_1Mb\tMb_0.1_1Mb\tN_1_5Mb\tMb_1_5Mb\tN_5_10Mb\tMb_5_10Mb\tN_gt10Mb\tMb_gt10Mb\n')
    for sid in sorted(meta):
        d=D[sid]; vals=[sid,meta[sid],d['n'],d['bp']/1e6,d['bp']/auto_len if auto_len else float('nan'),d['c1'],d['bp1']/1e6,d['c2'],d['bp2']/1e6,d['c3'],d['bp3']/1e6,d['c4'],d['bp4']/1e6]
        o.write('\t'.join(map(str,vals))+'\n')
