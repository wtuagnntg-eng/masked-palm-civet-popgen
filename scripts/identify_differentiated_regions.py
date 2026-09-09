#!/usr/bin/env python3
import argparse,csv,math

def read_pi(path):
    d={}
    with open(path) as f:
        r=csv.DictReader(f,delimiter='\t')
        for x in r:
            key=(x['CHROM'],int(x['BIN_START']),int(x['BIN_END']))
            d[key]=(int(x['N_VARIANTS']),float(x['PI']))
    return d

def read_fst(path):
    d={}
    with open(path) as f:
        r=csv.DictReader(f,delimiter='\t')
        for x in r:
            # VCFtools columns: CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST MEAN_FST
            key=(x['CHROM'],int(x['BIN_START']),int(x['BIN_END']))
            try: val=float(x['WEIGHTED_FST'])
            except: continue
            d[key]=(int(x['N_VARIANTS']),val)
    return d

def q(v,p):
    v=sorted(v); n=len(v)
    if not n: return float('nan')
    z=(n-1)*p; lo=int(math.floor(z)); hi=int(math.ceil(z))
    return v[lo] if lo==hi else v[lo]+(v[hi]-v[lo])*(z-lo)

p=argparse.ArgumentParser()
for k in ['wild_pi','farm_pi','fst','windows_out','candidate_windows_out','regions_out','thresholds_out']:
    p.add_argument('--'+k.replace('_','-'),dest=k,required=True)
p.add_argument('--min-snps',type=int,default=10); p.add_argument('--quantile',type=float,default=.95)
a=p.parse_args()
w=read_pi(a.wild_pi); f=read_pi(a.farm_pi); s=read_fst(a.fst)
rows=[]
for key in sorted(set(w)&set(f)&set(s), key=lambda z:(z[0],z[1])):
    nw,pw=w[key]; nf,pf=f[key]; ns,fs=s[key]
    if min(nw,nf,ns)<a.min_snps or pf<=0 or not all(map(math.isfinite,[pw,pf,fs])): continue
    rows.append((*key,nw,nf,ns,pw,pf,fs,pw/pf))
fst_thr=q([r[-2] for r in rows],a.quantile); ratio_thr=q([r[-1] for r in rows],a.quantile)
with open(a.windows_out,'w') as o:
    o.write('CHROM\tSTART\tEND\tN_WILD_PI\tN_FARM_PI\tN_FST\tPI_WILD\tPI_FARM\tWEIGHTED_FST\tPI_WILD_DIV_FARM\n')
    for r in rows:o.write('\t'.join(map(str,r))+'\n')
cands=[r for r in rows if r[-2]>=fst_thr and r[-1]>=ratio_thr]
with open(a.candidate_windows_out,'w') as o:
    for r in cands: o.write(f'{r[0]}\t{r[1]-1}\t{r[2]}\t{r[-2]}\t{r[-1]}\n')
# Merge physically overlapping sliding windows.
regions=[]
for r in cands:
    chrom,start,end=r[0],r[1]-1,r[2]
    if regions and regions[-1][0]==chrom and start < regions[-1][2]:
        regions[-1][2]=max(regions[-1][2],end)
    else: regions.append([chrom,start,end])
with open(a.regions_out,'w') as o:
    for r in regions:o.write('\t'.join(map(str,r))+'\n')
with open(a.thresholds_out,'w') as o:
    o.write('metric\tthreshold\n')
    o.write(f'weighted_FST_q{a.quantile}\t{fst_thr}\n')
    o.write(f'piWild_piFarm_q{a.quantile}\t{ratio_thr}\n')
    o.write(f'n_matched_windows\t{len(rows)}\n')
    o.write(f'n_candidate_windows\t{len(cands)}\n')
    o.write(f'n_merged_regions\t{len(regions)}\n')
