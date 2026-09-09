#!/usr/bin/env python3
import argparse,csv,os
p=argparse.ArgumentParser()
p.add_argument('--counts',required=True); p.add_argument('--lengths',required=True)
p.add_argument('--transcript-to-gene',required=True); p.add_argument('--sample-metadata',required=True)
p.add_argument('--output-transcript',required=True); p.add_argument('--output-gene',required=True)
a=p.parse_args()
lengths={}
with open(a.lengths) as f:
    for r in csv.DictReader(f,delimiter='\t'): lengths[r['transcript_id']]=float(r['cds_length_bp'])
t2g={}
with open(a.transcript_to_gene) as f:
    for r in csv.DictReader(f,delimiter='\t'): t2g[r['transcript_id']]=r['gene_name']
# Map featureCounts BAM-column basenames to sample/tissue labels.
label_by_bam={}
with open(a.sample_metadata) as f:
    for r in csv.DictReader(f,delimiter='\t'):
        sid=r['sample_id']; tissue=r.get('tissue','') or sid
        label_by_bam[sid+'.sorted.bam']=tissue
with open(a.counts) as f:
    lines=[x for x in f if not x.startswith('#')]
reader=csv.reader(lines,delimiter='\t'); hdr=next(reader); raw_cols=hdr[6:]
sample_cols=[]
for col in raw_cols:
    base=os.path.basename(col)
    sample_cols.append(label_by_bam.get(base, os.path.splitext(os.path.splitext(base)[0])[0]))
if len(set(sample_cols)) != len(sample_cols):
    raise SystemExit('Sample/tissue labels are not unique after metadata mapping.')
records=[]
for row in reader:
    tr=row[0]
    if tr not in lengths: continue
    counts=[float(x) for x in row[6:]]
    if sum(counts) <= 10 or sum(counts)/len(counts) <= 1: continue
    records.append((tr,counts))
rpk=[]
for tr,counts in records:
    kb=lengths[tr]/1000.0
    rpk.append((tr,[c/kb for c in counts]))
scales=[sum(vals[j] for _,vals in rpk)/1e6 for j in range(len(sample_cols))]
tpm=[]
for tr,vals in rpk:
    tpm.append((tr,[vals[j]/scales[j] if scales[j]>0 else 0 for j in range(len(sample_cols))]))
with open(a.output_transcript,'w') as o:
    o.write('transcript_id\t'+'\t'.join(sample_cols)+'\n')
    for tr,vals in tpm:o.write(tr+'\t'+'\t'.join(f'{v:.10g}' for v in vals)+'\n')
genes={}
for tr,vals in tpm:
    g=t2g.get(tr)
    if not g: continue
    genes.setdefault(g,[]).append(vals)
with open(a.output_gene,'w') as o:
    o.write('gene\t'+'\t'.join(sample_cols)+'\n')
    for g in sorted(genes):
        arr=genes[g]; means=[sum(x[j] for x in arr)/len(arr) for j in range(len(sample_cols))]
        o.write(g+'\t'+'\t'.join(f'{v:.10g}' for v in means)+'\n')
