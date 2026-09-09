#!/usr/bin/env python3
import argparse, gzip

IUPAC={
 frozenset(('A','G')):'R', frozenset(('C','T')):'Y', frozenset(('G','C')):'S',
 frozenset(('A','T')):'W', frozenset(('G','T')):'K', frozenset(('A','C')):'M'
}

def op(path):
    return gzip.open(path,'rt') if path.endswith('.gz') else open(path)

p=argparse.ArgumentParser(description='Convert a biallelic SNP VCF to a variable-site IUPAC FASTA alignment.')
p.add_argument('-i','--vcf',required=True)
p.add_argument('-o','--fasta',required=True)
a=p.parse_args()
seqs=None; samples=[]
with op(a.vcf) as f:
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples=line.rstrip().split('\t')[9:]
            seqs=[[] for _ in samples]
            continue
        x=line.rstrip().split('\t')
        ref=x[3].upper(); alts=x[4].upper().split(',')
        if len(ref)!=1 or len(alts)!=1 or len(alts[0])!=1: continue
        alleles=[ref,alts[0]]
        fmt=x[8].split(':');
        if 'GT' not in fmt: continue
        gi=fmt.index('GT')
        site=[]
        for s in x[9:]:
            fields=s.split(':'); gt=fields[gi] if gi < len(fields) else '.'
            if gt in ('.','./.','.|.'):
                site.append('N'); continue
            toks=gt.replace('|','/').split('/')
            try: bases=[alleles[int(z)] for z in toks if z!='.']
            except (ValueError,IndexError): bases=[]
            if len(bases)!=2: site.append('N')
            elif bases[0]==bases[1]: site.append(bases[0])
            else: site.append(IUPAC.get(frozenset(bases),'N'))
        # Keep variable SNP sites. A VCF is generally variant-only, but this check is explicit.
        observed=set(c for c in site if c in 'ACGT')
        hetero=any(c in 'RYSWKM' for c in site)
        if len(observed)<2 and not hetero: continue
        for i,c in enumerate(site): seqs[i].append(c)
with open(a.fasta,'w') as out:
    for s,seq in zip(samples,seqs):
        out.write(f'>{s}\n')
        st=''.join(seq)
        for i in range(0,len(st),80): out.write(st[i:i+80]+'\n')
