#!/usr/bin/env python
import argparse
import sys
from Bio import SeqUtils
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio import pairwise2
import primer3
import pandas
import json
import time


def sampleseq(args):
    f=open(args.json)
    lineage_acc={}
    lineage_seq={}
    lineage_subm_date={}
    for line in f:
        line_dict=json.loads(line)
        
        (acc,lineage,subm_date,seqlen)=(line_dict['covv_accession_id'],line_dict['covv_lineage'],line_dict['covv_subm_date'],line_dict['sequence_length'])
        subdate=lineage_subm_date.get(lineage,subm_date)
        if(seqlen>29000 and subm_date>=subdate):
            lineage_subm_date[lineage]=subm_date
            lineage_acc[lineage]=acc
            lineage_seq[lineage]=line_dict['sequence']
        
    lin_fasta=open("lineages.fasta",'w')
    lin_tsv=open("lineages.tsv",'w')
    for lineage,accession in lineage_acc.items():
        lin_fasta.write(f">{accession}\n")
        lin_fasta.write(lineage_seq[lineage].replace('\n',''))
        lin_fasta.write("\n")
        lin_tsv.write(f"{lineage}\t{accession}\t{lineage_subm_date[lineage]}\n")
    

def deltag(args):
    assay_data_df = pandas.read_excel(args.assays, sheet_name='Sheet1', usecols=['gene target', 'assay_name','country','type','sequence','position'])
    for fastarecord in SeqIO.parse(args.sequences, "fasta"):
        for index, row in assay_data_df.iterrows(): 
            if(row['type']=='F' or row['type']=='R'):   
                if(row['type']=='F'):
                    (sstart,send)=row['position'].strip().split('->')
                    query=row['sequence'].strip()
                    
                if(row['type']=='R'):
                    (send,sstart)=row['position'].strip().split('<-')
                    query=str(Seq(row['sequence'].strip()).reverse_complement())
                    
                search_start=int(sstart)-1001 if (int(sstart)-1001) >0 else 0
                search_end=int(send)+1000 if (int(send)+1000)<len(str(fastarecord.seq)) else len(str(fastarecord.seq))-1
                subject=str(fastarecord.seq)[search_start:search_end]
                #tic = time.perf_counter()
                mstart,mend=match(query,subject)
                #print(fastarecord.name,row['assay_name'],row['type'],query,subject[mstart-1:mend])
                if(mstart is not None):
                    #toc = time.perf_counter()
                    #print(f"matched in {toc - tic:0.4f} seconds")
                    
                    pbs = str(Seq(subject[mstart-21:mend+20]).reverse_complement())        
                    if(row['type']=='R'):
                        pbs = subject[mstart-21:mend+20]

                    thermal = primer3.bindings.calcEndStability(query,pbs)
                    print(fastarecord.name,row['assay_name'],row['type'],thermal.dg,sep="\t")
                else:
                    print(fastarecord.name,row['assay_name'],row['type'],'',sep="\t")                

def disvariants(args):
    print(args.reference,deltag.tsv)
    all=pandas.read_csv(deltag.tsv,sep="\t",header=None)
    ref=all[all[0]==args.reference]
    m=all.merge(ref,how='left',on=[1,2])
    m['dg_increase_pct']=m.apply(lambda row:-100*(row['3_x']-row['3_y'])/row['3_y'],axis=1)
    m.columns=['Sequence','Primer','Direction','dg_sequence','Ref','dg_ref','dg_increase_pct']
    result=m.loc[m['dg_increase_pct']>20,["Sequence","Primer","Direction","dg_sequence","dg_increase_pct"]]
    result.to_csv(sys.stdout,index=False,sep="\t")
    
def summary(args):
    print(args.reference,args.deltagtsv,args.lineagestsv)
    deltag=pandas.read_csv(args.deltagtsv,sep="\t",header=None)
    lineages=pandas.read_csv(args.lineagestsv,sep="\t",header=None)
    lineages.columns=['lineage','Sequence','submission_date']
    ref=deltag[deltag[0]==args.reference]
    m=deltag.merge(ref,how='left',on=[1,2])
    m['dg_increase_pct']=m.apply(lambda row:-100*(row['3_x']-row['3_y'])/row['3_y'],axis=1)
    m.columns=['Sequence','Primer','Direction','dg_sequence','Ref','dg_ref','dg_increase_pct']
    m2=m.merge(lineages,how='left',on=['Sequence'])
    m3=m2.loc[(m2['dg_increase_pct']>0) & (m2['dg_sequence']>0),['Primer','Direction','lineage']].groupby(['Primer','Direction']).agg(['count']).reset_index()
    m3['lineage_total']=lineages.lineage.count()
    m3.columns=['Primer','Direction','mismatch_variants','total_variants']
    m3['mismatch_pct']=round(m3['mismatch_variants']*100/m3['total_variants'],2)
    m3.to_csv('summary.csv',index=False,sep="\t")
    #m2[(m2['dg_increase_pct']>0) & (m2['dg_sequence']>0)]['Primer','Direction','lineage','Sequence','dg_sequence'].to_csv('details.csv',index=False,sep="\t")
    m2.loc[(m2['dg_increase_pct']>0) & (m2['dg_sequence']>0),['Primer','Direction','lineage','Sequence','dg_sequence']].to_csv('details.csv',index=False,sep="\t")

    
def match(query,subject):
    m=SeqUtils.nt_search(subject,query)
    matchstart=0
    mcount=len(m)-1
    mscore=0
    # if straightforward string match works
    if(mcount>0):
        matchstart=m[1]+1
        matchend=m[1]+len(query)
    # otherwise do local alignment with mismatches and gaps
    else:
        aln = pairwise2.align.localms(query,subject, 2, -1, -0.5, -.1,one_alignment_only=True)
        if(len(aln)>0):
            matchstart=aln[0].start+1
            matchend=aln[0].end
            mscore=aln[0].score
        else:
            return (None,None)
    if(mcount>0 or mscore>25):
        return(matchstart,matchend)
    else:
        return (None,None)
    
if __name__=='__main__':
    ap=argparse.ArgumentParser("assaym")
    sp = ap.add_subparsers(help='assaym commands')

    sampleseq_parser=sp.add_parser('sampleseq',help='Sample Sequences by lineage')
    sampleseq_parser.add_argument('-j','--json',help='gisaid json file',required=True)
    sampleseq_parser.set_defaults(func=sampleseq)
    
    
    getlineagefasta_parser=sp.add_parser('lineagefasta',help='extract fasta records of the give lineages from GISAID')
    deltag_parser=sp.add_parser('deltag',help='computer deltaG for primers against variants')
    deltag_parser.add_argument('-s','--sequences',help='fasta file name',required=True)
    deltag_parser.add_argument('-a','--assays',help='primers xlsx file',required=True)
    deltag_parser.set_defaults(func=deltag)
    
    disvariants_parser=sp.add_parser('disvariants',help='Filter variants that have higher deltaG compared to reference')
    disvariants_parser.add_argument('-r','--reference',help='Reference sequence name',default='NC_045512')
    disvariants_parser.add_argument('-d','--deltagtsv',help='deltaG tsv file name',required=True)
    disvariants_parser.set_defaults(func=disvariants)
    
    summary_parser=sp.add_parser('summary',help='Filter variants that have higher deltaG compared to reference')
    summary_parser.add_argument('-r','--reference',help='Reference sequence name',default='NC_045512')
    summary_parser.add_argument('-d','--deltagtsv',help='deltaG tsv file name',required=True)
    summary_parser.add_argument('-l','--lineagestsv',help='lineages tsv file name',default='lineages.tsv')
    summary_parser.set_defaults(func=summary)
    
    ns=ap.parse_args(args=None if sys.argv[1:] else ['--help'])
    ns.func(ns)
