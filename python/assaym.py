import argparse
import sys
from Bio import SeqUtils
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio import pairwise2
import primer3
import pandas

def deltag(args):
    assay_data_df = pandas.read_excel(args.assays, sheet_name='Sheet1', usecols=['gene target', 'assay_name','country','type','sequence','position'])
    for fastarecord in SeqIO.parse(args.sequences, "fasta"):
        for index, row in assay_data_df.iterrows(): 
            if(row['type']=='F' or row['type']=='R'):   
                query=row['sequence'].strip()
                subject=str(fastarecord.seq)
                if(row['type']=='R'):
                    query=str(Seq(query).reverse_complement())
                
                mstart,mend=match(query,subject)
                pbs = str(Seq(subject[mstart-21:mend+20]).reverse_complement())        
                if(row['type']=='R'):
                    pbs = subject[mstart-21:mend+20]

                thermal = primer3.bindings.calcEndStability(query,pbs)
                print(fastarecord.name,row['assay_name'],row['type'],thermal.dg,sep="\t")
                

def disvariants(args):
    print(args.reference,args.tsv)
    all=pandas.read_csv(args.tsv,sep="\t",header=None)
    ref=all[all[0]==args.reference]
    m=all.merge(ref,how='left',on=[1,2])
    m['dg_increase_pct']=m.apply(lambda row:-100*(row['3_x']-row['3_y'])/row['3_y'],axis=1)
    m.columns=['Sequence','Primer','Direction','dg_sequence','Ref','dg_ref','dg_increase_pct']
    result=m.loc[m['dg_increase_pct']>20,["Sequence","Primer","Direction","dg_sequence","dg_increase_pct"]]
    result.to_csv(sys.stdout,index=False,sep="\t")
    
    
        
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
        matchstart=aln[0].start+1
        matchend=aln[0].end
        mscore=aln[0].score
    if(mcount>0 or mscore>25):
        return(matchstart,matchend)
    else:
        return None
    
if __name__=='__main__':
    ap=argparse.ArgumentParser("assaym")
    sp = ap.add_subparsers(help='assaym commands')

    flaglineages_parser=sp.add_parser('flaglineages',help='identify lineages that have mutations in primers')
    getlineagefasta_parser=sp.add_parser('lineagefasta',help='extract fasta records of the give lineages from GISAID')
    deltag_parser=sp.add_parser('deltag',help='computer deltaG for primers against variants')
    deltag_parser.add_argument('-s','--sequences',help='fasta file name',required=True)
    deltag_parser.add_argument('-a','--assays',help='primers xlsx file',required=True)
    deltag_parser.set_defaults(func=deltag)
    
    disvariants_parser=sp.add_parser('disvariants',help='Filter variants that have higher deltaG compared to reference')
    disvariants_parser.add_argument('-r','--reference',help='Reference sequence name',default='NC_045512')
    disvariants_parser.add_argument('-t','--tsv',help='deltaG tsv file name',required=True)
    disvariants_parser.set_defaults(func=disvariants)
    
    ns=ap.parse_args(args=None if sys.argv[1:] else ['--help'])
    ns.func(ns)
