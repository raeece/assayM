# assayM

## Requirements
    pip install biopython pandas openpyxl xlrd primer3-py plotly dash dash-bootstrap-components snakemake

## Files  
* sequences.fasta - fasta sequences to be screened for deltaG
* assays.xlx - primer xlsx file

### deltaG command outputs a tsv a file with delta G values for each primer,sequence pair
    python assaym.py deltag -s sequences.fasta  -a cdc_assays.xlsx > deltag.tsv

### disvariants command filters the tsv file to display only when deltaG exceeds 20% compared to reference
    python assaym.py disvariants -d deltag.tsv

## Snakemake file runs the complete pipeline downloading all sequences from gisaid
* downloads all sequences from gisaid
* samples representative sequences for each lineage
* run deltag on all assays against representative sequences
* output results as summary.csv and details.csv


## Run the following commands to run the complete pipeline takes ~1 hour on a laptop 
    export GISAIDUSER=gisaidusername
    export GISAIDPASS=gisaidpass
    snakemake --delete-all-output
    snakemake --cores 1

