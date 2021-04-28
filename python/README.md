# assayM

## Requirements
    pip install biopython pandas openpyxl xlrd primer3-py

## Files  
* sequences.fasta - fasta sequences to be screened for deltaG
* assays.xlx - primer xlsx file

### deltaG command outputs a tsv a file with delta G values for each primer,sequence pair
    python assaym.py deltag -s sequences.fasta  -a cdc_assays.xlsx > deltag.tsv

### disvariants command filters the tsv file to display only when deltaG exceeds 20% compared to reference
    python assaym.py disvariants -t deltag.tsv

