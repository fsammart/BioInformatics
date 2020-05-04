# BioInformatics

Project for the BioInformatics course at ITBA.

## Installation 
```
python3 -m pip install --user virtualenv
python3 -m virtualenv env
source env/bin/activate
pip install -r requirements.txt
```
## Run
### Exercise1 
```
source env/bin/activate
python3 main.py
```

### Excercise 2

#### Requisites:

-Install Blast Locally at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

#### Steps:

1. Download db at ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz 
2. Unzip it in the /excercises directory
3. Convert it to blast format using the makeblastdb command


##### For MACOS users
>Command to convert the DB in blast format:
> ```/usr/local/ncbi/blast/bin/makeblastdb -in swissprot -dbtype prot```

