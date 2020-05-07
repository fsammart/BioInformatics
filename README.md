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
```
source env/bin/activate
python3 main.py
```

### Exercise 2

#### Requisites:

- Install Blast Locally at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

#### Steps:

1. Download db at ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz 
2. Unzip it in the /archives/swissprot directory
3. Convert it to blast format using the makeblastdb command

##### For MACOS users
>Command to convert the DB in blast format:
> ```/usr/local/ncbi/blast/bin/makeblastdb -in swissprot -dbtype prot```

### Exercise 3

#### Requisites

- Install ClustalW from here:  http://www.clustal.org/omega/