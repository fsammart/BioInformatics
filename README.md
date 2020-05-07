# BioInformatics

Project for the BioInformatics course at ITBA.

## Installation 
```
sudo apt install python3
sudo apt install python3-pip
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
### Exercise 1
Place your gb files in the archives/gb_files folder

### Exercise 2
Place your fasta files in the archives/prot_sequences folder

#### Requisites:
- Install Blast Locally at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- For Linux run
```
sudo apt-get install ncbi-blast+
```
#### Steps:
1. Download db at ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz 
2. Unzip it in the archives/swissprot directory
3. Convert it to blast format using the makeblastdb command

##### For MACOS users
>Command to convert the DB in blast format:
> ```/usr/local/ncbi/blast/bin/makeblastdb -in swissprot -dbtype prot```

### Exercise 3
Place your msa.fasta files in the archives/msa folder

#### Requisites
- Install ClustalW from here:  http://www.clustal.org/omega/
- For Linux run
```
sudo apt-get install -y clustalw 
```