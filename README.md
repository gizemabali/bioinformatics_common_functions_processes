# Finding The Common Functions, Processes For Given Proteins

## Project Details

Find/list Gene Ontology (GO) annotations (molecular function, biological process and cellular component) and pathways for given proteins

Database : Uniprot<br>
Tools : Python 3.x<br>

##STEPS
- [x] Connect to uniprot.
- [x] Get uniprot access ids of proteins which are given with ensembl and refseq ids
- [x] Use refseq ids and ensembl ids together in the command line
- [x] Get go annotations and pathways of proteins
- [x] Report GO annotations/pathways common to all protein set members.
- [ ] Report the frequency of occurrence in the protein set for each GO annotation/pathway


##Usage
  - run
      ``` python uniprot_acc_id_finder.py --refseq 'refseq_id' --ensembl 'ensembl_id' ```
      
  - example command<br>
      ``` python uniprot_acc_id_finder.py --refseq 'NP_000326.2','NP_001159374.1' --ensembl 'ENSG00000167110' ```<br>
      or <br>
      ``` python uniprot_acc_id_finder.py ensembl 'ENSG00000229215','ENSG00000235657' ```

##Requirements
  - Internet connection


##Team Members 

  * [Gizem Abalı] (https://github.com/gizemabali)
  * [Ümran İşler] (https://github.com/UmranIsler)
  * [Nahide Ergün] (https://github.com/nahidErgun)

