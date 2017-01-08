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
- [x] Report the frequency of occurrence in the protein set for each GO annotation/pathway


##Usage
  - run<br>
      ``` python uniprot_acc_id_finder.py --refseq 'refseq_id'[ --ensembl 'ensembl_id'[ --uniprot--id 'uniprot_id']] ```
      
  - example command<br>
      ``` python uniprot_acc_id_finder.py --refseq 'NP_000608.1','NP_001139633.1','NP_037305.2' ```<br>
      or <br>
      ``` python uniprot_acc_id_finder.py --ensembl 'ENSMUSG00000029682' --refseq 'NP_001166492.1' --uniprot_id 'HYALP_MACFA' ```
      or <br>
      ``` python uniprot_acc_id_finder.py --refseq 'WP_005082954.1','NP_416893.1','WP_003898649.1' --uniprot_id 'MNTH_BACSU' ```
       or <br>
      ``` python uniprot_acc_id_finder.py --ensembl 'ENSMUSG00000023030' ```
        
##Requirements
  - Internet connection


##Team Members 

  * [Gizem Abalı] (https://github.com/gizemabali)
  * [Ümran İşler] (https://github.com/UmranIsler)
  * [Nahide Ergün] (https://github.com/nahidErgun)

