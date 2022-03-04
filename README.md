# CASE STUDY: THE ORIGIN OF HUMAN MALARIA
- BINP29 Project Spring 2022
- Author: Elena Aramendía

## Software Used
- GeneMark-ES Suite v4.62_lic
- BLASTP v2.12.0+
- proteinortho v6.0.33
- BUSCO v5.3.0
- Clustal Omega - 1.2.4
- RAxML v8.2.12
- Phylip consense v3.697.0
- Figtree v1.4.4

## Procedure
In this project we analyze different malaria causing parasites from the phylum *Apicomplexa*.
We analyze 8 species:
- 6 *Plasmodium*: *P. berghei, P.vivax, P. yoelii, P.knowlesi, P.cynomolgi, P.falciparum*
- *Haemoproteus tartakovsky*
- *Toxoplasma gondii* as outgroup

We will follow the following steps to analyze the evolutionary relationships between these organisms.

1. *Plasmodium* Genome -> gene prediction for the provided plasmodium genomes
2. *Haemoproteus* Genome -> remove contamination
3. Gene prediction for *Haemoproteus*
4. Identify orthologs
  - Proteinortho
  - BUSCO
6. Alignments
7. Phylogenetic Trees

## *Plasmodium* data
First I performed a gene prediction step for the *Plasmodium* files, using
GeneMark with the gmes_petap.pl script.
This step was performed in the server.
I performed the gene prediction for P. berghei, so I had to use a
--min-contig of 5000

```
mkdir 1_genepred
cd 1_genepred
gmes_petap.pl
nohup gmes_petap.pl --ES --sequence Plasmodium_berghei.genome --min_contig 5000 &
```

The output file was renamed as Plasmodium_berghei.gtf and copied to
/tmp/Prediction on the server. The rest of the prediction files were obtained
from the server.

## *Haemoproteus tartakovskyi* data

Since the input reads derive both from the bird and the parasite,
bird scaffolds should be removed. To remove these scaffolds we filter them
according to GC content. Since we know that the mean GC-content for the parasites
is aproximately 26%, I decided to filter for sequences with a maximum of 30%
content.

### Filter for GC-content and length
I filtered for GC content and length, keeping scaffolds with a maximum of 30% GC-content
and a minimum length of 3000 nts.
For this I used a python script -> clean_scaffold.py
```
mkdir 1.1_cleanHt
python clean_scaffold.py -i Haemoproteus_tartakovskyi.genome -len 3000 -gc 30 -o 1.1_cleanHt/H_tartakovsky_clean30gc.fa
```

In this step 547 scaffolds removed, we keep 1696 scaffolds.

### Gene prediction
Make a gene prediction for H. tartakovsky with this new filtered data.
```
cd 1_genepred
gmes_petap.pl --ES --sequence ../H_tartakovsky_clean30gc.fa --min_contig 5000
```

## Gff to fasta
To create fasta sequences from the gff file and the genome file I used the
gffParse.pl script.

```
mkdir 2_gtf_to_fasta
cd 2_gtf_to_fasta
gffParse.pl -i ../H_tartakovsky_clean30gc.fa -g ../1_genepred/Ht.gff -c -p -b Ht
```

### Remove scaffolds that have genes that are from avian origin
I used blastp.

```
mkdir 3_blast
cd 3_blast
blastp -query ../2_gtf_to_fasta/Ht.faa -db SwissProt -evalue 1e-10 -out Ht_Swissprot.blastp -num_threads 10
```

I used the provided script datParse.py to get the list of scaffolds that contains genes from birds.

```
python datParser.py Ht2_Swissprot.blastp ../2_gtf_to_fasta/Ht.fna taxonomy.dat uniprot_sprot.dat > scaffolds.txt
```

Remove scaffolds from fasta file -> remove_scaffolds.py

16 scaffolds were removed. We perform a new gene prediction with these.

```
mkdir 4_newgenepred
cd 4_newgenepred
gmes_petap.pl --ES --sequence ../3_blast/only_ht.fna --min_contig 5000 --cores 5
```

Convert to fasta files:

```
gffParse.pl -i ../3_blast/only_ht.fna -g only_ht.gtf -c -p -b only_ht
```
Get some data about the files:
- Genome size:
```
for file in *.genome; do echo $file; cat $file | grep -v '>' | tr -d '\n\s' | wc -c; done
```
- Genes:
```
for file in *.gtf; do echo $file; cat $file | grep 'CDS' | cut -f 9 | cut -d ';' -f 1 | sort -u | wc -l; done
```
- Genomic GC -> python script genomic_gc.py
```
for file in *.genome; do python genomic_gc.py $file ; done
```

| Species | Host | Genome size | Genes |  Genomic GC   |
| ------- | ---- | ----------- | ----- | ------------- |
| Plasmodium berghei | rodents | 17,954,629 | 7282 | 23.72% |
| Plasmodium cynomolgi | macaques | 26,181,343 | 5787 | 40.38% |
| Plasmodium falciparum | humans | 23,270,305 | 5207 | 19.36% |
| Plasmodium knowlesi | lemures | 23,462,332 | 4953 | 38.83% |
| Plasmodium vivax | humans | 27,007,701 | 5682 | 42.28% |
| Plasmodium yoelii | rodents | 22,222,369 | 4919 | 21.77% |
| Haemoproteus tartakovskyi | birds | 20,431,038 | 4108 | 28.56%|
| Toxoplasma gondii | humans | 128,105,889 | 4193 | 52.35% |

# Phylogenetic Trees

## Identify orthologs
Print the protein sequences in fasta format for each genome with the aid of the
gffParse.pl program downloaded earlier.

```
for file in ../1_genepred/*.gtf; do id=$(echo $file | cut -d '/' -f 3 | cut -d '.' -f 1);
gffParse.pl -i ../$id.genome -g ../1_genepred/$id.gtf -b $id -c -p ;done
```
###  Proteinortho
Install with conda:
```
mkdir 5_proteinortho
cd 5_proteinortho
conda install proteinortho
```

Run proteinortho:

```
nohup proteinortho6.pl ../2_gtf_to_fasta/{Plasmodium_faciparum,Plasmodium_cynomolgi,Plasmodium_berghei,Plasmodium_knowlesi,Plasmodium_vivax,Plasmodium_yoelii,Tg,only_ht}.faa &

#change names:
for file in *; do mv $file only_ht.${file#myproject.}; done
```


### BUSCO

I first installed busco and made directories for each of the species.

```
mkdir Ht
mkdir Pb
mkdir Pf
mkdir Pk
mkdir Pv
mkdir Py
mkdir Tg
conda install busco
```
Since the names of my files where not very consistent and I had few files, I ran busco individually for each file in the corresponding directory and gave the outputs shorter names.
Run for the Pb files:
```
busco -i ../2_gtf_to_fasta/Plasmodium_berghei.faa -o Pb -m prot -l apicomplexa
```

I checked the full_table_??.tsv files for each species, most of them had between 300-400
single-copy complete BUSCOs, except Toxoplasma, which had only 4 singlecopy complete BUSCOs.
Because of this I decided to work with the complete single BUSCOs from the other species
and, in Toxoplasma, with the single-copy ones and one copy of the duplicated BUSCOs.

To create one protein fasta file per BUSCO, with all eight species present with one sequence,
I performed the following commands:

```
# First I get the busco ids and in how many of the species they appear (counts)
# At first I discarded Tg so I add those genes in the next step

for file in *.tsv; do grep -v "^#" $file | awk '$2=="Complete" {print $1}' >> complete_busco_ids.txt; done
for file in *.tsv; do grep 'Complete' $file | cut -f1 >> complete_busco_ids.txt; done


 # Get Tg ids
 grep 'Complete' full_table_Tg.tsv | cut -f1 > Tg_busco_ids.txt
 grep 'Duplicated' full_table_Tg.tsv | cut -f1 >> Tg_busco_ids.txt
 sort -u Tg_busco_ids.txt > Tg_uniq_busco.txt
 for line in Tg_uniq_busco.txt; do grep -m 1 $line | cut -f1,3;done

 sort complete_dpTg_busco_ids.txt |uniq -c > Tgcomplete_busco_ids_with_counts.txt

 cat Tg_uniq_busco.txt >> complete_busco_ids.txt


 ## Get buscos in all 8 species
 awk '$1 > 7 {print $2}' Tgcomplete_busco_ids_with_counts.txt > Tgcommon_buscos.txt

 ## Get file with the gene ids for each species
 for sp in Ht Pb Pc Pf Pk Pv Py Tg; do cat Tgcommon_buscos.txt | while read line; do grep -m 1 $line full_table_${sp}.tsv | cut -f1,3 >> ${sp}_common_gene_ids.txt ; done; done

 ## Get fasta sequences for each species
 # I couldn't use a loop because of the file names; this is an example command for the Ht files
 # I ran the same command for the rest of the species

 cat Ht_common_gene_ids.txt | while read line; do busco=$(echo $line | cut -d ' ' -f1); gene=$(echo $line | cut -d ' ' -f2); grep -m 1 -A 1 ">$gene" ../2_gtf_to_fasta/only_ht.faa | sed "s/>${gene}\slength=/>Ht\t${gene}\t${busco}\tlength=/" >> Ht_common_buscos.faa ; done

 ##Get BUSCO file with one sequence for each species
 mkdir ../7_busco_to_fasta # new directory for these files
 for sp in Tg Ht Pb Pc Pf Pk Pv Py; do cat Tgcommon_buscos.txt | while read line; do grep -A 1 [[:space:]]$line ${sp}_common_buscos.faa >> ../7_busco_to_fasta/${line}_all.faa; done; done
```  
We get a final file with **178 genes**.


 # ALIGNMENTS AND TREES

First we perform the alignments using Clustal Omega.
```
mkdir 8.1_alignments
 cat ../../6_busco/common_buscos.txt | while read line; do clustalo -i ../../7_busco_to_fasta/${line}_all.faa -o ../../8.1_alignments/${line}_aligned.faa -v; done
```
Make a tree from each of the files, using Tg as an outgroup:
```
 cat ../../6_busco/common_buscos.txt | while read line; do raxmlHPC -s ../../8.1_alignments/${line}_aligned.faa -n ${line}.tre -o Tg -m PROTGAMMABLOSUM62 -p 12345; done
```

To merge all the generated trees, I used the software phylip consense to get a consensus tree.
First get all trees in one file, this file is input for phylip consense.
```
 cat RAxML_bestTree.*.tre > all_bestTrees.tre
 phylip consense all_bestTrees.tre
```

To plot it I used figtree:
```
 figtree -graphic PDF outtree Pf_outgroup.pdf
```
Next I concatenated all the alignment sequences and generated a sequences from this "superalignment",
to compare it to this consensus tree.
For this I wrote the python script concatenate_alignments.py.
```
python concatenate_alignments.py -i ../6_busco/common_buscos.txt -d . > concatenated_alignments.faa
```

Make the tree:
```
raxmlHPC -s concatenated_alignments -n concPk_tree.tre -o Pk -m PROTGAMMABLOSUM62 -p 12345
```
To plot it I used figtree:
```
figtree -graphic PDF RAxML_parsimonyTree.concPk_tree.tre concPk_outgroup.pdf
```
