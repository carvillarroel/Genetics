
## Primero creamos el ambiente de conda con las herramientas bioinformaticas a usar
```bash
mamba create -n phylogeny -c bioconda blast mafft trimal iqtree seqkit
```
## Para crear una base de datos local, primero buscamos todas las proteinas del orden de virus llamado Nudiviridae usando la API de NCBI
```bash
esearch -db protein -query "Nudiviridae [ORGN]"| efetch -format acc > accession_list.txt
```
## Las descargamos con un EFETCH
```bash
efetch -db protein -input accession_list.txt  -format fasta > nudivirus_proteins.fasta
```

## Creamos la base de datos local
```bash

makeblastdb -in nudivirus_proteins.fasta -dbtype prot
```

## Buscamos nuestra proteina en la base de datos, por defecto nos da un resultado en formato largo de texto:
```bash

blastp -query query.fasta -db nudivirus_proteins.fasta
```
## 


## Para tener un resultado que se pueda trabajar en Excel (tabulado) usamos el outfmt "6"
```bash

blastp -query query.fasta -db nudivirus_proteins.fasta -outfmt 6 
```

## De esta forma nos muestra en cada columna estos datos:

```bash

qseqid qlen qstart qend sseqid slen sstart send evalue bitscore
```

## Podemos elegir mas o menos datos que mostrar segun esta lista de opciones:
         qseqid means Query Seq-id
               qgi means Query GI
              qacc means Query accesion
           qaccver means Query accesion.version
              qlen means Query sequence length
            sseqid means Subject Seq-id
         sallseqid means All subject Seq-id(s), separated by a ';'
               sgi means Subject GI
            sallgi means All subject GIs
              sacc means Subject accession
           saccver means Subject accession.version
           sallacc means All subject accessions
              slen means Subject sequence length
            qstart means Start of alignment in query
              qend means End of alignment in query
            sstart means Start of alignment in subject
              send means End of alignment in subject
              qseq means Aligned part of query sequence
              sseq means Aligned part of subject sequence
            evalue means Expect value
          bitscore means Bit score
             score means Raw score
            length means Alignment length
            pident means Percentage of identical matches
            nident means Number of identical matches
          mismatch means Number of mismatches
          positive means Number of positive-scoring matches
           gapopen means Number of gap openings
              gaps means Total number of gaps
              ppos means Percentage of positive-scoring matches
            frames means Query and subject frames separated by a '/'
            qframe means Query frame
            sframe means Subject frame
              btop means Blast traceback operations (BTOP)
            staxid means Subject Taxonomy ID
          ssciname means Subject Scientific Name
          scomname means Subject Common Name
        sblastname means Subject Blast Name
         sskingdom means Subject Super Kingdom
           staxids means unique Subject Taxonomy ID(s), separated by a ';'
                         (in numerical order)
         sscinames means unique Subject Scientific Name(s), separated by a ';'
         scomnames means unique Subject Common Name(s), separated by a ';'
        sblastnames means unique Subject Blast Name(s), separated by a ';'
                         (in alphabetical order)
        sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
                         (in alphabetical order)
            stitle means Subject Title
        salltitles means All Subject Title(s), separated by a '<>'
           sstrand means Subject Strand
             qcovs means Query Coverage Per Subject
           qcovhsp means Query Coverage Per HSP
            qcovus means Query Coverage Per Unique Subject (blastn only)



## Trabajemos con esta lista ahora:
```bash

blastp -query query.fasta -db nudivirus_proteins.fasta -evalue 1e-50 -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send evalue bitscore pident stitle"
```

## Ahora si queremos rescatar los archivos fasta de los HITS de blast para hacer un alineamiento usamos SEQKIT
## Primero copiemos la columna con el "sseqid" en un nuevo archivo de texto que llamaremos dnapol_genes.txt
```bash

blastp -query query.fasta -db lefavirales_proteins.fasta -evalue 1e-50 -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send evalue bitscore pident stitle"| awk '{print $5}' > dnapol_list.txt
uniq dnapol_list.txt > dnapol_list_unique.txt
seqkit faidx -l dnapol_list_unique.txt  lefavirales_proteins.fasta > dnapol_proteins.fasta
```

## A este archivo fasta le agregamos el QUERY y para realizar la filogenia con MAFFT, TRIMAL, y IQTREE
```bash

mafft dnapol_proteins.fasta > dnapol_al.fasta
trimal -phylip -in dnapol_al.fasta -out dnapol_trimal.phy
iqtree2 -s dnapol_trimal.phy
```


