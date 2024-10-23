# Clase RNAseq - Alineamiento de reads al genoma de referencia

### Primero creamos el ambiente de conda con las herramientas que usaremos
```bash
mamba create -n mapping -c bioconda trim-galore fastqc star multiqc samtools subread
```
### Exploremos los archivos de secuenciación de illumina (.fastq .fastq.gz .fq.gz etc)
```bash
zless yeast_lowN_R1_1.fq.gz
```

### Visualizamos calidad de secuenciacion con FASTQC
```bash
fastqc -t 8 *.fq.gz
```


### Usamos MULTIQC para juntar todos los informes de FASTQC
```bash
multiqc .
```


### Limpieza de regiones de baja calidad de secuenciación utilizando TRIM_GALORE
```bash

for i in yeast_control_R1 yeast_control_R2 yeast_control_R3 yeast_lowN_R1 yeast_lowN_R2 yeast_lowN_R3;do trim_galore -j 8 -q 20 -length 30 --paired ${i}_1.fq.gz ${i}_2.fq.gz; done

```

### Preparamos el genoma de referencia para que pueda ser usado por STAR
```bash
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir genome_index --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
```


### Realizamos el alineamiento (mapeo) de lecturas al genoma de referencia usando STAR
### Se realiza para cada muestra (pueden usar un FOR LOOP)
```bash

STAR --runThreadN 8 --genomeDir genome_index/ --readFilesIn yeast_control_R1_1_val_1.fq.gz yeast_control_R1_2_val_2.fq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted  --outFileNamePrefix yeast_control_R1  --outTmpDir ~/temp2/

```


### Al tener listo los mapeos (6 en total), contamos las lecturas que caen en cada gen usando featureCounts

```bash

featureCounts -T 8 -p -t gene -F GTF -a Saccharomyces_cerevisiae.R64-1-1.57.gff3 -o counts.tsv *.bam

```

### El archivo counts.tsv contiene el número de lecturas que mapean a cada gen y por cada condición y por cada replica

### Para visualizar los alineamientos, podemos usar el genome browser JBROWSE, pero primero debemos ordernar e indexar nuestros archivos .BAM usando SAMTOOLS
```bash
samtools sort -O bam -o yeast_control_R1Aligned.out_sorted.bam  yeast_control_R1Aligned.out.bam
samtools index yeast_control_R1Aligned.out_sorted.bam
```


