# Trabajo práctico N ° 2 – Genómica

En este práctico se analizarán secuenciaciones de genoma completo de SARS-CoV-2 obtenidas con el método de amplicones y secuenciadas en un equipo Illumina.  
Se realizarán los análisis en LINUX usando la aplicación disponible para Windows 10 (Windows Linux Subsystem).

### ¿Primera vez con Linux?
 Si es tu primera vez usando la terminal de comandos de Linux, para hacer el análisis debes poner atención a los siguientes puntos:
1)	Siempre escribir el nombre de cada comando y archivo de forma exacta, respetando espacios, mayúsculas y minúsculas.
2)	Escribir comandos en una sola línea. Solo presionar Enter al terminar de escribir el comando
3)	Se puede corregir lo escrito moviendo el cursor con el teclado (no usen el mouse!)
4)	Se puede pegar texto usando el botón derecho del mouse

### PRIMEROS PASOS
En plomo se indica que deben escribir en la terminal
1.-  Al iniciar el práctico, verificar si están en la carpeta que contiene los archivos a utilizar (/mnt/d/genomica)

```bash

pwd
```

2.- Para ver que archivos hay en la carpeta, usar el comando ls
```bash
ls
```

3.- En la carpeta debe estar el genoma de referencia (Wuhan-1.fasta), archivos fastq de secuenciación masiva (extensión “.fastq.gz”),y un archivo con información de partidores (primers.bed).
### ANALISIS DE CALIDAD
4.- Comenzaremos el análisis con la muestra Sample1. Esta muestra se secuenció con lecturas pareadas (paired-end) por lo tanto encontramos dos archivos : sample1_1.fastq.gz y sample2_2.fastq.gz
5.- Haremos una visualización de la calidad de secuenciación con el programa FASTQC
```bash
fastqc sample1_1.fastqc sample2_2.fastqc
```
6.- Este análisis produce dos archivos html (uno para cada fastq), estos se pueden encontrar usando el Explorador de Windows, y abriendo los archivos html en Chrome o Firefox.
### MAPEO DE SECUENCIAS
 7.- Para mapear las secuencias fastq al genoma de referencia utilizaremos el programa minimap2
```bash
minimap2 -a -x  sr -o sample1.sam  Wuhan-1.fasta sample1_1.fastq.gz sample1_1.fastq.gz
```
* -a dice a minimap2 que queremos el output en formato SAM
* -x sr dice a minimap2 que queremos optimizar el programa para usar short-reads (porque nuestras reads vienen de Illumina)

8.- Este programa nos produce un archivo de texto con información de los mapeos del tipo SAM, que debe ser optimizado para el uso posterior por otros programas, por lo que lo transformamos a BAM usando SAMtools
```bash
samtools sort -o sample1_untrimmed.bam sample1.sam
```
* sort le decimos a samtools que queremos usar su función de ordenar el archivo y transformar a bam
* -o le decimos a samtools el nombre del archivo de salida
### CORRECCION A LOS PARTIDORES 
9.- Como las secuencias genómicas se obtuvieron al amplificar los 30 kb de RNA original con partidores multiplex de PCR, se debe corregir que nuestro llamado de variantes no sea afectado por la secuencia de estos partidores. Para esto usamos el programa IVAR
```bash
ivar trim -x 5 -e -i sample1_untrimmed.bam -b primers.bed -p sample1_trimmed.unsorted
```
* trim le dice a IVAR que queremos usar esa función del programa
* -x 5 le dice a IVAR que queremos cortar cualquier secuencia que comience al menos a 5 bases de un partidor 
* -e le dice a IVAR que incluya las lecturas que no contienen partidores
* -i le dice a IVAR el nombre del archivo que queremos analizar
* -b le dice a IVAR el nombre del archivo que contiene información de la posición de los partidores
* -p  le dice a IVAR que nombre llevará el archivo de salida  (al que IVAR añadirá un .bam)

10.- IVAR nos produce un nuevo archivo .BAM pero que no esta optimizado, por lo que usamos SAMtools nuevamente para optimizar:
```bash
samtools sort -o sample1_trimmed.sorted.bam sample1_trimmed.unsorted.bam
```

### LLAMADO DE VARIANTES
11.- Para descubrir los nucleótidos que difieren con el genoma de referencia (SNPs) debemos usar otra función del programa SAMtools llamada mpileup (ATENCION CON EL DOBLE GUION)
```bash
samtools mpileup -A -aa -d 0 -Q 0 --reference  Wuhan-1.fasta sample1_trimmed.sorted.bam > sample1_pileup.txt
```
* mpileup dice a SAMtools que queremos usar su función para analizar variantes
* -A dice a samtools no descartar reads “anomalas”
* -aa dice a samtools que produzca información de todos los sitios del genoma (no solo los variantes). Esto para saber que sitios del genoma no fueron cubiertos (DEPTH = 0) en caso de que no haya amplificado un amplicon
* -d 0 dice a samtools no tener un limite de profundidad 
* -Q 0 dice  a samtools que no haga filtros de calidad 
* --reference le decimos donde esta nuestro archivo de referencia
* “>” con este símbolo le decimos que el output lo deje en este archivo “sample1_pileup.txt”
12.- Podemos visualizar las primeras 5 líneas del output con el comando head
```bash
head -n 5  sample1_pileup.txt
```

13.- Para analizar las variantes, usaremos IVAR 
```bash
cat sample1_pileup.txt | ivar variants -r Wuhan-1.fasta -p sample1_variants.tsv -m 10
```
* “cat” sirve para abrir el archivo pileup 
* “|” este símbolo sirve para enviar el contenido del archivo pileup al programa IVAR
* -m 10 dice al programa IVAR que queremos solo variantes que estén respaldadas por al menos 10 lecturas (DEPTH 10)

### GENERACION DEL GENOMA CONSENSO

14.- Crearemos un archivo FASTA “corregido” con la información del llamado de variantes. Es decir incluiremos las nuevas variantes al genoma de referencia, pero también incluiremos sitios donde no se supo si había o no una variante (DEPTH <10), que es donde se debe agregar una N.
```bash
cat sample1_pileup.txt | ivar consensus -p sample1_consensus.fasta -m 10 -t 0.5 -n N
```
* “cat” sirve para abrir el archivo pileup 
* “|” este símbolo sirve para enviar el contenido del archivo pileup al programa IVAR
* -m 10 dice al programa IVAR que queremos solo variantes que estén respaldadas por al menos 10 lecturas (DEPTH 10)
* -t 0.5 dice al programa IVAR que queremos solo variantes que tengan una proporción mayor al 50%
* -n N dice a IVAR que cuando no tengamos suficiente información (DEPTH < 10 o t < 0.5) incluir una N (en vez de un nucleótido)
15.- Para imprimir el archivo fasta en pantalla usar cat
```bash
cat sample1_consensus.fasta
```

### ANALIZAR EN NEXTCLADE
Ahora que tenemos un archivo FASTA con la secuencia de nuestra muestra clínica de SARS-CoV-2, podemos analizar si corresponde a una variante usando el sitio nextclade
Ir a https://clades.nextstrain.org/
Subir archivo .fasta (revisar en el Explorador de Windows, disco D)
Run

### Bonus:
Visualizar mapeos en un explorador de genoma
Para esto necesitamos hacer un paso con samtools:
samtools index sample1_trimmed.sorted.bam

Abrir este archivo en el programa de Windows IGV (seleccionando el genoma de SARS-CoV-2)




