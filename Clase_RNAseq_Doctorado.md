Clase\_RNAseq\_Doctorado
================
Carlos Villarroel
3 de diciembre de 2019

``` r
library(tidyr)
library(DESeq2)
library(tximport)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(stringr)
library(ggplot2)
library(pheatmap)
```

``` r
samples <- list.files(full.names = T, pattern="salmon$")

files <- file.path(samples, "quant.sf")

names(files) <- str_replace(samples, "./", "") %>% str_replace(".salmon", "")
```

``` r
txi <- tximport(files, type="salmon", countsFromAbundance="lengthScaledTPM",txOut = T)
```

    ## reading in files with read_tsv

    ## 1 2 3 4

``` r
data <- txi$counts %>% round() %>% data.frame()

sampletype <- factor(c(rep("control",2), rep("treatment", 2)))

meta <- data.frame(sampletype, row.names = colnames(txi$counts))
```

``` r
ggplot(data)+ geom_histogram(aes(x = W1), stat = "bin", bins = 200) + xlab("Raw expression counts")+ylab("Number of genes")
```

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
dds <- DESeqDataSetFromMatrix(data, colData = meta, design = ~sampletype)
```

    ## converting counts to integer mode

``` r
mean_counts <- apply(data[,1:4], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data[,1:4], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 43 rows containing missing values (geom_point).

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
head(counts(dds))
```

    ##            W1     W3    WM1    WM5
    ## YAL001C  2037   2334   1099   1610
    ## YAL002W  1288   1971   3052   3735
    ## YAL003W 15717  14581  14655  12669
    ## YAL005C 80472 128779 380239 415869
    ## YAL007C   955    646    307    484
    ## YAL008W   743    521    195    301

``` r
par(mfrow=c(1,2))
boxplot(log2(counts(dds, normalized=FALSE)+1), main="Raw counts", col=rep(c("lightblue","lightgreen"), each=2), ylim = c(0, 15),las=2)

boxplot(log2(counts(dds, normalized=TRUE)+1), main="Normalized counts", col=rep(c("lightblue","lightgreen"), each=2), ylim = c(0, 15),las=2)
```

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup=c("sampletype"))
```

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
res <- results(dds)
sum(res$padj < 0.05, na.rm=TRUE) #display number of genes with an adjusted p value less than .05
```

    ## [1] 2649

``` r
plotMA(res,ylim=c(-2,2),alpha=0.05,main="MA-plot raw values")
```

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
DEGSUP<-which(res$padj<0.05&res$log2FoldChange > 0.5)

DEGSDOWN<-which(res$padj<0.05&res$log2FoldChange < -0.5)

pheatmap(assay(rld)[DEGSDOWN,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE,scale = "row",fontsize_row = 4,main="Hierarchical cluster",border_color = NA)
```

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
allgenes<-rownames(data)

upgenes<-rownames(data)[DEGSUP]
downgenes<-rownames(data)[DEGSDOWN]
ego <- enrichGO(gene = downgenes, 
                 universe = allgenes,
                 keyType = "ENSEMBL",
                 OrgDb = org.Sc.sgd.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.05, 
                 readable = F, minGSSize = 5)
dotplot(ego, showCategory=10)
```

    ## wrong orderBy parameter; set to default `orderBy = "x"`

![](Clase_RNAseq_Doctorado_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
