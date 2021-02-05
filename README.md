# fastman

## Description
An **R package** for fast and efficient visualizing of GWAS results using Q-Q and Manhattan plots directly from PLINK output files.
* **Fast**: Drastically reduces time in plot generation compared to qqman. On a typical imputed PLINK assoc file of 10 million SNPs, plotting time is reduced from 737s in qqman to 60s.
* **Efficient**: Optimized memory management
* **Versatile**: Can handle various inputs from p-values, logarithms of p-values to FST scores. Compatible plotting with other genome-wide population genetic parameters (e.g. FST, pi and D statistics)
* **Familiar**: Has a very similar set of input arguments and code structure compared to qqman.
* **Non-model friendly**: Additional support for results from genomes of non-model organisms (often with hundreds of contigs or many scaffolds), alphabetical and other ordering options.

## Functions:

### 1. fastman

#### Description
Creates a Manhattan plot directly from a PLINK assoc output (or any data frame with chromosome, position, and p-value).

#### Usage
```
fastman (m, chr = "CHR", bp = "BP", p = "P", snp, chrlabs, speedup=TRUE, logp = TRUE, col="matlab", maxP, sortchr=TRUE,
        bybp=FALSE, chrsubset, bprange, highlight, annotateHighlight=FALSE, annotatePval, annotateTop=TRUE, annotationWinMb,
        annotateN, annotationCol, annotationAngle=45, suggestiveline, genomewideline, cex=0.4, cex.axis=0.6,
        xlab, ylab, xlim, ylim, ...)
```

#### Parameters:
* **m**	= A data frame with columns "BP," "CHR," "P," and optionally, "SNP”.
* **chr**	= A string denoting the column name for the chromosome. Defaults to “CHR”, which corresponds to the PLINK –assoc command output. The column must be numeric.
* **bp** = A string denoting the column name for the chromosomal position. Defaults to “BP”, which corresponds to the PLINK –assoc command output. The column must be numeric.
* **p**	= A string denoting the column name for the p-values or scores for the SNP association tests. Defaults to “P”, which corresponds to the PLINK –assoc command output. The column must be numeric.
* **snp**	= A string denoting the column name for the SNP name (rs number). The column must be character.
* **chrlabs**	= A character vector of length equal to the number of chromosomes specifying the chromosome labels (e.g., c (1:22, "X", "Y", "MT")).
* **speedup** = A logical value; if TRUE, the function employs the faster method where input values above 99.8% are rounded to 3 digits, and the rest is rounded to 2 digits.
* **logp**	= A logical value; if TRUE, negative logarithms (base 10) of p-values is plotted.
* **col**	= A string indicating the color scheme of the plot. Defaults to “matlab”.
* **maxP**	= A numeric value indicating the maximum negative logarithm of p-value till which user wants to visualize.
* **sortchr**	= A logical value; if TRUE, the table is sorted by chromosome number before plotting. If not specified by user, the function takes default value TRUE.
* **bybp** = A logical value; if TRUE, the table is sorted by chromosome positions before plotting. If not specified by user, the function takes default value FALSE.
* **chrsubset** = The range of chromosome numbers to be plotted.
* **bprange**	= The range of chromosome positions to be plotted.
* **highlight**	= A character vector of SNPs in the dataset to highlight. These SNPs should all be in the dataset.
* **annotateHighlight**	= A logical value; if TRUE, annotates all highlighted SNPs in case more specific annotation instructions are not provided.
* **annotatePval**	= A numeric value, if set, SNPs with p-values below this will be annotated on the plot.
* **annotateTop**	= A logical value; If TRUE, only annotates the top hit on each chromosome that is below the annotatePval threshold.
* **annotationWinMb**	= A numeric value, if set, will determine the megabase window within which the top SNP will be highlighted.
* **annotateN**	= A numeric value, if set, this number of top SNPs will be annotated on the plot.
* **annotationCol** = A string indicating the color of annotation. Defaults to grey.
* **annotationAngle**	= The angle of annotation, defaults to 45 degree.
* **suggestiveline**	= The position to draw a "suggestive" line. Defaults to -log10(1e-5).
* **genomewideline**	= The position to draw a “genome-wide significant” line. Defaults to -log10(5e-8).
* **cex** = A a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 0.4.
* **cex.axis**	= The magnification to be used for axis annotation relative to the current setting of cex. Defaults to 0.6.
* **xlab**	= A label for the x axis, defaults to a description of x.
* **ylab**	= A label for the y axis, defaults to a description of y.
* **xlim**	= the x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a ‘reversed axis’. The default value, NULL, indicates that the range of the finite values to be plotted should be used.
* **ylim** = The y limits of the plot.

#### Value
A Manhattan Plot

### 2. fastqq

#### Description
Creates a quick quantile-quantile plot from GWAS outputs. On a typical imputed assoc file of 10 million SNPs, it reduces plotting time from 400s in qqman to 24s.

#### Usage
```
fastqq (p, speedup=TRUE, lambda=TRUE, fix_zero=TRUE, cex=0.6, cex.axis=0.9, xlab, ylab, ...)
```

#### Parameters:
* **p**	= A numeric vector of p-values.
* **speedup**	= A logical value; if TRUE, the function employs the faster method where inputs are rounded to 3 digits.
* **lambda** = A logical value; if TRUE, the genomic inflation factor lambda is shown on the top-left corner.
* **fix_zero**	= A logical value; if TRUE, zero input values are converted to half of the minimum observed input value; if FALSE, zero input values are removed.
* **cex**	= A numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 0.6.
* **cex.axis**	= The magnification to be used for axis annotation relative to the current setting of cex. Defaults to 0.9.
* **xlab** = A label for the x axis, defaults to a description of x.
* **ylab**	= A label for the y axis, defaults to a description of y.

#### Value
* A Q-Q Plot
* The genomic inflation factor lambda

## Examples

The **fastman** package includes functions for creating Manhattan plots and Q-Q plots from GWAS results. Let us first try using the package on a regular PLINK assoc output file. This is a GWAS summary statistics output file from a study by the authors (https://doi.org/10.1038/ncomms10815), available from the GWAS Central repository ( https://www.gwascentral.org/study/HGVST2597)

### Reading a regular PLINK assoc output dataset
```
m=read.table("beard.assoc.linear",header=TRUE,stringsAsFactors=FALSE,sep="\t")
```
This dataset has results for 8,776,423 SNPs on 22 chromosomes. Let us take a look at the data.
```
str(m)
```
```
'data.frame':	8776423 obs. of  9 variables:
 $ CHR  : chr  "1" "1" "1" "1" ...
 $ SNP  : chr  "rs58108140" "rs180734498" "rs116400033" "rs62637813" ...
 $ BP   : chr  "10583" "13302" "51479" "52058" ...
 $ A1   : chr  "A" "T" "A" "C" ...
 $ TEST : chr  "ADD" "ADD" "ADD" "ADD" ...
 $ NMISS: chr  "1551" "1911" "1552" "2348" ...
 $ BETA : chr  "0.1486" "0.07006" "0.02572" "0.1391" ...
 $ STAT : chr  "1.793" "0.7667" "0.3814" "2.035" ...
 $ P    : chr  "0.07314" "0.4433" "0.703" "0.04199" ...
```
```
head(m)
```
```
  CHR         SNP    BP A1 TEST NMISS      BETA     STAT       P
1   1  rs58108140 10583  A  ADD  1551    0.1486    1.793 0.07314
2   1 rs180734498 13302  T  ADD  1911   0.07006   0.7667  0.4433
3   1 rs116400033 51479  A  ADD  1552   0.02572   0.3814   0.703
4   1  rs62637813 52058  C  ADD  2348    0.1391    2.035 0.04199
5   1 rs201374420 52185  T  ADD  2421 -0.004824 -0.04033  0.9678
6   1 rs150021059 52238  T  ADD  2290   -0.1036  -0.8664  0.3864
```
```
tail(m)
```
```
        CHR         SNP       BP A1 TEST NMISS      BETA    STAT      P
8776418  22 rs144549712 51229855  A  ADD  2206   0.02949  0.5249 0.5997
8776419  22  rs62240042 51233300  T  ADD  1536 -0.009488 -0.2751 0.7833
8776420  22 rs200507571 51236013 AT  ADD  1847  0.003518 0.08454 0.9326
8776421  22   rs3896457 51237063  C  ADD  1802   -0.0187 -0.5628 0.5737
8776422  22 rs149733995 51238249  C  ADD  2268   0.06418  0.9966  0.319
8776423  22 rs181833046 51243297  T  ADD  2250    0.1005   1.066 0.2866
```
### Creating Manhattan Plots
From the above dataset, lets generate a basic Manhattan plot.
```
png("md1.png", width=10, height=6, units="in", res=300)
fastman(m)
dev.off()
```

![](https://github.com/kaustubhad/fastman/blob/main/md1.png)

We can change some basic graph parameters. Let us increase the y-axis limit (```ylim=```) to 8, reduce the point size (```cex=```) to 30% and reduce the font size of the axis labels (```cex.axis=```) to 40%. We can also change the colour palette (```col=```) and remove the suggestive (```suggestiveline=```) and genome-wide (```genomewideline=```) significance lines.
```
png("md2.png", width=10, height=6, units="in", res=300)
fastman(m, ylim = c(0,8), cex = 0.3, cex.axis = 0.4, col = "rainbow1", suggestiveline = FALSE, genomewideline = FALSE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/md2.png)

We can now look into the SNPs of a single chromosome.
```
png("md3.png", width=10, height=6, units="in", res=300)
fastman(m, chrsubset=1)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/md3.png)

Let us say we are interested in highlighting some particular 1000 SNPs in chromosome 1. We have the name of the required SNPs in a character vector called snp1.
```
str(snp1)
```
```
chr [1:1000] "rs3094315" "rs3094315" "rs3131972" "rs3131972" "rs12124819" ...
```
```
png("md4.png", width=10, height=6, units="in", res=300)
fastman(m, highlight = snp1)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/md4.png)

We can annotate SNPs based on their p-value. By default, among the SNPs that exceed the provided threshold, only the top SNP in every chromosome is annotated.
```
png("md5.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-3)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/md5.png)

We can override the default rule, and annotate all the SNPs beyond our specified threshold.
```
png("md6.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-3, annotateTop = FALSE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/md6.png)

We can annotate among highlighted SNPs as well. By default, only the top SNP in every chromosome will be highlighted.
```
png("md7.png", width=10, height=6, units="in", res=300)
fastman(m, highlight = snp1, annotateHighlight = TRUE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/md7.png)
