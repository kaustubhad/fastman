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

The **fastman** package includes functions for creating Manhattan plots and Q-Q plots from GWAS results. Let us first try using the package on a regular PLINK assoc output file.

### Reading a regular PLINK assoc output dataset
```
m=read.delim("kd2_only_dz.10.a.assoc.linear",header=TRUE,stringsAsFactors=FALSE,sep=" ")
```
This dataset has results for 1,222,628 SNPs on 23 chromosomes. Let us take a look at the data.
```
str(m)
```
```
data.frame':	1222628 obs. of  9 variables:
 $ CHR   : int  1 1 1 1 1 1 1 1 1 1 ...
 $ SNP   : chr  "rs3094315" "rs3094315" "rs3131972" "rs3131972" ...
 $ BP    : int  752566 752566 752721 752721 776546 776546 798959 798959 800007 800007 ...
 $ allele: int  1 3 1 3 1 3 1 3 1 3 ...
 $ afreq : num  0.785 0.215 0.237 0.763 0.959 0.041 0.366 0.634 0.025 0.975 ...
 $ fams  : int  19 19 22 22 5 5 27 27 4 4 ...
 $ beta  : num  2.87 -2.87 -4.77 4.77 -0.9 ...
 $ Z     : num  1.115 -1.115 -1.715 1.715 -0.706 ...
 $ P     : num  0.265 0.265 0.0864 0.0864 0.4802 ...
```
```
head(m)
```
```
  CHR        SNP     BP allele afreq fams   beta      Z        P
1   1  rs3094315 752566      1 0.785   19  2.867  1.115 0.264981
2   1  rs3094315 752566      3 0.215   19 -2.867 -1.115 0.264981
3   1  rs3131972 752721      1 0.237   22 -4.767 -1.715 0.086419
4   1  rs3131972 752721      3 0.763   22  4.767  1.715 0.086419
5   1 rs12124819 776546      1 0.959    5 -0.900 -0.706 0.480177
6   1 rs12124819 776546      3 0.041    5  0.900  0.706 0.480177
```
```
tail(m)
```
```
        CHR       SNP        BP allele afreq fams   beta      Z        P
1222623  23 rs5940540 154833182      1 0.954    5  0.307  0.309 0.757520
1222624  23 rs5940540 154833182      3 0.046    5 -0.307 -0.309 0.757520
1222625  23  rs553678 154892230      1 0.261   16 -4.372 -2.081 0.037399
1222626  23  rs553678 154892230      3 0.739   16  4.372  2.081 0.037399
1222627  23  rs669237 154916845      1 0.261   16 -3.910 -1.893 0.058363
1222628  23  rs669237 154916845      2 0.739   16  3.910  1.893 0.058363
```
Let us see the distribution of SNPs across chromosomes.
```
as.data.frame(table(m$CHR))
```
```
   Var1  Freq
1     1 98286
2     2 97472
3     3 80762
4     4 69960
5     5 71526
6     6 81240
7     7 65046
8     8 64048
9     9 57040
10   10 66356
11   11 61970
12   12 60638
13   13 46970
14   14 39746
15   15 37068
16   16 37952
17   17 33426
18   18 36380
19   19 24388
20   20 31024
21   21 17676
22   22 17392
23   23 26262
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
