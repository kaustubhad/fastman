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
```
m=read.delim("kd2_only_dz.10.a.assoc.linear",header=TRUE,stringsAsFactors=FALSE,sep=" ");
```
This dataset has results for 1,222,628 SNPs on chromosomes. Let us take a look at the data.
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
