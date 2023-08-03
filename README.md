# fastman

## Description
An **R package** for fast and efficient visualizing of GWAS results using Q-Q and Manhattan plots directly from PLINK output files.
* **Fast**: Drastically reduces time in plot generation compared to qqman. On a typical imputed PLINK assoc file of 10 million SNPs, plotting time is reduced from 737s in qqman to 60s.
* **Efficient**: Optimized memory management
* **Versatile**: Can handle various inputs from p-values, logarithms of p-values to FST scores. Compatible plotting with other genome-wide population genetic parameters (e.g. FST, pi and D statistics). Allows both-sided scores, e.g. scores with negative values.
* **Non-model friendly**: Additional support for results from genomes of non-model organisms (often with hundreds of contigs or many scaffolds), alphabetical and other ordering options.
* **Annotation and Highlight versatility**: Has a wide set of options to customize annotating and highlighting SNPs of interest.
* **Familiar**: Has a very similar set of input arguments and code structure compared to qqman.
* **Missing Value Handling**: Can handle missing values in the input dataframe.

## Installation
```
devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
```

## Functions:

### 1. fastman

#### Description
Creates a Manhattan plot directly from a PLINK assoc output (or any data frame with chromosome, position, and p-value).

#### Usage
```
fastman (m, chr = "CHR", bp = "BP", p = "P", snp, chrlabs, speedup=TRUE, logp = TRUE, scattermore = FALSE,
        col="matlab", maxP=14, sortchr=TRUE, bybp=FALSE, chrsubset, bprange, highlight,
        annotateHighlight=FALSE, annotatePval, colAbovePval=FALSE, col2="greys", annotateTop=TRUE,
        annotationWinMb, annotateN, annotationCol, annotationAngle=45, baseline=NULL, suggestiveline,
        genomewideline, cex=0.4, cex.text=0.4, cex.axis=0.6, scattermoresize = c(3000,1800), xlab, ylab,
        xlim, ylim, ...)
```

#### Parameters:
* **m**	= A data frame containing data for producing the manhattan plot. Has to contain a minimum of three columns: base pair position, chromosome ID, and P-value: defaults are "BP", "CHR", "P" following the Plink assoc file convention. And optionally some sort of ID, e.g. the SNP ID (default "SNP”) if annotations are needed. See explanations of the next four parameters if your column names are different, e.g. for non-model organisms you can use contid ID instead of chromosome ID.
* **chr**	= A string denoting the column name for the chromosome. Defaults to “CHR”, which corresponds to the PLINK –assoc command output. For non-model organisms, this could be the contig ID. In case your chromosome column is actually numeric but has been converted into string during the reading of data in R, you must pay close attention to the sorting order of chromosomes. If you still want the chromosomes to be sorted in an increasing order of chromosome number then you must convert your chromosome column to numeric before using the function. If your data is already sorted then you do not need to convert the column to numeric, you can specify ```sortchr = FALSE``` in your input arguments instead.
* **bp** = A string denoting the column name for the chromosomal position. Defaults to “BP”, which corresponds to the PLINK –assoc command output. The column must be numeric.
* **p**	= A string denoting the column name for the p-values or scores for the SNP association tests. Defaults to “P”, which corresponds to the PLINK –assoc command output. The column must be numeric. You can also provide alread-computed log of p-values, e.g. from published summary statistics.
* **snp**	= A string denoting the column name for the SNP name (rs number). The column must be character.
* **chrlabs**	= An optional character vector of length equal to the number of chromosomes, specifying the chromosome labels. e.g., you can provide ```c(1:22, "X", "Y", "MT")``` to convert the Plink numerical notation of 23=X, 24=Y, etc. This character vector is used to create the axis labels of the manhattan plot. So, you must sort the character vector in the order you want the chromosome labels to appear in the final plot. For example, if your input dataframe has chromosome numbers in a particular order you specifically want, and you have used the option ```sortchr = FALSE``` to preserve the order for your final plot, then your ```chrlabs``` vector should also have the same order of chromosomes.
* **speedup** = A logical value; if TRUE, the function employs the faster method where input values at the extreme 0.2% are rounded to 3 digits, and the rest is rounded to 2 digits. The default value of this parameter is TRUE.
* **logp**	= A logical value; if TRUE, negative logarithms (base 10) of p-values are plotted. In case the user wants to use FST score type data or logarithm of p-values directly, then logp must be stated to be FALSE, as the default value of this parameter is TRUE.
* **scattermore** = A logical value; if TRUE, uses ```scattermore``` package to speed up plot generation faster. In case the user wants to use this feature, the ``scattermore`` package needs to be installed and loaded before running the command. The default value of this parameter is FALSE.
* **col**	= A string indicating the color scheme of the plot. Defaults to “matlab”. There are various options available for user. See below for details.
* **maxP**	= A numeric value indicating the maximum y-value till which user wants to visualize. The default value of this parameter is 14. If the data has negative values then both sides are truncated till the absolute value of the parameter. The user can provide NULL as input if truncation is not required.
* **sortchr**	= A logical value; if TRUE, the table is sorted by chromosome number before plotting. If not specified by user, the function takes default value TRUE.
* **bybp** = A logical value; if TRUE, the y-values are plotted against chromosome positions. In this case the table is not sorted by chromosome number before plotting. If not specified by user, the function takes default value FALSE. This feature is useful especially for plots where the user might be interested in studying the association p-values across contigs.
* **chrsubset** = The subset of chromosome numbers to be plotted.
* **bprange**	= The range of chromosome positions to be plotted. In case the user wants to subset the X-axis by region, then this should be the parameter of choice, not xlim.
* **highlight**	= A character vector of SNPs in the dataset to highlight. These SNPs should all be in the dataset.
* **annotateHighlight**	= A logical value; if TRUE, annotates all highlighted SNPs in case more specific annotation instructions are not provided.
* **annotatePval**	= A numeric value, if set, SNPs with p-values below this will be annotated on the plot. In case of p-value, the user can provide either the p-value or the negative logarithm of p-value as input for this argument, whichever is convenient. In case of scores, the user can provide the score cutoff directly as input. 
* **colAbovePval**	= A logical value, if TRUE, will colour all hits above the specified p-value threshold using colour scheme chosen in col argument (default "matlab"), while the points below the threshold will be coloured using the colour scheme chosen in col2 argument below (default "greys"). Defaults to FALSE.
* **col2** = A string indicating the color scheme of the part of the plot below the specified p-value threshold. Defaults to “greys”. There are various options available for user. See below for details.
* **annotateTop**	= A logical value; If TRUE, only annotates the top hit on each chromosome that is below the annotatePval threshold. This is just a modifier, and it works only when used with either annotatePval, annotateHighlight or annotateN.
* **annotationWinMb**	= A numeric value, if set, will determine the megabase window within which the top SNP will be annotated. This is just a modifier, and it works only when used with either annotatePval, annotateHighlight or annotateN.
* **annotateN**	= A numeric value, if set, this number of top SNPs will be annotated on the plot.
* **annotationCol** = A string indicating the color of annotation or the column name containing the colour information for the individual rows. The user can provide a column as a part of the input data frame which contains annotation colour information corresponding to each individual SNPs and specify the column name in this parameter. If the user does not want annotation for some particular SNP, NA can be provided in this column corresponding to that SNP. In case the user provides just a string indicating the annotation colour instead of a column name, then the same colour will be used for annotating all the SNPs. Defaults to grey.
* **annotationAngle**	= The angle of annotation, defaults to 45 degree.
* **baseline**	= The position to draw a baseline in black. Defaults to NULL, as a typical manhattan plot already has a baseline at y=0. In case the data has a left tail (e.g. two-sided scoes) the user might want to provide a baseline position for reference. In case multiple baselines are required, the user can provide a vector of positions.
* **suggestiveline**	= The position to draw a GWAS "suggestive significance" line. Defaults to -log10(1e-5). In case multiple suggestive lines are required, the user can provide a vector of positions.
* **genomewideline**	= The position to draw a GWAS “genome-wide significance” line. Defaults to -log10(5e-8). In case multiple genome-wide significant lines are required, the user can provide a vector of positions.
* **cex** = A numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 0.4.
* **cex.text** = A numerical vector giving the amount by which annotation text should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 0.4.
* **cex.axis**	= The magnification to be used for axis annotation relative to the current setting of cex. Defaults to 0.6.
* **scattermoresize** = A 2-element integer vector to specify the size of scattermore plot in pixels. Applicable only if ```scattermore = TRUE```. User should make sure that the value of ```scattermoresize``` parameter matches with the desired size of final image output for best quality plots. Defaults to c(3000,1800).
* **xlab**	= A label for the x axis, defaults to a description of x.
* **ylab**	= A label for the y axis, defaults to a description of y.
* **xlim**	= The x limits of the plot. The user should refrain from changing xlim in order to subset x-axis by region. The better option is to specify this in the bprange arguent, as changing xlim might lead to improper scaling and spacing in the plot.
* **ylim** = The y limits of the plot. The user should refrain from changing ylim in order to truncate y-axis. The better option is to specify the same in the maxP arguent, as changing ylim might lead to improper scaling and spacing in the plot.

#### Value
A Manhattan Plot

### 2. fastqq

#### Description
Creates a quick quantile-quantile plot from GWAS outputs. On a typical imputed assoc file of 10 million SNPs, it reduces plotting time from 400s in qqman to 24s.

#### Usage
```
fastqq (p1, p2=NULL, colour, logtransform=TRUE, pairwisecompare=TRUE, speedup=TRUE, lambda=TRUE, maxP=14,
        fix_zero=TRUE, cex=0.6, cex.axis=0.9, xlab, ylab, ...)
```

#### Parameters:
* **p1**	= A numeric vector of p-values. If the user has a single set of p-values, then this is the input for that. In case the user wants to compare two sets of p-values, then this argument takes the first set of p-values as input.
* **p2**	= A numeric vector of p-values. In case the user wants to compare two sets of p-values, this argument takes the second set of p-values as input. Otherwise, this argument defaults to NULL.
* **logtransform** = A logical value; if TRUE, negative logarithms (base 10) of p-values are plotted. In case the user wants to use FST score type data or logarithm of p-values directly, then logtransform must be stated to be FALSE, as the default value of this parameter is TRUE.
* **pairwisecompare** = A logical value. If the two sets of p-values are provided in pairs, the user should state this argument as TRUE. Otherwise, both sets of p-values are sorted individually and then plotted against each other.
* **speedup**	= A logical value; if TRUE, the function employs the faster method where inputs are rounded to 3 digits.
* **lambda** = A logical value; if TRUE, the genomic inflation factor lambda is shown on the top-left corner.
* **maxP**	= A numeric value indicating the maximum negative lograrithm of p-value till which user wants to visualize. The default value of this parameter is 14. If the data has negative values then both sides are truncated till the absolute value of the parameter. The user can provide NULL as input if truncation is not required.
* **fix_zero**	= A logical value; if TRUE, zero input values are converted to half of the minimum observed input value; if FALSE, zero input values are removed.
* **cex**	= A numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 0.6.
* **cex.axis**	= The magnification to be used for axis annotation relative to the current setting of cex. Defaults to 0.9.
* **xlab** = A label for the x axis, defaults to a description of x.
* **ylab**	= A label for the y axis, defaults to a description of y.

#### Value
* A Q-Q Plot
* The genomic inflation factor lambda

### 3. fastman_gg

#### Description
Creates a Manhattan ggplot object directly from a PLINK assoc output (or any data frame with chromosome, position, and p-value). User needs to install and load ```ggplot2``` package to use this function.

#### Usage
```
fastman_gg (m, chr = "CHR", bp = "BP", p = "P", snp, chrlabs, speedup=TRUE, logp = TRUE,
        scattermore = FALSE, repel = FALSE, col="matlab", maxP=14, sortchr=TRUE, bybp=FALSE, chrsubset,
        bprange, highlight, annotateHighlight=FALSE, annotatePval, colAbovePval=FALSE, col2="greys",
        annotateTop=TRUE, annotationWinMb, annotateN, annotationCol, annotationAngle=45, baseline=NULL,
        suggestiveline, genomewideline, cex=0.4, cex.text=0.4, cex.axis=0.6,
        scattermoresize = c(3000,1800), xlab, ylab, xlim, ylim, ...)
```

#### Parameters:
* **m**	= A data frame containing data for producing the manhattan plot. Has to contain a minimum of three columns: base pair position, chromosome ID, and P-value: defaults are "BP", "CHR", "P" following the Plink assoc file convention. And optionally some sort of ID, e.g. the SNP ID (default "SNP”) if annotations are needed. See explanations of the next four parameters if your column names are different, e.g. for non-model organisms you can use contid ID instead of chromosome ID.
* **chr**	= A string denoting the column name for the chromosome. Defaults to “CHR”, which corresponds to the PLINK –assoc command output. For non-model organisms, this could be the contig ID. In case your chromosome column is actually numeric but has been converted into string during the reading of data in R, you must pay close attention to the sorting order of chromosomes. If you still want the chromosomes to be sorted in an increasing order of chromosome number then you must convert your chromosome column to numeric before using the function. If your data is already sorted then you do not need to convert the column to numeric, you can specify ```sortchr = FALSE``` in your input arguments instead.
* **bp** = A string denoting the column name for the chromosomal position. Defaults to “BP”, which corresponds to the PLINK –assoc command output. The column must be numeric.
* **p**	= A string denoting the column name for the p-values or scores for the SNP association tests. Defaults to “P”, which corresponds to the PLINK –assoc command output. The column must be numeric. You can also provide alread-computed log of p-values, e.g. from published summary statistics.
* **snp**	= A string denoting the column name for the SNP name (rs number). The column must be character.
* **chrlabs**	= An optional character vector of length equal to the number of chromosomes, specifying the chromosome labels. e.g., you can provide ```c(1:22, "X", "Y", "MT")``` to convert the Plink numerical notation of 23=X, 24=Y, etc. This character vector is used to create the axis labels of the manhattan plot. So, you must sort the character vector in the order you want the chromosome labels to appear in the final plot. For example, if your input dataframe has chromosome numbers in a particular order you specifically want, and you have used the option ```sortchr = FALSE``` to preserve the order for your final plot, then your ```chrlabs``` vector should also have the same order of chromosomes.
* **speedup** = A logical value; if TRUE, the function employs the faster method where input values at the extreme 0.2% are rounded to 3 digits, and the rest is rounded to 2 digits. The default value of this parameter is TRUE.
* **logp**	= A logical value; if TRUE, negative logarithms (base 10) of p-values are plotted. In case the user wants to use FST score type data or logarithm of p-values directly, then logp must be stated to be FALSE, as the default value of this parameter is TRUE.
* **scattermore** = A logical value; if TRUE, uses ```scattermore``` package to speed up plot generation faster. In case the user wants to use this feature, the ``scattermore`` package needs to be installed and loaded before running the command. The default value of this parameter is FALSE.
* **repel**     = A logical value; if TRUE, uses ```ggrepel```package to repel overlapping text labels. In case the user wants to use this feature, the ```ggrepel``` package needs to be installed and loaded before running the command. The default value of this parameter is FALSE.
* **col**	= A string indicating the color scheme of the plot. Defaults to “matlab”. There are various options available for user. See below for details.
* **maxP**	= A numeric value indicating the maximum y-value till which user wants to visualize. The default value of this parameter is 14. If the data has negative values then both sides are truncated till the absolute value of the parameter. The user can provide NULL as input if truncation is not required.
* **sortchr**	= A logical value; if TRUE, the table is sorted by chromosome number before plotting. If not specified by user, the function takes default value TRUE.
* **bybp** = A logical value; if TRUE, the y-values are plotted against chromosome positions. In this case the table is not sorted by chromosome number before plotting. If not specified by user, the function takes default value FALSE. This feature is useful especially for plots where the user might be interested in studying the association p-values across contigs.
* **chrsubset** = The subset of chromosome numbers to be plotted.
* **bprange**	= The range of chromosome positions to be plotted. In case the user wants to subset the X-axis by region, then this should be the parameter of choice, not xlim.
* **highlight**	= A character vector of SNPs in the dataset to highlight. These SNPs should all be in the dataset.
* **annotateHighlight**	= A logical value; if TRUE, annotates all highlighted SNPs in case more specific annotation instructions are not provided.
* **annotatePval**	= A numeric value, if set, SNPs with p-values below this will be annotated on the plot. In case of p-value, the user can provide either the p-value or the negative logarithm of p-value as input for this argument, whichever is convenient. In case of scores, the user can provide the score cutoff directly as input. 
* **colAbovePval**	= A logical value, if TRUE, will colour all hits above the specified p-value threshold using colour scheme chosen in col argument (default "matlab"), while the points below the threshold will be coloured using the colour scheme chosen in col2 argument below (default "greys"). Defaults to FALSE.
* **col2** = A string indicating the color scheme of the part of the plot below the specified p-value threshold. Defaults to “greys”. There are various options available for user. See below for details.
* **annotateTop**	= A logical value; If TRUE, only annotates the top hit on each chromosome that is below the annotatePval threshold. This is just a modifier, and it works only when used with either annotatePval, annotateHighlight or annotateN.
* **annotationWinMb**	= A numeric value, if set, will determine the megabase window within which the top SNP will be annotated. This is just a modifier, and it works only when used with either annotatePval, annotateHighlight or annotateN.
* **annotateN**	= A numeric value, if set, this number of top SNPs will be annotated on the plot.
* **annotationCol** = A string indicating the color of annotation or the column name containing the colour information for the individual rows. The user can provide a column as a part of the input data frame which contains annotation colour information corresponding to each individual SNPs and specify the column name in this parameter. If the user does not want annotation for some particular SNP, NA can be provided in this column corresponding to that SNP. In case the user provides just a string indicating the annotation colour instead of a column name, then the same colour will be used for annotating all the SNPs. Defaults to grey.
* **annotationAngle**	= The angle of annotation, defaults to 45 degree.
* **baseline**	= The position to draw a baseline in black. Defaults to NULL, as a typical manhattan plot already has a baseline at y=0. In case the data has a left tail (e.g. two-sided scoes) the user might want to provide a baseline position for reference. In case multiple baselines are required, the user can provide a vector of positions.
* **suggestiveline**	= The position to draw a GWAS "suggestive significance" line. Defaults to -log10(1e-5). In case multiple suggestive lines are required, the user can provide a vector of positions.
* **genomewideline**	= The position to draw a GWAS “genome-wide significance” line. Defaults to -log10(5e-8). In case multiple genome-wide significant lines are required, the user can provide a vector of positions.
* **cex** = A numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 0.9.
* **cex.text** = A numerical vector giving the amount by which annotation text should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Defaults to 1.8.
* **cex.axis**	= The magnification to be used for axis annotation relative to the current setting of cex. Defaults to 0.6.
* **scattermoresize** = A 2-element integer vector to specify the size of scattermore plot in pixels. Applicable only if ```scattermore = TRUE```. User should make sure that the value of ```scattermoresize``` parameter matches with the desired size of final image output for best quality plots. Defaults to c(3000,1800).
* **xlab**	= A label for the x axis, defaults to a description of x.
* **ylab**	= A label for the y axis, defaults to a description of y.
* **xlim**	= The x limits of the plot. The user should refrain from changing xlim in order to subset x-axis by region. The better option is to specify this in the bprange arguent, as changing xlim might lead to improper scaling and spacing in the plot.
* **ylim** = The y limits of the plot. The user should refrain from changing ylim in order to truncate y-axis. The better option is to specify the same in the maxP arguent, as changing ylim might lead to improper scaling and spacing in the plot.

#### Value
A ggplot object

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
 $ CHR  : num  1 1 1 1 1 1 1 1 1 1 ...
 $ SNP  : chr  "rs58108140" "rs180734498" "rs116400033" "rs62637813" ...
 $ BP   : num  10583 13302 51479 52058 52185 ...
 $ A1   : chr  "A" "T" "A" "C" ...
 $ TEST : chr  "ADD" "ADD" "ADD" "ADD" ...
 $ NMISS: num  1551 1911 1552 2348 2421 ...
 $ BETA : num  0.1486 0.07006 0.02572 0.1391 -0.00482 ...
 $ STAT : num  1.793 0.7667 0.3814 2.035 -0.0403 ...
 $ P    : num  0.0731 0.4433 0.703 0.042 0.9678 ...
```
```
head(m)
```
```
  CHR         SNP    BP A1 TEST NMISS      BETA     STAT       P
1   1  rs58108140 10583  A  ADD  1551  0.148600  1.79300 0.07314
2   1 rs180734498 13302  T  ADD  1911  0.070060  0.76670 0.44330
3   1 rs116400033 51479  A  ADD  1552  0.025720  0.38140 0.70300
4   1  rs62637813 52058  C  ADD  2348  0.139100  2.03500 0.04199
5   1 rs201374420 52185  T  ADD  2421 -0.004824 -0.04033 0.96780
6   1 rs150021059 52238  T  ADD  2290 -0.103600 -0.86640 0.38640
```
```
tail(m)
```
```
        CHR         SNP       BP A1 TEST NMISS      BETA     STAT      P
8776418  22 rs144549712 51229855  A  ADD  2206  0.029490  0.52490 0.5997
8776419  22  rs62240042 51233300  T  ADD  1536 -0.009488 -0.27510 0.7833
8776420  22 rs200507571 51236013 AT  ADD  1847  0.003518  0.08454 0.9326
8776421  22   rs3896457 51237063  C  ADD  1802 -0.018700 -0.56280 0.5737
8776422  22 rs149733995 51238249  C  ADD  2268  0.064180  0.99660 0.3190
8776423  22 rs181833046 51243297  T  ADD  2250  0.100500  1.06600 0.2866
```
### Creating Manhattan Plots
From the above dataset, lets generate a basic Manhattan plot.
```
png("md1.png", width=10, height=6, units="in", res=300)
fastman(m)
dev.off()
```

![](https://github.com/kaustubhad/fastman/blob/main/plots/md1.png)

Before moving into any further details, let us look into the plethora of options that this package provides us in a single plot. We will understand these options more rigorously as we go into the later parts of this exercise.

![](https://github.com/kaustubhad/fastman/blob/main/plots/test2.png)

Let us compare the time of plot generation with qqman. For this purpose, we are going to use a library ```tictoc``` which will record the run time.
```
library(tictoc)
tic(); png("md1.png", width=10, height=6, units="in", res=300); fastman(m); dev.off(); toc();
```

This will record the run time of fastman.
```
77.238 sec elapsed
```

Now, lets run the same code for qqman.
```
library(qqman)
tic(); png("md1a.png", width=10, height=6, units="in", res=300); manhattan(m); dev.off(); toc();
```

Lets find out the run time of qqman.
```
630.995 sec elapsed
```

As seen above, fastman reduces the run time drastically.

We can change some basic graph parameters. Let us reduce the point size (```cex=```) to 30% and reduce the font size of the axis labels (```cex.axis=```) to 50%. We can also change the colour palette (```col=```) and remove the suggestive (```suggestiveline=```) and genome-wide (```genomewideline=```) significance lines.
```
png("md2.png", width=10, height=6, units="in", res=300)
fastman(m, cex=0.3, cex.axis=0.5, col="rainbow1", suggestiveline=FALSE, genomewideline=FALSE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md2.png)

In this plot, if we want to show p-values till 1E-10, we can set ```maxP=10```, instead of changing the ```ylim```.
```
png("md2a.png", width=10, height=6, units="in", res=300)
fastman(m, cex=0.3, cex.axis=0.5, col="rainbow1", suggestiveline=FALSE, genomewideline=FALSE, maxP=10)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md2a.png)

We can now look into the SNPs of a single chromosome.
```
png("md3.png", width=10, height=6, units="in", res=300)
fastman(m, chrsubset=3)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md3.png)

Let us say we are interested in highlighting some particular 10000 SNPs in chromosome 1. We have the name of the required SNPs in a character vector called snp1.
```
str(snp1)
```
```
chr [1:10000] "rs12618998" "rs201530365" "rs13001505" "rs6742479" "rs6742805" ...
```
```
png("md4.png", width=10, height=6, units="in", res=300)
fastman(m, highlight = snp1)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md4.png)

We can also annotate the highlighted SNPs. By default, only the top SNP in every chromosome will be annotated.
```
png("md7.png", width=10, height=6, units="in", res=300)
fastman(m, highlight = snp1, annotateHighlight = TRUE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md7.png)

Let us discuss the annotations in details. We have seven input parameters with respect to annotations: ```annotateHighlight```, ```annotatePval```, ```annotateTop```, ```annotationWinMb```, ```annotateN```, ```annotationCol``` and ```annotationAngle```. We have already discussed about ```annotateHighlight```. Among the rest, ```annotatePval``` and ```annotateN``` are annotation criteria, while the other four parameters are annotation modifiers.

As stated above, we can annotate our plot in only two ways, either by providing a p-value threshold or by providing the required number of top SNPs. The p-value threshold can be provided using the ```annotatePval``` criterion. The user can provide either p-value or the negative logarithm as an input for this parameter, whichever is convenient. For example, both 1E-7 and 7 are acceptable input values.
```
png("md5.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-7)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md5.png)

We can see that, among the SNPs that exceed the provided threshold, only the top SNP in every chromosome is annotated. We have observed the same for ```annotateHighlight``` as well. We will come to that later, when we discuss the annotation modifiers.

For now, let us explore the other annotation criterion ```annotateN```. Lets say we want to annotate the top 20 SNPs of the data.
```
png("md8.png", width=10, height=6, units="in", res=300)
fastman(m, annotateN=20)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md8.png)

Lets go back to the observation we made while using ```annotateHighlight``` and ```annotatePval```. By default, only the top SNP in every chromosome is annotated. We can override this default rule using the annotation modifier ```annotateTop```. This modifier has a default value ```TRUE``` unless mentioned otherwise by the user. Lets consider the example where we are annotating the SNPs beyond a specified p-value threshold. By setting ```annotateTop=FALSE```, we can annotate all the SNPs beyond our threshold instead of only the top SNP per chromosome.
```
png("md6.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-7, annotateTop=FALSE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md6.png)

Now we move on to the next annotation modifier ```annotationWinMb```. Instead of annotating the top SNP in every chromosome, if we want to annotate the top SNP within our chosen megabase window, then this is the modifier for us. For this example, lets first reduce our p-value threshold to 1E-5.
```
png("md9.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md9.png)

Now, lets specify a 5 megabase window.
```
png("md10.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, annotationWinMb=5)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md10.png)

Among the remaining two modifiers, ```annotationCol``` sets the annotation colour and ```annotationAngle``` specifies the angle of annotation. The default colour is "gray50" and the default angle is 45 degrees. Let us illustrate with an example.
```
png("md11.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, annotationCol="red", annotationAngle= 60)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md11.png)

As stated already, the default annotation colour in our function is "gray50". However, if the user chooses the colour palette "greys" for the plot, then the default annotation colour changes to "green4".
```
png("md18.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, col = "greys")
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md18.png)

We have the option to provide a column indicating the annotation colour intended for each SNP as a part of the input dataframe. In that case, we just need to specify the name of column in the ```annotationCol``` argument. Let us show this in the example below. In our example, the input dataframe has a separate column named ANNCOL for intended annotation colours for each SNP. Using the column we have annotated the first 4 chromosomes with green, chromosome 5-10 with red and the rest with blue as can be seen in the plot below.
```
png("md17.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, annotationCol="ANNCOL")
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md17.png)

In case we do not want to annotate any particular set of SNPs we can just provide NA values in the ```annotationCol``` column for those SNPs. In our example, we have another column ANNCOL2 in our input dataframe that has NA values corresponding to chromosomes 3,4,5 and 6. The rest of the values of ANNCOL are exactly same as that of the column ANNCOL. We have provided NA values in this column because we do not want to annotate SNPs in chromosomes 3,4,5 and 6. Let us see what the plot looks like when we run the command using this column for ```annoationCol``` argument.
```
png("md19.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, annotationCol="ANNCOL2")
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md19.png)

We can choose to colour only the points above our specified p-value threshold. The rest of the plot will become gray by default.
```
png("md12.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, colAbovePval=TRUE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md12.png)

In the above plot, we can change the colour of the region below our specified threshold as well.
```
png("md13.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-5, colAbovePval=TRUE, col2="greens")
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md13.png)

As stated previously, our package also supports plotting of results from genomes of non-model organisms (often with hundreds of contigs or many scaffolds). This is an incremental feature, as the qqman package does not support direct plotting of results from non-model organisms (Chromosome column for input dataset needs to be numeric for qqman).
```
m=read.table("plink.assoc.fisher",header=TRUE,stringsAsFactors=FALSE,sep= '')
```
This dataset has results for 57,712 SNPs on 760 chromosomes. Let us look at the data.
```
str(m)
```
```
'data.frame':	57712 obs. of  9 variables:
 $ CHR: chr  "NC_015762.1" "NC_015762.1" "NC_015762.1" "NC_015762.1" ...
 $ SNP: chr  "1139:88:-" "1139:45:-" "1139:28:-" "2518:86:+" ...
 $ BP : int  32821 32864 32881 58664 68491 83326 83332 92777 92937 92945 ...
 $ A1 : chr  "T" "C" "G" "T" ...
 $ F_A: num  0 0.5217 0 0.0222 0.5217 ...
 $ F_U: num  0.0454 0.375 0.0454 0.0465 0.3636 ...
 $ A2 : chr  "C" "T" "A" "C" ...
 $ P  : num  0.0551 0.0528 0.0551 0.436 0.0366 ...
 $ OR : num  0 1.818 0 0.466 1.909 ...
```
```
head(m)
```
```
          CHR       SNP    BP A1     F_A     F_U A2       P     OR
1 NC_015762.1 1139:88:- 32821  T 0.00000 0.04545  C 0.05513 0.0000
2 NC_015762.1 1139:45:- 32864  C 0.52170 0.37500  T 0.05276 1.8180
3 NC_015762.1 1139:28:- 32881  G 0.00000 0.04545  A 0.05513 0.0000
4 NC_015762.1 2518:86:+ 58664  T 0.02222 0.04651  C 0.43600 0.4659
5 NC_015762.1 3702:96:- 68491  C 0.52170 0.36360  T 0.03664 1.9090
6 NC_015762.1 4414:33:- 83326  A 0.50000 0.31820  G 0.01556 2.1430
```
```
tail(m)
```
```
                 CHR           SNP  BP A1     F_A     F_U A2      P     OR
57707 NW_003570988.1 18013690:37:- 547  G 0.06818 0.08333  C 0.7782 0.8049
57708 NW_003570988.1 18013681:25:- 550  A 0.00000 0.01190  G 0.4941 0.0000
57709 NW_003570988.1 18013687:31:- 550  A 0.01282 0.01389  G 1.0000 0.9221
57710 NW_003570988.1  18013670:8:- 556  T 0.12820 0.19230  A 0.3830 0.6176
57711 NW_003570988.1  18013666:3:- 557  A 0.05128 0.08571  G 0.5178 0.5766
57712 NW_003570988.1 18013689:26:- 557  A 0.02632 0.01282  G 0.6176 2.0810
```
The dataset contains P-values (column P) and two scores (F_A and F_U). From the dataset, lets generate a basic Manhattan plot for P-values.
```
png("md14.png", width=10, height=6, units="in", res=300)
fastman(m)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md14.png)

Now let us plot the Manhattan plots for one of the scores. We must note that the scores do not need log transformations and so we must specify that while running.
```
png("md15.png", width=10, height=6, units="in", res=300)
fastman(m, p = "F_A", logp = FALSE)
dev.off()
```
![](https://github.com/kaustubhad/fastman/blob/main/plots/md15.png)

The colour schemes available for this package are provided below.

![](https://github.com/kaustubhad/fastman/blob/main/colours/1.matlab.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/2.matlab2.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/3.Set1.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/4.Dark2.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/5.Set2.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/6.DarkSet2.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/7.reds.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/8.greens.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/9.blues.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/10.purples.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/11.greys.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/12.rgbs.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/13.all.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/14.rainbow1.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/15.rainbow2.png)
![](https://github.com/kaustubhad/fastman/blob/main/colours/16.rainbow3.png)


