\name{fastman}
\alias{fastman}
\title{Create a Manhattan plot}
\description{
Creates a Manhattan plot directly from a PLINK assoc output (or any data frame with chromosome, position, and p-value). See the package vignette for details:

\code{vignette("fastman")}
}

\usage{
\code{fastman <- function(m, chr = "CHR", bp = "BP", p = "P", snp, chrlabs, speedup=TRUE, logp = TRUE,
col="matlab", maxP=14, sortchr=TRUE, bybp=FALSE, chrsubset, bprange, highlight, annotateHighlight=FALSE,
annotatePval, colAbovePval=FALSE, col2="greys", annotateTop=TRUE, annotationWinMb, annotateN,
annotationCol, annotationAngle=45, baseline=NULL, suggestiveline, genomewideline, cex=0.4, cex.text=0.4,
cex.axis=0.6, xlab, ylab, xlim, ylim, ...)}
}

\examples{
## A basic Manhattan plot
\code{png("md1.png", width=10, height=6, units="in", res=300)
fastman(m)
dev.off()}

## A Manhattan plot annotating SNPs with -logP values higher than 7
\code{png("md5.png", width=10, height=6, units="in", res=300)
fastman(m, annotatePval=1E-7)
dev.off()}
}

\author{
Soumya Subhra Paria (\email{soumya.paria@open.ac.uk}), Sarthok Rasique Rahman (\email{sarthok@gmail.com}) and Kaustubh Adhikari (\email{kaustubh.adhikari@open.ac.uk})
}
