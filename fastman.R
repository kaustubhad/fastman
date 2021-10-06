fastman <- function(m, chr = "CHR", bp = "BP", p = "P", snp, chrlabs, speedup=TRUE, logp = TRUE, col="matlab", maxP=14, sortchr=TRUE, bybp=FALSE, chrsubset, bprange,
                           highlight, annotateHighlight=FALSE, annotatePval, colAbovePval=FALSE, col2="greys", annotateTop=TRUE, annotationWinMb, annotateN, annotationCol, annotationAngle=45, 
                           baseline=NULL, suggestiveline, genomewideline, cex=0.4, cex.text=0.4, cex.axis=0.6, xlab, ylab, xlim, ylim, ...) {

# use: source("fastman.R");
# example: tic(); fastman(m); toc();
# on a typical imputed assoc file, 10 million snps reduced to 170K, plotting time reduced to 737s in qqman to 60s.

# part 0: initialize parameters --------------------------------------------------------------------------------------------------------------------------------------------
if (missing(suggestiveline)) { suggestiveline=NULL; if (logp) { suggestiveline=-log10(1e-05); } } # if suggestiveline is not provided by user then 5 is set as default value
if (missing(genomewideline)) { genomewideline=NULL; if (logp) { genomewideline=-log10(5e-08); } } # if genomewideline is not provided by user then 8 is set as default value

if (missing(annotationCol)) { annotationCol="gray50"; if (col=="greys") { annotationCol="#6bbb00"; } } # if annotationCol is not provided by user then grey is set as default
if (!(annotationCol %in% colnames(m))) { m$annotationCol <- annotationCol; annotationCol <- "annotationCol"; } # if annotationCol is not a column in the data set then create a column with all entries taking the value of annotationCol provided
if (annotationAngle<0|annotationAngle>90) { warning("Text orientation angle should be between 0 & 90 degrees. reverting to default."); annotationAngle=45; } # if annotationAngle is not within 0 and 90 then it is set to 45

if (!missing(annotatePval)) { if (!is.numeric(annotatePval)) { stop("Non-numeric value supplied to annotatePval."); } } # if non-numeric annotatePval is provided by user then the code gives error message

if (!missing(annotationWinMb)) {
  if (!is.numeric(annotationWinMb)) { stop("Non-numeric value supplied to annotationWinMb."); } # if non-numeric annotationWinMb is provided by user then the code gives error message
  if (annotationWinMb<0) { stop("Negative value supplied to annotationWinMb."); } # if negative annotationWinMb is provided by user then the code gives error message
}

if (!missing(annotateN)) {
  if(!is.numeric(annotateN)) { stop("Non-positive value supplied to annotateN."); } # if non-numeric annotationN is provided by user then the code gives error message
  if (annotateN<1) { stop("Too small value supplied to annotateN."); } # if less than 1 input for annotationN is provided by user then the code gives error message
}

if (!missing(annotatePval)&logp) { if (annotatePval<=0) {stop("Non-positive value supplied to annotatePval with log transformation selected."); }; if (annotatePval<1) {annotatePval=-log10(annotatePval); }; } # if non-positive annotatePval is provided by user then the code gives error message

if (!missing(highlight)) { highlight=highlight[!is.na(highlight)]; }

if (nrow(m)==0) { stop("Input data has 0 rows, cannot be plotted!"); } # check if data has actual rows

if (!(chr %in% colnames(m))) { stop(paste("Column", chr, "not found!")); } # check if chr column is provided
if (!(bp %in% colnames(m))) { stop(paste("Column", bp, "not found!")); } # check if bp column is provided
if (!(p %in% colnames(m))) { stop(paste("Column", p, "not found!")); } # check if p column is provided

if (!missing(chrsubset)) { if (!is.numeric(chrsubset)) { stop("Non-numeric value supplied to chrsubset."); } } # if non-numeric chrsubset is provided by user then code gives error message

if (!missing(bprange)) { # range of BP provided, to subset the data
  if (!is.numeric(bprange)) { stop("Supplied BP range contains non-numeric values."); } # if non-numeric bprange is provided by user then code gives error message
  if (length(bprange)!=2) { stop("Supplied BP range must contain 2 values, start and end of range."); } # if both start and end of bprange are not provided by user then code gives error message
  if (bprange[2]<bprange[1]) { stop("Supplied BP range end point cannot be less than the start point."); } # if end point is less than the start point then code gives error message
  bybp=TRUE; # if BP range supplied, then plotting must be by position
}

# check if snp column is necessary
if (!missing(highlight) | !missing(annotatePval) | !missing(annotateN)) { showsnp=TRUE; } else { showsnp=FALSE; } # create an indicator showsnp to check whether snp column is needed
if (showsnp) { # snp column needed
  if (missing(snp)) { snp="SNP"; } # if snp column is not specified by user then the column named SNP in data is selected
  if (!(snp %in% colnames(m))) { stop(paste("Column", snp, "not found!")); } # if user doesn't specify snp column and there is no column called SNP in data then code gives error message
}

zerocount <- 5

# set color palettes; no gray in ggplot palettes
if (length(col)==1) {
  if (col=="matlab") { col=c("#D95319","#E4A100","#7E2F8E","#5EA500","#0095D4","#A2142F","#0C53AA"); }
  else if (col=="matlab2") { col=c("#A2142F","#0095D4","#7E2F8E","#5EA500","#D95319","#0C53AA","#E4A100"); }
  else if (col=="Set1") { col=c("#e41a1c","#377eb8","#f0c001","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf"); }
  else if (col=="Dark2") { col=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d"); }
  else if (col=="Set2") { col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#d1ad79"); }
  else if (col=="DarkSet2") { col=c("#1b9e77","#ffd92f","#d95f02","#a6d854","#7570b3","#d1ad79","#e7298a","#66c2a5","#e6ab02","#8da0cb","#a6761d","#fc8d62","#66a61e","#e78ac3"); }
  else if (col=="reds") { col=c("#e31a1c","#fd8f71"); }
  else if (col=="greens") { col=c("#33a02c","#97d52b"); }
  else if (col=="blues") { col=c("#1f78b4","#40a5eb"); }
  else if (col=="purples") { col=c("#984ea3","#f781dd"); }
  else if (col=="greys") { col=c("gray35","gray60"); }
  else if (col=="rgbs") { col=c("#e31a1c","#fd8f71","#33a02c","#97d52b","#1f78b4","#40a5eb"); }
  else if (col=="all") { col=c("#D95319","#E4A100","#7E2F8E","#5EA500","#0095D4","#A2142F","#97d52b","#0C53AA","#e6ab02","#40a5eb","#d95f02","#7570b3","#e7298a","#66a61e","#fd8f71","#a6761d","#1b9e77"); }
  else if (col=="rainbow1") { col=c("#9b001e","#e31a1c","#ff6726","#ff9249","#c58b00","#f0c001","#97d52b","#099100","#055600","#009476","#00abea","#000cff","#0c53aa","#7570b3","#7e2f8e","#de00ff","#f781bf","#e7298a"); }
  else if (col=="rainbow2") { col=c("#e31a1c","#ff7f00","#f0c001","#5EA500","#40a5eb","#0C53AA","#984ea3"); }
  else if (col=="rainbow3") { col=c("#e31a1c","#40a5eb","#7E2F8E","#ff7f00","#5EA500","#f0c001","#0C53AA"); }
  else { warning("inappropriate col color vector/tag provided. resorting to default"); col=c("#D95319","#E4A100","#7E2F8E","#5EA500","#0095D4","#A2142F","#0C53AA"); } # use matlab as default
}
if (length(col2)==1) {
  if (col2=="greys") { col2=c("gray35","gray60"); }
  else if (col2=="matlab2") { col2=c("#A2142F","#0095D4","#7E2F8E","#5EA500","#D95319","#0C53AA","#E4A100"); }
  else if (col2=="Set1") { col2=c("#e41a1c","#377eb8","#f0c001","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf"); }
  else if (col2=="Dark2") { col2=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d"); }
  else if (col2=="Set2") { col2=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#d1ad79"); }
  else if (col2=="DarkSet2") { col2=c("#1b9e77","#ffd92f","#d95f02","#a6d854","#7570b3","#d1ad79","#e7298a","#66c2a5","#e6ab02","#8da0cb","#a6761d","#fc8d62","#66a61e","#e78ac3"); }
  else if (col2=="reds") { col2=c("#e31a1c","#fd8f71"); }
  else if (col2=="greens") { col2=c("#33a02c","#97d52b"); }
  else if (col2=="blues") { col2=c("#1f78b4","#40a5eb"); }
  else if (col2=="purples") { col2=c("#984ea3","#f781dd"); }
  else if (col2=="matlab") { col2=c("#D95319","#E4A100","#7E2F8E","#5EA500","#0095D4","#A2142F","#0C53AA"); }
  else if (col2=="rgbs") { col2=c("#e31a1c","#fd8f71","#33a02c","#97d52b","#1f78b4","#40a5eb"); }
  else if (col2=="all") { col2=c("#D95319","#E4A100","#7E2F8E","#5EA500","#0095D4","#A2142F","#97d52b","#0C53AA","#e6ab02","#40a5eb","#d95f02","#7570b3","#e7298a","#66a61e","#fd8f71","#a6761d","#1b9e77"); }
  else if (col2=="rainbow1") { col2=c("#9b001e","#e31a1c","#ff6726","#ff9249","#c58b00","#f0c001","#97d52b","#099100","#055600","#009476","#00abea","#000cff","#0c53aa","#7570b3","#7e2f8e","#de00ff","#f781bf","#e7298a"); }
  else if (col2=="rainbow2") { col2=c("#e31a1c","#ff7f00","#f0c001","#5EA500","#40a5eb","#0C53AA","#984ea3"); }
  else if (col2=="rainbow3") { col2=c("#e31a1c","#40a5eb","#7E2F8E","#ff7f00","#5EA500","#f0c001","#0C53AA"); }
  else { warning("inappropriate col2 color vector/tag provided. resorting to default"); col2=c("gray35","gray60"); } # use greys as default
}


# part 1: calculate logP, etc. --------------------------------------------------------------------------------------------------------------------------------------------
if (showsnp) { m=m[,c(chr,bp,p,snp,annotationCol)]; colnames(m)=c("CHR","BP","P","SNP","annotationCol"); } # snp column needed
else { m=m[,c(chr,bp,p)]; colnames(m)=c("CHR","BP","P"); } # only chr, bp and p columns are kept, rest are dropped, snp column is kept only if needed
f=complete.cases(m[,1:3]); m=m[f,];

# check if BP and P columns are numeric
if (!is.numeric(m$BP)) { stop("BP column has non-numeric values"); } # if BP column has non-numeric values then code gives error message
if (!is.numeric(m$P)) { stop("P column has non-numeric values"); } # if P column has non-numeric values then code gives error message

if (!missing(chrsubset)) { f=(m$CHR %in% chrsubset); m=m[f,,drop=F]; } # subset by given chrsubset value(s)

if (nrow(m)==0) { stop("Input data has no usable rows, cannot be plotted!"); }  # check if data has actual usable rows

unc=unique(m$CHR); # list of unique chromosome numbers
if (sortchr) { unc=sort(unc); } # chromosome numbers might not be ordered in the file, so sort by chr first

numc=length(unc); # number of unique chromosomes

if (!missing(bprange)&numc>1) { stop("Only one chromosome can be used with BP range, otherwise it is ambiguous. Pass a single chromosome ID via chrsubset"); }

if (!missing(bprange)) { f=(m$BP>=bprange[1])&(m$BP<=bprange[2]); m=m[f,,drop=F]; } # subset by given BP range

if (nrow(m)==0) { stop("Input data has no usable rows, cannot be plotted!"); }  # check if data has actual usable rows

# if single chromosome, show mb position instead of scaled position
if (numc==1) { bybp=TRUE; }

# if to be sorted by position, no need to sort by chr name/number
if (bybp) { sortchr=FALSE; }

if (sortchr) { f=order(m$CHR); m=m[f,]; } # chromosome numbers might not be ordered in the file, so sort by chr first
# if (bybp) { f=order(m$BP); m=m[f,]; } # sort by bp. is not necessary for plotting, but for nice colors

m$BP=as.double(m$BP)/1E6; # convert to MB position

if (logp) { m$logP=-log10(m$P); } else { m$logP=m$P; } # if log transformation is required then log base 10 of p-values is calculated, else the given column is passed as logP

# if chrlabs provided, check length
if (numc>1&!missing(chrlabs)) { if (length(chrlabs)!=numc) { warning("number of chromosome labels provided do not match number of unique chromosomes. reverting to values in CHR column."); chrlabs=unc; } }
else { chrlabs=unc; }
  

  
# part 2: prepare plotting positions etc. --------------------------------------------------------------------------------------------------------------------------------------------
if (bybp) { # show mb position instead of scaled position. scale positions to make suitable for rounding. applies if single chromosome
  fac1c=23/(max(m$BP)-min(m$BP)); m$BPn=m$BP*fac1c; # calculating the scaling factor
  if (!missing(annotationWinMb)) { annotationWinMb=annotationWinMb*fac1c; } # scaling the mb window
}
else { # typical plot, combine positions across chromosomes
  cmat=matrix(NA,nrow=numc,ncol=4); # a matrix cmat is created with number of rows as number of unique chr and 3 columns
  for (i in 1:numc) { ch=unc[i]; f=m$CHR==ch; w=m$BP[f]; cmat[i,1]=min(w); cmat[i,2]=max(w); cmat[i,3]=cmat[i,2]-cmat[i,1]; m$BP[f]=w-cmat[i,1]; cmat[i,4]=median(diff(sort(w))); }
  cmat=as.data.frame(cmat); colnames(cmat)=c("min","max","width","medgap"); # each row of cmat corresponds to each chr, 4 columns of cmat are calculated as min BP, max BP, width, and median gap
  maxgap=max(cmat$medgap,na.rm=TRUE);
  cmat$base=0*cmat$min; cmat$midp=0*cmat$min;
  i=1; cmat$base[i]=0; cmat$midp[i]=cmat$width[i]/2; # base and midpoint are calculated for each chr with respect to width calculated previously
  if (numc>1) { for (i in 2:numc) { cmat$base[i]=cmat$base[i-1]+cmat$width[i-1]+maxgap; cmat$midp[i]=cmat$base[i]+(cmat$width[i]/2); } }
  fac=numc/cmat$midp[numc]; # calculating the scaling factor
  m$BP=fac*m$BP; cmat$basef=fac*cmat$base; cmat$midpf=fac*cmat$midp;
  m$BPn=m$BP;
  if (numc>1) { for (i in 2:numc) { ch=unc[i]; f=m$CHR==ch; w=m$BP[f]; m$BPn[f]=w+cmat$basef[i]; } }
  if (!missing(annotationWinMb)) { annotationWinMb=annotationWinMb*fac; } # scaling the mb window
}
m$C=m$CHR;
for (i in 1:numc) { ch=unc[i]; f=m$CHR==ch; m$C[f]=i; }




# part 2 B: annotations -----------------------------------
if (showsnp) { # pick selected snps, using either of 3 flags
  m=m[,c("C","BPn","logP","SNP","annotationCol")];
  if (!missing(highlight)) { # highlight listed snps
    f=(m$SNP %in% highlight); msnp=m[f,,drop=F]; # msnp will contain the SNP information for annotation, in this case msnp contains information of only SNPs to be highlighted
	if (annotateHighlight&annotateTop) { # annotate the top highlighted SNP per chromosome
	  m1=NULL; sunc=unique(msnp$C);
	  for (i in sunc) { m2=msnp[msnp$C==i,]; f=m2$logP==max(m2$logP); m1=rbind(m1,m2[f,]); } # loop over chr and identify the max logP value for every chr
	  msnp=m1; # SNP information of only max logP value per chr is stored in msnp
	  rm(m1,m2);
	}
  }
  else if (!missing(annotatePval)) { # highlight snps with pvalue beyond a threshold; change based on logp=TRUE or not
    f=(m$logP>=annotatePval); zerocount=sum(f); msnp=m[f,,drop=F]; # msnp will contain the SNP information for annotation, in this case msnp contains information of SNPs with p-value beyond threshold
  }
  else if (!missing(annotateN)) { # highlight top N snps
    if (!missing(annotationWinMb)) { # within each chr, only one snp within annotationWinMb window. annotationWinMb value adjusted by factor above
	  m0=NULL; sunc=unique(msnp$C);
	  for (i in sunc) { # loop over chr
	    m2=msnp[msnp$C==i,]; f=order(m2$logP,decreasing=TRUE); m2=m2[f,]; m1=m2[1,,drop=F];
		if (nrow(m2)>1) { for (j in 2:nrow(m2)) { f=abs(m1$BPn-m2$BPn[j]); if (min(f)>annotationWinMb) { m1=rbind(m1,m2[j,,drop=F]); } } }
		m0=rbind(m0,m1); # identify top SNP within annotationWinMb window for each chr
	  }
	  msnp=m0; # store information of only the top SNP within annotationWinMb window for each chr
	  rm(m0,m1,m2);
	  k=nrow(msnp)-annotateN+1; k=sort(msnp$logP,partial=k)[k]; # sorting the p-value column to identify the k-th highest value
	  f=(msnp$logP>=k); msnp=msnp[f,,drop=F]; # msnp will contain the SNP information for annotation, in this case msnp contains information of top N SNPs
    }
	else { # only highlight top N snps
	  k=nrow(m)-annotateN+1; k=sort(m$logP,partial=k)[k]; # sorting the p-value column to identify the k-th highest value
	  f=(m$logP>=k); msnp=m[f,,drop=F]; # msnp will contain the SNP information for annotation, in this case msnp contains information of top N SNPs
	}
  }
}

if (!missing(annotationWinMb)) { annotateTop=FALSE; } # if annotating by window, then top snp flag not applicable

if (missing(highlight)&!missing(annotatePval)&zerocount>0) { # not using highlight list, but just annotatePval instead. check annotateTop and annotationWinMb options
  if (annotateTop) { # if annotateTop = TRUE (default), only top snps in each chr is shown. else all top snps are shown
    m1=NULL; sunc=unique(msnp$C);
    for (i in sunc) { m2=msnp[msnp$C==i,]; f=m2$logP==max(m2$logP); m1=rbind(m1,m2[f,]); } # loop over chr and identify the max logP value for every chr
    msnp=m1; # SNP information of only max logP value per chr is stored in msnp
    rm(m1,m2);
  }
  else if (!missing(annotationWinMb)) { # within each chr, only one snp within annotationWinMb window. annotationWinMb value adjusted by factor above
    m0=NULL; sunc=unique(msnp$C);
    for (i in sunc) { # loop over chr
      m2=msnp[msnp$C==i,]; f=order(m2$logP,decreasing=TRUE); m2=m2[f,]; m1=m2[1,,drop=F];
      if (nrow(m2)>1) { for (j in 2:nrow(m2)) { f=abs(m1$BPn-m2$BPn[j]); if (min(f)>annotationWinMb) { m1=rbind(m1,m2[j,,drop=F]); } } }
      m0=rbind(m0,m1); # identify top SNP within annotationWinMb window for each chr
    }
    msnp=m0; # store information of only the top SNP within annotationWinMb window for each chr
    rm(m0,m1,m2);
  }
}

if (!is.null(maxP)) { # if user has not specified that there should not be any maxP truncation
  if (showsnp) { # if annotation is required
    f=msnp$logP>=maxP; msnp$logP[f]=maxP; # logP values above maxP are truncated to maxP in msnp
	f=msnp$logP<=-maxP; msnp$logP[f]=-maxP; # logP values below - maxP are truncated to -maxP in msnp
  }
  f=m$logP<=-maxP; m$logP[f]=-maxP; # logP values below -maxP are truncated to -maxP in m
  f=m$logP>=maxP; m$logP[f]=maxP; # logP values above maxP are truncated to maxP in m
}
else { # if user has specified that there should not be any maxP truncation
  if (showsnp) { # if annotation is required
    f=msnp$logP==Inf; msnp$logP[f]=sort(msnp$logP)[length(msnp$logP)-1]; # Inf values are truncated to highest finite value
	f=msnp$logP==-Inf; msnp$logP[f]=sort(msnp$logP)[2]; # -Inf values are truncated to lowest finite value
  }
  f=m$logP==Inf; m$logP[f]=sort(m$logP)[length(m$logP)-1]; # Inf values are truncated to highest finite value
  f=m$logP==-Inf; m$logP[f]=sort(m$logP)[2]; # -Inf values are truncated to lowest finite value
}

if (showsnp) { # if annotation is required
  msnp$SNP[is.na(msnp$annotationCol)] <- NA; # if user has specified NA as colour of a SNP then the annotation text is marked as NA
}


# part 2 C: finally keep 3 columns only, as msnp with snp names is separate now -----------------------------------
m=m[,c("C","BPn","logP")];

# temporarily scale y axis to make more suitable for plotting
facy=10/(max(m$logP)-min(m$logP)); m$logP=m$logP*facy;



# part 3: reduce size for fast plotting --------------------------------------------------------------------------------------------------------------------------------------------
if (speedup) { # fast method: below 0.2% and above 99.8% round to 3 digits, rest round to 2 digits
  quants=c(0,0.002,0.5,0.998,1); quants=quantile(m$logP,quants); # minimum, 0.2th percentile, median, 99.8th percentile and maximum are being calculated
  right=(quants[5] - quants[4])/(quants[4] - quants[3]); # measure of significance of right tail
  left=(quants[1] - quants[2])/(quants[2] - quants[3]); # measure of significance of left tail
  if (nrow(m)<1E5) { #if there are lest than 100k rows then full data is rounded to 3 digits
    digs=3; ms=m; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,];
  }
  else { # round lower and upper parts separately
    if (right>0.1) { # significant right tail
      if (left>0.1) { # significant left tail
	    outer=m$logP>quants[4]|m$logP<quants[2]; inner=!outer;
		digs=2; ms=m[inner,]; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,]; # inner part is rounded to 2 digits
		digs=3; m1=m[outer,]; m1$logP=round(m1$logP,digits=digs); m1$BPn=round(m1$BPn,digits=digs); f=duplicated(m1); m1=m1[!f,]; ms=rbind(ms,m1); # outer part is rounded to 3 digits
	  }
	  else { # insignificant left tail
	    outer=m$logP>quants[4]; inner=!outer;
		digs=2; ms=m[inner,]; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,]; # inner part is rounded to 2 digits
		digs=3; m1=m[outer,]; m1$logP=round(m1$logP,digits=digs); m1$BPn=round(m1$BPn,digits=digs); f=duplicated(m1); m1=m1[!f,]; ms=rbind(ms,m1); # outer part is rounded to 3 digits
	  }
    }
    else { # insignificant right tail
      if (left>0.1) { # significant left tail
	    outer=m$logP<quants[2]; inner=!outer;
		digs=2; ms=m[inner,]; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,]; # inner part is rounded to 2 digits
		digs=3; m1=m[outer,]; m1$logP=round(m1$logP,digits=digs); m1$BPn=round(m1$BPn,digits=digs); f=duplicated(m1); m1=m1[!f,]; ms=rbind(ms,m1); # outer part is rounded to 3 digits
	  }
	  else { # insignificant left tail
	    digs=3; ms=m; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,]; # as there is no significant tail full data is rounded to 3 digits
	  }
    }
    rm(m1);
  }
}
else # full data without any reduction
{ ms=m; }

# print(paste(nrow(m),"snps reduced to",nrow(ms)));

# remove variables to reduce memory hog before plotting
rm(m);

# if showing mb position instead of scaled position, rescale positions after rounding (e.g. if single chromosome)
if (bybp) {
  ms$BPn=ms$BPn/fac1c;
  if (showsnp) { msnp$BPn=msnp$BPn/fac1c; }
}

ms$logP=ms$logP/facy;



# part 4: plot --------------------------------------------------------------------------------------------------------------------------------------------
palette(col); # setting the colour palette

ybnd=c(floor(min(c(min(ms$logP),0))),ceiling(max(ms$logP))); # setting the y axis boundaries

# set x axis boundaries, depending on whether single chromosome or not
fac=0.015*(max(ms$BPn)-min(ms$BPn));
if (numc==1) { xbnd=c(min(ms$BPn)-fac,max(ms$BPn)+fac); } else { xbnd=c(-fac,max(ms$BPn)+fac); }

if (showsnp&zerocount>0) { # if annotation is required
	parmai=par("mai"); # measuring the margin size in inches
	plotx=par("pin")[1]; ploty=par("pin")[2]; # storing the plot width in inches within variable plotx and plot height in inches within variable ploty
	xwidth=xbnd[2]-xbnd[1]; ywidth=ybnd[2]-ybnd[1]; # calculating the inherent width and height of the plot

	m2=msnp;
	strgap=strheight("A",units="inches",cex=cex.text);
	m2$strwidth=strwidth(m2$SNP,units="inches",cex=cex.text)+strheight(m2$SNP,units="inches",cex=cex.text); # measuring the dimensions of the SNP texts in inches
	m2$width=m2$strwidth*cospi(annotationAngle/180); # measuring the x component of the SNP text in inches
	m2$height=m2$strwidth*sinpi(annotationAngle/180); # measuring the y component of the SNP text in inches
	m2$xmax=m2$BPn+fac/7+m2$width*xwidth/plotx; # measuring the maximum width required for SNP texts in plot X units
	m2$ymax=m2$logP+fac/7+m2$height*ywidth/ploty; # measuring the maximum height required for SNP text in plot Y units

	extrax=max( (max(m2$xmax)-xbnd[2])*plotx/xwidth - 0.8*parmai[4] , 0); # measuring the extra width required from margin in inches
	extray=max( (max(m2$ymax)-ybnd[2])*ploty/ywidth - 0.8*parmai[3] , 0); # measuring the extra height required from margin in inches

	parmai[3]=parmai[3]+extray; parmai[4]=parmai[4]+extrax; # adding the extra width and height
	par(mai=parmai); # setting the final margin
}

if (!missing(xlim)) { xbnd=xlim; }
if (!missing(ylim)) { ybnd=ylim; }

xlbl="Chromosome";
if (numc==1) { xlbl=paste(xlbl,unc[1],"(Mb)",sep=" "); } # if single chromosome, add chromosome number to X label
if (bybp&numc>1) { xlbl="Position (Mb)"; }
if (logp) { ylbl=expression(-log[10](italic(p))); } else { ylbl=p; } # adapt Y label to whether log transformed or not
if (!missing(xlab)) { xlbl=xlab; }
if (!missing(ylab)) { ylbl=ylab; }

par(xpd = NA)
if (colAbovePval) { # Colour all hits above p-value threshold, while points below are col2
  if (bybp) { # show mb position instead of scaled position
    palette(col2); # setting colour palette to col2 for hits below threshold
	plot(ms$BPn,ms$logP,pch=20,col=ms$C,cex=cex,cex.axis=cex.axis,las=1,xaxt="n",bty="n",xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
	axis(side=1,pos=0,las=1,cex.axis=cex.axis,col=NA,col.ticks="black");
	f=ms$logP>=annotatePval; ms=ms[f,]; # selecting hits above threshold
	palette(col); # setting colour palette to user defined / default for hits above threshold
	points(ms$BPn,ms$logP,col=ms$C,pch=20,cex=cex); # plotting hits above threshold
  }
  else { # typical plot
    palette(col2) # setting colour palette to col2 for hits below threshold
	plot(ms$BPn,ms$logP,pch=20,col=ms$C,cex=cex,cex.axis=cex.axis,las=1,xaxt="n",bty="n",xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
	axis(side=1,at=cmat$midpf,labels=chrlabs,pos=0,las=1,cex.axis=cex.axis,col=NA,col.ticks="black");
	f=ms$logP>=annotatePval; ms=ms[f,]; # selecting hits above threshold
	palette(col); # setting colour palette to user defined / default for hits above threshold
	points(ms$BPn,ms$logP,col=ms$C,pch=20,cex=cex); # plotting hits above threshold
  }
}
else { # one colour scheme for the entire plot
  if (bybp) { # show mb position instead of scaled position
    plot(ms$BPn,ms$logP,pch=20,col=ms$C,cex=cex,cex.axis=cex.axis,las=1,xaxt="n",bty="n",xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
	axis(side=1,pos=0,las=1,cex.axis=cex.axis,col=NA,col.ticks="black");
  }
  else { # typical plot
    plot(ms$BPn,ms$logP,pch=20,col=ms$C,cex=cex,cex.axis=cex.axis,las=1,xaxt="n",bty="n",xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
	axis(side=1,at=cmat$midpf,labels=chrlabs,pos=0,las=1,cex.axis=cex.axis,col=NA,col.ticks="black");
  }
}

if (is.numeric(baseline)) { abline(h=baseline,col="black", xpd=FALSE); } # plotting baseline
if (is.numeric(suggestiveline)) { abline(h=suggestiveline,col="blue", xpd=FALSE); } # plotting suggestiveline
if (is.numeric(genomewideline)) { abline(h=genomewideline,col="red", xpd=FALSE); } # plotting genomewideline


# part 5: plot highlights / annotations if necessary --------------------------------------------------------------------------------------------------------------------------------------------

adj=c(0,0);
if (annotationAngle==0) { adj=c(0,0.5); }

if (!missing(highlight)) { # show highlighted snps
  points(msnp$BPn,msnp$logP,col=msnp$annotationCol,pch=20,cex=cex); # plotting highlighted points
  if (annotateHighlight) { # annotate some/all of listed snps
    if (!missing(annotatePval)) { f=msnp$logP>=annotatePval; if (zerocount>0) { m2=msnp[f,]; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=m2$annotationCol,cex=cex.text,...); } } # only snps above threshold
    else { m2=msnp; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=m2$annotationCol,cex=cex.text,...); } # all snps
  }
}

if (!missing(annotateN)) { # not using highlight list, but annotatePval or annotateN instead
  m2=msnp; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=m2$annotationCol,cex=cex.text,...); # annotating the selected snps
}

if (!missing(annotatePval)&(zerocount>0)) { # not using highlight list, but annotatePval or annotateN instead
  m2=msnp; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=m2$annotationCol,cex=cex.text,...); # annotating the selected snps
}

}