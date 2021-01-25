fastman <- function(m, chr = "CHR", bp = "BP", p = "P", snp, chrlabs, speedup=TRUE, logp = TRUE, col="matlab", maxP, sortchr=TRUE, bybp=FALSE, chrsubset, bprange,
                           highlight, annotateHighlight=FALSE, annotatePval, annotateTop=TRUE, annotationWinMb, annotateN, annotationCol, annotationAngle=45, 
                           suggestiveline, genomewideline, cex=0.4, cex.axis=0.6, xlab, ylab, xlim, ylim, ...) {

# use: source("fastman.R");
# example: tic(); fastman(m); toc();
# on a typical imputed assoc file, 10 million snps reduced to 170K, plotting time reduced to 737s in qqman to 60s.

# part 0: initialize parameters --------------------------------------------------------------------------------------------------------------------------------------------
if (missing(maxP)) { maxP=NULL; if (logp) { maxP=14; } }
if (missing(suggestiveline)) { suggestiveline=NULL; if (logp) { suggestiveline=-log10(1e-05); } }
if (missing(genomewideline)) { genomewideline=NULL; if (logp) { genomewideline=-log10(5e-08); } }

if (missing(annotationCol)) { annotationCol="gray50"; if (col=="greys") { annotationCol="#6bbb00"; } }
if (annotationAngle<0|annotationAngle>90) { warning("Text orientation angle should be between 0 & 90 degrees. reverting to default."); annotationAngle=45; }

if (!missing(annotatePval)) { if (!is.numeric(annotatePval)) { stop("Non-numeric value supplied to annotatePval."); } }

if (!missing(annotationWinMb)) {
  if (!is.numeric(annotationWinMb)) { stop("Non-numeric value supplied to annotationWinMb."); }
  if (annotationWinMb<0) { stop("Negative value supplied to annotationWinMb."); }
}

if (!missing(annotateN)) {
  if(!is.numeric(annotateN)) { stop("Non-positive value supplied to annotateN."); }
  if (annotateN<1) { stop("Too small value supplied to annotateN."); }
}

if (!missing(annotatePval)&logp) { if (annotatePval<=0) {stop("Non-positive value supplied to annotatePval with log transformation selected."); }; annotatePval=-log10(annotatePval); }
if (!missing(highlight)) { highlight=highlight[!is.na(highlight)]; }

if (nrow(m)==0) { stop("Input data has 0 rows, cannot be plotted!"); } # check if data has actual rows

if (!(chr %in% colnames(m))) { stop(paste("Column", chr, "not found!")); }
if (!(bp %in% colnames(m))) { stop(paste("Column", bp, "not found!")); }
if (!(p %in% colnames(m))) { stop(paste("Column", p, "not found!")); }

if (!missing(chrsubset)) { if (!is.numeric(chrsubset)) { stop("Non-numeric value supplied to chrsubset."); } } # subset of chr values provided, to subset the data

if (!missing(bprange)) { # range of BP provided, to subset the data
  if (!is.numeric(bprange)) { stop("Supplied BP range contains non-numeric values."); }
  if (length(bprange)!=2) { stop("Supplied BP range must contain 2 values, start and end of range."); }
  if (bprange[2]<bprange[1]) { stop("Supplied BP range end point cannot be less than the start point."); }
  bybp=TRUE; # if BP range supplied, then plotting must be by position
}

# check if snp column is necessary
if (!missing(highlight) | !missing(annotatePval) | !missing(annotateN)) { showsnp=TRUE; } else { showsnp=FALSE; }
if (showsnp) { # snp column needed
  if (missing(snp)) { snp="SNP"; }
  if (!(snp %in% colnames(m))) { stop(paste("Column", snp, "not found!")); }
}


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



# part 1: calculate logP, etc. --------------------------------------------------------------------------------------------------------------------------------------------
if (showsnp) { m=m[,c(chr,bp,p,snp)]; colnames(m)=c("CHR","BP","P","SNP"); } # snp column needed
else { m=m[,c(chr,bp,p)]; colnames(m)=c("CHR","BP","P"); }
f=complete.cases(m[,1:3]); m=m[f,];

# check if BP and P columns are numeric
if (!is.numeric(m$BP)) { stop("BP column has non-numeric values"); }
if (!is.numeric(m$P)) { stop("P column has non-numeric values"); }

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

if (logp) { m$logP=-log10(m$P); } else { m$logP=m$P; }

# if chrlabs provided, check length
if (numc>1&!missing(chrlabs)) { if (length(chrlabs)!=numc) { warning("number of chromosome labels provided do not match number of unique chromosomes. reverting to values in CHR column."); chrlabs=unc; } }
else { chrlabs=unc; }
  

  
# part 2: prepare plotting positions etc. --------------------------------------------------------------------------------------------------------------------------------------------
if (bybp) { # show mb position instead of scaled position. scale positions to make suitable for rounding. applies if single chromosome
  fac1c=23/(max(m$BP)-min(m$BP)); m$BPn=m$BP*fac1c;
  if (!missing(annotationWinMb)) { annotationWinMb=annotationWinMb*fac1c; }
}
else { # typical plot, combine positions across chromosomes
  cmat=matrix(NA,nrow=numc,ncol=3);
  for (i in 1:numc) { ch=unc[i]; f=m$CHR==ch; w=m$BP[f]; cmat[i,1]=min(w); cmat[i,2]=max(w); cmat[i,3]=cmat[i,2]-cmat[i,1]; m$BP[f]=w-cmat[i,1]; }
  cmat=as.data.frame(cmat); colnames(cmat)=c("min","max","width");
  cmat$base=0*cmat$min; cmat$midp=0*cmat$min;
  i=1; cmat$base[i]=0; cmat$midp[i]=cmat$width[i]/2;
  if (numc>1) { for (i in 2:numc) { cmat$base[i]=cmat$base[i-1]+cmat$width[i-1]; cmat$midp[i]=cmat$base[i]+(cmat$width[i]/2); } }
  fac=numc/cmat$midp[numc];
  m$BP=fac*m$BP; cmat$basef=fac*cmat$base; cmat$midpf=fac*cmat$midp;
  m$BPn=m$BP;
  if (numc>1) { for (i in 2:numc) { ch=unc[i]; f=m$CHR==ch; w=m$BP[f]; m$BPn[f]=w+cmat$basef[i]; } }
  if (!missing(annotationWinMb)) { annotationWinMb=annotationWinMb*fac; }
}
m$C=m$CHR;
for (i in 1:numc) { ch=unc[i]; f=m$CHR==ch; m$C[f]=i; }

if (!is.null(maxP)) {
  maxP=maxP-0.1;
  f=m$logP>=maxP; m$logP[f]=maxP;
}



# part 2 B: annotations -----------------------------------
if (showsnp) { # pick selected snps, using either of 3 flags
  m=m[,c("C","BPn","logP","SNP")];
  if (!missing(highlight)) { # highlight listed snps
    f=(m$SNP %in% highlight); msnp=m[f,,drop=F];
  }
  else if (!missing(annotatePval)) { # highlight snps with pvalue beyond a threshold; change based on logp=TRUE or not
    f=(m$logP>=annotatePval); msnp=m[f,,drop=F];
  }
  else if (!missing(annotateN)) { # highlight top N snps
    k=nrow(m)-annotateN+1; k=sort(m$logP,partial=k)[k];
    f=(m$logP>=k); msnp=m[f,,drop=F];
  }
}

if (!missing(annotationWinMb)) { annotateTop=FALSE; } # if annotating by window, then top snp flag not applicable

if (missing(highlight)&!missing(annotatePval)) { # not using highlight list, but just annotatePval instead. check annotateTop and annotationWinMb options
  if (annotateTop) { # if annotateTop = TRUE (default), only top snps in each chr is shown. else all top snps are shown
    m1=NULL; sunc=unique(msnp$C);
    for (i in sunc) { m2=msnp[msnp$C==i,]; f=m2$logP==max(m2$logP); m1=rbind(m1,m2[f,]); }
    msnp=m1;
    rm(m1,m2);
  }
  else if (!missing(annotationWinMb)) { # within each chr, only one snp within annotationWinMb window. annotationWinMb value adjusted by factor above
    m0=NULL; sunc=unique(msnp$C);
    for (i in sunc) {
      m2=msnp[msnp$C==i,]; f=order(m2$logP,decreasing=TRUE); m2=m2[f,]; m1=m2[1,,drop=F];
      if (nrow(m2)>1) { for (j in 2:nrow(m2)) { f=abs(m1$BPn-m2$BPn[j]); if (min(f)>annotationWinMb) { m1=rbind(m1,m2[j,,drop=F]); } } }
      m0=rbind(m0,m1);
    }
    msnp=m0;
    rm(m0,m1,m2);
  }
}


# part 2 C: finally keep 3 columns only, as msnp with snp names is separate now -----------------------------------
m=m[,c("C","BPn","logP")];

# temporarily scale y axis to make more suitable for plotting
facy=10/(max(m$logP)-min(m$logP)); m$logP=m$logP*facy;



# part 3: reduce size for fast plotting --------------------------------------------------------------------------------------------------------------------------------------------
if (speedup) { # fast method: below 99.8% round to 2 digits, 99.8% above round to 3 digits
  quants=0.998; quants=quantile(m$logP,quants);
  f1=m$logP<=quants; f2=!f1;
  gap=max(m$logP)-min(m$logP);
  if (nrow(m)<1E5|sum(f2)==0) { # full data to be rounded to 3 digits
    digs=3; ms=m; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,];
  }
  else { # round lower and upper parts separately
    digs=2; ms=m[f1,]; ms$logP=round(ms$logP,digits=digs); ms$BPn=round(ms$BPn,digits=digs); f=duplicated(ms); ms=ms[!f,];
    digs=3; m1=m[f2,]; m1$logP=round(m1$logP,digits=digs); m1$BPn=round(m1$BPn,digits=digs); f=duplicated(m1); m1=m1[!f,]; ms=rbind(ms,m1);
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
palette(col);

ybnd=c(floor(min(c(min(ms$logP),0))),ceiling(max(ms$logP)));

# set x axis boundaries, depending on whether single chromosome or not
fac=0.015*(max(ms$BPn)-min(ms$BPn));
if (numc==1) { xbnd=c(min(ms$BPn)-fac,max(ms$BPn)+fac); } else { xbnd=c(-fac,max(ms$BPn)+fac); }

if (!missing(xlim)) { xbnd=xlim; }
if (!missing(ylim)) { ybnd=ylim; }

if (annotateHighlight|!missing(annotatePval)|!missing(annotateN)) { xbnd[2]=xbnd[2]+(xbnd[2]-xbnd[1])*0.1*cex*cos(annotationAngle*pi/180); }
if (annotateHighlight|!missing(annotatePval)|!missing(annotateN)) { ybnd[2]=ybnd[2]+(ybnd[2]-ybnd[1])*0.3*cex*sin(annotationAngle*pi/180); }

xlbl="Chromosome";
if (numc==1) { xlbl=paste(xlbl,unc[1],"(Mb)",sep=" "); } # if single chromosome, add chromosome number to X label
if (bybp&numc>1) { xlbl="Position (Mb)"; }
if (logp) { ylbl=expression(-log[10](italic(p))); } else { ylbl=expression(italic(P)); } # adapt Y label to whether log transformed or not
if (!missing(xlab)) { xlbl=xlab; }
if (!missing(ylab)) { ylbl=ylab; }

if (bybp) { # show mb position instead of scaled position
  plot(ms$BPn,ms$logP,pch=20,col=ms$C,cex=cex,cex.axis=cex.axis,las=1,xaxt="n",bty="n",xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
  axis(side=1,pos=0,las=1,cex.axis=cex.axis,col=NA,col.ticks="black");
}
else { # typical plot
  plot(ms$BPn,ms$logP,pch=20,col=ms$C,cex=cex,cex.axis=cex.axis,las=1,xaxt="n",bty="n",xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
  axis(side=1,at=cmat$midpf,labels=chrlabs,pos=0,las=1,cex.axis=cex.axis,col=NA,col.ticks="black");
}

if (is.numeric(suggestiveline)) { abline(h=suggestiveline,col="blue"); }
if (is.numeric(genomewideline)) { abline(h=genomewideline,col="red"); }



# part 5: plot highlights / annotations if necessary --------------------------------------------------------------------------------------------------------------------------------------------

adj=c(0,0);
if (annotationAngle==0) { adj=c(0,0.5); }
# if (annotationAngle==90) { adj=c(0.5,0); }

if (!missing(highlight)) { # show highlighted snps
  points(msnp$BPn,msnp$logP,col=annotationCol,pch=20,cex=cex);
  if (annotateHighlight) { # annotate some/all of listed snps
    if (!missing(annotatePval)) { f=msnp$logP>=annotatePval; m2=msnp[f,]; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=annotationCol,cex=cex); } # only snps above threshold
    else if (annotateTop) { f=msnp$logP==max(msnp$logP); m2=msnp[f,]; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=annotationCol,cex=cex); } # only top snp
    else { m2=msnp; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=annotationCol,cex=cex); } # all snps
  }
}

if (showsnp&missing(highlight)) { # not using highlight list, but annotatePval or annotateN instead
  points(msnp$BPn,msnp$logP,col=annotationCol,pch=20,cex=cex*0.7);
  m2=msnp; text(x=m2$BPn+fac/7,y=m2$logP+fac/7,labels=m2$SNP,adj=adj,srt=annotationAngle,col=annotationCol,cex=cex);
}


}