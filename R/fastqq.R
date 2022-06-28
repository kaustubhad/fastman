fastqq <- function(p1, p2=NULL, logtransform=TRUE, speedup=TRUE, lambda=TRUE, maxP=14, fix_zero=TRUE, cex=0.6, cex.axis=0.9, xlab, ylab, ...) {

# use: source("fastqq.R")
# example: lambda=fastqq(p1)
# on a typical imputed assoc file, 10 million snps reduced to 8K, plotting time reduced from 400s in qqman to 24s.
# if speedup is TRUE, rounds everything to 3 digits
# if lambda is TRUE, the genomic inflation factor lambda is shown on the top-left corner.
# regardless of the display option lambda, its value is returned as function output.
# if fix_zero = TRUE, p1=0 cases are converted to minimum p-value/2. otherwise they are removed.

# part 1: check values
if (!is.numeric(p1)) { stop("p1 must be numeric."); } #  check input

# use valid values only
f=!is.na(p1) & !is.na(p1) & !is.nan(p1) & !is.na(p1) & !is.null(p1) & is.finite(p1);
p1=p1[f];

if (any(p1<0)) { warning("negative p-values found in p1. will be excluded."); }
if (any(p1==0)) { warning("some p-values in p1 equal to zero. check fix_zero behavior."); }
if (any(p1>1)) { warning("some p-values > 1 in p1. will be excluded."); }

f=(p1 <= 1) & (p1>= 0);
p1=p1[f];


# fix p1=zero cases, if necessary by setting them minimum
f=(p1==0);
if (any(f)) {
if (fix_zero) { mx=min(p1[!f]); p1[f]=mx; }
else { p1=p1[!f]; }
}

# part 2: calcualte genomic inflation factor lambda

lmb=qchisq(median(p1),1,lower.tail=F)/0.4549364;

# part 3: if there are two sets of p-values, repeat the process for the other p-value, else sort and create expected distribution

if (length(p2)>1) # if user opts to compare two sets of p-values
{
  if (length(p2)!=length(p1)) { stop("p1 and p2 must have same length."); } #  check input
  else {
    if (!is.numeric(p2)) { stop("p2 must be numeric."); } #  check input
    
    f=!is.na(p2) & !is.na(p2) & !is.nan(p2) & !is.na(p2) & !is.null(p2) & is.finite(p2);
    p2=p2[f];
    
    if (any(p2<0)) { warning("negative p-values found in p2. will be excluded."); }
    if (any(p2==0)) { warning("some p-values in p2 equal to zero. check fix_zero behavior."); }
    if (any(p2>1)) { warning("some p-values > 1 in p2. will be excluded."); }
    
    f=(p2 <= 1) & (p2>= 0);
    p2=p2[f];
    
    f=(p2==0);
    if (any(f)) {
      if (fix_zero) { mx=min(p2[!f]); p2[f]=mx; }
      else { p2=p2[!f]; }
    }
    
    lmb=qchisq(median(p1),1,lower.tail=F)/qchisq(median(p2),1,lower.tail=F);
    
    if (logtransform)
    {
      p2=-log10(sort(p2));
      p1=-log10(sort(p1));
    }
    else
    {
      p2=sort(p2);
      p1=sort(p1);
    }
  }
  
}
else #  sort and create expected distribution
{
  if (logtransform)
  {
    p2=-log10(ppoints(length(p1)));
    p1=-log10(sort(p1));
  }
  else
  {
    p2=ppoints(length(p1));
    p1=sort(p1);
  }
}

# part 4: reduce size for fast plotting
m=as.data.frame(cbind(p1,p2));

if (speedup) # rounds round everything to 3 digits
{
digs=3; ms=m; ms$p1=round(ms$p,digits=digs); ms$p2=round(ms$p2,digits=digs); f=duplicated(ms); ms=ms[!f,];
# print(paste(nrow(m),"snps reduced to",nrow(ms)));
}
else # full data without any reduction
{ ms=m; }

# remove variables to reduce memory hog before plotting
rm(m,p1,p2);

# part 5: plot
xlbl=expression(Expected ~ ~-log[10](italic(p1)));
ylbl=expression(Observed ~ ~-log[10](italic(p1)));
if (!missing(xlab)) { xlbl=xlab; } # if xlab provided, use that instead
if (!missing(ylab)) { ylbl=ylab; } # if ylab provided, use that instead
fac=0.4; xbnd=c(0,max(ms$p2)+fac); ybnd=c(0,max(ms$p)+fac); # add a bit of space in top corner

plot(ms$p2,ms$p1,pch=20,cex=cex,cex.axis=cex.axis,las=1,xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
abline(0,1,col = "red");
if (lambda) { text(x=max(ms$p2)*0.05,y=max(ms$p1),labels=bquote(lambda == .(round(lmb,digits=4))),adj=c(0,1)); } # show labmda unless requested

#end
return(lmb);
}
