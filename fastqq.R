fastqq <- function(p, q=1, logtransform=TRUE, compare=FALSE, speedup=TRUE, lambda=TRUE, fix_zero=TRUE, cex=0.6, cex.axis=0.9, xlab, ylab, ...) {

# use: source("fastqq.R")
# example: lambda=fastqq(p)
# on a typical imputed assoc file, 10 million snps reduced to 8K, plotting time reduced from 400s in qqman to 24s.
# if speedup is TRUE, rounds everything to 3 digits
# if lambda is TRUE, the genomic inflation factor lambda is shown on the top-left corner.
# regardless of the display option lambda, its value is returned as function output.
# if fix_zero = TRUE, p=0 cases are converted to minimum p-value/2. otherwise they are removed.

# part 1: check values
if (!is.numeric(p)) { stop("Input must be numeric."); } #  check input

# use valid values only
f=!is.na(p) & !is.na(p) & !is.nan(p) & !is.na(p) & !is.null(p) & is.finite(p);
p=p[f];

if (any(p<0)) { warning("negative p-values found. will be excluded."); }
if (any(p==0)) { warning("some p-values equal to zero. check fix_zero behavior."); }
if (any(p>1)) { warning("some p-values > 1. will be excluded."); }

f=(p <= 1) & (p>= 0);
p=p[f];


# fix p=zero cases, if necessary by setting them minimum
f=(p==0);
if (any(f)) {
if (fix_zero) { mx=min(p[!f]); p[f]=mx; }
else { p=p[!f]; }
}


# part 2: calcualte genomic inflation factor lambda
lmb=qchisq(median(p),1,lower.tail=F)/0.4549364;

# part 3: if compare, repeat the process for the other p-value, else sort and create expected distribution

if (compare) # if user opts to compare two sets of p-values
{
if (!is.numeric(q)) { stop("Input must be numeric."); } #  check input

f=!is.na(q) & !is.na(q) & !is.nan(q) & !is.na(q) & !is.null(q) & is.finite(q);
q=q[f];

if (any(q<0)) { warning("negative p-values found. will be excluded."); }
if (any(q==0)) { warning("some p-values equal to zero. check fix_zero behavior."); }
if (any(q>1)) { warning("some p-values > 1. will be excluded."); }

f=(q <= 1) & (q>= 0);
q=q[f];

f=(q==0);
if (any(f)) {
if (fix_zero) { mx=min(q[!f]); q[f]=mx; }
else { q=q[!f]; }
}

lmb=qchisq(median(p),1,lower.tail=F)/qchisq(median(q),1,lower.tail=F);

if (logtransform)
{
q=-log10(sort(q));
p=-log10(sort(p));
}
else
{
q=sort(q);
p=sort(p);
}
}
else #  sort and create expected distribution
{
if (logtransform)
{
q=-log10(ppoints(length(p)));
p=-log10(sort(p));
}
else
{
q=ppoints(length(p));
p=sort(p);
}
}

# part 4: reduce size for fast plotting
m=as.data.frame(cbind(p,q));

if (speedup) # rounds round everything to 3 digits
{
digs=3; ms=m; ms$p=round(ms$p,digits=digs); ms$q=round(ms$q,digits=digs); f=duplicated(ms); ms=ms[!f,];
# print(paste(nrow(m),"snps reduced to",nrow(ms)));
}
else # full data without any reduction
{ ms=m; }

# remove variables to reduce memory hog before plotting
rm(m,p,q);

# part 5: plot
xlbl=expression(Expected ~ ~-log[10](italic(p)));
ylbl=expression(Observed ~ ~-log[10](italic(p)));
if (!missing(xlab)) { xlbl=xlab; } # if xlab provided, use that instead
if (!missing(ylab)) { ylbl=ylab; } # if ylab provided, use that instead
fac=0.4; xbnd=c(0,max(ms$q)+fac); ybnd=c(0,max(ms$p)+fac); # add a bit of space in top corner

plot(ms$q,ms$p,pch=20,cex=cex,cex.axis=cex.axis,las=1,xaxs="i",yaxs="i",xlim=xbnd,ylim=ybnd,xlab=xlbl,ylab=ylbl,...);
abline(0,1,col = "red");
if (lambda) { text(x=max(ms$q)*0.05,y=max(ms$p),labels=bquote(lambda == .(round(lmb,digits=4))),adj=c(0,1)); } # show labmda unless requested

#end
return(lmb);
}
