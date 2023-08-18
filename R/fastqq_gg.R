fastqq_gg <- function(p1, p2=NULL, colour, logtransform=TRUE, pairwisecompare=TRUE, speedup=TRUE, lambda=TRUE, maxP=14, fix_zero=TRUE, size=0.8, cex.axis=0.9, xlab, ylab, ...) {
  
  # use: source("fastqq.R")
  # example: lambda=fastqq(p1)
  # on a typical imputed assoc file, 10 million snps reduced to 8K, plotting time reduced from 400s in qqman to 24s.
  # if speedup is TRUE, rounds everything to 3 digits
  # if lambda is TRUE, the genomic inflation factor lambda is shown on the top-left corner.
  # regardless of the display option lambda, its value is returned as function output.
  # if fix_zero = TRUE, p1=0 cases are converted to minimum p-value/2. otherwise they are removed.
  
  # part 1: check values
  if (!is.numeric(p1)) { stop("p1 must be numeric."); } #  check input
  
  if (missing(colour)) { colour='black' ;}
  
  # use valid values only
  f=!is.na(p1) & !is.na(p1) & !is.nan(p1) & !is.na(p1) & !is.null(p1) & is.finite(p1);
  p1=p1[f];
  
  if (logtransform) {
    if (any(p1<0)) { warning("negative p-values found in p1. will be truncated to 0 first and then converted to minimum p-value/2."); }
    if (any(p1==0)) { warning("some p-values in p1 equal to zero. check fix_zero behavior."); }
    if (any(p1>1)) { warning("some p-values > 1 in p1. will be truncated to 1."); }
    
    f=p1 < 0;
    p1[f]=0;
    f=p1 > 1;
    p1[f]=1;
    
    # fix p1=zero cases, if necessary by setting them minimum
    f=(p1==0);
    if (any(f)) {
      if (fix_zero) { mx=min(p1[!f]); p1[f]=mx; }
      else { p1=p1[!f]; }
    }
    
  }
  
  
  # part 2: calculate genomic inflation factor lambda
  
  lmb=qchisq(median(p1),1,lower.tail=F)/0.4549364;
  
  # part 3: if there are two sets of p-values, repeat the process for the other p-value, else sort and create expected distribution
  
  if (length(p2)>1) # if user opts to compare two sets of p-values
  {
    if (length(p2)!=length(p1)) { stop("p1 and p2 must have same length."); } #  check input
    else {
      if (!is.numeric(p2)) { stop("p2 must be numeric."); } #  check input
      
      f=!is.na(p2) & !is.na(p2) & !is.nan(p2) & !is.na(p2) & !is.null(p2) & is.finite(p2);
      p2=p2[f];
      
      if (logtransform) {
        if (any(p2<0)) { warning("negative p-values found in p1. will be truncated to 0 first and then converted to minimum p-value/2."); }
        if (any(p2==0)) { warning("some p-values in p2 equal to zero. check fix_zero behavior."); }
        if (any(p2>1)) { warning("some p-values > 1 in p1. will be truncated to 1."); }
        
        f=(p2 < 0);
        p2[f]=0;
        f=(p2 > 1);
        p2[f]=1;
        
        f=(p2==0);
        if (any(f)) {
          if (fix_zero) { mx=min(p2[!f]); p2[f]=mx; }
          else { p2=p2[!f]; }
        }
      }
      
      
      lmb=qchisq(median(p1),1,lower.tail=F)/qchisq(median(p2),1,lower.tail=F);
      
      if (pairwisecompare)
      {
        if (logtransform)
        {
          p2=-log10(p2);
          p1=-log10(p1);
          xlbl=expression(-log[10](italic(p2)));
          ylbl=expression(-log[10](italic(p1)));
        }
        else
        {
          xlbl=expression(italic(p2));
          ylbl=expression(italic(p1));
        }
      }
      else
      {
        if (logtransform)
        {
          p2=-log10(sort(p2));
          p1=-log10(sort(p1));
          xlbl=expression(-log[10](italic(p2)));
          ylbl=expression(-log[10](italic(p1)));
        }
        else
        {
          p2=sort(p2);
          p1=sort(p1);
          xlbl=expression(italic(p2));
          ylbl=expression(italic(p1));
        }
      }
      
    }
    
  }
  else #  sort and create expected distribution
  {
    if (logtransform)
    {
      p2=-log10(ppoints(length(p1)));
      p1=-log10(sort(p1));
      xlbl=expression(Expected ~ ~-log[10](italic(p1)));
      ylbl=expression(Observed ~ ~-log[10](italic(p1)));
    }
    else
    {
      p2=ppoints(length(p1));
      p1=sort(p1);
      xlbl=expression(Expected ~ ~italic(p1));
      ylbl=expression(Observed ~ ~italic(p1));
    }
  }
  
  
  if (!is.null(maxP)) { # if user has not specified that there should not be any maxP truncation
    f=p1<=-maxP; p1[f]=-maxP; # logP values below -maxP are truncated to -maxP
    f=p1>=maxP; p1[f]=maxP; # logP values above maxP are truncated to maxP
    f=p2<=-maxP; p2[f]=-maxP; # logP values below -maxP are truncated to -maxP
    f=p2>=maxP; p2[f]=maxP; # logP values above maxP are truncated to maxP
  }
  else { # if user has specified that there should not be any maxP truncation
    f=p1==Inf; p1[f]=sort(p1)[length(p1)-1]; # Inf values are truncated to highest finite value
    f=p1==-Inf; p1[f]=sort(p1)[2]; # -Inf values are truncated to lowest finite value
    f=p2==Inf; p2[f]=sort(p2)[length(p2)-1]; # Inf values are truncated to highest finite value
    f=p2==-Inf; p2[f]=sort(p2)[2]; # -Inf values are truncated to lowest finite value
  }
  
  # part 4: reduce size for fast plotting
  if (length(colour)>1) # if user provides a colour vector
  {
    if (length(colour)!=length(p1)) { stop("p-value and colour must have same length."); } #  check input
    else {
      m=data.frame(p1,p2,colour);
    }
    
    if (speedup) # rounds round everything to 3 digits
    {
      digs=3; m$p1=round(m$p1,digits=digs); m$p2=round(m$p2,digits=digs); f=duplicated(m); m=m[!f,];
      # print('paste(nrow(m),"snps reduced to",nrow(m)));
    }
    
    # remove variables to reduce memory hog before plotting
    rm(p1,p2,colour);
    
    # part 5: plot
    if (!missing(xlab)) { xlbl=xlab; } # if xlab provided, use that instead
    if (!missing(ylab)) { ylbl=ylab; } # if ylab provided, use that instead
    fac=0.4; xbnd=c(0,max(m$p2)+fac); ybnd=c(0,max(m$p1)+fac); # add a bit of space in top corner
    
    p=ggplot() + geom_point(aes(x=m$p2, y=m$p1), color=colour, size=size) +
      scale_x_continuous(limits=xbnd, name=xlbl) +
      scale_y_continuous(name=ylbl, limits=ybnd) +
      theme_bw() + theme(legend.position = "none") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_abline(intercept=0, slope=1, color="red")
    
    lambda_text <- paste("\u03BB", "=", round(lmb, 4))
    
    if (lambda) { p=p + geom_text(aes(x = max(m$p2) * 0.05, y = max(m$p1), label = lambda_text), hjust = 0, vjust = 1) } # show lambda unless requested
    
    #end
    return(lmb);
  }
  else
  {
    m=as.data.frame(cbind(p1,p2));
    
    if (speedup) # rounds round everything to 3 digits
    {
      digs=3; m$p1=round(m$p1,digits=digs); m$p2=round(m$p2,digits=digs); f=duplicated(m); m=m[!f,];
      # print('paste(nrow(m),"snps reduced to",nrow(m)));
    }
    
    
    # remove variables to reduce memory hog before plotting
    rm(p1,p2);
    
    # part 5: plot
    if (!missing(xlab)) { xlbl=xlab; } # if xlab provided, use that instead
    if (!missing(ylab)) { ylbl=ylab; } # if ylab provided, use that instead
    fac=0.4; xbnd=c(0,max(m$p2)+fac); ybnd=c(0,max(m$p1)+fac); # add a bit of space in top corner
    
    p=ggplot() + geom_point(aes(x=m$p2, y=m$p1), color=colour, size=size) +
      scale_x_continuous(limits=xbnd, name=xlbl) +
      scale_y_continuous(name=ylbl, limits=ybnd) +
      theme_bw() + theme(legend.position = "none") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_abline(intercept=0, slope=1, color="red")
    
    lambda_text <- paste("\u03BB", "=", round(lmb, 4))
    
    if (lambda) { p=p + geom_text(aes(x = max(m$p2) * 0.05, y = max(m$p1), label = lambda_text), hjust = 0, vjust = 1) } # show lambda unless requested
    
    #end
    return(list(p, lmb))
  }
  
}
