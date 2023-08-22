gene.annotate.closest <- function(data, build, chr = "CHR", bp = "BP", sep="|", border=0) {

# find the gene closest to a snp in the data
# data could be the data frame from an gwas result assoc file
# genelist should follow the Plink gene-list style, with 4 columns: chromosome, start, end, gene name
# the user should ensure that the data and gene list use the same build
# if a snp falls within multiple genes, the separator sep is used to concatenate the gene names into a single string
# border (in BP units) is used to allow snps that are not within a gene but within a distance 'border' to a gene to be considered as sitting within that gene

# part 1: define the core function to annotate genes
gene_closest <- function(pos,gc,sep) {
	dleft=gc$start-pos;
	dleft=pmax(dleft,0*dleft);
	dright=pos-gc$end;
	dright=pmax(dright,0*dright);
	dmin=pmax(dleft,dright);
	f=which(dmin==min(dmin));
	if (any(f)) { return(paste0(gc$gene[f],collapse=sep)); } else { return(NA); }
}

# part 2: read genelist from provided build
if (is.numeric(build)) { build <- ifelse(build > 30, build - 18, build); genelist_name <- paste("hg", build, sep = ""); } # if numeric build input then create genelist_name from build number
else { genelist_name <- build; } # if alphanumeric build input then build input is taken as genelist name

if (!exists(genelist_name)) {
  stop("Invalid build") # check whether genelist name exists
}

genelist <- get(genelist_name)

# part 3: select columns and prepare data
colnames(genelist)=c("chr","start","end","gene");

# adjust boundary of genes if necessary
genelist$start=genelist$start-border;
genelist$end=genelist$end+border;

data=as.data.frame(data);

if (!(chr %in% colnames(data))) { stop(paste("Column", chr, "not found!")); }
if (!(bp %in% colnames(data))) { stop(paste("Column", bp, "not found!")); }

m=data[,c(chr,bp)]; colnames(m)=c("CHR","BP");
m$gene=NA;

# make list of chromosomes that are relevant
unc=intersect(m$CHR,genelist$chr);

# part 3: main loop to annotate
for (i in unc) {
  f=m$CHR==i; mpc=m$BP[f]; gc=genelist[genelist$chr==i,];
  pg=lapply(mpc,gene_closest,gc=gc,sep=sep);
  m$gene[f]=unlist(pg);
}

data$gene=m$gene;
return(data);
}
