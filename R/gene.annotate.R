gene.annotate <- function(data, genelist, chr = "CHR", bp = "BP", sep="|", border=0) {

# data could be the data frame from an gwas result assoc file
# genelist should follow the Plink gene-list style, with 4 columns: chromosome, start, end, gene name
# the user should ensure that the data and gene list use the same build
# if a snp falls within multiple genes, the separator sep is used to concatenate the gene names into a single string
# border (in BP units) is used to allow snps that are not within a gene but within a distance 'border' to a gene to be annotated with that gene

# part 1: define the core function to annotate genes
gene_member <- function(pos,gc,sep) { f=(gc$start<=pos)&(gc$end>=pos); if (any(f)) { return(paste0(gc$gene[f],collapse=sep)); } else { return(NA); } }

# part 2: select columns and prepare data
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
  pg=lapply(mpc,gene_member,gc=gc,sep=sep);
  m$gene[f]=pg;
}

data$gene=unlist(m$gene);
return(data);
}
