#!/usr/bin/env Rscript
library(argparser, quietly=TRUE)
library(edgeR, quietly=TRUE)
library(qvalue, quietly=TRUE)
library(plyr, quietly=TRUE)

WriteTable <- function(data, filename, index.name) { 
    datafile <- gzfile(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))

    writeLines(paste0(dim(data), collapse="\t"), con=datafile, sep="\n")
    writeLines(paste0(header, collapse="\t"), con=datafile, sep='\n')
    write.table(format(data, digits=6), datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

p <- arg_parser("Run DE analysis using edgeR.")
p <- add_argument(p, "expression_matrix", help="Matrix feature counts across for all groups.")
p <- add_argument(p, "group_list", help="File listing column group integers in expression matrix (1/2's).")
p <- add_argument(p, "prefix", help="Prefix for output file.") 
p <- add_argument(p, "--fdr", default=0.05, help="") 
p <- add_argument(p, "--output_dir", short="-o", default='.', help="") 

argv <- parse_args(p)

if(!dir.exists(argv$output_dir)){
    dir.create(argv$output_dir)
}

groups = scan(argv$group_list, numeric())
design = model.matrix(~groups) 

df <- read.table(argv$expression_matrix, sep='\t', skip=1, header=TRUE, row.names=1)

y <- DGEList(counts=df, group=groups)

keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]

y <- estimateDisp(y, design)

fit <- glmQLFit(y)
qlf <- glmQLFTest(fit, coef=2)

qvalobj <- qvalue(qlf$table$PValue)

cat(" * Proportion of significant phenotypes (1-pi0): ", round(1-qvalobj$pi0, 2), "\n", sep="")
cat(" * eGenes @ FDR ", argv$fdr, ":	", sum(qvalobj$qvalues < argv$fdr), "\n", sep="")

res.df <- data.frame(intercept=qlf$coefficients[,1], coefficient=qlf$coefficients[,2], ave_logCPM=qlf$AveLogCPM, pvalue=qlf$table$PValue, logFC=qlf$table$logFC, qvalue=qvalobj$qvalues)

WriteTable(res.df, file.path(argv$output_dir, paste0(argv$prefix, '.summary.txt.gz')), 'gene_id') 

egenes.df <- res.df[res.df$qvalue<argv$fdr,]
egenes.df <- egenes.df[order(egenes.df$qvalue),]

WriteTable(egenes.df, file.path(argv$output_dir, paste0(argv$prefix, '.eGenes.txt.gz')), 'gene_id') 
