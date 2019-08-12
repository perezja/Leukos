library(argparser, quietly=TRUE)
library(dplyr)

p <- arg_parser("PCA of sequencing count features.")
p <- add_argument(p, "count_matrix", type="character", help="Count matrix")
p <- add_argument(p, "prefix", type="character", help="Outfile prefix")
p <- add_argument(p, "--output_dir", short="-o", type="character", default=".", help="")

args <- parse_args(p)

if(!dir.exists(args$output_dir)){ dir.create(args$output_dir)}

outfile <- file.path(args$output_dir, paste0(args$prefix, '.pca.txt')) 

count_df <- read.table(args$count_matrix, header=TRUE, sep='\t')

rownames(count_df) <- count_df[[1]]

count_df <- count_df[,2:ncol(count_df)]

count_matrix <- as.matrix(count_df)

pca_obj = prcomp(count_matrix)
pcomps = pca_obj$rotation

pcomps <- as.data.frame(pcomps) %>% tibble::rownames_to_column(., var='sample')

outfile = file.path(args$output_dir, paste0(args$prefix, '.pca.txt'))
write.table(pcomps, outfile, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)

