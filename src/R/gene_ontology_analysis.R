library(argparser, quietly=TRUE)
library(topGO, quietly=TRUE)
library(genefilter, quietly=TRUE)
library(plyr, quietly=TRUE)

WriteTable <- function(data, filename, index.name) { 
    datafile <- gzfile(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))

    writeLines(paste0(dim(data), collapse="\t"), con=datafile, sep="\n")
    writeLines(paste0(header, collapse="\t"), con=datafile, sep='\n')
    write.table(format(data, digits=6), datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

p <- arg_parser("Gene ontology analysis using topGO")
p <- add_argument(p, "diff_feature_file", type="character", help="Custom results file from differential expression.")
p <- add_argument(p, "ensembl2go", type="character", help="")
p <- add_argument(p, "prefix", type="character", help="Outfile prefix")
p <- add_argument(p, '--atac', flag=FALSE, help='Keep note that feature level is peaks.')
p <- add_argument(p, '--peak_to_gene', short='-p', help='A peak to gene annotation map.')
p <- add_argument(p, "--output_dir", short="-o", type="character", default=".", help="")

args <- parse_args(p)

GEAnalysis <- function(res.df){

    ## gene set enrichment

    # get background gene set with similar average expression 
    
    sig_idx <- which(res.df$qvalue < 0.05)
    
    gfobj <- genefinder(as.matrix(subset(res.df, select=ave_logCPM)), sig_idx, 1000, method='manhattan')
    
    bg_idx <- as.vector(sapply(gfobj, function(x) x$indices))
    
    # remove sig genes if any from background
    
    bg_idx <- setdiff(unique(bg_idx), sig_idx)
    
    bg_genes <- rownames(res.df)[bg_idx] 
    sig_genes <- rownames(res.df)[sig_idx]
    gene_universe <- c(bg_genes, sig_genes)
    
    alg <- factor(as.integer(gene_universe %in% sig_genes))
    names(alg) <- gene_universe
    
    # cellular component (CC), biological processes (BP), and molecular function (MF)("Writing output to : '", argv$output_dir)
    
    onts = c( "MF", "BP", "CC" )
    tab = as.list(onts)
    names(tab) = onts
    
    for(i in 1:3){
#        tgd <- new( "topGOdata", ontology=onts[i], allGenes=alg, nodeSize=10, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        tgd <- new( "topGOdata", ontology=onts[i], allGenes=alg, nodeSize=10, annot = annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
    
        resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
        resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    
        tab[[i]] <- GenTable(tgd, Fisher.elim = resultTopGO.elim, 
            Fisher.classic = resultTopGO.classic,
            orderBy = "Fisher.classic" , topNodes = 10)
    }
    
    topGOResults.df <- rbind.fill(tab)
    rownames(topGOResults.df) <- topGOResults.df$GO.ID
    topGOResults.df$GO.ID <- NULL 
    
    WriteTable(topGOResults.df, file=file.path(args$output_dir, paste0(args$prefix, '.topGOResults.txt.gz')), 'GO.ID')

}

res.df = read.table(file=args$diff_feature_file, skip=1, sep='\t', header=TRUE, row.names=1)
GEAnalysis(res.df)

