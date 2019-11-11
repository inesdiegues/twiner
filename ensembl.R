# ENSEMBL IDS THAT ARE PROTEIN CODING 

# source("http://bioconductor.org/biocLite.R")
# biocLite(c('ensembldb', 'EnsDb.Hsapiens.v86'))

suppressPackageStartupMessages({
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(futile.logger)
  .Last.value <- flog.layout(layout.format('~m'))
})

edb <- EnsDb.Hsapiens.v86

ensembl.protein.coding <- genes(edb,
                                filter = list(GenebiotypeFilter('protein_coding')),
                                columns = c('gene_id', 'gene_name'))
{
  flog.info('         Granges: %d', nrow(ensembl.protein.coding@elementMetadata))
  flog.info('Metadata columns: %d', ncol(ensembl.protein.coding@elementMetadata))
}

# CDDS - PROJECT THAT IDENTIFIES WHICH ARE THE MEANINGFUL PROTEINS IN HUMANS 

ccds <- read.table(url('ftp://ftp.ncbi.nih.gov/pub/CCDS/current_human/CCDS.current.txt'),
                   sep = '\t',
                   header = T,
                   comment.char = "|", 
                   stringsAsFactors = FALSE)
flog.info('Size of ccds: %d x %d', nrow(ccds), ncol(ccds))

ensembl.genes         <- sort(unique(ensembl.protein.coding@elementMetadata$gene_name))
ensembl.genes.ensg.id <- sort(ensembl.protein.coding@elementMetadata$gene_id)

ccds.genes <- sort(unique(ccds$gene))
ccds.extra.genes <- sort(ccds.genes[(!ccds.genes %in% ensembl.genes)])
ccds.extra.genes.ensg.id <- genes(edb, filter = list(GenenameFilter(ccds.extra.genes)),
                                  columns = c('gene_id', 'gene_name'))

ensg.id <- sort(unique(c(ensembl.protein.coding@elementMetadata$gene_id, ensembl.genes.ensg.id)))

all.genes <- sort(unique(c(ensembl.genes, ccds.extra.genes)))
