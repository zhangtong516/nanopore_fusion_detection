library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(txdbmaker)
library(optparse)


option_list <- list(
  make_option(c("--samplename"), type = "character", default = NULL),
  make_option(c("--jaffaresults"), type = "character", default = NULL), 
  make_option(c("--gtf"), type = "character", default = NULL),
)
parser <- OptionParser(option_list = option_list)

opt <- parse_args(parser)

## load jaffa results
jaffa_results.csv <- opt$jaffaresults
jaffa_results_table <- fread(jaffa_results.csv)
fusionGeneNames <- jaffa_results_table[, get("fusion genes")]


# Step1:  generate fiusion fasta file 
## make reference
txdb <- txdbmaker::makeTxDbFromGFF(opt$gtf) 
exonsByGene <- exonsBy(txdb, 'gene')
exonsByTx <- exonsBy(txdb,"tx", use.names = TRUE)


## preapre gene/tx name conversion

gtf_dt = fread(opt$gtf)
genes <- gtf[V3 == "gene"]
genes[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', V9)]
genes[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', V9)]
gene_table <- unique(genes[, .(gene_id, gene_name)])


# ensemblAnnotations.genes <- read.delim(file = 'Homo_sapiens.GRCh38.91.annotations-genes.txt',header=TRUE)
# ensemblAnnotations.genes <- data.table(ensemblAnnotations.genes, keep.rownames = TRUE)
# ensemblAnnotations.txs <- read.delim(file = 'Homo_sapiens.GRCh38.91.annotations-transcripts.txt',header=TRUE)
# ensemblAnnotations.txs <- data.table(ensemblAnnotations.txs, keep.rownames = TRUE)


## Load transcript sequence information 
geneSeq <- readDNAStringSet(file='Homo_sapiens.GRCh38.91.primary_assembly.fa')
listNames <- unlist(lapply(strsplit(names(geneSeq)," "),'[[',1))

## generate fusion sequences 

fusionGeneSequence = rbindlist(lapply(fusionGeneNames, function(s){    
    genevec <- unlist(strsplit(s, ":"))
    fusionseq = paste(unlist(lapply(seq_along(genevec), function(g){
        geneid <- gene_table[gene_name == genevec[g], gene_id] 
        if(length(geneid) == 0) {return("") }
        ## in case of multiple  IDs to same gene symbol, use the first 1 -- unlikely to happen 
        if(length(geneid) > 1 ) {geneid = geneid[1]}
        tmp_range <- exonsByGene[geneid][[1]]
        seqlevelsStyle(tmp_range) = "NCBI"

        seq_pos <- min(start(tmp_range)):max(end(tmp_range))
        seqChar <- geneSeq[[match(as.character(unique(seqnames(tmp_range))), listNames)]][seq_pos]
        if(as.character(unique(strand(tmp_range))) == "-"){
            seqChar <- Biostrings::reverseComplement(seqChar)
        }
        as.character(seqChar)
    })), collapse = "")
    
    return(data.table(fusion_gene = s, fusion_seq = fusionseq))
})) 

fasta_sequences = Biostrings::DNAStringSet(fusionGeneSequence$fusion_seq)
names(fasta_sequences) = Biostrings::fusionGeneSequence$fusion_gene
Biostrings::writeXStringSet(fasta_sequences, file = paste0(opt$samplename, "_fusionGene_sequences.fa" )) #TODO: define name of the file 


# step2: generate fusion annotation
fusionTx <- getFusionTxRange(fusionGeneNames, ensemblAnnotations.genes, exonsByGene, ensemblAnnotations.txs, exonsByTx,jaffa_results_table)
anno.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs_corrected.gtf"
fusionAnnotations <- getFusionAnnotations(fusionGeneNames, anno.file,ensemblAnnotations.genes,ensemblAnnotations.txs,fusionTx)

## Re-align reads to fusion chromosomes
full_annotations <- bambu:::prepareAnnotationsFromGTF(anno.file)
filtered_annotations <- full_annotations[ensemblAnnotations.txs[hgnc_symbol %in% unlist(strsplit(fusionGeneNames, ":"))]$ensembl_transcript_id]
gr <- unlist(filtered_annotations)
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr))
write.table(df, file="fusion.bed", quote=F, sep="\t", row.names=F, col.names=F)


