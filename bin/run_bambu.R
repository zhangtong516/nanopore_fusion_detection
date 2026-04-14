library(bambu)
library(GenomicAlignments)
suppressWarnings(suppressMessages(require(optparse)))

load("bambu_functions.R")

option_list <- list(
  make_option(c("--samplename"), type = "character", default = NULL), 
  make_option(c("--reads"), type = "character", default = NULL),
  make_option(c("--rc_out_dir"), type = "character", default = "readClassSaveDir"),
  make_option(c("--annotations"), type = "character", default = NULL),
  make_option(c("--genome"), type = "character", default = NULL),
  make_option(c("--fusion_mode"), type = "logical", default = TRUE),
  make_option(c("--stranded"), type = "logical", default = FALSE),
  make_option(c("--ncore"), type = "integer", default = 8),
  make_option(c("--yield_size"), type = "integer", default = 1000000),
  make_option(c("--verbose"), type = "logical", default = FALSE),
)
parser <- OptionParser(option_list = option_list)

opt <- parse_args(parser)

if (is.null(opt$samplename) || is.null(opt$reads) || is.null(opt$annotations) || is.null(opt$genome)) {
  stop("Missing required arguments: --samplename, --reads, --annotations, --genome")
}
dir.create(opt$rc_out_dir, showWarnings = FALSE, recursive = TRUE)
ann <- opt$annotations
if (grepl("\\\\.rds$", ann, ignore.case = TRUE)) {
  annotations_obj <- readRDS(ann)
} else {
  annotations_obj <- prepareAnnotations(ann)
}

rcSaveDir = paste0(opt$samplename, "_readClassSaveDir/") 
out_rds = paste0(opt$samplename, "_seFusion.rds")


## run BAMBU 
seFusion <- bambu(
  reads = opt$reads,
  rcOutDir = opt$rc_out_dir,
  annotations = annotations_obj,
  genome = opt$genome,
  fusionMode = opt$fusion_mode,
  stranded = opt$stranded,
  ncore = opt$ncore,
  yieldSize = opt$yield_size,
  verbose = opt$verbose
)
saveRDS(seFusion, file = opt$out_rds)


# visualization of fusion for coverage
reads <-  readGAlignments(file=opt$reads, 
        param=ScanBamParam(flag=scanBamFlag(
            isSecondaryAlignment=FALSE, 
            isDuplicate = FALSE), 
        what=c("qual", "flag","mapq")),
        use.names=T) 

readsFiltered <- reads[mcols(reads)$flag %in% c(0,16)]

## add read count column to BAMBU output 
annotations <- rowRanges(seFusion)
mcols(annotations)$readCount <- NULL
mcols(annotations)$readCount <- rowSums(assays(seFusion)$CPM)

break.pointList <- getFusionBreakpoints(
    fusionGeneNames, 
    ensemblAnnotations.genes, 
    exonByGene, 
    jaffa_results_table)

fusionGene <- getFusionGeneRange(
    fusionGeneNames, 
    ensemblAnnotations.genes,
    exonsByGene, 
    jaffa_results_table)

prime5Gene <- GetPrimeGene(
    fusionGeneNames, 
    ensemblAnnotations.genes, 
    exonsByGene, 
    jaffa_results_table, 
    prime = 5)
prime3Gene <- GetPrimeGene(
    fusionGeneNames, 
    ensemblAnnotations.genes, 
    exonsByGene, 
    jaffa_results_table, 
    prime = 3)

gene_type <- defineGeneType(
    fusionGeneNames, 
    annotations, 
    fusionGene,
    prime3Gene,
    prime5Gene)

## plot with annotations on top and read coverage in the bottom.
library(ggbio)
#gr.tx <- annotations[as.character(unlist(unique(seqnames(annotations))))]
gene_tmp <- fusionGeneNames[1]
break.point <- break.pointList[gene_tmp][[1]]
# strand.info <- as.character(unique(strand(fusionAnnotations[gene_tmp][[1]])))
fusionGene <- GRangesList(fusionGene)
names(fusionGene) <- fusionGeneNames

p1 <- ggbio::autoplot(fusionGene[gene_tmp], aes( col = as.factor(score), fill = as.factor(score)), group.selfish = TRUE)+
      #guides(col = FALSE, fill = FALSE)+
      geom_vline(xintercept = break.point, col = "red", alpha = 0.3)+
      # geom_vline(xintercept = validated_breakpoint, col = "blue", alpha = 0.8)+
      scale_color_brewer(type = "qual", guide = FALSE)+scale_fill_brewer(type = "qual", labels = c("5' Gene","3' Gene"), name = "Gene type")
    
    
p2 <- ggbio::autoplot(GRangesList(fusionTx[[gene_tmp]]), aes(col = as.factor(score), fill = as.factor(score)), group.selfish = TRUE)+geom_vline(xintercept = break.point, col = "red", alpha = 0.3)+scale_color_brewer(type = "qual", guide = FALSE)+scale_fill_brewer(type = "qual", labels = c("5' Gene","3' Gene"), name = "Gene type")
    
tmp <- readsFiltered[as.character(seqnames(readsFiltered))==gene_tmp]
p3 <- ggbio::autoplot(tmp, method = "raw", stat = "coverage", aes(col = coverage, fill = coverage))+scale_fill_distiller(type = "seq",  direction = -1,n.breaks = 3)+scale_color_distiller(type = "seq", direction = -1,n.breaks = 3)

    
txvec_id <- unlist(lapply(3:1, function(i){
      gene_pos <- c("5PrimeGene","3PrimeGene","fusion_gene")[i]
      txvec_id <- which(as.character(unique(seqnames(annotations)))==gene_tmp&(mcols(annotations)$readCount>=10)&names(annotations) %in% gene_type[fusion_gene == gene_tmp &(gene_type == gene_pos)]$tx_name)
      tmp_anno <- annotations[txvec_id]
      sort.id <- order(min(start(tmp_anno)),max(end(tmp_anno)))
      return(txvec_id[sort.id])
}))
    
    
tmp_anno <- annotations[txvec_id] # 
tmp_anno <- tmp_anno[!is.na(mcols(tmp_anno)$readCount)]
p4 <- ggbio::autoplot(tmp_anno,aes(col = log10(readCount), fill = log10(readCount)), group.selfish = TRUE)+scale_fill_distiller(type = "seq",  direction = 1)+scale_color_distiller(type = "seq", direction = 1)
    
    
tracks(Gene = p1,Transcript = p2,  Coverage = p3,Count = p4,
                heights = c(1,3,3,6))+theme_classic()