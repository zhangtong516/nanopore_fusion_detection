library(devtools)
install_github("GoekeLab/bambu", ref = "restore_fusion_function")
library(dplyr)
library(bambu)
getFusionAnnotations <- function(fusionGeneNames, anno.file,ensemblAnnotations.genes,ensemblAnnotations.txs,fusionTx){
    annotationRanges <- prepareAnnotations(anno.file)
    fusionAnnotation <- lapply(fusionGeneNames, function(s){
    print(s)
    genevec <- unlist(strsplit(s, ":"))
    tmp_range <- do.call("c",lapply(genevec, function(g){
    geneid <- ensemblAnnotations.genes[hgnc_symbol==g]$ensembl_gene_id
    txid <- ensemblAnnotations.txs[ensembl_gene_id==geneid]$ensembl_transcript_id
    print(length(txid))
    tmp <- do.call("c",lapply(txid, function(t){
      print(t)
      tmp <- annotationRanges[t]
      fusionTmp <- fusionTx[[s]][[t]]
      seqlevels(tmp, pruning.mode = "tidy") <- as.character(unique(seqnames(tmp)))
      seqlevels(tmp) <- as.character(unique(seqnames(fusionTmp)))                                 
      start(tmp[[1]]) <- start(fusionTmp)
      end(tmp[[1]]) <- end(fusionTmp)
      # if negative strand
      if(unique(strand(tmp[[1]]))=="-"){
        tmp[[1]]$exon_rank <- max(tmp[[1]]$exon_rank)-tmp[[1]]$exon_rank+1
        tmp[[1]]$exon_endRank <- max(tmp[[1]]$exon_endRank)-tmp[[1]]$exon_endRank+1
      }
      strand(tmp[[1]]) <- strand(fusionTmp)
      return(tmp)
    }))
    return(tmp)
  }))
  return(tmp_range)
})
keep.id <- which(!sapply(fusionAnnotation, is.null))
fusionAnnotation <- fusionAnnotation[keep.id]
fusionAnnotation <- do.call("c", fusionAnnotation)
return(fusionAnnotation)
}

getFusionTxRange <- function(fusionGeneNames, ensemblAnnotations.genes, exonsByGene, ensemblAnnotations.txs, exonsByTx, jaffa_results_table){
    fusionTx <- lapply(seq_along(fusionGeneNames), function(s){
  print(s)
  tmp <- fusionGeneNames[s]
  genevec <- unlist(strsplit(tmp, ":"))
  geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[1]]$ensembl_gene_id
  tmp_range <- exonsByGene[geneid][[1]]
  min_start <- min(start(tmp_range)) 
  length_g1 <- max(end(tmp_range)) - min_start + 1
  tmp_range <- lapply(seq_along(genevec), function(g){
    geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[g]]$ensembl_gene_id
    tmp_range <- exonsByGene[geneid][[1]]
    min_start <- min(start(tmp_range)) 
    length_g_tmp <- max(end(tmp_range)) 
    strand.info <- as.character(unique(strand(tmp_range)))
    txid <- ensemblAnnotations.txs[(ensembl_gene_id == geneid)]$ensembl_transcript_id
    new_range <- lapply(txid, function(t){
      tmp_range <- exonsByTx[t][[1]]
      if(g==1){
        scoreNum = ifelse(strand.info=="-", length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base1+1,
                          jaffa_results_table[fusion.genes == tmp]$base1-min_start+1)
      }else{
        scoreNum = ifelse(strand.info=="-",length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base2+1+length_g1,
                          jaffa_results_table[fusion.genes == tmp]$base2-min_start+1+length_g1)
      }
      if(strand.info=="-"){
        new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                             ranges = IRanges(start = length_g_tmp - end(tmp_range) + 1, 
                                              end = length_g_tmp - start(tmp_range) + 1 ),
                             strand = Rle("+",length(tmp_range)),
                             score = scoreNum)
      }else{
        new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                             ranges = IRanges(start = start(tmp_range) - min_start+1, 
                                              end = end(tmp_range) - min_start+1),
                             strand = Rle(strand.info,length(tmp_range)),
                             score = scoreNum)
      }
      if(g == 2){
        new_range <- GenomicRanges::shift(new_range, shift = length_g1)
        
      }
      return(new_range)  
    })
    names(new_range) <- txid
    return(new_range)
  })
  
  return(do.call("c",tmp_range))
})
names(fusionTx) <- fusionGeneNames
return(fusionTx)
}

getFusionGeneRange <- function(fusionGeneNames, ensemblAnnotations.genes,exonsByGene, jaffa_results_table){
  fusionGene <- lapply(seq_along(fusionGeneNames), function(s){
  print(s)
  tmp <- fusionGeneNames[s]
  genevec <- unlist(strsplit(tmp, ":"))
  geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[1]]$ensembl_gene_id
  tmp_range <- exonsByGene[geneid][[1]]
  min_start <- min(start(tmp_range)) 
  length_g1 <- max(end(tmp_range)) - min_start + 1
  tmp_range <- lapply(seq_along(genevec), function(g){
    geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[g]]$ensembl_gene_id
    tmp_range <- exonsByGene[geneid][[1]]
    min_start <- min(start(tmp_range)) 
    length_g_tmp <- max(end(tmp_range)) 
    strand.info <- as.character(unique(strand(tmp_range)))
    if(g==1){
      scoreNum = ifelse(strand.info=="-", length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base1+1,
                        jaffa_results_table[fusion.genes == tmp]$base1-min_start+1)
    }else{
      scoreNum = ifelse(strand.info=="-",length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base2+1+length_g1,
                        jaffa_results_table[fusion.genes == tmp]$base2-min_start+1+length_g1)
    }
    if(strand.info=="-"){
      new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                           ranges = IRanges(start = rev(length_g_tmp - end(tmp_range) + 1 ), 
                                            end = rev(length_g_tmp - start(tmp_range) + 1 )),
                           strand = Rle("+",length(tmp_range)), # change strand information to +
                           score = scoreNum)
    }else{
      new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                           ranges = IRanges(start = start(tmp_range) - min_start+1, 
                                            end = end(tmp_range) - min_start+1),
                           strand = Rle(strand.info,length(tmp_range)),
                           score = scoreNum)
    }
    if(g == 2){
      new_range <- GenomicRanges::shift(new_range, shift = length_g1)
      
    }
    
    return(new_range)
    
  })
  
  return(do.call("c",tmp_range))
})
  return(fusionGene)
}

# Get Fusion Break Points 
getFusionBreakpoints <- function(fusionGeneNames, ensemblAnnotations.genes, exonByGene, jaffa_results_table){
    break.pointList <- lapply(seq_along(fusionGeneNames), function(s){
  print(s)
  tmp <- fusionGeneNames[s]
  genevec <- unlist(strsplit(tmp, ":"))
  geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[1]]$ensembl_gene_id
  tmp_range <- exonsByGene[geneid][[1]]
  min_start <- min(start(tmp_range)) 
  length_g1 <- max(end(tmp_range)) - min_start + 1
  scoreNum <- unlist(lapply(seq_along(genevec), function(g){
      geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[g]]$ensembl_gene_id
    tmp_range <- exonsByGene[geneid][[1]]
    min_start <- min(start(tmp_range)) 
    length_g_tmp <- max(end(tmp_range))
    strand.info <- as.character(unique(strand(tmp_range)))
    
    if(g==1){
      if(strand.info=="-"){
        scoreNum <- length_g_tmp-unique(jaffa_results_table[fusion.genes == tmp]$base1)+1
      }else{
        scoreNum <-  unique(jaffa_results_table[fusion.genes == tmp]$base1)-min_start+1
      }
      
    }else{
      if(strand.info=="-"){
        scoreNum <- length_g_tmp-unique(jaffa_results_table[fusion.genes == tmp]$base2)+1+length_g1
      }else{
        scoreNum <- unique(jaffa_results_table[fusion.genes == tmp]$base2)-min_start+1+length_g1
      }
      }
    }))
    return(scoreNum)
  })
  names(break.pointList) <- fusionGeneNames
  return(break.pointList)
}


GetPrimeGene <- function(fusionGeneNames,ensemblAnnotations.genes,exonsByGene,jaffa_results_table,prime = 5){
primeGene <- lapply(seq_along(fusionGeneNames), function(s){
  print(s)
  tmp <- fusionGeneNames[s]
  genevec <- unlist(strsplit(tmp, ":"))
  geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[1]]$ensembl_gene_id
  tmp_range <- exonsByGene[geneid][[1]]
  min_start <- min(start(tmp_range)) 
  length_g1 <- max(end(tmp_range)) - min_start + 1
  tmp_range <- lapply(which(c(5,3)==prime), function(g){
    
      geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[g]]$ensembl_gene_id
    tmp_range <- exonsByGene[geneid][[1]]
    min_start <- min(start(tmp_range)) 
    length_g_tmp <- max(end(tmp_range)) 
    strand.info <- as.character(unique(strand(tmp_range)))
    if(g==1){
      scoreNum = ifelse(strand.info=="-", length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base1+1,
                        jaffa_results_table[fusion.genes == tmp]$base1-min_start+1)
    }else{
      scoreNum = ifelse(strand.info=="-",length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base2+1+length_g1,
                        jaffa_results_table[fusion.genes == tmp]$base2-min_start+1+length_g1)
    }
    if(strand.info=="-"){
      new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                           ranges = IRanges(start = rev(length_g_tmp - end(tmp_range) + 1 ), 
                                            end = rev(length_g_tmp - start(tmp_range) + 1 )),
                           strand = Rle("+",length(tmp_range)),
                           score = scoreNum)
    }else{
      new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                           ranges = IRanges(start = start(tmp_range) - min_start+1, 
                                            end = end(tmp_range) - min_start+1),
                           strand = Rle(strand.info,length(tmp_range)),
                           score = scoreNum)
    }
    if(g == 2){
      new_range <- GenomicRanges::shift(new_range, shift = length_g1)
      
    }
    
    return(new_range)
    
  })
  
  return(do.call("c",tmp_range))
})
    return(primeGene)
}

defineGeneType <- function(fusionGeneNames, annotations, fusionGene,prime3Gene, prime5Gene){
    gene_type <- do.call("rbind",lapply(seq_along(fusionGeneNames), function(s){
  
  print(s)
  genevec <- unlist(strsplit(fusionGeneNames[s], ":"))
  fusionHits <- unique(c(which(overlapsAny(annotations, range(GRangesList(fusionGene[s])), type = "within", ignore.strand = TRUE)),
                         queryHits(findOverlaps(annotations, range(GRangesList(fusionGene[s])), type = "within", ignore.strand = TRUE))))
  match_type <- data.table(fusion_gene = fusionGeneNames[s],
                           gene_type = "fusion_gene",
                           tx_name = names(annotations)[fusionHits])
  prime5Hits <- unique(c(which(overlapsAny(annotations, range(GRangesList(prime5Gene[s])), type = "within", ignore.strand = TRUE)),
                         queryHits(findOverlaps(annotations, range(GRangesList(prime5Gene[s])), type = "within", ignore.strand = TRUE))))
  prime3Hits <- unique(c(which(overlapsAny(annotations, range(GRangesList(prime3Gene[s])), type = "within", ignore.strand = TRUE)),
                         queryHits(findOverlaps(annotations, range(GRangesList(prime3Gene[s])), type = "within", ignore.strand = TRUE))))
  match_type[tx_name %in% names(annotations)[prime5Hits], gene_type := "5PrimeGene"]
  match_type[tx_name %in% names(annotations)[prime3Hits], gene_type := "3PrimeGene"]
  return(match_type)
}))
    return(gene_type)
}
