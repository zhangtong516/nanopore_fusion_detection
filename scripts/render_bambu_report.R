#!/usr/bin/env Rscript
# Render bambu.Rmd report with specified parameters
# Usage: Rscript render_bambu_report.R --sample-name <name> --jaffa-results <path> --gtf <path> --ref-genome <path> --output-dir <dir>

library(optparse)

option_list <- list(
  make_option(c("--samplename"), type = "character", default = "test_sample",
              help = "Sample name for the analysis [default: %default]"),
  make_option(c("--bamfile"), type = "character", default = "alignment.bam",
              help = "Path to BAM file from dorado aligner [default: %default]"),
  make_option(c("--jaffaresults"), type = "character", default = "jaffaresults.tsv",
              help = "Path to JAFFA results file [default: %default]"),
  make_option(c("--gtf"), type = "character", default = "grch38.gtf",
              help = "Path to GTF annotation file [default: %default]"),
  make_option(c("--refgenome"), type = "character", default = "grch38.fasta",
              help = "Path to reference genome FASTA file [default: %default]"),
  make_option(c("--rmdfile"), type = "character", default = "bambu.Rmd",
              help = "Path to the Rmarkdown file to render [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


cat("Rendering bambu report with the following parameters:\n")
cat("  Sample name:", opt$`samplename`, "\n")
cat("  BAM file:", opt$`bamfile`, "\n")
cat("  JAFFA results:", opt$`jaffaresults`, "\n")
cat("  GTF file:", opt$gtf, "\n")
cat("  Reference genome:", opt$`refgenome`, "\n")
cat("  Rmarkdown file:", opt$`rmdfile`, "\n")

# Render the Rmarkdown file
tryCatch({
  # Copy the Rmarkdown file to current directory if it's not already there
  rmd_source <- opt$rmdfile
  rmd_local <- basename(rmd_source)
  
  if (!file.exists(rmd_local) && file.exists(rmd_source)) {
    file.copy(rmd_source, rmd_local, overwrite = TRUE)
    cat("Copied Rmarkdown file to:", rmd_local, "\n")
  }
  
  # Render using the local copy
  rmarkdown::render(
    rmd_local,
    params = list(
      samplename = opt$samplename,
      bam = opt$bamfile,
      jaffaresults = opt$jaffaresults,
      gtf = opt$gtf,
      refgenomefasta = opt$refgenome
    ),
    output_file = paste0(opt$samplename, "_bambu.html"),
    quiet = FALSE
  )
  cat("\n✓ Report rendered successfully to:", paste0(opt$samplename, "_bambu.html"), "\n")
}, error = function(e) {
  cat("\n✗ Error rendering report:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})
