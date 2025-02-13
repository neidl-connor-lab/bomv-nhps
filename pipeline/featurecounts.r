#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Rsubread))
options(stringsAsFactors=FALSE)

# helper functions
err <- function(msg) {
  print(paste("[ERR]", msg))
  q(save="no", status=1)
}

## read arguments and check inputs ---------------------------------------------
args <- ArgumentParser()
args$add_argument("-f", "--feature", help="GTF quantification feature (transcript)", default="transcript")
args$add_argument("-m",  "--meta", help="GTF meta-feature (gene_name)", default="gene_name")
args$add_argument("-p", "--paired", help="paired-end libraries (FALSE)", default=FALSE, action="store_true")
args$add_argument("-a", "--gtf", help="GTF annotation file", required=TRUE)
args$add_argument("-o", "--ofile", help="output TSV", required=TRUE)
args$add_argument("bams", nargs="+", help="aligned BAM file(s)")
args <- args$parse_args()

# GTF
if(!file.exists(args$gtf)) {
  err(paste("Invalid GTF:", args$gtf))
}

# create output directory if it doesn't exist
if(!dir.exists(dirname(args$ofile))) {
  dir.create(dirname(args$ofile), recursive=TRUE)
}

# check all input files
for(i in args$bams) {
  if(!file.exists(i)) {
    err(paste("Invalid BAM:", i))
  }
}
rm(i)

## run featureCounts and format ------------------------------------------------
df <- featureCounts(args$bams, 
                    annot.ext=args$gtf, 
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType=args$feature, 
                    GTF.attrType="gene_id", 
                    GTF.attrType.extra=args$meta,
                    allowMultiOverlap=TRUE, 
                    countMultiMappingReads=TRUE,
		    ignoreDup=TRUE,
                    isPairedEnd=args$paired)

# extract counts matrix and update column names
counts <- df$counts %>%
          as.data.frame() %>%
          rownames_to_column("GeneID")
colnames(counts) <- colnames(counts) %>% 
                    str_remove("\\.bam")

# if using an additional meta-feature (not gene_id), 
# add in any missing values with GeneIDs
anot <- df$annotation[, args$meta]
replacements <- which(is.na(anot))
gids <- df$annotation[replacements, "GeneID"]
anot[replacements] <- gids
# replace names in counts matrix using left merge
counts <- data.frame(Gene=anot,
                    GeneID=df$annotation$GeneID) %>%
          left_join(counts, by="GeneID") %>%
          select(-GeneID)
rm(anot, gids, replacements)

# write counts matrix
write.table(counts, file=args$ofile, sep="\t", row.names=FALSE)

# write log
colnames(df$stat)[-1] <- colnames(counts)[-1]
write.table(df$stat, file=paste0(args$ofile, ".summary"), sep="\t", row.names=FALSE)

# done; print environment
sessionInfo()

