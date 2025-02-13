#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

## process arguments -----------------------------------------------------------
args <- ArgumentParser()
args$add_argument("-d", "--dir", help="input and output directory", 
                  required=TRUE)
args$add_argument("-a", "--annotation", help="viral GTF",
                  required=TRUE)
args$add_argument("-t", "--threshold", help="minimum aligned read depth",
                  default=100, type="integer")
args <- args$parse_args()

## read inputs -----------------------------------------------------------------
# annotation file
if(file.exists(args$annotation)) { 
  anno <- read.delim(args$annotation,
                     comment.char="#",
                     col.names=c("Segment", "Source", "Feature", "Start", 
                                 "Stop", "Score", "Strand", "Phase", 
                                 "Info")) %>%
          filter(Feature=="CDS") %>%
          select(Segment, Start, Stop)
} else {
  print("[ERR] Annotation file does not exist.")
  q(save="no", status=1)
}

# check input directory
if(!dir.exists(args$dir)) {
  print("[ERR] Input directory does not exist.")
  q(save="no", status=1)
}

# load coverage
covs <- paste0(args$dir, "/coverage.tsv")
if(file.exists(covs)) {
  covs <- read.delim(covs, col.names=c("Segment", "Position", "Depth"))
} else {
  print("[ERR] Coverage file does not exist.")
  q(save="no", status=1)
}

# load SNVs and collapse to 1 SNV per position
snvs <- paste0(args$dir, "/snvs.vcf")
if(file.exists(snvs)) {
  snvs <- snvs %>%
          read.delim(comment.char="#",
                     col.names=c("Segment", "Position", "ID", "NT.ref", 
                                 "NT.alt", "QUAL", "FILTER", "INFO"),
                     colClasses=c(Segment="character", 
                                  Position="numeric",
                                  NT.ref="character",
                                  NT.alt="character")) %>%
          mutate(Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[^;$]+"))) %>%
          filter(FILTER=="PASS",
                 Frequency > 0.001) %>% # LOD for Illumina
          group_by(Segment, Position) %>%
          summarise(Frequency=sum(Frequency),
                    .groups="drop")
} else {
  print("[ERR] VCF file does not exist.")
  q(save="no", status=1)
  
}

## build profile and plot ------------------------------------------------------
# add SNV frequencies to coverage data
prof <- left_join(covs, snvs, by=c("Segment", "Position")) %>%
        # add in missing SNV 0%
        replace_na(list(Frequency=0))

# coverage profile
p1 <- prof %>%
      ggplot() +
      geom_rect(data=anno, 
                aes(xmin=Start, 
                    xmax=Stop, 
                    ymin=0.1, 
                    ymax=max(c(covs$Depth, args$threshold))),
                fill="lightgrey") +
      geom_hline(yintercept=args$threshold, linetype=2, col="darkgrey") +
      geom_line(aes(Position, Depth+0.1)) +
      scale_y_continuous(limits=c(0.1, NA), trans="log10") +
      facet_wrap(~Segment, nrow=1, scales="free_x") +
      labs(x="Nucleotide position",
           y="Aligned read depth",
           title="Genome coverage")

# SNV profile
p2 <- prof %>%
      ggplot() +
      geom_rect(data=anno, 
                aes(xmin=Start, xmax=Stop, ymin=0, ymax=1),
                fill="lightgrey") +
      geom_hline(yintercept=0.5, linetype=2, col="darkgrey") +
      geom_line(aes(Position, Frequency)) +
      scale_y_continuous(limits=c(0, 1.05), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
      facet_wrap(~Segment, nrow=1, scales="free_x") +
      labs(x="Nucleotide position",
           y="Aligned read depth",
           title="SNV profile")

# plot together and save
p0 <- cowplot::plot_grid(p1, p2, ncol=1, labels="AUTO")
ggsave(paste0(args$dir, "/profile.png"),
       plot=p0,
       units="cm", width=20, height=15)

## FIN -------------------------------------------------------------------------
q(save="no", status=0)
