#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables

# helper functions
read.coverage <- function(id) {
  # init output data frame with zeros
  coverage <- data.frame(ID=id,
                         Position=1:19041,
                         Depth=0)
  
  # if the coverage file exists, read it
  fname <- paste0("data/genomics/", id, "/coverage.tsv")
  if(file.exists(fname)) {
    coverage <- read.delim(fname, 
                           col.names=c("ID", "Position", "Depth")) %>%
                mutate(ID=id)
  }
  
  # return the data frame
  return(coverage)
}

read.snvs <- function(id) {
  # init output data frame with nothing
  snvs <- matrix(nrow=0, ncol=6) %>%
          as.data.frame()
  colnames(snvs) <- c("ID", "Position", "NT.ID",
                      "NT.ref", "NT.alt", "Frequency")
  
  # if the coverage file exists, read it
  fname <- paste0("data/genomics/", id, "/snvs.vcf")
  if(file.exists(fname)) {
    snvs <- read.delim(fname, 
                       comment.char="#",
                       col.names=c("ID", "Position", "NT.ID", 
                                   "NT.ref", "NT.alt", "QUAL", 
                                   "FILTER", "INFO")) %>%
            filter(FILTER=="PASS") %>%
            mutate(ID=id,
                   NT.ID=paste0(Position, "-", NT.ref, "-", NT.alt),
                   Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[^;$]+"))) %>%
            select(-QUAL, -FILTER, -INFO)
  }
  # return the data frame
  return(snvs)
}


## load data -------------------------------------------------------------------
# sample sheet
meta <- read.csv("samplesheet.csv", na.strings="")

## coverage --------------------------------------------------------------------
# read in coverage files
covs <- meta$ID %>%
        lapply(function(i) {
          paste0("data/genomics/", i, "/coverage.tsv") %>%
            read.delim(col.names=c("ID", "NT.position", "Depth")) %>%
            mutate(ID=i)
        }) %>%
        do.call(rbind, .)

# median coverage
covs %>%
  group_by(ID) %>%
  summarise(MedCov=median(Depth),
            .groups="drop") %>%
  filter(MedCov >= 10) %>%
  left_join(meta, by="ID")

# filter to 7 DPI and plot coverage
supA <- unique(meta$NHP) %>%
        lapply(function(i) {
          meta %>%
            filter(DPI==7, NHP==i) %>%
            left_join(covs, by="ID") %>%
            ggplot(aes(NT.position, Depth+1)) +
            geom_line() +
            geom_hline(yintercept=10, linetype=2, col="lightgrey") +
            scale_y_continuous(limits=c(1, max(covs$Depth)),
                               transform="log10") +
            labs(x="Nucleotide position",
                 y="Aligned read depth",
                 title=paste("BOMV NHP", i)) +
            theme(axis.title=element_blank())
        })
x <- cowplot::get_plot_component(supA[[1]] + theme(axis.title.x=element_text()),
                                 "xlab-b")
y <- cowplot::get_plot_component(supA[[1]] + theme(axis.title.y=element_text()),
                                 "ylab-l")
supA <- cowplot::plot_grid(plotlist=supA, ncol=1)
supA <- cowplot::plot_grid(supA, x, ncol=1, rel_heights=c(25, 1))
supA <- cowplot::plot_grid(y, supA, nrow=1, rel_widths=c(1, 25))
supA
ggsave("analysis/coverage.png",
       units="in", width=3.75, height=6)

# clean up
rm(x, y)

## variation profile -----------------------------------------------------------
snvs <- meta %>%
        filter(DPI==7,
               NHP != "B") %>%
        select(ID) %>%
        unlist() %>%
        lapply(function(i) {
          paste0("data/genomics/", i, "/snvs.vcf") %>%
            read.delim(comment.char="#",
                       col.names=c("ID", "NT.position", "NT.ID", 
                                   "NT.ref", "NT.alt", "QUAL", 
                                   "FILTER", "INFO")) %>%
            filter(FILTER=="PASS") %>%
            mutate(ID=i,
                   Frequency=as.numeric(str_extract(INFO, 
                                                    "(?<=AF=)[^;$]+"))) %>%
            select(ID, NT.position, Frequency)
        }) %>%
        do.call(rbind, .)

# join with coverage to fill in the gaps
supB <- covs %>%
        filter(ID %in% snvs$ID) %>%
        left_join(snvs, by=c("ID", "NT.position")) %>%
        replace_na(list(Frequency=0)) %>%
        left_join(meta, by="ID")
supB <- unique(supB$NHP) %>%
        lapply(function(i) {
          supB %>%
            filter(NHP==i) %>%
            ggplot(aes(NT.position, Frequency)) +
            geom_line() +
            geom_hline(yintercept=0.5, linetype=2, col="lightgrey") +
            ylim(0, 1) +
            labs(x="Nucleotide position",
                 y="SNV frequency",
                 title=paste("BOMV NHP", i)) +
            theme(axis.title=element_blank())
        })
x <- cowplot::get_plot_component(supB[[1]] + theme(axis.title.x=element_text()),
                                 "xlab-b")
y <- cowplot::get_plot_component(supB[[1]] + theme(axis.title.y=element_text()),
                                 "ylab-l")
supB <- cowplot::plot_grid(plotlist=supB, ncol=1)
supB <- cowplot::plot_grid(supB, x, ncol=1, rel_heights=c(25, 1))
supB <- cowplot::plot_grid(y, supB, nrow=1, rel_widths=c(1, 25))
supB
ggsave("analysis/snvs.png",
       units="in", width=3.75, height=6)

# clean up
rm(x, y)

## nonsynonymous mutations > 10% -----------------------------------------------
meta %>%
  filter(DPI==7,
         NHP != "B") %>%
  select(ID) %>%
  unlist() %>%
  lapply(function(i) {
    paste0("data/genomics/", i, "/snvs.csv") %>%
      read.csv() %>%
      mutate(ID=i)
  }) %>%
  do.call(rbind, .) %>%
  filter(Frequency >= 0.1,
         Type=="Nonsynonymous") %>%
  mutate(AA.ID=paste0(Gene.name, " ", AA.ref, AA.pos, AA.alt),
         NT.ID=paste(NT.position, NT.ref, ">", NT.alt)) %>%
  left_join(meta, by="ID") %>%
  select(NHP, NT.ID, AA.ID, Frequency) %>%
  arrange(NT.ID) %>%
  write.csv("analysis/snvs.csv",
            row.names=FALSE)
