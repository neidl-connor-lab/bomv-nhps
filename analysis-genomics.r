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
meta %>%
  filter(DPI==7) %>%
  left_join(covs, by="ID") %>%
  ggplot(aes(NT.position, Depth+1)) +
  geom_line() +
  geom_hline(yintercept=10, linetype=2, col="lightgrey") +
  scale_y_log10() +
  facet_wrap(~NHP, ncol=1) +
  labs(x="Nucleotide position",
       y="Aligned read depth")
ggsave("analysis/coverage.png",
       units="cm", width=8, height=10)

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

# join with coverage to fill in the gaps, then plot
meta %>%
  filter(DPI==7, 
         NHP != "B") %>%
  left_join(covs, by="ID") %>%
  left_join(snvs, by=c("ID", "NT.position")) %>%
  replace_na(list(Frequency=0)) %>%
  ggplot(aes(NT.position, Frequency)) +
  geom_line() +
  geom_hline(yintercept=0.5, linetype=2, col="lightgrey") +
  ylim(0, 1) +
  facet_wrap(~NHP, ncol=1) +
  labs(x="Nucleotide position",
       y="SNV frequency")
ggsave("analysis/snvs.png",
       units="cm", width=8, height=10)

## nonsynonymous mutations -----------------------------------------------------
# nonsynonymous mutations > 10%
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
  write.csv("analysis/snvs.csv",
            row.names=FALSE)
