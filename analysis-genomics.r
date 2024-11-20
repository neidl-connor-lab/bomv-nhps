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
                         Position=1:19043,
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
meta <- read.csv("analysis/samplesheet.csv", na.strings="") %>%
        # TEMPORARILY REMOVING STOCK
        filter(DPI != "Stock")

# coverage
covs <- lapply(meta$ID, read.coverage) %>%
        do.call(rbind, .)

# check median coverage
medcov <- covs %>%
          group_by(ID) %>% 
          summarize(Depth=median(Depth),
                    .groups="drop") %>% 
          left_join(meta, by="ID")
medcov %>%
  ggplot(aes(as.integer(DPI), Depth+1, col=NHP)) +
  geom_hline(yintercept=100, linetype=2) +
  geom_line() +
  geom_point(size=3) +
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(breaks=c(-8, 1, 3, 5, 7, 10, 14, 21, 28)) +
  scale_y_continuous(limits=c(1, 1e3), trans="log10") +
  labs(x="Days postchallenge",
       y="Median BOMV aligned read depth") +
  theme(legend.position=c(0.1, 0.35), 
        legend.box.background=element_rect(color="black"))
ggsave("analysis/median-coverage.png", 
       units="cm", width=10, height=8)

# subset to samples with >100x median coverage only
medcov <- medcov[medcov$Depth > 100, ]$ID
meta <- filter(meta, ID %in% medcov)
covs <- filter(covs, ID %in% medcov)
rm(medcov)

# load SNVs
snvs <- lapply(meta$ID, read.snvs) %>%
        do.call(rbind, .)

# build profile in coverage
covs <- snvs %>%
        # collapse by position and ID
        group_by(ID, Position) %>%
        summarise(Frequency=sum(Frequency),
                  .groups="drop") %>%
        right_join(covs, by=c("ID", "Position")) %>%
        # fill in frequency zeros
        replace_na(list(Frequency=0)) %>%
        # add in metadata
        left_join(meta, by="ID")

## plot profiles ---------------------------------------------------------------
# coverage
p1 <- covs %>%
      ggplot(aes(Position, Depth+1)) +
      geom_hline(yintercept=100, linetype=2, col="lightgrey") +
      geom_line() +
      scale_y_log10() +
      facet_wrap(~NHP, ncol=1) +
      labs(x="Nucleotide position",
           y="Aligned read depth")

# SNVs
p2 <- covs %>%
      ggplot(aes(Position, Frequency)) +
      geom_hline(yintercept=0.5, linetype=2, col="lightgrey") +
      geom_line() +
      facet_wrap(~NHP, ncol=1) +
      labs(x="Nucleotide position",
           y="SNV frequency")

# plot together
cowplot::plot_grid(p1, p2, labels="AUTO")
ggsave("analysis/profiles.png",
       units="cm", width=16, height=8)

# SNV table
snvs %>%
  filter(Frequency > 0.5) %>%
  left_join(meta, by="ID") %>%
  select(NHP, DPI, NT.ID, Frequency) %>%
  write.csv("analysis/consensus-snvs.csv",
            row.names=FALSE)
