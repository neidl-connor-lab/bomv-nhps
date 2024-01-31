#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(Rsamtools))
options(stringsAsFactors=FALSE)
theme_set(theme_classic() + 
            theme(strip.text=element_text(color="white"),
                  strip.background=element_rect(fill="black"),
                  legend.box.background=element_rect(color="black", 
                                                     linewidth=1)))

# colors
cols.grp <- c("Stock"="black", "7 DPI"="darkgrey")
cols.prd <- c(sGP="black", GP="#e41a1c", ssGP="#377eb8", Other="#4daf4a")

# helper variables
orfs <- read.csv("pipeline/annotations/bomv-MF319186.1.csv") %>%
        mutate(Gene=factor(Gene, levels=Gene))

# helper functions
annot.coding <- function(i, anno) {
  # format row as its own data.frame since reading in with apply
  i <- data.frame(SNV=as.character(i["SNV"]),
                  Position=as.integer(i["Position"]),
                  Ref=as.character(i["Ref"]),
                  Alt=as.character(i["Alt"]))
  
  # set variables
  pos <- i$Position
  alt <- i$Alt
  
  # get ORF (thankfully no overlapping ORFs)
  orf <- anno %>%
         filter(End >= pos,
                Start <= pos)
  # extract & format nucleotide sequence 
  nt.seq <- orf[1, "NucSeq"] %>%
            str_split("") %>%
            unlist()
  # re-number nucleotides
  names(nt.seq) <- as.character(seq(orf$Start, orf$End, by=1))
  # extract & format AA sequence 
  aa.seq <- orf[1, "AASeq"] %>%
            str_split("") %>%
            unlist()
  
  # check for deletion, which requires special handling
  if(nchar(i$Ref) > 1) {
    # how many nucletides deleted?
    del <- nchar(i$Ref)-1
    # get positions to be deleted
    del <- seq(pos+1, pos+del, by=1) %>%
           as.character()
    # use that to get indices of positions to keep
    del <- which(!(names(nt.seq) %in% del))
    # update nt.seq
    nt.seq <- nt.seq[del]
  }
  
  # substitute nucleotide, translate, and format
  nt.seq[as.character(pos)] <- alt
  nt.seq <- nt.seq %>%
            paste(collapse="") %>%
            DNAString() %>%
            translate() %>%
            as.character() %>%
            strsplit("") %>%
            unlist()
  # compare to aaseq
  aa.mut <- which(aa.seq != nt.seq)
  
  # if no aa.mut, the SNV is synonymous! 
  # format output
  if(length(aa.mut)==0) {
    i$Type <- "Synonymous"
    i$Gene <- orf$Gene
    i$Label <- paste(orf$Gene, "synonymous", format(i$Position, big.mark=","),
                     i$Ref, ">", i$Alt)
  } else { # if AA change, the SNV is nonsynonymous
    # if insertion, there will be many AA changes, so just take the 1st one
    aa.mut <- aa.mut[1]
    # build return data.frame
    # set up as nonsynonymous, insertion, or deletion
    if(nchar(i$Ref) > 1) {
      i$Type <- "Deletion"
    } else if(nchar(alt) > 1) {
      i$Type <- "Insertion"
    } else {
      i$Type <- "Nonsynonymous"
    }
    i$Gene <- orf$Gene
    i$AA.ref <- aa.seq[aa.mut]
    i$AA.alt <- nt.seq[aa.mut]
    i$AA.position <- aa.mut
    i$Label <- paste0(orf$Gene, " ", aa.seq[aa.mut], aa.mut, nt.seq[aa.mut])
  }
  
  # return annotated data.frame
  return(i)
}
load.spanning <- function(id, seq="MF319186.1", a1=6910, a7=6916, bookend=7) {
  # define range
  site.bookend <- GRanges(seqnames=seq, 
                          ranges=IRanges(start=a1-bookend,
                                         end=a7+bookend))
  # load from bam
  paste0("data/snvs/", id, ".bam") %>%
    scanBam(param=ScanBamParam(which=site.bookend,
                               what=scanBamWhat())) %>%
    # outputs a list of 1 item, so extract it
    magrittr::extract2(1) %>%
    # list of vectors, so coerce to data.frame
    as.data.frame() %>%
    # make columns defining what bases each aligned read spans
    mutate(Start=pos,
           End=pos+qwidth-1) %>%
    # only keep reads that span stutter+bookends
    # don't removed based on fwd/reverse -- an alignment is an aligment
    filter(Start <= a1-bookend,
           End >= a7+bookend) %>%
    # summarize duplicate reads
    group_by(seq, Start, End) %>%
    summarise(Reads=n(),
              .groups="drop") %>%
    # add ID
    mutate(ID=id)
}
annot.spanning <- function(i, a1=6910, a7=6916, bookend=7) {
  # format input as data.frame (receiving as vector from apply)
  i <- data.frame(Sequence=i["seq"],
                  Start=as.integer(i["Start"]),
                  End=as.integer(i["End"]))
  
  # get nucleotide sequence and label with positions
  nt.seq <- i$Sequence %>%
            strsplit("") %>%
            unlist()
  names(nt.seq) <- seq(i$Start, i$End, by=1)
  
  # extract stutter plus bookends
  pos <- seq(a1-bookend, a7+bookend, by=1)
  i$Span <- nt.seq[as.character(pos)] %>%
            paste(collapse="")
  i$Stutter <- str_extract(i$Span, "(?<=AAT).+(?=..TCA)")
  
  # return
  return(i)
}

## load data -------------------------------------------------------------------
# metadata
meta <- read.csv("data/metadata.csv") %>%
        filter(DPI %in% c("Stock", "7")) %>%
        mutate(Group=factor(DPI, levels=c("Stock", "7"),
                            labels=names(cols.grp)),
               Label=factor(NHP, levels=c("Stock", "A", "B", "C", "D"),
                            labels=c("Stock", "NHP A", "NHP B", 
                                     "NHP C", "NHP D"))) %>%
        select(ID, Label, Group, Genome.copies)

# load and format SNVs
# TMP!!! LOAD STOCK SEPARATELY AS PLACEHOLDER
snvs <- read.delim("data/snvs/sample18.vcf",
                   comment.char="#",
                   col.names=c("Segment", "Position", "ID", "Ref", "Alt", 
                               "QUAL", "FILTER", "INFO")) %>%
        filter(FILTER=="PASS") %>%
        mutate(SNV=paste0(Position, "-", Ref, "-", Alt),
               Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[0-9\\.]+")),
               InDel=str_detect(INFO, "INDEL"),
               ID="sample01") %>%
        select(ID, SNV, Position, Frequency, InDel, Ref, Alt)

snvs <- meta$ID[-1] %>%
        lapply(function(i) {
          paste0("data/snvs/", i, ".vcf") %>%
            read.delim(comment.char="#",
                       col.names=c("Segment", "Position", "ID", "Ref", "Alt", 
                                   "QUAL", "FILTER", "INFO")) %>%
            filter(FILTER=="PASS") %>%
            mutate(SNV=paste0(Position, "-", Ref, "-", Alt),
                   Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[0-9\\.]+")),
                   InDel=str_detect(INFO, "INDEL"),
                   ID=i) 
        }) %>%
        do.call(rbind, .) %>%
        select(ID, SNV, Position, Frequency, InDel, Ref, Alt) %>%
        # remove SNVs below illumina threshold (frequency<0.001)
        filter(Frequency > 0.001) %>% 
        ### TMP!! ADD IN FAKE STOCK
        rbind(snvs)

# add in coverage to make profiles
# collapse SNVs to single SNV/nucleotide/sample
x <- snvs %>%
     group_by(ID, Position) %>%
     summarise(Frequency=sum(Frequency),
               .groups="drop")
# load coverage
# TMP!!! LOAD STOCK SEPARATELY AS PLACEHOLDER
prof <- read.delim("data/snvs/sample18.tsv",
                   col.names=c("Segment", "Position", "Depth")) %>%
        mutate(ID="sample01") %>%
        select(ID, Position, Depth)

prof <- meta$ID[-1] %>%
        lapply(function(i) {
          paste0("data/snvs/", i, ".tsv") %>%
            read.delim(col.names=c("Segment", "Position", "Depth")) %>%
            mutate(ID=i)
        }) %>%
        do.call(rbind, .) %>%
        select(ID, Position, Depth) %>%
        # TMP!!!! ADD IN STOCK
        rbind(prof) %>%
        # BACK TO NORMAL
        # add in SNV info
        left_join(x, by=c("ID", "Position")) %>%
        # fill in missing values
        replace_na(list(Frequency=0)) %>%
        # add in metadata
        left_join(meta, by="ID")
rm(x)

# save profile table
prof %>%
  select(ID, Position, Depth, Frequency) %>%
  write.csv("analysis/profile.csv", row.names=FALSE)

## coverage --------------------------------------------------------------------
prof %>%
  ggplot() +
  geom_rect(data=orfs, 
            aes(xmin=Start, xmax=End, 
                ymin=1, ymax=1e4,
                fill=Gene),
            alpha=0.5) +
  scale_fill_brewer(palette="Set3") +
  geom_line(aes(Position, Depth+1)) +
  geom_hline(yintercept=10, linetype=2, col="red") +
  scale_y_continuous(limits=c(1, 1e4), breaks=c(10, 1000), trans="log10") +
  facet_wrap(~Label, ncol=1) +
  labs(x="Nucleotide position",
       y="Aligned read depth")
ggsave("analysis/coverage.png",
       units="cm", width=12, height=12)

# median coverage: remove sample 19
prof %>%
  group_by(Label, ID) %>%
  summarise(MedianDepth=median(Depth),
            .groups="drop")
prof <- filter(prof, ID!="sample19")
snvs <- filter(snvs, ID!="sample19")

## SNV profile -----------------------------------------------------------------
prof %>%
  ggplot() +
  geom_rect(data=orfs, 
            aes(xmin=Start, xmax=End, 
                ymin=0, ymax=1,
                fill=Gene),
            alpha=0.5) +
  scale_fill_brewer(palette="Set3") +
  geom_line(aes(Position, Frequency)) +
  geom_hline(yintercept=0.5, linetype=2, col="red") +
  facet_wrap(~Label, ncol=1) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  labs(x="Nucleotide position",
       y="SNV frequency")
ggsave("analysis/snvs-trace.png",
       units="cm", width=12, height=10)

# clean up
rm(prof)

## annotate SNVs ---------------------------------------------------------------
# get set of unique SNVs
unique.set <- snvs %>%
              select(SNV, Position, Ref, Alt) %>%
              distinct()

# get set of all nucleotides within ORFs
coding <- apply(orfs, 1, function(i){ 
            seq(as.integer(i["Start"]), 
                as.integer(i["End"]), 
                by=1)
          }) %>%
          unlist()

# subset noncoding SNVs
noncod <- unique.set %>%
          filter(!(Position %in% coding)) %>%
          mutate(Type="Noncoding",
                 Label=paste("Noncoding", format(Position, big.mark=","),
                             Ref, ">", Alt))
rm(coding)
unique.set <- filter(unique.set, !(SNV %in% noncod$SNV))

# subset coding SNVs that are *not* in GP after stutter (6915-8060)
# and are not deletions (insertions have the same handling as substitutions)
coding <- unique.set %>%
          filter(!(Position %in% 6915:8060)) %>%
          apply(1, annot.coding, anno=orfs) %>%
          plyr::join_all(type="full") 

# subset coding SNVs that *are* in GP after stutter (6915-8060)
# make a "fake" GP ORF range that starts 1 position earlier so that the 
# number of nucleotides is correct
x <- orfs
x[x$Gene=="GP", "Start"] <- x[x$Gene=="GP", "Start"]-1
coding <- unique.set %>%
          filter(Position %in% 6915:8060) %>%
          apply(1, annot.coding, anno=x) %>%
          plyr::join_all(type="full") %>%
          # add in the rest of the coding mutations
          rbind(coding)
rm(x)

# combine coding and noncoding, then add to SNV frequency table
snvs <- full_join(noncod, coding) %>%
        right_join(snvs) %>%
        # rename nucleotide change columns
        dplyr::rename(NT.ref=Ref,
                      NT.alt=Alt,
                      NT.position=Position) %>%
        # re-order columns
        select(ID, SNV, Label, Frequency, NT.position, NT.ref, NT.alt, 
               Type, Gene, AA.position, AA.ref, AA.alt) %>% 
        # order by position
        arrange(ID, NT.position)
# save full SNV table
write.csv(snvs, "analysis/snvs.csv", row.names=FALSE, na="")

# clean up
rm(coding, noncod, unique.set)

# major SNVs
major.snvs <- snvs %>% 
              filter(Frequency > 0.5) %>%
              group_by(SNV, Label, NT.position) %>%
              summarise(Count=n(),
                        .groups="drop") %>%
              dplyr::arrange(NT.position) %>%
              select(Label, SNV, Count)
major.snvs 
snvs %>%
  filter(SNV %in% major.snvs$SNV) %>%
  dplyr::rename(AA.ID=Label) %>%
  mutate(AA.ID=factor(AA.ID, levels=major.snvs$Label)) %>%
  left_join(meta, by="ID") %>%
  ggplot(aes(Label, Frequency, fill=Group)) +
  geom_col(col="black") +
  scale_fill_manual(values=cols.grp) +
  geom_text(aes(label=format(Frequency, digits=2)), nudge_y=0.1) +
  facet_wrap(~AA.ID) +
  scale_y_continuous(limits=c(0, 1.2), breaks=c(0, 0.5, 1)) +
  labs(x=element_blank(),
       y="SNV frequency",
       fill="Source") 
ggsave("analysis/snvs-major.png",
       units="cm", width=13, height=10)

# hilighter plot
snvs %>%
  filter(SNV %in% major.snvs$SNV) %>%
  dplyr::rename(Name=Label) %>%
  left_join(meta, by="ID") %>%
  select(Label, Name, NT.position) %>%
  # set y-axis heights
  mutate(Label=droplevels(Label),
         y=as.integer(Label),
         Name=factor(Name, levels=major.snvs$Label)) %>%
  ggplot() +
  # lines per genome
  geom_segment(aes(x=1, y=y, xend=19043, yend=y), linewidth=1) +
  # vertical lines per mutation
  geom_segment(aes(x=NT.position, y=y-0.2, 
                   xend=NT.position, yend=y+0.2, 
                   col=Name),
               linewidth=1) +
  scale_color_brewer(palette="Dark2") +
  geom_text(aes(label=Label, x=19043/2, y=y-0.5)) +
  labs(x=element_blank(),
       y=element_blank(),
       col="SNV") +
  scale_y_reverse() +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
ggsave("analysis/hilighter.png",
       units="cm", width=12, height=8)
rm(major.snvs)

## GP stutter site distribution ------------------------------------------------
#TMP TMP TMP! LOADING STOCK PLACEHOLDER
spanning.reads <- load.spanning("sample18") %>%
                  mutate(ID="sample01")
spanning.reads <- meta$ID[-1] %>%
                  lapply(load.spanning) %>% 
                  do.call(rbind, .) %>%
                  rbind(spanning.reads)

# get unique reads only
unique.reads <- spanning.reads %>%
                select(seq, Start, End) %>%
                distinct()

# extract stutter site and annotate to get frame shift and final protein product
annot.span <- unique.reads %>%
              apply(1, annot.spanning) %>%
              do.call(rbind, .) %>%
              mutate(Frame=nchar(Stutter)-7, # indels relative to 7U
                     Frameshift=Frame %% 3, # frameshift relative to 7U
                     Product=factor(Frameshift, levels=c(0, 1, 2, -1), 
                                    labels=names(cols.prd))) %>% 
              replace_na(list(Product="Other"))

# link extracted stutter site to the number of reads per sequence per sample
spanning.reads <- spanning.reads %>%
                  dplyr::rename(Sequence=seq) %>%
                  left_join(annot.span, by=c("Sequence", "Start", "End")) 

# calculate percent per span (per sample) and remove those below
# Illumina sequencing threshold (≤0.1%)
spanning.reads <- spanning.reads %>%
                  group_by(ID, Span) %>%
                  summarise(Reads=sum(Reads),
                            .groups="drop") %>%
                  group_by(ID) %>%
                  mutate(Percent=100*Reads/sum(Reads)) %>%
                  ungroup() %>%
                  filter(Percent > 0.1) %>%
                  select(ID, Span) %>%
                  left_join(spanning.reads, by=c("ID", "Span")) %>%
                  # format for saving
                  select(ID, Sequence, Span, Stutter, Frame, Frameshift,
                         Product, Reads) %>%
                  arrange(ID, desc(Reads))
write.csv(spanning.reads, "analysis/gp-spanning.csv",
          row.names=FALSE, na="")

# summarize to get reads per product per sample
gp.prods <- spanning.reads %>%
            group_by(ID, Product) %>%
            summarise(Reads=sum(Reads),
                      .groups="drop") %>%
            # calculate percentages
            group_by(ID) %>%
            mutate(Percent=100*Reads/sum(Reads)) %>%
            ungroup() %>%
            left_join(meta, by="ID")

# how many representative reads per sample?
gp.prods %>%
  group_by(ID) %>%
  summarise(Reads=sum(Reads),
            .groups="drop") %>%
  left_join(meta, by="ID") %>%
  select(Label, Reads)

# plot by relative abundance
gp.prods %>%
  mutate(Product=factor(Product, levels=rev(names(cols.prd)))) %>%
  ggplot(aes(Label, Percent, fill=Product)) +
  geom_col(col="black", position="Stack") +
  geom_text(data=filter(gp.prods, Product=="sGP"),
            aes(label=paste0(round(Percent), "%"), y=50),
            col="white") +
  scale_fill_manual(values=cols.prd) +
  # add asterisk to NHP D since it has few reads
  geom_text(aes(x="NHP D", y=105, label="*"), 
            col="red", size=10) +
  labs(x=element_blank(),
       y="Relative abundance (%)")
ggsave("analysis/gp-spanning.png",
       units="cm", width=10, height=10)
