#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Biostrings))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper functions
mutate.seq <- function(snvinfo, nt.seq) {
  # parse snv info
  nt.ref <- as.character(snvinfo["NT.ref"])
  nt.alt <- as.character(snvinfo["NT.alt"])
  nt.pos <- as.integer(snvinfo["NT.position"])
  
  # is the nt.ref multiple characters? If so, it's a deletion and
  # requires special handling
  if(nchar(nt.ref) > 1) { # deletion
    # how many nucleotides deleted?
    nt.ref <- nchar(nt.ref)-1
    # what positions are these?
    nt.pos <- (nt.pos+1):(nt.pos+1)+nt.ref
    # set these positions to NA, then filter NAs
    nt.seq[nt.seq$NT.position %in% nt.pos, "Nucleotide"] <- NA
    nt.seq <- na.omit(nt.seq)
    
  } else { # substitution or insertion
    nt.seq[nt.seq$NT.position==nt.pos, "Nucleotide"] <- nt.alt
  }
  
  # collapse the nucleotide sequence and compare to the AA sequence
  nt.seq$Nucleotide %>%
    paste(collapse="") %>%
    DNAString()
}

## process arguments and inputs ------------------------------------------------
args <- ArgumentParser()
args$add_argument("-d", "--dir", help="input and output directory", 
                  required=TRUE)
args$add_argument("-a", "--annotation", help="viral GTF",
                  required=TRUE)
args$add_argument("-g", "--genome", help="viral genome FASTA",
                  required=TRUE)
args <- args$parse_args()

# check inputs
# annotation file
if(!file.exists(args$annotation)) {
  print("[ERR] Annotation file does not exist.")
  q(save="no", status=1)
}
# genome file
if(!file.exists(args$genome)) {
  print("[ERR] Viral genome FASTA file does not exist.")
  q(save="no", status=1)
}

# input directory
if(!dir.exists(args$dir)) {
  print("[ERR] Input directory does not exist.")
  q(save="no", status=1)
}

# SNV file
snvs <- paste0(args$dir, "/snvs.vcf")
if(!file.exists(snvs)) {
  print("[ERR] SNV VCF does not exist.")
  q(save="no", status=1)
}

## load data -------------------------------------------------------------------
# load SNV VCF
snvs <- read.delim(snvs,
                   comment.char="#",
                   col.names=c("Segment", "NT.position", "ID", "NT.ref", 
                               "NT.alt", "QUAL", "FILTER", "INFO"),
                   colClasses=c(Segment="character",
                                NT.position="integer",
                                NT.ref="character",
                                NT.alt="character")) %>%
        # extract frequency, tag indels, and construct GRange keys
        mutate(Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[^;$]+")),
               InDel=str_detect(INFO, "INDEL"),
               NT.GRange=paste0(Segment, ":", NT.position)) %>%
        # filter with LoFreq filter and Illumina error
        filter(FILTER=="PASS",
               Frequency > 0.001) %>%
        select(-ID, -QUAL, -FILTER, -INFO)

# check: if there aren't any SNVs in the file, quit without an error
if(nrow(snvs)==0) {
  print("[MSG] SNV VCF is empty. Exiting...")
  q(save="no", status=0)
}

# load annotation
anno <- args$annotation %>%
        read.delim(header=FALSE,
                   comment.char="#",
                   col.names=c("Segment", "Source", "Type", "Start", "Stop",
                               "Score", "Strand", "Phase", "Info")) %>%
        # subset to CDS 
        filter(Type=="CDS") %>%
        # pull out gene name, protein ID, and build GRange string
        mutate(Gene.ID=str_extract(Info, "(?<=gene_id )[^;$]+"),
               Gene.name=str_extract(Info, "(?<=gene_name )[^;$]+"),
               Protein.ID=str_extract(Info, "(?<=protein_id )[^;$]+"),
               Product=str_extract(Info, "(?<=product )[^;$]+"),
               NT.GRange=paste0(Segment, ":", Start, "-", Stop, ":+"),
               AA.GRange=paste0(Segment, ":", Start, "-", Stop, ":", Strand))

# load genome and parse segment names to match annotation
genome <- Biostrings::readDNAStringSet(args$genome)
names(genome) <- str_extract(names(genome), "^[^ ]+")

# build nucleotide CDS GRanges object
cds <- GRanges(anno$AA.GRange,
               Protein.ID=anno$Protein.ID,
               Product=anno$Product,
               Gene.ID=anno$Gene.ID,
               Gene.name=anno$Gene.name)


## annotate SNVs ---------------------------------------------------------------
# split into coding vs. noncoding SNVs
coding <- (GRanges(snvs$NT.GRange) %within% cds)
noncoding <- snvs[!coding, ] %>%
             mutate(Type="Noncoding")

# annotate the SNVs that fall within each protein ID
snvs <- cds$Protein.ID %>%
        unique() %>% # there may be multiple CDS per protein
        lapply(function(i) {
          # subset CDS
          c <- cds[cds$Protein.ID==i]
          # subset SNVs
          s <- snvs[(GRanges(snvs$NT.GRange) %within% c), ]
          
          # Sanity check! Are there any SNVs within this protein? If not, return.
          if(dim(s)[1]==0) {
            return()
          }
          
          # To load the NT sequence, update CDS so the strand is always "+".
          # This ensures the nucleotides will be in genome coordinate order.
          nt.ref <- c
          strand(nt.ref) <- "+"
          # load the NT sequence for this protein and format as a char vector
          nt.ref <- BSgenome::getSeq(genome, nt.ref) %>%
                    as.character() %>%
                    strsplit("") %>%
                    unlist()
          # Make a vector of NT positions, stitching together CDS entries.
          # Strand doesn't matter for getting the NT position range, so we
          # can use "c" here. 
          x <- ranges(c) %>% 
               as.data.frame() %>% 
               apply(1, function(i) { 
                 seq(from=i["start"], to=i["end"]) 
                }) %>% 
               unlist()
          # combine NTs and positions into a data frame
          nt.ref <- data.frame(Nucleotide=nt.ref,
                                  NT.position=x)
          
          # Format the reference NT sequence so we can translate it to get a
          # reference AA sequence
          aa.ref <- nt.ref$Nucleotide %>%
                    paste(collapse="") %>%
                    DNAString()
          # Before we translate the reference NT sequence to get the reference
          # AA sequence, we need to make the reverse complement for CDS on the
          # "-" strand
          x <- unique(strand(c)) # using "unique" here to collapse multiple CDS
          if(x == "-") {
            aa.ref <- reverseComplement(aa.ref)
          }
          # Now we can translate and format for easy comparison later
          aa.ref <- aa.ref %>%
                    translate() %>%
                    as.character() %>%
                    strsplit("") %>%
                    unlist()
          
          # loop through each SNV and get the mutated AA sequence
          aa.alt <- s %>%
                    apply(1, mutate.seq, nt.seq=nt.ref) %>%
                    # collapse the list into a DNAStringSet object
                    DNAStringSet()
          # if the protein is on the "-" strand, take the reverse complement
          if(x == "-") {
            aa.alt <- reverseComplement(aa.alt)
          }
          # translate the mutated proteins
          aa.alt <- aa.alt %>%
                    translate() %>%
                    # don't show INDEL warnings where CDS isn't divisible by 3
                    suppressWarnings() %>% 
                    # format for comparison
                    as.character() %>%
                    strsplit("")
          
          # compare each mutated sequence to the AA reference
          s$AA.pos <- aa.alt %>%
                      lapply(function(j){ 
                        # make sure AA.alt and AA.ref are the same length
                        # this is only a concern with INDELs
                        x <- min(length(aa.ref), length(j))
                        j <- j[1:x]
                        a <- aa.ref[1:x]
                        # return only the first difference, if any, to control 
                        # for frame shifts
                        which(j != a)[1] 
                      }) %>%
                      unlist()
          
          # fill in the AA.ref
          s$AA.ref <- aa.ref[s$AA.pos]
          
          # fill in the AA.alt by looping through the AA.alt and SNV data frame
          s$AA.alt <- 1:dim(s)[1] %>%
                      lapply(function(j) { 
                        aseq <- aa.alt[[j]]
                        apos <- s[j, "AA.pos"]
                        return(aseq[apos])
                      }) %>% 
                      unlist()
          
          # SNVs with "NA" in the AA.pos are synonymous
          s$Type <- "Nonsynonymous"
          s[is.na(s$AA.pos), "Type"] <- "Synonymous"
          
          # return the updated SNV data frame with the protein ID
          s$Protein.ID <- i
          return(s)
        }) %>%
        # now that we have all the annotations, bind the rows together
        do.call(rbind, .) 

# use the protein IDs to link gene info
anno %>%
  select(Gene.ID, Gene.name, Protein.ID, Product) %>%
  distinct() %>% # keep only unique sets of gene id, gene name, and protein ID
  right_join(snvs, by="Protein.ID") %>% # add coding SNVs
  # add in the noncoding SNVs
  full_join(noncoding, by=c("Segment", "NT.position", "NT.ref", "NT.alt", 
                            "Frequency", "InDel", "NT.GRange", "Type")) %>%
  # sort by segment and position
  arrange(Segment, NT.position) %>%
  # re-arrange columns
  select(Segment, NT.position, NT.ref, NT.alt, Frequency, Type, InDel,
         Gene.ID, Gene.name, Protein.ID, Product, AA.pos, AA.ref, AA.alt) %>%
  # write to CSV
  write.csv(paste0(args$dir, "/snvs.csv"), 
            na="", row.names=FALSE)
