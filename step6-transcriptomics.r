#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
options(stringsAsFactors=FALSE)
theme_set(theme_classic() + 
          theme(strip.text=element_text(color="white"),
                strip.background=element_rect(fill="black"),
                legend.box.background=element_rect(color="black", 
                                                   linewidth=1)))

# colors
cols.dpi <- c("-8"="white", "1"="#ffffb2", "3"="#fecc5c", 
              "5"="#fd8d3c", "7"="#e31a1c", 
              "10"="#bd0026", "15"="#800026", 
              "21"="#525252", "28"="#969696")
cols.reg <- c(Up="#e41a1c", None="darkgrey", Down="#377eb8")

# other global variables
markers <- c(IFNB1="Interferon β (IFNB1)", OAS1="ISGs (OAS1)", 
             CXCL10="CXCL10/IP-10", IL6="IL-6", 
             ELANE="Neutrophil elasetase (ELANE)", 
             S100A8="Calprotectin (S100A8)",
             CCL7="CCL7/MCP3",CXCL8="CXCL8/IL-8")

# helper functions
calculate.pca <- function(v) {
  # run PCA
  x <- assay(v) %>%
         t() %>%
         prcomp()
  # pull out PC % variance explained
  y <- summary(x)
  y <- y$importance["Proportion of Variance", c("PC1", "PC2")]
  y <- round(100*y)
  y <- c(paste0("PC1 (", y[1], "%)"),
         paste0("PC2 (", y[2], "%)"))
  # format PCA
  x <- x$x %>%
       as.data.frame() %>%
       rownames_to_column("ID") %>%
       select(ID, PC1, PC2) %>%
       left_join(as.data.frame(colData(v)), by="ID")
  # return PCA info
  list(PCA=x,
       PCs=y)
}

## inputs ----------------------------------------------------------------------
# metadata (sans stock)
meta <- read.csv("data/metadata.csv") %>%
        filter(Group!="Stock") %>%
        mutate(DPI=factor(DPI, levels=names(cols.dpi)))
rownames(meta) <- meta$ID

# counts 
counts <- read.delim("data/featurecounts/rhesus.tsv",
                     comment.char="#", na.strings="") %>%
          dplyr::rename(Gene=Geneid) %>%
          select(Chr, Gene, ends_with("bam")) %>%
          # remove unnamed genes
          na.omit()
# consolidate viral genes
bomv <- counts %>%
        filter(Chr=="MF319186.1") %>%
        select(-Chr, -Gene) %>%
        colSums() %>%
        data.frame(BOMV=.) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Gene")
counts <- counts %>%
          filter(Chr!="MF319186.1") %>%
          select(-Chr) %>%
          full_join(bomv, by=colnames(bomv)) %>%
          column_to_rownames("Gene")
rm(bomv)
# update column names
colnames(counts) <- colnames(counts) %>%
                    str_extract("sample[0-9]+")

# all samples pass qc (see multiqc for rhesus)
# align column
x <- intersect(meta$ID, colnames(counts))
meta <- meta[x, ]
counts <- counts[ , x]
rm(x)

# remove genes with <10 counts total
x <- counts %>%
     rowSums() %>%
     data.frame(Counts=.) %>% 
     rownames_to_column("Gene") %>%
     filter(Counts > 10) %>%
     select(Gene) %>%
     unlist()
counts <- counts[x, ]
rm(x)

# save raw counts
write.csv(counts, "analysis/counts-raw.csv")

## sampling & PCA --------------------------------------------------------------
# make DESeq object then calculate PCA
dds <- DESeqDataSetFromMatrix(counts, meta, ~DPI)
pca <- dds %>%
       varianceStabilizingTransformation(blind=FALSE) %>%
       calculate.pca()

# PCA by DPI
pca$PCA %>%
  ggplot(aes(PC1, PC2, fill=DPI)) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_point(pch=21, size=5) +
  scale_fill_manual(values=cols.dpi) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       fill="Timepoint") 
ggsave("analysis/pca.png",
       units="cm", width=15, height=10)

# save PCA
pca$PCA %>%
  select(ID, DPI, Group, PC1, PC2) %>%
  write.csv("analysis/pca.csv", row.names=FALSE)
rm(pca)

## DE by DPI -------------------------------------------------------------------
# run DE and get results
dds <- DESeq(dds, quiet=TRUE)
res <- lapply(c(1, 5, 7, 15, 21, 28), 
            function(i) {
              lfcShrink(dds, coef=paste0("DPI_", i, "_vs_.8"), 
                        type="apeglm", quiet=TRUE) %>%
                as.data.frame() %>%
                rownames_to_column("Gene") %>%
                mutate(DPI=as.integer(i))
            }) %>%
       do.call(rbind, .) %>%
       mutate(psig=(!is.na(padj) & padj < 0.05),
             fsig=(abs(log2FoldChange) > 2),
             Significant=(psig & fsig),
             Regulation="None") %>%
       select(-psig, -fsig)
res[res$Significant & res$log2FoldChange < 0, "Regulation"] <- "Down"
res[res$Significant & res$log2FoldChange > 0, "Regulation"] <- "Up"
res$Regulation <- factor(res$Regulation, levels=names(cols.reg))

# write results
write.csv(res, "analysis/deseq.csv", row.names=FALSE, na="")

# format results for IPA
res %>%
  mutate(DPI=factor(DPI, levels=c(1, 5, 7, 15, 21, 28),
                    labels=c("DPI01", "DPI05", "DPI07", 
                             "DPI15", "DPI21", "DPI28"))) %>%
  reshape2::melt(id.vars=c("Gene", "DPI"),
                 measure.vars=c("padj", "log2FoldChange")) %>%
  mutate(variable=factor(variable, levels=c("padj", "log2FoldChange"),
                         labels=c("padj", "lfc")),
         label=paste0(DPI, ".", variable)) %>%
  reshape2::dcast(Gene ~ label, value.var="value") %>%
  write.csv("analysis/ipa-input.csv", row.names=FALSE)

# total DE genes
totals <- res %>% 
          filter(Significant) %>%
          group_by(DPI, Regulation) %>%
          summarise(Genes=n(),
                    .groups="drop") %>%
          # fill in zeros
          right_join(expand.grid(Regulation=c("Down", "Up"),
                                 DPI=c(0, unique(res$DPI))),
                     by=c("Regulation", "DPI")) %>%
          replace_na(list(Genes=0))
totals %>%
  reshape2::dcast(DPI ~ Regulation, value.var="Genes") %>%
  write.csv("analysis/degenes.csv", row.names=FALSE)
totals %>%
  ggplot(aes(DPI, Genes, group=Regulation)) +
  geom_line(aes(linetype=Regulation)) +
  geom_point(aes(fill=Regulation), pch=21, size=5) +
  scale_fill_manual(values=cols.reg) +
  scale_linetype_manual(values=c(Up=1, Down=2)) +
  labs(x="Days postchallenge",
       y="Significantly DE genes") +
  theme(legend.position=c(0.8, 0.8))
ggsave("analysis/degenes.png",
       units="cm", width=15, height=10)
rm(totals)

# volcano plots: read in genes to label
x <- read.csv("analysis/degene-shortlist.csv") %>%
      mutate(Label=Gene) %>%
      full_join(res, by=c("Gene", "DPI")) %>%
      mutate(DPI=factor(DPI, levels=unique(DPI),
                        labels=paste(unique(DPI), "DPI")),
             y=-log10(padj),
             labmod=log2FoldChange/abs(log2FoldChange))
x %>%
  ggplot(aes(log2FoldChange, y)) +
  geom_point(aes(size=Regulation, fill=Regulation), 
             pch=21, na.rm=TRUE, alpha=0.8) +
  scale_size_manual(values=c(Down=1.5, None=1, Up=1.5)) +
  scale_fill_manual(values=cols.reg) +
  ggrepel::geom_text_repel(aes(label=Label),
                           nudge_x=x$labmod*5,
                           #nudge_y=2,
                           force = 50,
                           size=2, max.time=20, na.rm=TRUE) +
  facet_wrap(~DPI, nrow=1) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value (-log10)") +
  theme(legend.position=c(0.08, 0.8))
ggsave("analysis/volcano.png", 
       units="cm", width=20, height=10)

## marker genes ----------------------------------------------------------------
names(markers) %>%
  lapply(function(i) {
    plotCounts(dds, i, intgroup="DPI", returnData=TRUE) %>%
      rownames_to_column("ID") %>%
      mutate(Gene=i)
  }) %>%
  do.call(rbind, .) %>%
  mutate(Label=factor(Gene, levels=names(markers),
                      labels=markers),
         DPI=as.integer(as.character(DPI))) %>%
  left_join(select(res, Gene, DPI, Regulation),
            by=c("Gene", "DPI")) %>% 
  left_join(select(meta, ID, NHP),
            by="ID") %>%
  replace_na(list(Regulation="None")) %>%
  ggplot(aes(DPI, count)) +
  geom_line(aes(linetype=NHP)) +
  geom_point(aes(fill=Regulation), pch=21, size=2) +
  scale_fill_manual(values=cols.reg) +
  scale_y_log10() +
  scale_x_continuous(limits=c(-10, 30)) +
  facet_wrap(~Label, ncol=2) +
  labs(x="Days post-challenge",
       y="Normalized counts")
ggsave("analysis/markers.png",
       units="cm", width=20, height=20)
