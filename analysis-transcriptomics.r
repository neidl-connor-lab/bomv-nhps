#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables
cols.batch <- c(A="black", B="darkgrey")
cols.dpi <- c("-8"="white", "1"="#ffffcc", "3"="#fed976",
              "5"="#fd8d3c", "7"="#fc4e2a", "10"="#e31a1c",
              "15"="#bd0026", "21"="darkgrey", "28"="lightgrey")
cols.sympt <- c(Baseline="white", Asymptomatic="#ffffcc",
                Symptomatic="#e31a1c", Recovered="lightgrey")
cols.reg <- c(Down="#377eb8", None="lightgrey", Up="#e41a1c")
vgenes <- c("NP", "VP35", "VP40", "GP", "VP30", "VP24", "L")

# helper functions
read.counts <- function(counts.matrix) {
  read.delim(counts.matrix) %>%
    # collapse genes by name where possible
    reshape2::melt(id.vars="Gene",
                   variable.name="Sample",
                   value.name="Counts") %>%
    group_by(Gene, Sample) %>%
    summarise(Counts=sum(Counts),
              .groups="drop") %>%
    reshape2::dcast(Gene ~ Sample, value.var="Counts") %>%
    # remove any rRNA since it should have been depleted
    filter(!str_detect(Gene, "rRNA")) %>%
    # format as matrix
    column_to_rownames("Gene") %>%
    as.matrix()
}
filter.lowcounts <- function(counts.matrix, threshold=10) {
  # calculate total gene counts
  x <- rowSums(counts.matrix)
  # filter to get a vector of keep/discard
  x <- (x > threshold)
  # filter data set and return
  counts.matrix[x, ]
}
run.pca <- function(dds) {
  # run PCA
  p <- dds %>%
       varianceStabilizingTransformation(blind=FALSE) %>%
       assay() %>%
       t() %>%
       prcomp()
  # pull out PCs
  x <- summary(p)$importance["Proportion of Variance", c("PC1", "PC2")]
  x <- round(100* x)
  x <- paste0(c("PC1 (", "PC2 ("), x, "%)")
  # format PCA matrix
  p <- p$x %>%
       as.data.frame() %>%
       rownames_to_column("ID") %>%
       select(ID, PC1, PC2) %>%
       left_join(as.data.frame(colData(dds)))
  # return list with PCA matrix and PCs
  list(PCA=p,
       PCs=x)
}
run.results <- function(dds, name, label, sig.p=0.05, sig.lfc=2) {
  # calculate results 
  r <- results(dds, name=name, cooksCutoff=FALSE)
  # shrink lfc
  r <- lfcShrink(dds, coef=name, res=r, type="apeglm", quiet=TRUE)
  # format results
  r <- r %>%
       as.data.frame() %>%
       rownames_to_column("Gene") %>%
       dplyr::rename(log2FC=log2FoldChange) %>%
       mutate(psig=(!is.na(padj) & padj < sig.p),
              fsig=(abs(log2FC) > sig.lfc),
              Significant=(psig & fsig),
              Condition=label)
  # add regulation column
  r$Regulation <- "None"
  r[r$Significant & r$log2FC > 0, "Regulation"] <- "Up"
  r[r$Significant & r$log2FC < 0, "Regulation"] <- "Down"
  r$Regulation <- factor(r$Regulation, levels=c("Up", "None", "Down"))
  # return results
  select(r, Gene, Condition, log2FC, padj, Significant, Regulation)
}


## cyno or rhesus annotation? --------------------------------------------------
cyno <- read.delim("analysis/quantification-cyno.tsv") %>%
        reshape2::melt(id.vars="Status",
                       variable.name="ID",
                       value.name="Reads") %>%
        mutate(Species="Cyno")
rhesus <- read.delim("analysis/quantification-rhesus.tsv") %>%
          reshape2::melt(id.vars="Status", 
                         variable.name="ID",
                         value.name="Reads") %>%
          mutate(Species="Rhesus")

# compare quantification
cyno %>%
  full_join(rhesus, by=colnames(cyno)) %>%
      reshape2::dcast(ID + Species ~ Status, value.var="Reads") %>%
      mutate(TotalReads=Assigned+Unassigned_Unmapped+Unassigned_NoFeatures,
             Percent=100*Assigned/TotalReads) %>%
      select(ID, Species, Assigned, Percent) %>%
  ggplot(aes(Species, Percent)) +
  geom_hline(yintercept=75, col="red", linetype=2) +
  geom_violin(fill="darkgrey", alpha=0.8) +
  geom_jitter(height=0, width=0.1) +
  ggpubr::stat_compare_means(comparisons=list(c("Cyno", "Rhesus")), 
                             paired=TRUE, label="p.signif") +
  ylim(0, 100) +
  labs(x=element_blank(),
       y="Read quantification (%)")
ggsave("analysis/quantification.png",
       units="cm", width=10, height=10)
rm(cyno, rhesus)

## inputs ----------------------------------------------------------------------
# metadata
meta <- read.csv("analysis/samplesheet.csv") %>%
        filter(ID != "sample01") %>%
        mutate(DPI=as.integer(DPI),
               Group=factor(Group, levels=unique(Group)),
               Condition=factor(DPI, levels=unique(DPI)))
rownames(meta) <- meta$ID

# counts 
cmat <- read.counts("data/featurecounts/cyno.tsv") 
mode(cmat) <- "integer"

# align samples
x <- intersect(rownames(meta), colnames(cmat))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

## PCA without batch correction ------------------------------------------------
pca <- cmat %>%
       filter.lowcounts() %>%
       DESeqDataSetFromMatrix(meta, ~Batch) %>%
       run.pca()
# plot by DPI
pca$PCA %>%
  ggplot(aes(PC1, PC2, fill=Batch)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(pch=21, size=5) +
  scale_fill_manual(values=cols.batch) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       title="PCA before batch regression") +
  theme(legend.box.background=element_rect(color="black"),
        legend.position=c(0.9, 0.8))
ggsave("analysis/pca-raw.png",
       units="cm", height=10, width=10)
rm(pca)

# removing genes only expressed in 1 batch
# what genes would we keep if only had batch A samples?
batchA <- meta[meta$Batch=="A", "ID"]
batchA <- cmat[, batchA] %>% 
          filter.lowcounts() %>%
          rownames()
# what genes would we keep if only had batch B samples?
batchB <- meta[meta$Batch=="B", "ID"]
batchB <- cmat[, batchB] %>%
          filter.lowcounts() %>%
          rownames()
# find intersection and subset counts
keepgenes <- intersect(batchA, batchB)
cmat <- cmat[keepgenes, ]
# clean up
rm(keepgenes, batchA, batchB)

## batch correction and PCA to check -------------------------------------------
# batch correction with Com-Bat-seq; can't include DPI due to confounding
cmat <- sva::ComBat_seq(cmat, meta$Batch, group=meta$Group,
                        covar_mod=select(meta, NHP), full_mod=TRUE)
cmat %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  write.csv("analysis/counts-batchcorrected.csv",
            row.names=FALSE)

# repeat PCA
pca <- cmat %>%
       filter.lowcounts() %>%
       DESeqDataSetFromMatrix(meta, ~Condition) %>%
       run.pca()
# color by batch
pca$PCA %>%
  ggplot(aes(PC1, PC2, fill=Batch)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(pch=21, size=5) +
  scale_fill_manual(values=cols.batch) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       title="PCA after batch regression") +
  theme(legend.box.background=element_rect(color="black"),
        legend.position=c(0.9, 0.8))
ggsave("analysis/pca-regressed.png",
       units="cm", width=10, height=10)
# color by DPI
pca$PCA %>%
  ggplot(aes(PC1, PC2, fill=Condition)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(pch=21, size=5) +
  scale_fill_manual(values=cols.dpi) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       fill="DPI")
ggsave("analysis/pca-dpi.png",
       units="cm", width=15, height=10)
# color by symptoms
pca$PCA %>%
  ggplot(aes(PC1, PC2, fill=Group)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(pch=21, size=5) +
  scale_fill_manual(values=cols.sympt) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       fill="DPI")
ggsave("analysis/pca-symptoms.png",
       units="cm", width=15, height=10)
rm(pca)

## DE analysis: DPI ------------------------------------------------------------
# run DESeq and calculate results
desq <- cmat %>%
        filter.lowcounts() %>%
        DESeqDataSetFromMatrix(meta, ~Condition) %>%
        DESeq(quiet=TRUE)
rmat <- list(run.results(desq, "Condition_1_vs_.8", 1),
             run.results(desq, "Condition_3_vs_.8", 3),
             run.results(desq, "Condition_5_vs_.8", 5),
             run.results(desq, "Condition_7_vs_.8", 7),
             run.results(desq, "Condition_10_vs_.8", 10),
             run.results(desq, "Condition_15_vs_.8", 15),
             run.results(desq, "Condition_21_vs_.8", 21),
             run.results(desq, "Condition_28_vs_.8", 28)) %>%
        do.call(rbind, .)

# save results, plus make IPA input matrix
write.csv(rmat, "analysis/results-dpi.csv", na="", row.names=FALSE)
rmat %>%
  reshape2::melt(id.vars=c("Gene", "Condition"),
                 measure.vars=c("log2FC", "padj")) %>%
  mutate(Condition=paste0("dpi", Condition, ".", variable)) %>%
  reshape2::dcast(Gene ~ Condition, value.var="value") %>%
  write.csv("analysis/ipa-input-dpi.csv", na="", row.names=FALSE)

# "squashed" volcano plot
r <- rmat %>%
     na.omit() %>%
     mutate(Condition=factor(paste(Condition, "DPI"),
                             levels=paste(unique(Condition), "DPI")))
r[r$padj < 1e-10, "padj"] <- 1e-10
topgenes <- r %>%
            filter(Significant) %>%
            group_by(Condition) %>%
            top_n(n=10, wt=abs(log2FC)) %>%
            ungroup()
r %>%
  ggplot(aes(log2FC, -log10(padj), size=Regulation, col=Regulation)) +
  geom_point(alpha=0.8) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           col="black", size=2, max.overlaps=100) +
  scale_size_manual(values=c(Down=1, None=0.5, Up=1)) +
  scale_color_manual(values=cols.reg) +
  scale_y_continuous(limits=c(NA, 11),
                     breaks=c(1, 5, 10),
                     labels=c("1", "1e-5", "<1e-10")) +
  scale_x_continuous(limits=c(-10, 20),
                     breaks=c(-10, 0, 10, 20)) +
  facet_wrap(~Condition, nrow=2) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value (-log10)")
ggsave("analysis/volcano-dpi.png",
       units="cm", width=20, height=10)
rm(r, topgenes)

# total degenes
rmat %>%
  filter(Significant) %>%
  group_by(Regulation, Condition) %>%
  summarise(Genes=n(),
            .groups="drop") %>%
  # add in any missing values
  right_join(expand.grid(Condition=as.numeric(names(cols.dpi)),
                         Regulation=c("Down", "Up")),
             by=c("Regulation", "Condition")) %>%
  replace_na(list(Genes=0)) %>%
  # plot it
  ggplot(aes(Condition, Genes, group=Regulation, fill=Regulation)) +
  geom_line(aes(linetype=Regulation)) +
  geom_point(pch=21, size=3) +
  scale_fill_manual(values=cols.reg) +
  labs(x="Days postinfection",
       y="Significantly DE genes") +
  theme(legend.box.background=element_rect(color="black"),
        legend.position=c(0.8, 0.8))
ggsave("analysis/degenes-dpi.png",
       units="cm", width=10, height=8)

# define up-regulated gene modules
# get list of all up-regulated genes
gene.mod <- rmat[rmat$Regulation=="Up", "Gene"]
gene.mod <- rmat %>%
            filter(Gene %in% gene.mod,
                   !(Gene %in% vgenes)) %>%
            # format DPI
            mutate(DPI=paste0(Condition, ".DPI")) %>%
            # expand to matrix with DPI as columns and lfc as values
            reshape2::dcast(Gene ~ DPI, value.var="log2FC") %>%
            column_to_rownames("Gene") %>%
            as.matrix()

# how many clusters? Using 3
lapply(1:10, function(i) {
  data.frame(k=i,
             tot.withinss=kmeans(gene.mod, i, iter.max=1e9)$tot.withinss)
}) %>%
  do.call(rbind, .) %>%
  ggplot(aes(k, tot.withinss)) +
  geom_line() +
  geom_vline(xintercept=3, col="red", linetype=2) +
  labs(x="# clusters",
       y="Total within sum of squares",
       title="K-means elbow plot")

# add cluster to data set
gene.mod <- kmeans(gene.mod, 3, iter.max=1e9)$cluster
gene.mod <- data.frame(Gene=names(gene.mod),
                       Cluster=gene.mod) %>%
            left_join(filter(rmat, Gene %in% .$Gene),
                      by="Gene") %>%
            select(Gene, Cluster, Condition, log2FC)
# add in 0 DPI
gene.mod <- gene.mod %>%
            select(Gene, Cluster) %>%
            distinct() %>%
            mutate(Condition=-8,
                   log2FC=0) %>%
            rbind(gene.mod)
# format gene module names
# first get total genes per module
x <- gene.mod %>%
     select(Gene, Cluster) %>%
     distinct() %>%
     group_by(Cluster) %>%
     summarise(Genes=n(),
               .groups="drop")
anchor.genes <- c("Early"="IFI44",
                  "Late"="LCN2",
                  "Both"="OAS1")
gene.mod <- gene.mod %>%
            filter(Gene %in% anchor.genes) %>%
            group_by(Cluster, Gene) %>%
            summarise(.groups="drop") %>%
            left_join(x, by="Cluster") %>%
            mutate(Gene=factor(Gene, 
                              levels=anchor.genes, 
                              labels=names(anchor.genes)),
                  Label=paste0(Gene, " (n=", Genes, ")")) %>%
            select(Cluster, Label) %>%
            right_join(gene.mod, by="Cluster")
  
# get some examples per cluster
gene.examp <- c("OAS1", "CXCL10", "SOCS1",
                "IFI44", "CCL23", "CCL8",
                "LCN2", "MMP8")
gene.examp <- gene.mod %>%
              filter(Gene %in% gene.examp)
    
# group by cluster and plot
gene.mod %>%
  group_by(Condition, Label) %>%
  summarise(StDev=sd(log2FC),
            log2FC=median(log2FC),
            .groups="drop") %>%
  ggplot(aes(Condition)) +
  geom_ribbon(aes(ymin=log2FC-StDev,
                  ymax=log2FC+StDev),
              fill="lightgrey",
              alpha=0.5) +
  geom_line(data=gene.examp, aes(y=log2FC, col=Gene)) +
  scale_color_brewer(palette="Paired") +
  facet_wrap(~Label, ncol=1) +
  labs(x="Days postinfection",
       y="Fold change (log2)")
ggsave("analysis/modules-dpi.png", 
       units="cm", width=10, height=10)

# clean up
rm(x, gene.mod, gene.examp, anchor.genes, rmat, desq)

## DE analysis: symptoms -------------------------------------------------------
# run DESeq and calculate results
desq <- cmat %>%
        filter.lowcounts() %>%
        DESeqDataSetFromMatrix(meta, ~Group) %>%
        DESeq(quiet=TRUE)
rmat <- list(run.results(desq, "Group_Asymptomatic_vs_Baseline", "Asymptomatic"),
             run.results(desq, "Group_Symptomatic_vs_Baseline", "Symptomatic"),
             run.results(desq, "Group_Recovered_vs_Baseline", "Recovered")) %>%
        do.call(rbind, .) %>%
        mutate(Condition=factor(Condition, levels=levels(meta$Group)))

# save results, plus make IPA input matrix
write.csv(rmat, "analysis/results-symptoms.csv", na="", row.names=FALSE)
rmat %>%
  reshape2::melt(id.vars=c("Gene", "Condition"),
                 measure.vars=c("log2FC", "padj")) %>%
  mutate(Condition=paste0("dpi", Condition, ".", variable)) %>%
  reshape2::dcast(Gene ~ Condition, value.var="value") %>%
  write.csv("analysis/ipa-input-symptoms.csv", na="", row.names=FALSE)

# "squashed" volcano plot
r <- na.omit(rmat)
r[r$padj < 1e-10, "padj"] <- 1e-10
topgenes <- r %>%
            filter(Significant) %>%
            group_by(Condition) %>%
            top_n(n=15, wt=abs(log2FC)) %>%
            ungroup()
r %>%
  ggplot(aes(log2FC, -log10(padj), size=Regulation, col=Regulation)) +
  geom_point(alpha=0.8) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           col="black", size=2, max.overlaps=100) +
  scale_size_manual(values=c(Down=1, None=0.5, Up=1)) +
  scale_color_manual(values=cols.reg) +
  scale_y_continuous(limits=c(NA, 11),
                     breaks=c(1, 5, 10),
                     labels=c("1", "1e-5", "<1e-10")) +
  scale_x_continuous(limits=c(-10, 20),
                     breaks=c(-10, 0, 10, 20)) +
  facet_wrap(~Condition, nrow=1) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value (-log10)")
ggsave("analysis/volcano-symptoms.png",
       units="cm", width=20, height=10)
rm(r, topgenes)

# total degenes
rmat %>%
  filter(Significant) %>%
  group_by(Regulation, Condition) %>%
  summarise(Genes=n(),
            .groups="drop") %>%
  # add in any missing values
  right_join(expand.grid(Condition=levels(meta$Group),
                         Regulation=c("Down", "Up")),
             by=c("Regulation", "Condition")) %>%
  replace_na(list(Genes=0)) %>%
  # plot it
  ggplot(aes(Condition, Genes, group=Regulation, fill=Regulation)) +
  geom_line(aes(linetype=Regulation)) +
  geom_point(pch=21, size=3) +
  scale_fill_manual(values=cols.reg) +
  labs(x=element_blank(),
       y="Significantly DE genes") +
  theme(legend.box.background=element_rect(color="black"),
        legend.position=c(0.2, 0.8))
ggsave("analysis/degenes-symptoms.png",
       units="cm", width=10, height=8)
