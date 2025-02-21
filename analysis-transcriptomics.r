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
scalefact <- 1.5

# helper functions
read.counts <- function(counts.matrix, 
                        viral=vgenes) {
  counts.matrix <- read.delim(counts.matrix) %>%
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
  # pull out viral counts and sum per sample; format as matrix
  vcounts <- colSums(counts.matrix[viral, ])
  vcounts <- matrix(vcounts,
                    nrow=1,
                    dimnames=list("BOMV", names(vcounts)))
  # remove viral counts from the original matrix
  counts.matrix <- counts.matrix[!(rownames(counts.matrix) %in% viral), ]
  # add back in viral counts under "BOMV" and return
  rbind(counts.matrix, vcounts)
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

## metadata --------------------------------------------------------------------
meta <- read.csv("samplesheet.csv") %>%
        mutate(DPI=as.integer(DPI),
               Group=factor(Group, levels=unique(Group)),
               Condition=factor(DPI, levels=unique(DPI)))
rownames(meta) <- meta$ID

## quantification QC -----------------------------------------------------------
quanti <- read.delim("data/transcriptomics/counts.tsv.summary") %>%
          reshape2::melt(id.vars="Status",
                         variable.name="ID",
                         value.name="Reads") %>%
          filter(Reads!=0) %>%
          left_join(meta, by="ID") %>%
          mutate(Status=factor(Status, 
                               levels=c("Unassigned_Unmapped",
                                        "Unassigned_NoFeatures", 
                                        "Unassigned_Duplicate",
                                        "Assigned"),
                               labels=c("Unmapped (unassigned)",
                                        "No features (unassigned)",
                                        "PCR duplicates (unassigned)",
                                        "Assigned")),
                 DPI=factor(DPI, levels=names(cols.dpi),
                            labels=paste(names(cols.dpi), "DPI")))
quanti <- quanti$DPI %>%
          unique() %>%
          lapply(function(i) {
            quanti %>%
              filter(DPI==i) %>%
              ggplot(aes(Reads, NHP)) +
              geom_col(aes(fill=Status), col="black") +
              geom_vline(xintercept=1e7, linetype=2, col="lightgrey") +
              scale_fill_manual(values=c("Unmapped (unassigned)"="white",
                                         "No features (unassigned)"="lightgrey",
                                         "PCR duplicates (unassigned)"="darkgrey",
                                         "Assigned"="black")) +
              scale_y_discrete(limits=rev(unique(meta$NHP))) +
              scale_x_continuous(expand=c(0, 0),
                                 limits=c(0, 9.1e7),
                                 breaks=c(0, 1e7, 3e7, 5e7, 7e7, 9e7),
                                 labels=c("0", "1e7", "3e7", "5e7", "7e7", "9e7")) +
              labs(x="Reads",
                   y="BOMV NHPs",
                   fill=element_blank(),
                   title=i) +
              theme(legend.position="none",
                    axis.title.y=element_blank(),
                    axis.title.x=element_blank())
          })
# pull out legend
leg <- cowplot::get_plot_component(quanti[[1]] + theme(legend.position="bottom"),
                                   "guide-box-bottom")
y <- cowplot::get_plot_component(quanti[[1]] + theme(axis.title.y=element_text()),
                                 "ylab-l")
x <- cowplot::get_plot_component(quanti[[1]] + theme(axis.title.x=element_text()),
                                 "xlab-b")
quanti <- cowplot::plot_grid(plotlist=quanti, ncol=3)
quanti <- cowplot::plot_grid(y, quanti, nrow=1, rel_widths=c(1, 50))
supA <- cowplot::plot_grid(quanti, x, leg, 
                           ncol=1, rel_heights=c(15, 1, 1))
supA
ggsave("analysis/quantification.png", scale=scalefact,
       units="in", width=7.5, height=3)
rm(quanti, x, y, leg)

## counts ----------------------------------------------------------------------
cmat <- read.counts("data/transcriptomics/counts.tsv") 
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
supB <- pca$PCA %>%
        ggplot(aes(PC1, PC2, fill=Batch)) +
        geom_hline(yintercept=0, linetype=3) +
        geom_vline(xintercept=0, linetype=3) +
        geom_point(pch=21, size=5) +
        scale_fill_manual(values=cols.batch) +
        labs(x=pca$PCs[1],
             y=pca$PCs[2],
             title="PCA before batch regression") +
        theme(legend.position=c(0.9, 0.8))
supB
ggsave("analysis/pca-raw.png", scale=scalefact,
       units="in", height=2, width=2.33)
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
supC <- pca$PCA %>%
        ggplot(aes(PC1, PC2, fill=Batch)) +
        geom_hline(yintercept=0, linetype=3) +
        geom_vline(xintercept=0, linetype=3) +
        geom_point(pch=21, size=5) +
        scale_fill_manual(values=cols.batch) +
        labs(x=pca$PCs[1],
             y=pca$PCs[2],
             title="PCA after batch regression") +
        theme(legend.position=c(0.9, 0.85))
supC
ggsave("analysis/pca-regressed.png", scale=scalefact,
       units="in", width=2.33, height=2)

# color by DPI
supD <- pca$PCA %>%
        ggplot(aes(PC1, PC2, fill=Condition)) +
        geom_hline(yintercept=0, linetype=3) +
        geom_vline(xintercept=0, linetype=3) +
        geom_point(pch=21, size=5) +
        scale_fill_manual(values=cols.dpi) +
        labs(x=pca$PCs[1],
             y=pca$PCs[2],
             fill="DPI",
             title="PCA after batch regression")
supD
ggsave("analysis/pca-dpi.png", scale=scalefact,
       units="in", width=2.84, height=2)
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
            filter(Significant,
                   !str_detect(Gene, "^ENS")) %>%
            group_by(Condition) %>%
            top_n(n=10, wt=abs(log2FC)) %>%
            ungroup()
vplt <- r$Condition %>%
        unique() %>%
        lapply(function(i) {
          t <- filter(topgenes, Condition==i)
          r %>%
            filter(Condition==i) %>%
            ggplot(aes(log2FC, -log10(padj), size=Regulation, col=Regulation)) +
            geom_point(alpha=0.8) +
            ggrepel::geom_text_repel(data=t, aes(label=Gene), 
                                     col="black", size=3, max.overlaps=100) +
            scale_size_manual(values=c(Down=1, None=0.5, Up=1)) +
            scale_color_manual(values=cols.reg) +
            scale_y_continuous(limits=c(NA, 11),
                               breaks=c(1, 5, 10),
                               labels=c("1", "1e-5", "<1e-10")) +
            scale_x_continuous(limits=c(-10, 20),
                               breaks=c(-10, 0, 10, 20)) +
            labs(x="Fold change (log2)",
                 y="FDR-adjusted p-value (-log10)",
                 title=i) +
            theme(legend.position="none",
                  axis.title.y=element_blank(),
                  axis.title.x=element_blank())
        })
# get legend and axis labels
x <- cowplot::get_plot_component(vplt[[1]] + theme(axis.title.x=element_text()),
                                 "xlab-b")
y <- cowplot::get_plot_component(vplt[[1]] + theme(axis.title.y=element_text()),
                                 "ylab-l")
leg <- cowplot::get_plot_component(vplt[[1]] + theme(legend.position="right"),
                                   "guide-box-right")
# assemble figures
vplt <- cowplot::plot_grid(plotlist=vplt, ncol=4)
vplt <- cowplot::plot_grid(y, vplt, leg, nrow=1, rel_widths=c(1, 50, 4))
supE <- cowplot::plot_grid(vplt, x, ncol=1, rel_heights=c(15, 1))
supE
ggsave("analysis/volcano-dpi.png", scale=scalefact, 
       units="in", width=7.5, height=3)
rm(r, topgenes, vplt, x, y, leg)

# total degenes
figA <- rmat %>%
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
        scale_x_continuous(breaks=c(-8, 1, 3, 5, 7, 10, 15, 21, 28)) +
        labs(x="Days postinfection",
             y="Significantly DE genes") +
        theme(legend.position=c(0.8, 0.8))
figA
ggsave("analysis/degenes-dpi.png", scale=scalefact,
       units="in", width=3.75, height=2)

# define up-regulated gene modules
# get list of all up-regulated genes
gene.mod <- rmat$Gene[rmat$Regulation=="Up"]
gene.mod <- rmat %>%
            filter(Gene %in% gene.mod) %>%
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
                  "Prolonged"="OAS1")
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
glist <- c(OAS1="#6a51a3", CXCL10="#9e9ac8", SOCS1="#810f7c",
           IFI44="#990000", CCL23="#fc8d59", CCL8="#ef6548",
           LCN2="#1f78b4", MMP8="#034e7b")
gene.examp <- gene.mod %>%
              filter(Gene %in% names(glist))
    
# group by cluster and plot
modplt <- gene.mod$Label %>%
          unique() %>%
          sort() %>%
          lapply(function(i) {
            gene.mod %>%
              filter(Label==i) %>%
              group_by(Condition) %>%
              summarise(StDev=sd(log2FC),
                        log2FC=median(log2FC),
                        .groups="drop") %>%
              ggplot(aes(Condition)) +
              geom_ribbon(aes(ymin=log2FC-StDev,
                              ymax=log2FC+StDev),
                          fill="lightgrey",
                          alpha=0.5) +
              geom_line(data=filter(gene.examp, Label==i), 
                        aes(y=log2FC, col=Gene)) +
              scale_color_manual(values=glist) +
              scale_x_continuous(limits=c(-8, 30),
                                 breaks=c(-8, 1, 3, 5, 7, 10, 15, 21, 28)) +
              labs(x="Days postinfection",
                   y="Fold change (log2)",
                   col=element_blank(),
                   title=i) +
              theme(axis.title=element_blank(),
                    legend.position=c(0.9, 0.9))
          })
# extract axis labels
x <- cowplot::get_plot_component(modplt[[1]] + theme(axis.title.x=element_text()),
                                 "xlab-b")
y <- cowplot::get_plot_component(modplt[[1]] + theme(axis.title.y=element_text()),
                                 "ylab-l")
# assemble figure
modplt <- cowplot::plot_grid(plotlist=modplt, ncol=1)
modplt <- cowplot::plot_grid(modplt, x, ncol=1, rel_heights=c(25, 1))
figB <- cowplot::plot_grid(y, modplt, nrow=1, rel_widths=c(1, 20))
figB
ggsave("analysis/modules-dpi.png", scale=scalefact,
       units="in", width=3.75, height=4)

# clean up
rm(x, y, modplt, glist, gene.mod, gene.examp, anchor.genes, rmat, desq)

## IPA heatmap -----------------------------------------------------------------
# canonical and functional
zcol <- circlize::colorRamp2(breaks=c(-3, 0, 5), 
                             colors=c("#377eb8", "white", "#e41a1c"))
figC <- read.csv("analysis/ipa-output.csv") %>%
        filter(Type != "Upstream") %>%
        select(-Type) %>%
        column_to_rownames("Pathway") %>% 
        as.matrix() %>% 
        ComplexHeatmap::Heatmap(name="z-score", 
                                heatmap_legend_param=list(direction="horizontal",
                                                          title_position="lefttop"),
                                col=zcol, na_col="white",
                                border=TRUE,
                                cluster_columns=FALSE, 
                                cluster_rows=TRUE,
                                row_names_gp=grid::gpar(fontsize=8),
                                column_title="Days postinfection",
                                column_title_side="bottom",
                                column_labels=c(1, 3, 5, 7, 10, 15, 21, 28),
                                column_names_centered=TRUE,
                                column_names_rot=0) %>%
        ComplexHeatmap::draw(heatmap_legend_side="bottom") %>%
        grid::grid.grabExpr()

# upstream regulators
zcol <- circlize::colorRamp2(breaks=c(-4, 0, 9), 
                             colors=c("#377eb8", "white", "#e41a1c"))
figD <- read.csv("analysis/ipa-output.csv") %>%
        filter(Type == "Upstream") %>%
        select(-Type) %>%
        column_to_rownames("Pathway") %>% 
        as.matrix() %>% 
        ComplexHeatmap::Heatmap(name="z-score", 
                                col=zcol, na_col="white",
                                border=TRUE,
                                cluster_columns=FALSE, 
                                cluster_rows=TRUE,
                                row_names_gp=grid::gpar(fontsize=8),
                                column_title="Days postinfection",
                                column_title_side="bottom",
                                column_labels=c(1, 3, 5, 7, 10, 15, 21, 28),
                                column_names_centered=TRUE,
                                column_names_rot=0) %>%
        ComplexHeatmap::draw() %>%
        grid::grid.grabExpr()

rm(zcol)

## knit together the figure and supplemental figure ----------------------------
# main figure
x <- cowplot::plot_grid(figA, figC, ncol=1, labels=c("A", "C"), rel_heights=c(1, 1.6))
y <- cowplot::plot_grid(figB, figD, ncol=1, labels=c("B", "D"), rel_heights=c(2, 1))
cowplot::plot_grid(x, y, nrow=1)
ggsave("analysis/figure5.png",
       units="in", width=7.5, height=8)
rm(x, y)

# supplemental figure
x <- cowplot::plot_grid(supB, supC, supD, 
                        nrow=1, labels=c("B", "C", "D"))
cowplot::plot_grid(supA, x, supE, 
                   ncol=1, labels=c("A", NA, "E"),
                   rel_heights=c(2, 1, 2))
ggsave("analysis/supplemental7.png", scale=1.5, 
       units="in", width=7.5, height=10)

## fin -------------------------------------------------------------------------
sessionInfo()
