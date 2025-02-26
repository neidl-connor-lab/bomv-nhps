#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggsurvfit))
theme_set(theme_classic())
options(stringsAsFactors=FALSE)

# helper variables
cols.virus <- c("EBOV (historical)"="#e41a1c", BOMV="black")
cols.dpi <- c("0"="lightblue", "15"="#3182bd", "28"="darkblue")
lines.nhps <- c(A=1, B=2, C=3, D=4, E=5, "F"=6, G=7, H=8)
shapes.nhps <- c(A=21, B=22, C=23, D=24, E=25, "F"=3, G=4, H=8)
sample.days <- c(0, 1, 3, 5, 7, 10, 15, 21, 28)
scalefact <- 1.5

# read in data file
sheets <- readxl::excel_sheets("data.xlsx")
data <- sheets %>%
        lapply(function(i) {
          readxl::read_excel("data.xlsx", sheet=i)
        })
names(data) <- sheets
rm(sheets)

## survival KM curve -----------------------------------------------------------
surv <- data$survival %>%
        mutate(Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                            labels=c("BOMV", "EBOV (historical)"))) %>%
        survfit2(Surv(DPI, Censoring) ~ Virus, data=.) %>%
        ggsurvfit(theme=theme_classic(), 
                  size=0.5) +
        scale_color_manual(values=cols.virus) +
        labs(x="Days postinfection",
             y="Survival",
             col="Virus") +
        theme(legend.position=c(0.8, 0.4)) +
        scale_x_continuous(limits=c(0, 28),
                           breaks=sample.days) + 
        theme(legend.title=element_text())
surv
ggsave("figures/survival.png", scale=scalefact,
       units="in", width=3.75, height=2)

## clinical illness ------------------------------------------------------------
clin <- data$clinical %>%
        reshape2::melt(id.vars=c("NHP", "DPI"),
                       variable.name="Parameter") %>%
        ggplot(aes(DPI, NHP, fill=value)) +
        geom_tile(col="black") +
        scale_fill_manual(values=c("TRUE"="black", "FALSE"="white")) +
        facet_wrap(~Parameter, ncol=1) +
        scale_y_discrete(limits=c("D", "C", "B", "A")) +
        scale_x_continuous(expand=c(0, 0),
                           breaks=sample.days) +
        labs(x="Days postinfection",
            y="BOMV NHPs") +
        theme(legend.position="none",
              axis.line.y.left=element_blank(),
              strip.background=element_blank(),
              strip.text=element_text(size=10,
                                      hjust=0,
                                      margin=margin(l=0)))
clin
ggsave("figures/clinical.png", scale=scalefact,
       units="in", width=3.75, height=6)

## viremia: viral genomes and replicating virus --------------------------------
viremia <- list()
viremia$full <- data$viremia %>%
                # add pseudo-count and transform to log10 scale
                mutate(`PFU/mL`=log10(`PFU/mL`+1), 
                       `Genomes/mL`=log10(`Genomes/mL`+1),
                       Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                                    labels=c("BOMV", "EBOV (historical)")))

# format data: distributions for EBOV and spaghetti for BOMV
viremia$bomv <- viremia$full %>%
                filter(Virus=="BOMV")
viremia$ebov <- viremia$full %>%
                filter(Virus=="EBOV (historical)") %>%
                group_by(DPI, Virus) %>%
                summarise(StDevGenomes=sd(`Genomes/mL`, na.rm=TRUE),
                          `Genomes/mL`=median(`Genomes/mL`, na.rm=TRUE),
                          StDevPFU=sd(`PFU/mL`, na.rm=TRUE),
                          `PFU/mL`=median(`PFU/mL`, na.rm=TRUE),
                          .groups="drop") %>%
                # add upper and lower error bar bounds
                mutate(UpperGenomes=`Genomes/mL`+StDevGenomes,
                       LowerGenomes=`Genomes/mL`-StDevGenomes,
                       UpperPFU=`PFU/mL`+StDevPFU,
                       LowerPFU=`PFU/mL`-StDevPFU)
# can't have negative values, so cap lower bounds at zero
viremia$ebov[viremia$ebov < 0] <- 0

# stats: get significance by DPI (3, 5, & 7)
sig <- list()
x <- filter(viremia$full, !(NHP %in% c("E", "F", "G", "H")))
sig$pfu <- ggpubr::compare_means(`PFU/mL` ~ Virus, 
                                 group.by="DPI", 
                                 data=x) %>%
           # add y-coordinate using the upper error bar
           left_join(viremia$ebov, by="DPI") %>%
           mutate(y.position=UpperPFU+0.2) %>%
           select(DPI, p.signif, y.position)
sig$pcr <- ggpubr::compare_means(`Genomes/mL` ~ Virus, 
                                 group.by="DPI", 
                                 data=x) %>%
           # add y-coordinate using the upper error bar
           left_join(viremia$ebov, by="DPI") %>%
           mutate(y.position=UpperGenomes+0.2) %>%
           select(DPI, p.signif, y.position)
rm(x)

# plot viral genomes
pcr <- viremia$ebov %>%
      ggplot(aes(DPI, `Genomes/mL`)) +
      # plot EBOV with error bars
      geom_line(aes(col=Virus)) +
      geom_errorbar(aes(ymax=UpperGenomes, ymin=LowerGenomes, col=Virus),
                    width=0.5, linewidth=0.25) +
      geom_point(aes(fill=Virus), pch=21, size=2) +
      # plot individual BOMV 
      geom_line(data=viremia$bomv, aes(col=Virus, linetype=NHP)) +
      geom_point(data=viremia$bomv, aes(fill=Virus, shape=NHP), size=2) +
      scale_fill_manual(values=cols.virus) +
      scale_color_manual(values=cols.virus) +
      scale_shape_manual(values=shapes.nhps) +
      scale_linetype_manual(values=lines.nhps) +
      # add p-values
      geom_text(data=sig$pcr, aes(label=p.signif, y=y.position), size=5) +
      # format the axes
      labs(x="Days post infection",
           y="GEq/mL",
           fill="Virus",
           col="Virus",
           shape="BOMV NHPs",
           linetype="BOMV NHPs") + 
      scale_y_continuous(limits=c(0, 13), breaks=c(0, 4, 8, 12),
                         labels=c("LOD", "1e4", "1e8", "1e12")) +
      scale_x_continuous(limits=c(NA, 28), breaks=sample.days) +
      guides(fill=guide_legend(ncol=2,
                               override.aes=list(pch=21)),
             shape=guide_legend(ncol=2,
                                override.aes=list(fill=cols.virus["BOMV"]))) +
      theme(legend.position=c(0.8, 0.55))
pcr
ggsave("figures/pcr.png", scale=scalefact,
       units="in", width=3.75, height=2)

# plot replicating virus
pfu <- viremia$ebov %>% 
      ggplot(aes(DPI, `PFU/mL`)) +
      # plot EBOV with error bars
      geom_line(aes(col=Virus)) +
      geom_errorbar(aes(ymax=UpperPFU, ymin=LowerPFU, col=Virus),
                    width=0.5, linewidth=0.25) +
      geom_point(pch=21, aes(fill=Virus), size=2) +
      # plot individual BOMV 
      geom_line(data=viremia$bomv, aes(col=Virus, linetype=NHP)) +
      geom_point(data=viremia$bomv, aes(fill=Virus, shape=NHP), size=2) +
      scale_fill_manual(values=cols.virus) +
      scale_color_manual(values=cols.virus) +
      scale_shape_manual(values=shapes.nhps) +
      scale_linetype_manual(values=lines.nhps) +
      # add p-values
      geom_text(data=sig$pfu, aes(label=p.signif, y=y.position), size=5) +
      # format the axes
      labs(x="Days post infection",
           y="PFU/mL",
           fill="Virus",
           col="Virus",
           shape="BOMV NHPs",
           linetype="BOMV NHPs") +
      scale_y_continuous(limits=c(0, 13), breaks=c(0, 4, 8, 12),
                         labels=c("LOD", "1e4", "1e8", "1e12")) +
      scale_x_continuous(limits=c(NA, 28), breaks=sample.days) +
      guides(fill=guide_legend(ncol=2,
                               override.aes=list(pch=21)),
             shape=guide_legend(ncol=2,
                                override.aes=list(fill=cols.virus["BOMV"]))) +
      theme(legend.position=c(0.8, 0.55))
pfu
ggsave("figures/viremia.png", scale=scalefact,
       units="in", width=3.75, height=2)

# assemble figure 1
f1 <- cowplot::plot_grid(surv,
                         pcr + guides(color="none", fill="none"),
                         pfu + theme(legend.position="none"),
                         ncol=1, labels=c("A", "C", "D"))
cowplot::plot_grid(f1, clin, labels=c(NA, "B"), nrow=1)
ggsave("figures/figure1.png", units="in", width=7.5, height=6)

# clean up
rm(f1, surv, clin, pcr, pfu, sig, viremia)

## hematology ------------------------------------------------------------------
hema <- list()
hema$full <- data$hematology %>%
             filter(DPI >= 0) %>%
             select(-starts_with("pct")) %>%
             reshape2::melt(id.vars=c("NHP", "DPI", "Virus"),
                            variable.name="Analyte") %>%
             mutate(Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                                 labels=c("BOMV", "EBOV (historical)")))

# split into BOMV and EBOV distribution
hema$bomv <- filter(hema$full, Virus=="BOMV")
hema$ebov <- hema$full %>%
             filter(Virus=="EBOV (historical)") %>%
             group_by(DPI, Virus, Analyte) %>%
             summarise(StDev=sd(value),
                       value=median(value),
                       .groups="drop") %>%
             mutate(Upper=value+StDev,
                    Lower=value-StDev)
# concentrations can't be < 0 so set the cap
hema$ebov[hema$ebov < 0] <- 0

# loop over each analyte and plot
analytes <- unique(hema$full$Analyte)
hema <- analytes %>%
        lapply(function(i) {
          # subset data
          bomv <- filter(hema$bomv, Analyte==i)
          ebov <- filter(hema$ebov, Analyte==i)
          
          # format y-axis label using units
          x <- data$units$Units[data$units$Analyte==i]
          
          # plot it
          ebov %>%
            ggplot(aes(DPI, value)) +
            geom_line(aes(col=Virus)) +
            geom_errorbar(aes(col=Virus, ymax=Upper, ymin=Lower),
                          width=1, linewidth=0.25) +
            geom_point(aes(fill=Virus), pch=21) +
            # add in BOMV
            geom_line(data=bomv, aes(col=Virus, linetype=NHP)) +
            geom_point(data=bomv, aes(fill=Virus, shape=NHP)) +
            scale_color_manual(values=cols.virus) +
            scale_fill_manual(values=cols.virus) +
            scale_shape_manual(values=shapes.nhps) +
            scale_linetype_manual(values=lines.nhps) +
            scale_x_continuous(limits=c(NA, 28), breaks=sample.days) +
            labs(x="Days postinfection",
                 y=x,
                 col="Virus",
                 fill="Virus",
                 shape="BOMV NHPs",
                 linetype="BOMV NHPs",
                 title=i) +
            theme(legend.position="none", title=element_text(size=7)) +
            guides(fill=guide_legend(nrow=1, 
                                     override.aes=list(pch=21)),
                   shape=guide_legend(nrow=1, 
                                      override.aes=list(fill=cols.virus["BOMV"])))
        })
names(hema) <- analytes

# extract legend
leg <- cowplot::get_plot_component(hema[[1]] + theme(legend.position="bottom"), 
                                   "guide-box-bottom")
# plot panels
x <- cowplot::plot_grid(plotlist=hema, ncol=3, labels="AUTO")
cowplot::plot_grid(x, leg, ncol=1, rel_heights=c(20, 1))
ggsave("figures/supplemental2.png",
       units="in", width=7.5, height=10)

# plot selected analytes with legend
leg <- cowplot::get_plot_component(hema[[1]] + 
                                     guides(fill=guide_legend(ncol=2),
                                            shape=guide_legend(ncol=2, ,
                                                               override.aes=list(fill=cols.virus["BOMV"]))) +
                                     theme(legend.position="right"), 
                                   "guide-box-right")
f2AB <- cowplot::plot_grid(hema$Platelets, hema$Lymphocytes, leg, 
                           nrow=1, labels=c("A", "B"))

# clean up
rm(x, leg, hema, analytes)

## clinical chemistry ----------------------------------------------------------
chem <- list()
chem$full <- data$chemistry %>%
             filter(DPI >= 0) %>%
             reshape2::melt(id.vars=c("NHP", "DPI", "Virus"),
                            variable.name="Analyte") %>%
             na.omit() %>%
             # log10 for plotting concentration
             mutate(value=log10(value),
                    Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                                 labels=c("BOMV", "EBOV (historical)")))

# split into BOMV and EBOV distribution
chem$bomv <- filter(chem$full, Virus=="BOMV")
chem$ebov <- chem$full %>%
             filter(Virus=="EBOV (historical)") %>%
             group_by(DPI, Virus, Analyte) %>%
             summarise(StDev=sd(value),
                       value=median(value),
                       .groups="drop") %>%
             mutate(Upper=value+StDev,
                    Lower=value-StDev)

# loop over each analyte and plot
analytes <- unique(chem$full$Analyte)
chem <- analytes %>%
        lapply(function(i) {
          # subset data
          bomv <- filter(chem$bomv, Analyte==i)
          ebov <- filter(chem$ebov, Analyte==i)
          
          # format y-axis label using units
          x <- data$units[data$units$Analyte==i, "Units"]
          x <- paste(x, "(log10)")
          
          # plot it
          ebov %>%
            ggplot(aes(DPI, value)) +
            geom_line(aes(col=Virus)) +
            geom_errorbar(aes(col=Virus, ymax=Upper, ymin=Lower),
                          width=1, linewidth=0.25) +
            geom_point(aes(fill=Virus), pch=21) +
            # add in BOMV
            geom_line(data=bomv, aes(col=Virus, linetype=NHP)) +
            geom_point(data=bomv, aes(fill=Virus, shape=NHP)) +
            scale_color_manual(values=cols.virus) +
            scale_fill_manual(values=cols.virus) +
            scale_shape_manual(values=shapes.nhps) +
            scale_linetype_manual(values=lines.nhps) +
            scale_x_continuous(limits=c(NA, 28), breaks=sample.days) +
            labs(x="Days postinfection",
                 y=x,
                 col="Virus",
                 fill="Virus",
                 shape="BOMV NHPs",
                 linetype="BOMV NHPs",
                 title=i) +
            theme(legend.position="none", 
                  title=element_text(size=7)) +
            guides(fill=guide_legend(nrow=1, 
                                     override.aes=list(pch=21)),
                   shape=guide_legend(nrow=1, 
                                      override.aes=list(fill=cols.virus["BOMV"])))
        })
names(chem) <- analytes

# extract legend
leg <- cowplot::get_plot_component(chem[[1]] + theme(legend.position="bottom"), 
                                   "guide-box-bottom")
# plot panels
x <- cowplot::plot_grid(plotlist=chem, ncol=3, labels="AUTO")
cowplot::plot_grid(x, leg, ncol=1, rel_heights=c(20, 1))
ggsave("figures/supplemental3.png",
       units="in", width=7.5, height=8.5)

# plot selected analytes with selected hematology
f2CDE <- c("Aspartate transferase (AST)",
           "Blood urea nitrogen (BUN)",
           "C-reactive protein (CRP)")
f2CDE <- cowplot::plot_grid(plotlist=chem[f2CDE], 
                            nrow=1, labels=c("C", "D", "E"))
cowplot::plot_grid(f2AB, f2CDE, ncol=1)
ggsave("figures/figure2.png",
       units="in", width=7.5, height=4.5)

# clean up
rm(x, leg, chem, f2AB, f2CDE, analytes)

## coag legendplex panels ------------------------------------------------------
coag <- data$coag %>%
        reshape2::melt(id.vars=c("NHP", "DPI", "Replicate"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        # average the technical duplicates
        group_by(NHP, DPI, Analyte) %>%
        summarise(Concentration=mean(Concentration),
                  .groups="drop") 

# pull out 0 DPI and calculate fold change
dpi0 <- coag %>%
        filter(DPI==0) %>%
        select(-DPI) %>%
        rename(DPI0=Concentration)
coag <- coag %>%
        left_join(dpi0, by=c("NHP", "Analyte")) %>%
        mutate(log2fc=log2(Concentration/DPI0))
rm(dpi0)

# plot log2FC over time for each analyte
analytes <- unique(coag$Analyte)
coag <- analytes %>%
        lapply(function(i) {
          coag %>%
            filter(Analyte==i) %>%
            ggplot(aes(DPI, log2fc)) +
            geom_hline(yintercept=0, linetype=2, col="lightgrey") +
            geom_line(aes(linetype=NHP)) +
            geom_point(aes(shape=NHP), fill=cols.virus["BOMV"]) +
            scale_shape_manual(values=shapes.nhps) +
            scale_linetype_manual(values=lines.nhps) +
            scale_y_continuous(limits=c(-4.5, 6),
                               breaks=c(-4, -2, 0, 2, 4, 6)) +
            scale_x_continuous(limits=c(0, 28), breaks=sample.days) +
            labs(x="Days postinfection",
                 y="Fold change (log2)",
                 shape="BOMV NHPs",
                 linetype="BOMV NHPs",
                 title=i) +
            theme(legend.position="none",
                  axis.title=element_text(size=8))
        })
names(coag) <- analytes

# pull out legend
leg <- cowplot::get_plot_component(coag[[1]] + theme(legend.position="bottom"),
                                   "guide-box-bottom")
x <- cowplot::plot_grid(plotlist=coag, labels="AUTO", ncol=3)
cowplot::plot_grid(x, leg, ncol=1, rel_heights=c(10, 1))
ggsave("figures/supplemental4.png", 
       units="in", width=7.5, height=8.5)

# plot selected markers
analytes <- c("D-Dimer", "tPA", "PAI-1")
f3ABC <- cowplot::plot_grid(plotlist=coag[analytes], labels="AUTO", nrow=1)

# clean up
rm(coag, analytes, x, leg)

## inflammation ----------------------------------------------------------------
inflam <- data$inflammation %>%
          reshape2::melt(id.vars=c("NHP", "DPI"),
                         variable.name="Analyte") %>%
          # collapse technical duplicates
          group_by(NHP, DPI, Analyte) %>%
          summarise(value=mean(value),
                    .groups="drop") %>%
          # add a pseudo-count for log2FC handling
          mutate(value=value+1,
                 NHP=factor(NHP))

# calculate log2fc
dpi0 <- inflam %>%
        filter(DPI == 0) %>%
        select(-DPI) %>%
        rename(DPI0=value)
inflam <- inflam %>%
          left_join(dpi0, by=c("NHP", "Analyte")) %>%
          mutate(log2fc=log2(value/DPI0))
rm(dpi0)

# plot log2FC over time for each analyte
analytes <- unique(inflam$Analyte)
inflam <- analytes %>%
          lapply(function(i) {
            inflam %>%
              filter(Analyte==i) %>%
              ggplot(aes(DPI, log2fc)) +
              geom_hline(yintercept=0, linetype=2, col="lightgrey") +
              geom_line(aes(linetype=NHP)) +
              geom_point(aes(shape=NHP), fill=cols.virus["BOMV"]) +
              scale_shape_manual(values=shapes.nhps) +
              scale_linetype_manual(values=lines.nhps) +
              scale_y_continuous(limits=c(-6.2, 11),
                                 breaks=c(-6, -3, 0, 3, 6, 9)) +
              scale_x_continuous(limits=c(0, 28), breaks=sample.days) +
              labs(x="Days postinfection",
                   y="Fold change (log2)",
                   shape="BOMV NHPs",
                   linetype="BOMV NHPs",
                   title=i) +
              theme(legend.position="none", 
                    axis.title=element_text(size=8))
          })
names(inflam) <- analytes

# pull out legend
leg <- cowplot::get_plot_component(inflam[[1]] + theme(legend.position="bottom"),
                                   "guide-box-bottom")
x <- cowplot::plot_grid(plotlist=inflam, 
                        labels="AUTO", ncol=4)
cowplot::plot_grid(x, leg, ncol=1, rel_heights=c(25, 1))
ggsave("figures/supplemental5.png", 
       units="in", width=7.5, height=10)

# plot selected markers with coag selected markers for figure 3
analytes <- c("IFNg", "IL-1ra", "IL-6")
f3DEF <- cowplot::plot_grid(plotlist=inflam[analytes], 
                            labels=c("D", "E", "F"), nrow=1)
cowplot::plot_grid(f3ABC, f3DEF, leg, ncol=1, rel_heights=c(10, 10, 1))
ggsave("figures/figure3.png", 
       units="in", width=7.5, height=5)

# clean up
rm(inflam, analytes, x, leg, f3ABC, f3DEF)

## ELISAs ----------------------------------------------------------------------
elisas <- data$elisa %>%
          reshape2::melt(id.vars=c("NHP", "DPI"),
                         variable.name="Analyte",
                         value.name="Titer")
elisas <- elisas$Analyte %>%
          unique() %>%
          lapply(function(i) {
            # plot on a log scale using a pseudo-count
            elisas %>%
              filter(Analyte==i) %>%
              ggplot(aes(DPI, Titer+1)) +
              geom_line(aes(linetype=NHP),
                        col=cols.virus["BOMV"]) +
              geom_point(aes(shape=NHP), 
                         fill=cols.virus["BOMV"], 
                         size=2) +
              scale_shape_manual(values=shapes.nhps) +
              scale_linetype_manual(values=lines.nhps) +
              scale_y_continuous(limits=c(NA, 1e5),
                                 breaks=c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5),
                                 labels=c("<LOD", "10", "100", "1000", 
                                          "10,000", "100,000"),
                                 transform="log10") +
              scale_x_continuous(limits=c(NA, 28), 
                                 breaks=c(0, 5, 7, 10, 15, 28)) +
              labs(x="Days postinfection",
                   y="Endpoint titer",
                   shape="BOMV NHPs",
                   linetype="BOMV NHPs",
                   title=paste("Anti-BOMV", i)) +
              theme(legend.position="none")
          })

# pull out legend
leg <- cowplot::get_plot_component(elisas[[1]] + theme(legend.position="right"),
                                   "guide-box-right")
f4AB <- cowplot::plot_grid(plotlist=elisas, nrow=1, labels="AUTO")
f4AB <- cowplot::plot_grid(f4AB, leg, nrow=1, rel_widths=c(7, 1))
f4AB

# clean up
rm(elisas, leg)

## PRNTs -----------------------------------------------------------------------
prnt <- data$prnt %>%
        mutate(DPI=factor(DPI, levels=c(0, 15, 28))) %>%
        # collapse replicates
        group_by(NHP, DPI, Virus, Dilution) %>%
        summarise(Plaques=mean(Plaques), 
                  .groups="drop")

# pull out virus control and fill in gaps
vc <- prnt %>%
      filter(NHP=="Control") %>%
      select(Virus, Plaques) %>%
      rename(Control=Plaques)
prnt <- prnt %>%
        filter(NHP!="Control") %>%
        left_join(vc, by="Virus") %>%
        mutate(Reduction=Control-Plaques,
               PercentReduction=100*Reduction/Control) 

# cap the reduction at 0%
prnt[prnt$PercentReduction < 0, "PercentReduction"] <- 0
rm(vc)

# plot homologous neutralization curves
curves <- prnt$NHP %>%
          unique() %>%
          lapply(function(i) {
            prnt %>%
              filter(Virus=="BOMV",
                     NHP==i) %>%
              ggplot(aes(Dilution, PercentReduction)) +
              geom_hline(yintercept=50, linetype=2, col="lightgrey") +
              geom_line(aes(group=DPI, col=DPI, linetype=NHP)) +
              geom_point(aes(fill=DPI, shape=NHP), size=3) +
              scale_color_manual(values=cols.dpi, guide="none") +
              scale_fill_manual(values=cols.dpi) +
              scale_linetype_manual(values=lines.nhps, guide="none") +
              scale_shape_manual(values=shapes.nhps, guide="none") +
              scale_x_continuous(limits=c(10, 640),
                                 breaks=c(10, 20, 40, 80, 160, 320, 640),
                                 transform="log2") +
              ylim(0, 100) +
              labs(x="Dilution factor",
                   y="Reduction (%)",
                   fill="DPI",
                   title=paste("BOMV NHP", i)) +
              guides(fill=guide_legend(override.aes=list(pch=21, size=3))) +
              theme(legend.position="none")
          })

# pull out legend and plot figure 4
leg <- cowplot::get_plot_component(curves[[1]] + theme(legend.position="right"),
                                   "guide-box-right")
f4CDEF <- cowplot::plot_grid(plotlist=curves, ncol=2, 
                             labels=c("C", "D", "E", "F"))
f4CDEF <- cowplot::plot_grid(f4CDEF, leg, nrow=1, rel_widths=c(10, 1))
cowplot::plot_grid(f4AB, f4CDEF, ncol=1, rel_heights=c(1, 2))
ggsave("figures/figure4.png",
       units="in", width=7.5, height=6)

# plot heterologous neutralization curves
curves <- prnt %>%
          filter(Virus!="BOMV") %>%
          select(NHP, Virus) %>%
          arrange(Virus) %>%
          distinct() %>%
          apply(1, function(i) {
            prnt %>%
              filter(Virus==i["Virus"],
                     NHP==i["NHP"]) %>%
              ggplot(aes(Dilution, PercentReduction)) +
              geom_hline(yintercept=50, linetype=2, col="lightgrey") +
              geom_line(aes(group=DPI, col=DPI, linetype=NHP)) +
              geom_point(aes(fill=DPI, shape=NHP), size=3) +
              scale_color_manual(values=cols.dpi, guide="none") +
              scale_fill_manual(values=cols.dpi) +
              scale_linetype_manual(values=lines.nhps, guide="none") +
              scale_shape_manual(values=shapes.nhps, guide="none") +
              scale_x_continuous(limits=c(10, 640),
                                 breaks=c(10, 20, 40, 80, 160, 320, 640),
                                 transform="log2") +
              ylim(0, 100) +
              labs(x="Dilution factor",
                   y="Reduction (%)",
                   fill="DPI",
                   title=paste("BOMV NHP", i["NHP"], "vs.", i["Virus"])) +
              guides(fill=guide_legend(override.aes=list(pch=21, size=3))) +
              theme(legend.position="none")
          })
# pull out legend
leg <- cowplot::get_plot_component(curves[[1]] + theme(legend.position="bottom"),
                                   "guide-box-bottom")
curves <- cowplot::plot_grid(plotlist=curves, ncol=2, 
                             labels="AUTO")
cowplot::plot_grid(curves, leg, 
                   ncol=1, rel_heights=c(15, 1))
ggsave("figures/supplemental6.png",
       units="in", width=7.5, height=8)

# clean up
rm(prnt, leg, curves)


## done! -----------------------------------------------------------------------
sessionInfo()
