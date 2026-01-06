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

# read in data file
sheets <- readxl::excel_sheets("data/data.xlsx")
data <- sheets %>%
  lapply(function(i) {
    readxl::read_excel("data/data.xlsx", sheet=i)
  })
names(data) <- sheets
rm(sheets)

# init figure variables
plotvals <- list()
fig1 <- list()
fig2 <- list() 
fig3 <- list()

## survival KM curve -----------------------------------------------------------
plotvals$`Figure 1a` <- data$survival %>%
                        mutate(Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                                            labels=c("BOMV", "EBOV (historical)")))
fig1$a <- plotvals$`Figure 1a` %>%
          survfit2(Surv(DPI, Censoring) ~ Virus, data=.) %>%
          ggsurvfit(theme=theme_classic(), 
                    linewidth=0.5) +
          scale_color_manual(values=cols.virus) +
          labs(x="Days postinfection",
               y="Survival",
               col="Virus") +
          theme(legend.position=c(0.8, 0.4)) +
          scale_x_continuous(limits=c(0, 28),
                             breaks=sample.days) + 
          theme(legend.title=element_text())

## clinical illness ------------------------------------------------------------
plotvals$`Figure 1b` <- data$clinical %>%
                        reshape2::melt(id.vars=c("NHP", "DPI"),
                                       variable.name="Parameter")
fig1$b <- plotvals$`Figure 1b` %>%
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

## viremia: viral genomes and replicating virus --------------------------------
viremia <- list()
viremia$full <- data$viremia %>%
                # add pseudo-count and transform to log10 scale
                mutate(`PFU/mL`=log10(`PFU/mL`+1), 
                       `Genomes/mL`=log10(`Genomes/mL`+1),
                       Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                                    labels=c("BOMV", "EBOV (historical)")))

# save values
plotvals$`Figure 1c` <- viremia$full %>%
                        select(NHP, DPI, Virus, `Genomes/mL`) %>%
                        rename(`log10 GEq/mL`=`Genomes/mL`)
plotvals$`Figure 1d` <- viremia$full %>%
                        select(NHP, DPI, Virus, `PFU/mL`) %>%
                        rename(`log10 PFU/mL`=`PFU/mL`)

# format data for plotting: distributions for EBOV and spaghetti for BOMV
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
fig1$c <- viremia$ebov %>%
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

# plot replicating virus
fig1$d <- viremia$ebov %>% 
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

# clean up
rm(sig, viremia)

## hematology ------------------------------------------------------------------
hema <- list()
hema$full <- data$hematology %>%
             filter(DPI >= 0) %>%
             select(-starts_with("pct")) %>%
             reshape2::melt(id.vars=c("NHP", "DPI", "Virus"),
                            variable.name="Analyte") %>%
             mutate(Virus=factor(Virus, levels=c("BOMV", "EBOV"),
                                 labels=c("BOMV", "EBOV (historical)")))

# compare lymphopenia and thrombocytopenia; report mean and st. dev.
# (minimum values lymphocytes & platelets)
hema$full %>%
  filter(Analyte %in% c("Lymphocytes", "Platelets")) %>%
  group_by(NHP, Virus, Analyte) %>%
  summarise(Concentration=min(value, na.rm=TRUE),
            .groups="drop") %>%
  ggpubr::compare_means(Concentration ~ Virus, data=., 
                        group.by="Analyte", p.adjust.method="fdr")
hema$full %>%
  filter(Analyte %in% c("Lymphocytes", "Platelets")) %>%
  group_by(Virus, Analyte) %>%
  summarise(Mean=min(value, na.rm=TRUE),
            StDev=sd(value, na.rm=TRUE),
            .groups="drop")

# save values
plotvals$`Supplemental 2` <- hema$full %>%
                             reshape2::dcast(NHP + DPI + Virus ~ Analyte,
                                             value.var="value")
plotvals$`Figure 2a` <- plotvals$`Supplemental 2` %>%
                        select(NHP, DPI, Virus, Platelets)
plotvals$`Figure 2b` <- plotvals$`Supplemental 2` %>%
                        select(NHP, DPI, Virus, Lymphocytes)

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

# supplemental figure 2 has all hematology analytes
sup2 <- hema
names(sup2) <- letters[1:length(sup2)]

# figure 2 a-b are platelets and lymphocytes, respectively
fig2$a <- hema$Platelets
fig2$b <- hema$Lymphocytes

# clean up
rm(hema, analytes)

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

# save values
plotvals$`Supplemental 3` <- chem$full %>%
                             reshape2::dcast(NHP + DPI + Virus ~ Analyte,
                                             value.var="value")
plotvals$`Figure 2c` <- plotvals$`Supplemental 3` %>%
                        select(NHP, DPI, Virus, `Aspartate transferase (AST)`)
plotvals$`Figure 2d` <- plotvals$`Supplemental 3` %>%
                        select(NHP, DPI, Virus, `Blood urea nitrogen (BUN)`)
plotvals$`Figure 2e` <- plotvals$`Supplemental 3` %>%
                        select(NHP, DPI, Virus, `C-reactive protein (CRP)`)

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

# supplemental figure 3 has all chemistry analytes
sup3 <- chem
names(sup3) <- letters[1:length(sup3)]

# figure 2 c-e are selected analytes
fig2$c <- chem$`Aspartate transferase (AST)`
fig2$d <- chem$`Blood urea nitrogen (BUN)`
fig2$e <- chem$`C-reactive protein (CRP)`

# clean up
rm(chem, analytes)

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

# save values
plotvals$`Supplemental 4` <- coag %>%
                             reshape2::dcast(NHP + DPI ~ Analyte,
                                             value.var="log2fc")
plotvals$`Figure 3a` <- plotvals$`Supplemental 4` %>%
                        select(NHP, DPI, `D-Dimer`)
plotvals$`Figure 3b` <- plotvals$`Supplemental 4` %>%
                        select(NHP, DPI, tPA)
plotvals$`Figure 3c` <- plotvals$`Supplemental 4` %>%
                        select(NHP, DPI, `PAI-1`)

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

# supplemental figure 4 has all coagulation analytes
sup4 <- coag
names(sup4) <- letters[1:length(sup4)]

# figure 3 a-c are selected analytes
fig3$a <- coag$`D-Dimer`
fig3$b <- coag$tPA
fig3$c <- coag$`PAI-1`

# clean up
rm(coag, analytes)

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

# save values
plotvals$`Supplemental 5` <- inflam %>%
                             reshape2::dcast(NHP + DPI ~ Analyte,
                                             value.var="log2fc")
plotvals$`Figure 3d` <- plotvals$`Supplemental 5` %>%
                        select(NHP, DPI, IFNg)
plotvals$`Figure 3e` <- plotvals$`Supplemental 5` %>%
                        select(NHP, DPI, `IL-1ra`)
plotvals$`Figure 3f` <- plotvals$`Supplemental 5` %>%
                        select(NHP, DPI, `IL-6`)

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

# supplemental figure 5 has all inflammation analytes
sup5 <- inflam
names(sup5) <- letters[1:length(sup5)]

# figure 3 e-g are selected analytes
fig3$e <- inflam$IFNg
fig3$f <- inflam$`IL-1ra`
fig3$g <- inflam$`IL-6`

# clean up
rm(inflam, analytes)

## ELISAs ----------------------------------------------------------------------
# save values
plotvals$`Figure 5a` <- data$elisa %>%
                        filter(Antigen=="BOMV",
                               Immunoglobulin=="IgM")
plotvals$`Figure 5b` <- data$elisa %>%
                        filter(Antigen=="BOMV",
                               Immunoglobulin=="IgG")

# get range for BOMV 28 DPI titers
data$elisa %>%
  filter(DPI==28,
         Antigen=="BOMV") %>%
  group_by(Immunoglobulin) %>%
  summarise(Minimum=min(`Endpoint titer`),
            Maximum=max(`Endpoint titer`),
            .groups="drop")

# plot ELISAs
fig5 <- data$elisa$Immunoglobulin %>%
        unique() %>%
        lapply(function(i) {
          # plot on a log scale using a pseudo-count
          data$elisa %>%
            filter(Immunoglobulin==i,
                   Antigen=="BOMV") %>%
            ggplot(aes(DPI, `Endpoint titer`+1)) +
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
names(fig5) <- letters[1:length(fig5)]

# plot other antigens
elisa <- filter(data$elisa, Antigen != "BOMV")
plotvals$`Figure 5d` <- elisa
fig5$d <- elisa %>%
          mutate(DPI=as.factor(DPI)) %>%
          group_by(Antigen, DPI) %>%
          summarize(StDev=sd(`Endpoint titer`),
                    `Endpoint titer`=mean(`Endpoint titer`),
                    .groups="drop") %>%
          ggplot(aes(Antigen, `Endpoint titer`+1, group=DPI)) +
          geom_col(aes(fill=DPI), col="black", position="dodge") +
          geom_errorbar(aes(ymin=`Endpoint titer`-StDev, 
                            ymax=`Endpoint titer`+StDev),
                        position=position_dodge(0.9), width=0.5) +
          geom_point(data=elisa,
                     aes(shape=NHP),
                     fill="black",
                     position=position_jitterdodge(jitter.width=0.5, 
                                                   jitter.height=0)) +
          scale_shape_manual(values=shapes.nhps) +
          scale_fill_manual(values=c("0"="white", "28"="grey80")) +
          scale_y_continuous(limits=c(0.9, 1e4),
                             breaks=c(1e0, 1e1, 1e2, 1e3, 1e4),
                             labels=c("<LOD", "10", "100", "1000", "10,000"),
                             transform="log10") +
          labs(x="ELISA antigen",
               y="Endpoint titer",
               title="Anti-ebolavirus IgG") +
          guides(shape="none")

# clean up
rm(elisa)

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

# what are PRNT50 at 28 DPI?
prnt %>% 
  filter(DPI==28,
         PercentReduction >= 50) %>%
  group_by(NHP) %>%
  top_n(n=1, wt=Dilution) %>%
  select(NHP, Dilution, PercentReduction)

# save percent BOMV neutralization
plotvals$`Figure 5c` <- prnt %>%
                        filter(Virus=="BOMV", Dilution==10) %>%
                        select(NHP, DPI, PercentReduction)

# plot 1:10 percent BOMV neutralization
fig5$c <- prnt %>%
          filter(Virus=="BOMV",
                 Dilution==10) %>%
          ggplot(aes(DPI, PercentReduction)) +
          geom_boxplot(fill="grey80", outliers=FALSE) +
          geom_jitter(aes(shape=NHP), fill="black", width=0.2, height=0) +
          ggpubr::stat_compare_means(comparisons=list(c("0", "15"),
                                                      c("0", "28"),
                                                      c("15", "28")),
                                     label="p.signif", 
                                     step.increase=0.25) +
          scale_shape_manual(values=shapes.nhps) +
          scale_y_continuous(limits=c(0, NA),
                             breaks=c(0, 25, 50, 75, 100),
                             expand=expansion(mult = c(0, 0.1))) +
          labs(x="Days postinfection",
               y="Reduction (%)",
               shape="BOMV NHPs",
               title="BOMV neutralization")

# save neutralization curves
plotvals$`Supplemental 6` <- select(prnt,
                                    NHP, DPI, Virus, 
                                    Dilution, PercentReduction)

# plot neutralization curves
sup6 <- prnt %>%
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
                 title=paste("NHP", i["NHP"], "vs.", i["Virus"])) +
            guides(fill=guide_legend(override.aes=list(pch=21, size=3))) +
            theme(legend.position="none")
        })
names(sup6) <- letters[1:length(sup6)]

# clean up
rm(prnt)

## save plotting data file -----------------------------------------------------
openxlsx::write.xlsx(plotvals, "analysis/plotting.xlsx")

## assemble figures ------------------------------------------------------------
# figure 1
x <- cowplot::plot_grid(fig1$a,
                        fig1$c + guides(color="none", fill="none"),
                        fig1$d + theme(legend.position="none"),
                        ncol=1, labels=c("a", "c", "d"))
fig1 <- cowplot::plot_grid(x, fig1$b, labels=c(NA, "b"), nrow=1)
ggsave("analysis/figure1.pdf", fig1, units="in", width=7.5, height=6)

# figure 2
# plot selected analytes with legend
x <- cowplot::get_plot_component(fig2$a + 
                                 guides(fill=guide_legend(ncol=2),
                                        shape=guide_legend(ncol=2, ,
                                                           override.aes=list(fill=cols.virus["BOMV"]))) +
                                 theme(legend.position="right"), 
                                 "guide-box-right")
fig2 <- cowplot::plot_grid(fig2$a, fig2$b, x,
                           fig2$c, fig2$d, fig2$e,
                           nrow=2, labels=c("a", "b", "", "c", "d", "e"))
ggsave("analysis/figure2.pdf", fig2, units="in", width=7.5, height=4.5)

# figure 3
x <- cowplot::get_plot_component(fig3$a + theme(legend.position="bottom"),
                                  "guide-box-bottom")
fig3 <- cowplot::plot_grid(plotlist=fig3, labels=names(fig3), nrow=2)
fig3 <- cowplot::plot_grid(fig3, x, ncol=1, rel_heights=c(20, 1))
ggsave("analysis/figure3.pdf", fig3, units="in", width=7.5, height=5)

# figure 5: need to re-order by name
fig5 <- fig5[order(names(fig5))]
fig5 <- cowplot::plot_grid(plotlist=fig5, labels=names(fig5))
ggsave("analysis/figure5.pdf", fig5, units="in", width=7.5, height=5)

# supplemental 2
x <- cowplot::get_plot_component(sup2$a + theme(legend.position="bottom"), 
                                 "guide-box-bottom")
sup2 <- cowplot::plot_grid(plotlist=sup2, ncol=3, labels=names(sup2))
sup2 <- cowplot::plot_grid(sup2, x, ncol=1, rel_heights=c(20, 1))
ggsave("analysis/supplemental2.png", sup2, units="in", width=7.5, height=10)

# supplemental 3
x <- cowplot::get_plot_component(sup3$a + theme(legend.position="bottom"), 
                                 "guide-box-bottom")
sup3 <- cowplot::plot_grid(plotlist=sup3, ncol=3, labels=names(sup3))
sup3 <- cowplot::plot_grid(sup3, x, ncol=1, rel_heights=c(20, 1))
ggsave("analysis/supplemental3.png", sup3, units="in", width=7.5, height=8.5)

# supplemental 4
x <- cowplot::get_plot_component(sup4$a + theme(legend.position="bottom"),
                                 "guide-box-bottom")
sup4 <- cowplot::plot_grid(plotlist=sup4, labels=names(sup4), ncol=3)
sup4 <- cowplot::plot_grid(sup4, x, ncol=1, rel_heights=c(20, 1))
ggsave("analysis/supplemental4.png", sup4, units="in", width=7.5, height=8.5)

# supplemental 5
x <- cowplot::get_plot_component(sup5$a + theme(legend.position="bottom"),
                                 "guide-box-bottom")
sup5 <- cowplot::plot_grid(plotlist=sup5, labels=names(sup5), ncol=4)
sup5 <- cowplot::plot_grid(sup5, x, ncol=1, rel_heights=c(25, 1))
ggsave("analysis/supplemental5.png", sup5, units="in", width=7.5, height=10)

# supplemental 6
x <- cowplot::get_plot_component(sup6$a + theme(legend.position="bottom"),
                                 "guide-box-bottom")
sup6 <- cowplot::plot_grid(plotlist=sup6, ncol=3, labels=names(sup6))
sup6 <- cowplot::plot_grid(sup6, x, ncol=1, rel_heights=c(15, 1))
ggsave("analysis/supplemental6.png",
       units="in", width=7.5, height=8)

## done! -----------------------------------------------------------------------
sessionInfo()
