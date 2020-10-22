library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(scales)
library(grid)
library(gtable)
library(gridExtra)
library(directlabels)

# From https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2 ####
# Reversed log transformation for 5' Upstream regions
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

dir_with_rscript = "VGP"
getwd()
setwd(dir_with_rscript)
################################################################
########### [1. GC content of all VGP species ] ################
################################################################
summary_df <- read.csv("./input/Fig7/clade_GC_content/forR.gc.clade.tsv", sep="\t")
summary_df$GC_mean = summary_df$GC_mean*100
summary_df$GC_se = summary_df$GC_se * 100


vertebrate_clade <- c("bird", "mammal", "reptile", "skate", "amphibian", "fish")
exon_order <- c("UTR5", "FirstCDS",  "InternalCDS", "LastCDS", "UTR3")
intron_order <- c("5'UTR-intron","First-intron","Internal-intron", "Last-intron","3'UTR-intron")

up_30Kb_order <- c()
for (i in c(1:300)) {
    up_30Kb_order <- c(up_30Kb_order, paste0("five-prime-30kb-", 100*(i), "bp"))
}

down_30kb_order <- c()
for (i in c(1:300)) {
    down_30kb_order <- c(down_30kb_order, paste0("three-prime-30kb-", 100*(i), "bp"))
}

###############################################################
##### [ 1. GC content of all VGP species- 1.Genic region ] ####
###############################################################

# (1-1) Exon
exon_df <-summary_df %>%
    filter(position %in% exon_order)
exon_df$numeric_position <- as.numeric(factor(exon_df$position, level=exon_order))
exon_df$clade <- factor(exon_df$clade, level=vertebrate_clade)

# (1-2) Intron
intron_df <-summary_df %>%
    filter(position %in% intron_order)
intron_df$numeric_position <- as.numeric(factor(intron_df$position, level=intron_order))
intron_df$clade <- factor(intron_df$clade, level=vertebrate_clade)
intron_df$numeric_position <- ifelse(intron_df$numeric_position == 1, 1.5, # 5'UTR Intron
                                     ifelse(intron_df$numeric_position == 2, 2.5, # First Intron
                                            ifelse(intron_df$numeric_position == 3, 3, # Internal Intron
                                                   ifelse(intron_df$numeric_position == 4, 3.5, # Last Intron
                                                          ifelse(intron_df$numeric_position == 5,4.5, 100))))) # 3'UTR Intron

# (1-3) Merge exon and intron dataframes
exon_df$type <- "exon"
intron_df$type <- "intron"
genic_df <- rbind(exon_df, intron_df)
genic_df$type <- factor(genic_df$type, level=c("exon","intron"))


###############################################################
#### [1. GC content of all VGP species- 2. Up/downstream ] ####
###############################################################

# (2-1) Upstream 

up_30kb_df <- summary_df %>%
    filter(position %in% up_30Kb_order)
up_30kb_df$numeric_position <- as.numeric(factor(up_30kb_df$position, level=up_30Kb_order))
up_30kb_df$clade <- factor(up_30kb_df$clade, level=vertebrate_clade)
up_30kb_labels <- c("-0","-10","-20","-30  ")

# (2-2) Downstream 

down_30kb_df <- summary_df %>%
    filter(position %in% down_30kb_order)
down_30kb_df$numeric_position <- as.numeric(factor(down_30kb_df$position, level=down_30kb_order))
down_30kb_df$clade <- factor(down_30kb_df$clade, level=vertebrate_clade)
down_30kb_labels <- c("+0","+10","+20", "  +30")

###############################################################
#### [ 1.GC content of all VGP species- 3. Drawing plots ] ####
###############################################################
# (3-1) gene body plot
gene_body_plot <- ggplot(genic_df, aes(x =numeric_position, y = GC_mean, color=clade)) +
    geom_vline(xintercept = 1.5, linetype="dashed", color = "grey", size=0.1) +
    geom_vline(xintercept = 2.5, linetype="dashed", color = "grey", size=0.1) +
    geom_vline(xintercept = 3.5, linetype="dashed", color = "grey", size=0.1) +
    geom_vline(xintercept = 4.5, linetype="dashed", color = "grey", size=0.1) +
    geom_vline(xintercept = 1,   linetype="solid", color = "#5e5e5e", size=0.1) +
    geom_vline(xintercept = 2,   linetype="solid", color = "#5e5e5e", size=0.1) +
    geom_vline(xintercept = 3,   linetype="solid", color = "#5e5e5e", size=0.1) +
    geom_vline(xintercept = 4,   linetype="solid", color = "#5e5e5e", size=0.1) +
    geom_vline(xintercept = 5,   linetype="solid", color = "#5e5e5e", size=0.1) +
    geom_line(aes(group=interaction(type,clade), linetype = type), lwd=0.5) +
    geom_point(aes(group=interaction(type,clade), shape = type), size=0.5) +
    geom_errorbar(aes(ymin=GC_mean-GC_se,ymax=GC_mean+GC_se, color=clade, group=interaction(type,clade), linetype = type), width=0.1, alpha=0.5) +
    scale_x_continuous(name="Intron", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'UTR", "First   ", "Internal", "   Last", "3'UTR"),
                       sec.axis=sec_axis(~. + 0, name="Exon", breaks = c(1,2,3,4,5), labels=c("5'UTR", "First", "Internal", "Last", "3'UTR"))) + 
    scale_y_continuous(breaks = c(40, 50, 60, 70))+
    scale_shape_manual(values=c(19,1), labels=c("Exon or upstream/downstream sequence", "Intron"))+
    scale_linetype(labels=c("Exon or upstream/downstream sequence", "Intron"))+
    scale_color_discrete(labels=c("Bird", "Mammal", "Reptile", "Skate", "Amphibian", "Fish"))+
    coord_cartesian(ylim=c(35,75), xlim=c(0.75,5.25), expand=FALSE)+#, expand=FALSE)+#coord_cartesian(ylim=c(30,75), expand=FALSE)+
    labs(color="Clade", title = "Gene body", linetype="Type", shape="Type") + 
    theme_pubr() +
    theme(plot.title = element_text(size=8,face='bold'),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(size=0.1),
          legend.direction = "horizontal",
          legend.position = "bottom",
          #legend.key.width = unit(1.5,"cm"),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.text.x =  element_text(size=6,hjust=0.5,vjust=0.4),
          plot.margin = margin(t = 0, r = 0.5, b = 8, l = 0, unit = "pt"),
          axis.ticks = element_line(size=0.1),
          #axis.ticks.length = unit(0.05, "cm"),
          axis.title = element_text(size=6, face='bold')) +     
    guides(linetype=guide_legend(nrow = 2),
           shape=guide_legend(nrow = 2),
           color = guide_legend(nrow = 2))
gene_body_plot

# (3-2) upstream 30Kbp plot
up30_plot <- ggplot(up_30kb_df, aes(x =numeric_position, y = GC_mean, color=clade)) +
    geom_line(aes(group=clade), lwd=0.5) +
    geom_point(aes(group=clade), size=0.5) +
    geom_errorbar(aes(ymin=GC_mean-GC_se,ymax=GC_mean+GC_se, color=clade), width=0.06, alpha=0.3, size=0.5) +
    scale_x_continuous(breaks=c(1,100,200,300), labels=up_30kb_labels, name="(Kbp)", trans=reverselog_trans(10),
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    scale_y_continuous(breaks = c(40, 50, 60, 70))+
    coord_cartesian(ylim=c(35,75), expand=TRUE)+
    # labs(color="Clade", title="Upstream", x="(Kbp)") + 
    theme_pubr()+
    theme(plot.title = element_text(size=8,face='bold'),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(size=0.1),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.width = unit(1.5,"cm"),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.text.x =  element_text(size=6,hjust=0.5),
          plot.margin = margin(t = 0, r = 0.5, b = 8, l = 0, unit = "pt"),
          axis.ticks = element_line(size=0.1),
          #axis.ticks.length = unit(0.05, "cm"),
          axis.title = element_text(size=6, face='bold'),
          axis.title.x = element_text(size=6, face='bold', hjust=1)) +
    guides(color = guide_legend(nrow = 1))
up30_plot

# (3-3) downstream 30Kbp plot
down30_plot <- ggplot(down_30kb_df, aes(x =numeric_position, y = GC_mean, color=clade)) +
    geom_line(aes(group=clade), lwd=0.5) +
    geom_point(aes(group=clade), size=0.5) +
    geom_errorbar(aes(ymin=GC_mean-GC_se,ymax=GC_mean+GC_se, color=clade), width=0.06, alpha=0.3, size=0.5) +
    scale_x_log10(breaks=c(1,100,200,300), labels=down_30kb_labels, name="(Kbp)",#Downstream (Kbp)",
                  sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    scale_y_continuous(breaks = c(40, 50, 60, 70))+
    coord_cartesian(ylim=c(35,75), expand=TRUE)+
    # labs(color="Clade", title="Downstream", x="(Kbp)") + 
    theme_pubr()+
    theme(plot.title = element_text(size=8,face='bold'),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(size=0.1),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.width = unit(1.5,"cm"),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.text.x =  element_text(size=6,hjust=0.5),
          plot.margin = margin(t = 0, r = 4, b = 8, l = 0, unit = "pt"),
          axis.ticks = element_line(size=0.1),
          #axis.ticks.length = unit(0.05, "cm"),
          axis.title = element_text(size=6, face='bold'),
          axis.title.x = element_text(size=6, face='bold', hjust=1)) +
    guides(color = guide_legend(nrow = 1))
down30_plot

# (3-4) combine all plots into one
## Note: need to change plot size
legend = gtable_filter(ggplotGrob(gene_body_plot), "guide-box")
clade_gc_plot <- arrangeGrob(up30_plot + theme(legend.position="none", axis.title.y = element_blank()),
                             gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()), 
                             down30_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                             legend,
                             layout_matrix = rbind(c(1,2,3),
                                                   c(1,2,3),
                                                   c(1,2,3),
                                                   c(1,2,3),
                                                   c(1,2,3),
                                                   c(1,2,3),
                                                   c(4,4,4)),
                             left = textGrob("GC content (%)", rot = 90, vjust = 1, gp = gpar(cex = 0.5, fontface = "bold"))) 
#legend, 
#nrow=2,heights=c(10, 1))
out_png_file <- paste0("./output/Fig7_gc_content_by_clades.png")
ggsave(out_png_file, clade_gc_plot, width=5, height=2.5,units = "in")

###############################################################
################ [ 2. GC & missing ratio ] ####################
###############################################################

genic_order <- c("UTR5", "5'UTR-intron", "FirstCDS","First-intron",  "InternalCDS", "Internal-intron","LastCDS","Last-intron", "UTR3", "3'UTR-intron")

up_2Kb_order_c <- c()
for (i in c(1:20)) {
    up_2Kb_order_c <- c(up_2Kb_order_c, paste0("five-prime-", 100*(21-i), "bp"))
}

down_2kb_order_c <- c()
for (i in c(1:20)) {
    down_2kb_order_c <- c(down_2kb_order_c, paste0("three-prime-", 100*(i), "bp"))
}

###############################################################
####### [ 2. GC & missing ratio - 1. Preparing data ] #########
###############################################################

zebrafinch_aln_ratio_df <- read.csv("./input/missing/bTaeGut1.forR.gc_and_aln.tsv", sep="\t")
platypus_aln_ratio_df <- read.csv("./input/missing/mOrnAna1.forR.gc_and_aln.tsv", sep="\t")
anna_aln_ratio_df <- read.csv("./input/missing/bCalAnn1.forR.gc_and_aln.tsv", sep="\t")
climbingperch_aln_ratio_df <- read.csv("./input/missing/fAnaTes1.2.forR.gc_and_aln.tsv", sep="\t")

zebrafinch_aln_ratio_df$species <- "Zebra finch"
platypus_aln_ratio_df$species <- "Platypus"
anna_aln_ratio_df$species <- "Anna's hummingbird"
climbingperch_aln_ratio_df$species <- "Climbing perch"

###############################################################
####### [ 2. GC & missing ratio - 2. Drawing plots ] #########
###############################################################

# version 1. original layout: 
#    [Title] Species name
##   [Subtitle] Upstream,Gene body,Downstream
###  [X axis title] Exon,Intron, (kbp) (Upside and downside of gene body plot)
missing_plot_maker_original_layout <- function(species_df, species_s) {
    # (1) gene body
    # (1-1) Exon
    aln_exon_df <-species_df %>%
        filter(position %in% exon_order)
    aln_exon_df$num_position <- as.numeric(factor(aln_exon_df$position, level=exon_order))
    # (1-2) Intron
    aln_intron_df <-species_df %>%
        filter(position %in% intron_order)
    aln_intron_df$num_position <- as.numeric(factor(aln_intron_df$position, level=intron_order))
    aln_intron_df$num_position <- ifelse(aln_intron_df$num_position == 1, 1.5, 
                                         ifelse(aln_intron_df$num_position == 2, 2.5,
                                                ifelse(aln_intron_df$num_position == 3, 3,
                                                       ifelse(aln_intron_df$num_position == 4, 3.5,
                                                              ifelse(aln_intron_df$num_position == 5,4.5, 100)))))
    aln_intron_df$num_position <- as.numeric(aln_intron_df$num_position)
    aln_exon_df$num_position <- as.numeric(aln_exon_df$num_position)
    aln_exon_df$type <- "exon"
    aln_intron_df$type <- "intron"
    aln_genic_df<- rbind(aln_exon_df, aln_intron_df)
    aln_genic_df$type <- factor(aln_genic_df$type, level=c("exon","intron"))
    
    # (2-1) Upstream 
    aln_up_2kb_df <- species_df %>%
        filter(position %in% up_2Kb_order_c)
    aln_up_2kb_df$num_position <- as.numeric(factor(aln_up_2kb_df$position, level=up_2Kb_order_c))
    up_2kb_labels <- c("-2","-1","-0")
    
    # (2-2) Downstream 
    aln_down_2kb_df <- species_df %>%
        filter(position %in% down_2kb_order_c)
    aln_down_2kb_df$num_position <- as.numeric(factor(aln_down_2kb_df$position, level=down_2kb_order_c))
    down_2kb_labels <- c("+0","+1","+2")
    colnames(aln_genic_df) <- c("position", "GC content", "Ratio of\nmissing sequences", "species", "num_position","type")
    colnames(aln_up_2kb_df) <- c("position", "GC content", "Ratio of\nmissing sequences", "species", "num_position")
    colnames(aln_down_2kb_df) <- c("position", "GC content", "Ratio of\nmissing sequences", "species", "num_position")
    
    exon_and_intron_melt_df <- melt(aln_genic_df,id.vars=c("position","num_position","type"),
                                    measure.vars=c("GC content","Ratio of\nmissing sequences"))
    up_2kb_melt_df<- melt(aln_up_2kb_df,id.vars=c("position","num_position"),
                          measure.vars=c("GC content","Ratio of\nmissing sequences"))
    down_2kb_melt_df<- melt(aln_down_2kb_df,id.vars=c("position","num_position"),
                            measure.vars=c("GC content","Ratio of\nmissing sequences"))
    
    aln_gene_body_plot <- ggplot(exon_and_intron_melt_df, aes(x =num_position, y = value, color=variable)) +
        geom_line(aes(linetype = type), lwd=0.1) +
        geom_point(aes(shape = type), size=0.01) +
        geom_vline(xintercept = 1.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 2.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 3.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 4.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 1, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 2, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 3, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 4, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 5, linetype="solid", color = "#5e5e5e", size=0.05) +
        scale_x_continuous(name="Intron", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'", "F", "I", "L", "3'"),
                           sec.axis=sec_axis(~. + 0, name="Exon", breaks = c(1,2,3,4,5), labels=c("5'", "F", "I", "L", "3'"))) + 
        scale_y_continuous(breaks = c(20,40,60))+
        scale_shape_manual(values=c(19,1), labels=c("Exon or upstream/downstream sequence", "Intron"))+
        scale_color_manual(values=c("red","black"),labels=c("GC content", "Ratio of missing sequences"))+
        scale_linetype_manual(values=c("solid","dotted"), labels=c("Exon or upstream/downstream sequence", "Intron"))+
        coord_cartesian(ylim=c(0,75), xlim=c(0.75,5.25), expand=FALSE)+#, expand=FALSE)+#coord_cartesian(ylim=c(30,75), expand=FALSE)+
        labs(color="Color", title = "Gene body", linetype="Type", shape="Type") + 
        theme_pubr() +
        theme(plot.title = element_text(size=6,face='bold'),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.background = element_blank(),
              legend.text = element_text(size=6),
              legend.title = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.text.x =  element_text(size=6,hjust=0.5,vjust=0.4),
              plot.margin = margin(t = 0, r = 0.5, b = 2, l = 0, unit = "pt"),
              axis.line = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              axis.title = element_text(size=6, face='bold')) +
        guides(linetype=guide_legend(nrow = 2),
               shape=guide_legend(nrow = 2),
               color = guide_legend(nrow = 2))
    
    aln_up2kb_plot <- ggplot(up_2kb_melt_df, aes(x =num_position, y = value, color=variable)) +
        geom_line(lwd=0.1) +
        #geom_dl(aes(label=paste0(variable,"\n"), color=variable), method = list("first.points",cex = 0.45, hjust=0 , vjust=0))+
        geom_point(size=0.01) +
        scale_x_continuous(breaks=c(1,10,20), labels=c("-2","-1","0"), name="(Kbp)", 
                           sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
        scale_y_continuous(breaks = c(20,40,60))+
        scale_color_manual(values=c("red","black"),labels=c("GC content","Ratio of missing sequences"))+
        coord_cartesian(ylim=c(0,75),expand=TRUE)+
        labs(color="Color", title = "Upstream", x="(Kbp)") + 
        theme_pubr() +
        theme(plot.title = element_text(size=6,face='bold'),
              strip.text = element_blank(),
              strip.background = element_blank(),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.text = element_text(size=6),
              legend.title = element_text(size=6),
              axis.line = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              axis.text.y = element_text(size=6),
              axis.text.x =  element_text(size=6),
              plot.margin = margin(t = 0, r = 0.5, b = 2, l = 0, unit = "pt"),
              axis.title = element_text(size=6, hjust=1,face='bold'))+
        guides(color = guide_legend(nrow = 1))
    
    
    aln_down2kb_plot <- ggplot(down_2kb_melt_df, aes(x =num_position, y = value, color=variable)) +
        geom_line(lwd=0.1) +
        geom_point(size=0.01) +
        scale_x_continuous(breaks=c(1,10,20), labels=c("0","+1","+2"), name="(Kbp)",
                           sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
        scale_y_continuous(breaks = c(20,40,60))+
        scale_color_manual(values=c("red","black"),labels=c("GC content", "Ratio of missing sequences"))+
        coord_cartesian(ylim=c(0,75),expand=TRUE)+
        labs(color="Color", title = "Downstream", x="(Kbp)") + 
        theme_pubr() +
        theme(plot.title = element_text(size=6,face='bold'),
              strip.text = element_blank(),
              strip.background = element_blank(),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.text = element_text(size=6),
              legend.title = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.text.x = element_text(size=6),
              axis.line = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              plot.margin = margin(t = 0, r = 4, b = 2, l = 0, unit = "pt"),
              axis.title = element_text(size=6, hjust=1,face='bold'))+
        guides(color = guide_legend(nrow = 1))
    
    legend = gtable_filter(ggplotGrob(aln_gene_body_plot), "guide-box")
    species_title = textGrob(species_s, vjust = 0.5, gp = gpar(fontface = "bold", cex = 0.6))
    gc_n_missing_plot <- arrangeGrob(aln_up2kb_plot + theme(legend.position="none", axis.title.y = element_blank()),
                                     aln_gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),#, axis.ticks.y = element_blank(), axis.line.y.left = element_blank()), 
                                     aln_down2kb_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),#, axis.ticks.y = element_blank(), axis.line.y.left = element_blank()),
                                     ncol = 3,
                                     left = textGrob("Percentage (%)", rot = 90, vjust = 0, gp = gpar(cex = 0.5, fontface = "bold")),
                                     top = species_title)
    #))
    #legend, 
    #nrow=2,heights=c(10, 1))
    out_png_file <- paste0("./output/", species_s, ".gc_n_missing.png")
    print(out_png_file)
    ggsave(out_png_file, gc_n_missing_plot, width=3, height=1.8,units = "in")
    if (species_s == "Climbing perch") {
        plot_list <- list(species_title,
                          aln_up2kb_plot + theme(legend.position="none", axis.title.y = element_blank()),
                          aln_gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                          aln_down2kb_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                          legend)
    } else {plot_list <- list(species_title,
                              aln_up2kb_plot + theme(legend.position="none", axis.title.y = element_blank()),
                              aln_gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                              aln_down2kb_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()))}
    return(plot_list)
}

finch_plot_list <- missing_plot_maker_original_layout(species_df=zebrafinch_aln_ratio_df, species_s="Zebra finch")
anna_plot_list <- missing_plot_maker_original_layout(species_df=anna_aln_ratio_df, species_s="Anna's hummingbird")
platypus_plot_list <- missing_plot_maker_original_layout(species_df=platypus_aln_ratio_df, species_s="Platypus")
perch_plot_list <- missing_plot_maker_original_layout(species_df=climbingperch_aln_ratio_df, species_s="Climbing perch")

all_species_plots <- c(finch_plot_list, anna_plot_list, platypus_plot_list, perch_plot_list)

all_species_gc_n_missing_plot <- arrangeGrob(grobs=all_species_plots, 
                                             layout_matrix = rbind(c(1,1,1,5,5,5),
                                                                   c(2,3,4,6,7,8),
                                                                   c(2,3,4,6,7,8),
                                                                   c(2,3,4,6,7,8),
                                                                   c(2,3,4,6,7,8),
                                                                   c(9,9,9,13,13,13), 
                                                                   c(10,11,12,14,15,16), 
                                                                   c(10,11,12,14,15,16),                                               
                                                                   c(10,11,12,14,15,16),
                                                                   c(10,11,12,14,15,16),
                                                                   c(17,17,17,17,17,17)),
                                             left = textGrob("Percentage (%)", rot = 90, vjust = 0, gp = gpar(cex = 0.5, fontface = "bold")))#,

missing_out_png_file <- paste0("./output/all_species.gc_n_missing.png")
print(missing_out_png_file)
ggsave(missing_out_png_file, all_species_gc_n_missing_plot, width=6, height=5,units = "in")

# version 2. modified layout: try to reduce the height of plots
#    [Title] Species name
##   [X axis title] Upstream (kbp),Gene body (Exon), Gene body (Intron), Downstream (kbp)
missing_plot_maker_modified_layout <- function(species_df, species_s) {
    # (1) gene body
    # (1-1) Exon
    aln_exon_df <-species_df %>%
        filter(position %in% exon_order)
    aln_exon_df$num_position <- as.numeric(factor(aln_exon_df$position, level=exon_order))
    # (1-2) Intron
    aln_intron_df <-species_df %>%
        filter(position %in% intron_order)
    aln_intron_df$num_position <- as.numeric(factor(aln_intron_df$position, level=intron_order))
    aln_intron_df$num_position <- ifelse(aln_intron_df$num_position == 1, 1.5, 
                                         ifelse(aln_intron_df$num_position == 2, 2.5,
                                                ifelse(aln_intron_df$num_position == 3, 3,
                                                       ifelse(aln_intron_df$num_position == 4, 3.5,
                                                              ifelse(aln_intron_df$num_position == 5,4.5, 100)))))
    aln_intron_df$num_position <- as.numeric(aln_intron_df$num_position)
    aln_exon_df$num_position <- as.numeric(aln_exon_df$num_position)
    aln_exon_df$type <- "exon"
    aln_intron_df$type <- "intron"
    aln_genic_df<- rbind(aln_exon_df, aln_intron_df)
    aln_genic_df$type <- factor(aln_genic_df$type, level=c("exon","intron"))
    
    # (2-1) Upstream 
    aln_up_2kb_df <- species_df %>%
        filter(position %in% up_2Kb_order_c)
    aln_up_2kb_df$num_position <- as.numeric(factor(aln_up_2kb_df$position, level=up_2Kb_order_c))
    up_2kb_labels <- c("-2","-1","-0")
    
    # (2-2) Downstream 
    aln_down_2kb_df <- species_df %>%
        filter(position %in% down_2kb_order_c)
    aln_down_2kb_df$num_position <- as.numeric(factor(aln_down_2kb_df$position, level=down_2kb_order_c))
    down_2kb_labels <- c("+0","+1","+2")
    colnames(aln_genic_df) <- c("position", "GC content", "Ratio of\nmissing sequences", "species", "num_position","type")
    colnames(aln_up_2kb_df) <- c("position", "GC content", "Ratio of\nmissing sequences", "species", "num_position")
    colnames(aln_down_2kb_df) <- c("position", "GC content", "Ratio of\nmissing sequences", "species", "num_position")
    
    exon_and_intron_melt_df <- melt(aln_genic_df,id.vars=c("position","num_position","type"),
                                    measure.vars=c("GC content","Ratio of\nmissing sequences"))
    up_2kb_melt_df<- melt(aln_up_2kb_df,id.vars=c("position","num_position"),
                          measure.vars=c("GC content","Ratio of\nmissing sequences"))
    down_2kb_melt_df<- melt(aln_down_2kb_df,id.vars=c("position","num_position"),
                            measure.vars=c("GC content","Ratio of\nmissing sequences"))
    
    aln_gene_body_plot <- ggplot(exon_and_intron_melt_df, aes(x =num_position, y = value, color=variable)) +
        geom_line(aes(linetype = type), lwd=0.1) +
        geom_point(aes(shape = type), size=0.01) +
        geom_vline(xintercept = 1.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 2.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 3.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 4.5, linetype="dashed", color = "grey", size=0.05) +
        geom_vline(xintercept = 1, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 2, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 3, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 4, linetype="solid", color = "#5e5e5e", size=0.05) +
        geom_vline(xintercept = 5, linetype="solid", color = "#5e5e5e", size=0.05) +
        scale_x_continuous(name="Gene body (intron)", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'", "F", "I", "L", "3'"),
                           sec.axis=sec_axis(~. + 0, name="Gene body (exon)", breaks = c(1,2,3,4,5), labels=c("5'", "F", "I", "L", "3'"))) + 
        scale_y_continuous(breaks = c(20,40,60))+
        scale_shape_manual(values=c(19,1), labels=c("Exon or upstream/downstream sequence", "Intron"))+
        scale_color_manual(values=c("red","black"),labels=c("GC content", "Ratio of missing sequences"))+
        scale_linetype_manual(values=c("solid","dotted"), labels=c("Exon or upstream/downstream sequence", "Intron"))+
        coord_cartesian(ylim=c(0,75), xlim=c(0.75,5.25), expand=FALSE)+#, expand=FALSE)+#coord_cartesian(ylim=c(30,75), expand=FALSE)+
        labs(color="Color", title = species_s, linetype="Type", shape="Type") + 
        theme_pubr() +
        theme(plot.title = element_text(size=6,face='bold',hjust=0.5),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.background = element_blank(),
              legend.text = element_text(size=6),
              legend.title = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.text.x =  element_text(size=6,hjust=0.5,vjust=0.4),
              plot.margin = margin(t = 0, r = 0.5, b = 3, l = 0, unit = "pt"),
              axis.line = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              axis.title = element_text(size=6, face='bold')) +
        guides(linetype=guide_legend(nrow = 2),
               shape=guide_legend(nrow = 2),
               color = guide_legend(nrow = 2))
    
    aln_up2kb_plot <- ggplot(up_2kb_melt_df, aes(x =num_position, y = value, color=variable)) +
        geom_line(lwd=0.1) +
        #geom_dl(aes(label=paste0(variable,"\n"), color=variable), method = list("first.points",cex = 0.45, hjust=0 , vjust=0))+
        geom_point(size=0.01) +
        scale_x_continuous(breaks=c(1,10,20), labels=c("-2","-1","0"), name="Upstream (Kbp)", 
                           sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
        scale_y_continuous(breaks = c(20,40,60))+
        scale_color_manual(values=c("red","black"),labels=c("GC content","Ratio of missing sequences"))+
        coord_cartesian(ylim=c(0,75),expand=TRUE)+
        labs(color="Color", title = "", x="Upstream (Kbp)") + 
        theme_pubr() +
        theme(plot.title = element_text(size=6,face='bold'),
              strip.text = element_blank(),
              strip.background = element_blank(),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.text = element_text(size=6),
              legend.title = element_text(size=6),
              axis.line = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              axis.text.y = element_text(size=6),
              axis.text.x =  element_text(size=6),
              plot.margin = margin(t = 0, r = 0.5, b = 3, l = 0, unit = "pt"),
              axis.title = element_text(size=6, hjust=0.5,face='bold'))+
        guides(color = guide_legend(nrow = 1))
    
    
    aln_down2kb_plot <- ggplot(down_2kb_melt_df, aes(x =num_position, y = value, color=variable)) +
        geom_line(lwd=0.1) +
        geom_point(size=0.01) +
        scale_x_continuous(breaks=c(1,10,20), labels=c("0","+1","+2"), name="Downstream (Kbp)",
                           sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
        scale_y_continuous(breaks = c(20,40,60))+
        scale_color_manual(values=c("red","black"),labels=c("GC content", "Ratio of missing sequences"))+
        coord_cartesian(ylim=c(0,75),expand=TRUE)+
        labs(color="Color", title = " ", x="Downstream (Kbp)") + 
        theme_pubr() +
        theme(plot.title = element_text(size=6,face='bold'),
              strip.text = element_blank(),
              strip.background = element_blank(),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.text = element_text(size=6),
              legend.title = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.text.x = element_text(size=6),
              axis.line = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              plot.margin = margin(t = 0, r = 4, b = 3, l = 0, unit = "pt"),
              axis.title = element_text(size=6, hjust=0.5,face='bold'))+
        guides(color = guide_legend(nrow = 1))
    
    legend = gtable_filter(ggplotGrob(aln_gene_body_plot), "guide-box")
    species_title = textGrob(species_s, hjust = 0.5, vjust = 0, gp = gpar(fontface = "bold", cex = 0.6))
    gc_n_missing_plot <- arrangeGrob(aln_up2kb_plot + theme(legend.position="none", axis.title.y = element_blank()),
                                     aln_gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),#, axis.ticks.y = element_blank(), axis.line.y.left = element_blank()), 
                                     aln_down2kb_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),#, axis.ticks.y = element_blank(), axis.line.y.left = element_blank()),
                                     ncol = 3,
                                     left = textGrob("Percentage (%)", rot = 90, vjust = 0, gp = gpar(cex = 0.5, fontface = "bold")))
    
    out_png_file <- paste0("./output/", species_s, ".modified.gc_n_missing.png")
    print(out_png_file)
    ggsave(out_png_file, gc_n_missing_plot, width=3, height=1.8,units = "in")
    if (species_s == "Climbing perch") {
        plot_list <- list(aln_up2kb_plot + theme(legend.position="none", axis.title.y = element_blank()),
                          aln_gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                          aln_down2kb_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                          legend)
    } else {plot_list <- list(aln_up2kb_plot + theme(legend.position="none", axis.title.y = element_blank()),
                              aln_gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                              aln_down2kb_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()))}
    return(plot_list)
}

finch_plot_list <- missing_plot_maker_modified_layout(species_df=zebrafinch_aln_ratio_df, species_s="Zebra finch")
anna_plot_list <- missing_plot_maker_modified_layout(species_df=anna_aln_ratio_df, species_s="Anna's hummingbird")
platypus_plot_list <- missing_plot_maker_modified_layout(species_df=platypus_aln_ratio_df, species_s="Platypus")
perch_plot_list <- missing_plot_maker_modified_layout(species_df=climbingperch_aln_ratio_df, species_s="Climbing perch")
all_species_plots <- c(finch_plot_list, anna_plot_list, platypus_plot_list, perch_plot_list)

all_species_gc_n_missing_plot <- arrangeGrob(grobs=all_species_plots, 
                                             layout_matrix = rbind(c(1,2,3,4,5,6),
                                                                   c(1,2,3,4,5,6),
                                                                   c(1,2,3,4,5,6),
                                                                   c(1,2,3,4,5,6),
                                                                   c(1,2,3,4,5,6),
                                                                   c(7,8,9,10,11,12), 
                                                                   c(7,8,9,10,11,12), 
                                                                   c(7,8,9,10,11,12),                                               
                                                                   c(7,8,9,10,11,12),
                                                                   c(7,8,9,10,11,12),
                                                                   c(13,13,13,13,13,13)),
                                             left = textGrob("Percentage (%)", rot = 90, vjust = 0, gp = gpar(cex = 0.5, fontface = "bold")))#,

missing_out_png_file <- paste0("./output/all_species.modified.gc_n_missing.png")
print(missing_out_png_file)
ggsave(missing_out_png_file, all_species_gc_n_missing_plot, width=6, height=4,units = "in")



###############################################################
#########[ 2-1. GC & missing ratio: Scatter plot ] ############
###############################################################

all_species_df <- rbind(zebrafinch_aln_ratio_df, platypus_aln_ratio_df, anna_aln_ratio_df, climbingperch_aln_ratio_df)
genic_order <- c("UTR5", "5'UTR-intron", "FirstCDS","First-intron",  "InternalCDS", "Internal-intron","LastCDS","Last-intron", "UTR3", "3'UTR-intron")
all_order = c(up_2Kb_order_c, genic_order, down_2kb_order_c)

gc_n_missing_df <- all_species_df %>%
    filter(position %in% all_order)
gc_n_missing_df$position <- factor(gc_n_missing_df$position, level=all_order)
gc_n_missing_df$general_type <- ifelse(gc_n_missing_df$position %in% up_2Kb_order_c, "upstream",
                                       ifelse(gc_n_missing_df$position %in% down_2kb_order_c, "downstream",
                                              ifelse(gc_n_missing_df$position %in% genic_order, "gene body", "")))
gc_n_missing_df$general_type <- factor(gc_n_missing_df$general_type,levels=c("upstream", "gene body", "downstream"))


one_panel_plot <- ggplot(gc_n_missing_df, aes(x=GC_mean, y=Not_aln_mean, shape=general_type, color=species)) + 
    geom_point(size=3) + 
    scale_shape_manual(values=c(16,3,2))+
    theme_pubr() + 
    labs(x="GC content (%)", y="Ratio of missing sequences (%)", shape="Position", color="Species") + 
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(), 
          text=element_text(size=6))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
one_panel_missing_png_file <- "./output/single_panel.modified.gc_n_missing.png"
print(one_panel_plot)
ggsave(one_panel_missing_png_file, one_panel_plot, width=4, height=5,units = "in")