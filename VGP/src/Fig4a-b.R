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

draw_gene_body_lines <- list(
    geom_vline(xintercept = 1.5, linetype="dashed", color = "grey",    size=0.1),
    geom_vline(xintercept = 2.5, linetype="dashed", color = "grey",    size=0.1),
    geom_vline(xintercept = 3.5, linetype="dashed", color = "grey",    size=0.1),
    geom_vline(xintercept = 4.5, linetype="dashed", color = "grey",    size=0.1),
    geom_vline(xintercept = 1,   linetype="solid",  color = "#5e5e5e", size=0.1),
    geom_vline(xintercept = 2,   linetype="solid",  color = "#5e5e5e", size=0.1),
    geom_vline(xintercept = 3,   linetype="solid",  color = "#5e5e5e", size=0.1),
    geom_vline(xintercept = 4,   linetype="solid",  color = "#5e5e5e", size=0.1),
    geom_vline(xintercept = 5,   linetype="solid",  color = "#5e5e5e", size=0.1),
    scale_x_continuous(name="Intron", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'UTR", "First", " ", "Last", "3'UTR"),
                       sec.axis=sec_axis(~. + 0, name="Exon", breaks = c(1,2,3,4,5), labels=c("5'UTR", "First", "Internal", "Last", "3'UTR  ")))
)

shared_theme <- function(my_col, y_breaks) {
    list(
    theme_pubr(),
    scale_color_manual(values = my_col),
    scale_y_continuous(breaks = y_breaks),
    theme(plot.title = element_text(size=6,face='bold', hjust = 0.5, margin=margin(0,0,0,0)),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(size=0.2),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size=5),
          legend.title = element_text(size=5),
          legend.margin=margin(c(0,0,0,0)),
          axis.text.y = element_text(size=5),
          axis.text.x =  element_text(size=5,hjust=0.5,vjust=0.4),
          plot.margin = margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"),
          axis.ticks = element_line(size=0.1),
          axis.title = element_text(size=5)),
    guides(linetype=guide_legend(nrow = 2, keyheight = unit(3, "mm")),
           shape=guide_legend(nrow = 2, keyheight = unit(3, "mm")),
           color = guide_legend(nrow = 2, keyheight = unit(3, "mm")))
    )
}

transform_intron_position <- function(column) {
    column =   ifelse(column == 1, 1.5, # 5'UTR Intron
                ifelse(column == 2, 2.5, # First Intron
                ifelse(column == 3, 3,   # Internal Intron
                ifelse(column == 4, 3.5, # Last Intron
                ifelse(column == 5, 4.5, 100))))) # 3'UTR Intron
    return(column)
}

## Styles
LINE_W = 0.2
G_POINT_STROKE = 0.5
G_POINT_SIZE = 0.7
S_POINT_STROKE = 0.3
S_POINT_SIZE = 0.3
setwd("VGP")
getwd()

############################################
########### [ Shared Variables ] ###########
############################################
exon_order <- c("UTR5", "FirstCDS",  "InternalCDS", "LastCDS", "UTR3")
intron_order <- c("5'UTR-intron","First-intron","Internal-intron", "Last-intron","3'UTR-intron")
genic_order <- c("UTR5", "5'UTR-intron", "FirstCDS","First-intron",  "InternalCDS", "Internal-intron","LastCDS","Last-intron", "UTR3", "3'UTR-intron")

###########################################################
########### [ Fig. 4a GC & Missing Ratio ] ################
###########################################################

up_2Kb_order_c <- c()
for (i in c(1:20)) {
    up_2Kb_order_c <- c(up_2Kb_order_c, paste0("five-prime-", 100*(21-i), "bp"))
}

down_2kb_order_c <- c()
for (i in c(1:20)) {
    down_2kb_order_c <- c(down_2kb_order_c, paste0("three-prime-", 100*(i), "bp"))
}

species_order=c("Anna's hummingbird", "Zebra finch", "Platypus", "Climbing perch")

###############################################
####### [ Fig. 4a 1. Preparing data ] #########
###############################################

zebrafinch_aln_ratio_df <- read.csv("input/Fig4/missing/bTaeGut1.forR.gc_and_aln.tsv", sep="\t")
platypus_aln_ratio_df <- read.csv("input/Fig4/missing/mOrnAna1.forR.gc_and_aln.tsv", sep="\t")
anna_aln_ratio_df <- read.csv("input/Fig4/missing/bCalAnn1.forR.gc_and_aln.tsv", sep="\t")
climbingperch_aln_ratio_df <- read.csv("input/Fig4/missing/fAnaTes1.2.forR.gc_and_aln.tsv", sep="\t")

anna_aln_ratio_df$species <- "Anna's hummingbird"
zebrafinch_aln_ratio_df$species <- "Zebra finch"
platypus_aln_ratio_df$species <- "Platypus"
climbingperch_aln_ratio_df$species <- "Climbing perch"

# Combine all 4 to one table
all_4 <- rbind(zebrafinch_aln_ratio_df, anna_aln_ratio_df, platypus_aln_ratio_df, climbingperch_aln_ratio_df)
all_4$species <- factor(all_4$species, level=species_order)

gc_mean <- cbind(all_4[, -which(names(all_4) %in% c("Not_aln_mean"))], "Found")
colnames(gc_mean) <- c("position", "mean", "species", "sequence")
not_aln_mean <- cbind(all_4[, -which(names(all_4) %in% c("GC_mean"))], "Missing")
colnames(not_aln_mean) <- c("position", "mean", "species", "sequence")
all_4 <- rbind(gc_mean, not_aln_mean)
all_4$sequence <- factor(all_4$sequence, level=c("Found","Missing"))

# All 4 in one plot
# (1) gene body
# (1-1) Exon
aln_exon_df <- all_4 %>% filter(position %in% exon_order)
aln_exon_df$num_position <- as.numeric(factor(aln_exon_df$position, level=exon_order))
aln_exon_df$type <- "Exon"

# (1-2) Intron
aln_intron_df <- all_4 %>%
    filter(position %in% intron_order)
aln_intron_df$num_position <- as.numeric(factor(aln_intron_df$position, level=intron_order))
aln_intron_df$num_position <- transform_intron_position(aln_intron_df$num_position)
aln_intron_df$num_position <- as.numeric(aln_intron_df$num_position)
aln_intron_df$type <- "Intron"

aln_genic_df<- rbind(aln_exon_df, aln_intron_df)
aln_genic_df$type <- factor(aln_genic_df$type, level=c("Exon","Intron"))

# (2-2) Upstream
aln_up_2kb_df <- all_4 %>%
    filter(position %in% up_2Kb_order_c)
aln_up_2kb_df$num_position <- as.numeric(factor(aln_up_2kb_df$position, level=up_2Kb_order_c))
up_2kb_labels <- c("-2","-1","-0")

head(aln_up_2kb_df)
tail(aln_up_2kb_df)

# (2-2) Downstream 
aln_down_2kb_df <- all_4 %>%
    filter(position %in% down_2kb_order_c)
aln_down_2kb_df$num_position <- as.numeric(factor(aln_down_2kb_df$position, level=down_2kb_order_c))
down_2kb_labels <- c("+0","+1","+2")

##############################################
####### [ Fig. 4a 2. Drawing plots ] #########
##############################################
my_col2=c("black", "#d7191c")   # Color by sequence - Found and Missing
y_breaks <- c(0, 10, 20, 30, 40, 50, 60, 70)

gene_body_plot <- ggplot(aln_genic_df, aes(x = num_position, y = mean, color=sequence)) +
    draw_gene_body_lines +
    geom_line( data = aln_intron_df, mapping = aes(y = mean, group = interaction(sequence, species)), linetype = "dashed", lwd=LINE_W) +
    geom_line( data = aln_exon_df, mapping = aes(y = mean, group = interaction(sequence, species)), linetype = "solid", lwd=LINE_W) +
    geom_point(data = aln_intron_df, mapping = aes(group=interaction(sequence, species), shape = species), size=G_POINT_SIZE, stroke = G_POINT_STROKE) +
    geom_point(data = aln_exon_df, mapping = aes(group=interaction(type,species), shape = species), size=G_POINT_SIZE, stroke = G_POINT_STROKE) +
    scale_shape_manual(values = c(2,3,4,8), name = "Species") +
    scale_x_continuous(name="Intron", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'UTR", "First", " ", "Last", "3'UTR"),
                       sec.axis=sec_axis(~. + 0, name="Exon", breaks = c(1,2,3,4,5), labels=c("5'UTR", "First", "Internal", "Last", "3'UTR  "))) +
    coord_cartesian(ylim=c(0, 75), xlim=c(0.75,5.25), expand=FALSE)+
    labs(title = "Gene Body", color = "Sequence") + 
    shared_theme(my_col2, y_breaks)
gene_body_plot

up_plot <- ggplot(aln_up_2kb_df, aes(x = num_position, y = mean, color=sequence)) +
    geom_line(aes(group=interaction(sequence,species)),  lwd=LINE_W, alpha=0.7) +
    geom_point(aes(group=interaction(sequence,species), shape = species), size=G_POINT_SIZE, stroke = G_POINT_STROKE) +
    scale_shape_manual(values=c(2,3,4,8), name = "Species") +
    scale_x_continuous(breaks=c(1,10,20), labels=up_2kb_labels, name=" ", 
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(0, 75), xlim=c(0.5, 20.5), expand=FALSE)+
    labs(color="Sequence", title="Upstream") + 
    shared_theme(my_col2, y_breaks)
up_plot

down_plot <- ggplot(aln_down_2kb_df, aes(x = num_position, y = mean, color=sequence)) +
    geom_line(aes(group=interaction(sequence,species)),  lwd=LINE_W, alpha=0.7) +
    geom_point(aes(group=interaction(sequence,species), shape = species), size=G_POINT_SIZE, stroke = G_POINT_STROKE) +
    scale_shape_manual(values=c(2,3,4,8), name = "Species")+
    scale_x_continuous(breaks=c(1,10,20), labels=down_2kb_labels, name="(kbp)",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(0, 75), xlim=c(0.5, 20.5), expand=FALSE)+
    labs(color="Sequence", title="Downstream") + 
    shared_theme(my_col2, y_breaks) +
    theme(axis.title.x = element_text(size=5, face='bold', hjust=1))
down_plot

gc_missing_plot <- arrangeGrob(up_plot + theme(legend.position="none", axis.title.y = element_blank()),
                                 gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                                 down_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()),
                                 ncol = 3,
                                 left = textGrob("GC Content (%)", rot = 90, vjust = 0, gp = gpar(cex = 0.5, fontface = "bold")))
legend1 = gtable_filter(ggplotGrob(down_plot), "guide-box")

#####################################################################
########### [Fig. 4b GC content of all VGP species ] ################
#####################################################################
summary_df <- read.csv("./input/Fig4/clade_GC_content/gc.clade.tsv", sep="\t")
head(summary_df)
summary_df$GC_mean  = summary_df$GC_mean*100
summary_df$GC_sd    = summary_df$GC_sd * 100

vertebrate_clade <- c("Bird", "Mammal", "Reptile", "Skate", "Amphibian", "Fish")

up_30Kb_order <- c()
for (i in c(1:300)) {
    up_30Kb_order <- c(up_30Kb_order, paste0("five-prime-30kb-", 100*(i), "bp"))
}

down_30kb_order <- c()
for (i in c(1:300)) {
    down_30kb_order <- c(down_30kb_order, paste0("three-prime-30kb-", 100*(i), "bp"))
}

#####################################
##### [ Fig. 4b 1.Genic region ] ####
#####################################

# (1-1) Exon
exon_df <-summary_df %>% filter(position %in% exon_order)
exon_df$numeric_position <- as.numeric(factor(exon_df$position, level=exon_order))
exon_df$clade <- factor(exon_df$clade, level=vertebrate_clade)

# (1-2) Intron
intron_df <-summary_df %>% filter(position %in% intron_order)
intron_df$numeric_position <- as.numeric(factor(intron_df$position, level=intron_order))
intron_df$clade <- factor(intron_df$clade, level=vertebrate_clade)
intron_df$numeric_position <- transform_intron_position(intron_df$numeric_position)

# (1-3) Merge exon and intron dataframes
exon_df$type <- "Exon"
intron_df$type <- "Intron"
genic_df <- rbind(exon_df, intron_df)
genic_df$type <- factor(genic_df$type, level=c("Exon","Intron"))


######################################
#### [ Fig. 4b 2. Up/downstream ] ####
######################################

# (2-1) Upstream 

up_30kb_df <- summary_df %>%
    filter(position %in% up_30Kb_order)
up_30kb_df$numeric_position <- as.numeric(factor(up_30kb_df$position, level=up_30Kb_order))
up_30kb_df$clade <- factor(up_30kb_df$clade, level=vertebrate_clade)
up_30kb_labels <- c("-0","-10"," -20","-30 ")

# (2-2) Downstream 
down_30kb_df <- summary_df %>%
    filter(position %in% down_30kb_order)
down_30kb_df$numeric_position <- as.numeric(factor(down_30kb_df$position, level=down_30kb_order))
down_30kb_df$clade <- factor(down_30kb_df$clade, level=vertebrate_clade)
down_30kb_labels <- c("+0","+10","+20  ", " +30")

head(genic_df)
head(up_30kb_df)
head(down_30kb_df)

######################################
#### [ Fig. 4b 3. Drawing plots ] ####
######################################
my_col=c("#d73027", "#fc8d59", "#fee090", "#762a83", "#1b7837", "#4575b4")
y_breaks = c(30, 40, 50, 60, 70)

# (3-1) gene body plot
gene_body_plot <- ggplot(genic_df, aes(x =numeric_position, y = GC_mean, color=clade)) +
    draw_gene_body_lines +
    geom_line(aes(group=interaction(type,clade), linetype = type), lwd=LINE_W) +
    geom_point(aes(group=interaction(type,clade), shape = type), size=S_POINT_SIZE) +
    geom_errorbar(aes(ymin = GC_mean - GC_sd, ymax = GC_mean + GC_sd, color = clade, group = interaction(type,clade), linetype = type), width=0.1, alpha=0.5, size = LINE_W) +
    scale_shape_manual(values=c(19,1), labels=c("Exon or Up/Downstream", "Intron"))+
    scale_linetype(labels=c("Exon or Up/Downstream", "Intron")) +
    coord_cartesian(ylim=c(25,75), xlim=c(0.75,5.25), expand=FALSE)+
    labs(color="Clade", title = "Gene Body", linetype="Type", shape="Type") + 
    shared_theme(my_col, y_breaks)
gene_body_plot

# (3-2) upstream 30Kbp plot
up_plot <- ggplot(up_30kb_df, aes(x = numeric_position, y = GC_mean, color=clade)) +
    geom_line(aes(group=clade), lwd=LINE_W) +
    geom_point(aes(group=clade), size=S_POINT_SIZE) +
    geom_errorbar(aes(ymin = GC_mean - GC_sd, ymax = GC_mean + GC_sd, color = clade), width = 0.05, alpha = 0.5, size = LINE_W) +
    scale_x_continuous(breaks = c(1,100,200,300), labels = up_30kb_labels, trans = reverselog_trans(10), name="",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(25,75), xlim = c(350, 0.8), expand=FALSE) +
    labs(color="Clade", title="Upstream") + 
    shared_theme(my_col, y_breaks)
up_plot

# (3-3) downstream 30Kbp plot
down_plot <- ggplot(down_30kb_df, aes(x = numeric_position, y = GC_mean, color=clade)) +
    geom_line(aes(group=clade), lwd=LINE_W) +
    geom_point(aes(group=clade), size=S_POINT_SIZE) +
    geom_errorbar(aes(ymin = GC_mean - GC_sd, ymax = GC_mean + GC_sd, color = clade), width = 0.05, alpha = 0.5, size = LINE_W) +
    scale_x_log10(breaks=c(1,100,200,300), labels = down_30kb_labels, name="(kbp)",
                  sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(25,75), xlim = c(0.8, 350), expand=FALSE) +
    labs(color="Clade", title="Downstream", x="(kbp)") + 
    shared_theme(my_col, y_breaks) +
    theme(axis.title.x = element_text(size=5, face='bold', hjust=1))
down_plot

# (3-4) combine all plots into one
## Note: need to change plot size
clade_gc_plot <- arrangeGrob(up_plot + theme(legend.position="none", axis.title.y = element_blank(), plot.title = element_blank()),
                             gene_body_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title = element_blank()),
                             down_plot + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title = element_blank()),
                             ncol = 3,
                             left = textGrob("GC Content (%)", rot = 90, vjust = 0, gp = gpar(cex = 0.5, fontface = "bold")))
#ggsave("output/Fig6_gc_content_by_clades_sd.png", clade_gc_plot, width=6, height=2, units = "in")
#ggsave("output/pub/Fig4b_gc_content_by_clades_sd.pdf", clade_gc_plot, width=120, height=45, units = "mm")

p <- arrangeGrob(gc_missing_plot, clade_gc_plot, nrow = 2)
ggsave("output/pub/Fig4.pdf", p, width=120, height=90, units = "mm", device=cairo_pdf)

legend2 = gtable_filter(ggplotGrob(gene_body_plot), "guide-box")

########################################
#### [ Fig. 4b 3. Drawing legends ] ####
########################################
labels_plot <- arrangeGrob(legend1)
ggsave("output/pub/Fig4_legend1.pdf", labels_plot, width=3.2, height=0.3, units = "in", device=cairo_pdf)
labels_plot <- arrangeGrob(legend2)
ggsave("output/pub/Fig4_legend2.pdf", labels_plot, width=3.5, height=0.3, units = "in", device=cairo_pdf)

