library(ggplot2)
library(scales)

# Read the data from the TSV file
data <- read.table("input/chm13v2.0_num_variants_accessible.tsv", header = TRUE)

# Convert the chromosome column to a factor, use order as levels
data$Chr <- factor(data$Chr, levels = c("1", "2", "3", "4", "5",
                                        "6", "7", "8", "9", "10",
                                        "11", "12", "13", "14", "15",
                                        "16", "17", "18", "19", "20",
                                        "21", "22", "X", "Y"))

# Variants are numeric, need to convert
data$Variants <- as.numeric(data$Variants)

max(data$Variants)

head(data)
# Create the bar plot
ggplot(data, aes(x = Chr, y = Variants, fill = Ref)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#CCCCCC", "#AAAAAA")) +
  theme_classic() +
  theme(
        axis.text = element_text(size = 5),
        legend.position = "none",
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous(expand = c(0, 0), labels = label_comma(),
                     limits = c(0, 100000))

ggsave(file = "output/num_variants_accessible.png",
       width = 70, height = 30, units = "mm")
ggsave(file = "output/num_variants_accessible.pdf",
       width = 70, height = 30, units = "mm")
