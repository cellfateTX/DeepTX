# import R package
options(warn=-1)
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(SingleCellExperiment)
library(ggthemes)
library(stringr)
library(scales)
setwd("D:/academic_relate_code_two/Nessie-main/DeepTXcopy4/DeepTX-main/plotsAndAnalysis/figureS10-12")
figure_path = "D:/academic_relate_code_two/Nessie-main/DeepTXcopy4/DeepTX-main/plotsAndAnalysis/figure/"


group1 <- rnorm(50, mean = 5)
group2 <- rnorm(50, mean = 7)
group3 <- rnorm(50, mean = 6)

data_deepTX <- read.csv("data/helling_result_deepTX.csv")
data_abc <- read.csv("data/hellinger_results_abc.csv")
data_txBurst <- read.csv("data/helling_txburst.csv")

data_abc_helling = data_abc$HellingerDistance
abc_clean <- data_abc_helling[!is.nan(data_abc_helling) & data_abc_helling != 0]



df <- data.frame(
  value = c(data_deepTX$hellinger_distance[1:100], abc_clean[1:100], data_txBurst$X0[1:100]),
  group = factor(rep(c("DeepTX", "txABC", "txBurst"), each = 100))
)
# 
# df <- data.frame(
#   value = c(group1, group2, group3),
#   group = factor(rep(c("DeepTX", "txABC", "txBurst"), each = 50))
# )



# 画箱线图
P1 = ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.01,  # 控制离群点大小
               size = 0.01          # 控制箱线图边框粗细
  ) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group", y = "Hellinger distance") +
  theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
  )
P1

val <- c(0.3, 3, 0.15)
names(val) <- c("DeepTX", "txABC", "txBurst")
bar_df <- data.frame(
  group = factor(names(val), levels = names(val)),
  value = val
)
P2 <- ggplot(bar_df, aes(x = group, y = value, fill = group)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_bw() +
  labs( x = "", y = "Inference efficiency (s/sample)") +
  theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
  )
P2
combined_plot <- ggarrange(P1, P2,ncol = 2, nrow = 1)
combined_plot
ggsave(str_glue(figure_path,'figureS10.pdf'), width = 4, height = 1.5, useDingbats = FALSE)


