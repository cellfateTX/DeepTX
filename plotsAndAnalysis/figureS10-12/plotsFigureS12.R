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

burst2 = read.csv(file = "data/burst_0.2.csv")
burst3 = read.csv(file = "data/burst_0.3.csv")
burst5 = read.csv(file = "data/burst_0.5.csv")

df <- data.frame(
  x1 = burst5$bs,
  y1 = burst3$bs,
  x2 = burst5$bf,
  y2 = burst3$bf,
  x3 = burst5$bs,
  y3 = burst2$bs,
  x4 = burst5$bf,
  y4 = burst2$bf
)


# 创建四个散点图
p1 <- ggplot(df, aes(x = x1, y = y1)) +
  geom_point(color = "steelblue", size =0.1) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1,size = 2) +
  # labs(x = "BS (seq_depth=0.5)", y = "BS (seq_depth=0.3)") +
  labs(
    x = expression(BS(alpha == 0.5)),
    y = expression(BS(alpha == 0.3))
  )+
  theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
        #         legend.position = "none",
        #         legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
  )

p2 <- ggplot(df, aes(x = x2, y = y2)) +
  geom_point(color = "darkgreen", size =0.1) +
  # stat_cor(method = "pearson", label.x = 1, label.y = 1,size = 2) +
  labs(
    x = expression(BF(alpha == 0.5)),
    y = expression(BF(alpha == 0.3))
  )+
  stat_cor(method = "pearson",size = 2) +
  
  theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
        #         legend.position = "none",
        #         legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
  )

p3 <- ggplot(df, aes(x = x3, y = y3)) +
  geom_point(color = "tomato", size =0.1) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1,size = 2) +
  # labs(x = "BS (seq_depth=0.5)", y = "BS (seq_depth=0.2)") +
  labs(
    x = expression(BS(alpha == 0.5)),
    y = expression(BS(alpha == 0.3))
  )+
  theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
        #         legend.position = "none",
        #         legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
  )

p4 <- ggplot(df, aes(x = x4, y = y4)) +
  geom_point(color = "purple", size =0.1) +
  stat_cor(method = "pearson",size = 2) +
  labs(
    x = expression(BF(alpha == 0.5)),
    y = expression(BF(alpha == 0.2))
  )+  theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
        #         legend.position = "none",
        #         legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
  )

# 拼接四个图
combined_plot <- ggarrange(p1, p2, p3, p4,
                           ncol = 4, nrow = 1)
combined_plot

ggsave(str_glue(figure_path,'figureS12.pdf'), width = 6.3, height = 1.5, useDingbats = FALSE)


