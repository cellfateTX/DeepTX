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

library(RColorBrewer)
colormap<- brewer.pal(9,"Blues")[2:6]
colors = c("#D2E0FF","#7BA4FF","#D2E0FF")
# colors = c("#89A1CF","#bac7e4","#F2F3F9")
# load the estimated data
color.feedback <- c( "#547CBE", "#3AA438","#E71419")
point.color = color.feedback[2]
density.color = "darkblue"
vline.color = "green"

trueSolution = read.csv(file = "data/posteriorDist/true_params.csv")

i=12
estimated_file_path = sprintf("data/posteriorDist/AllSolution_nn_%s.csv", i )
estimatedSolution = read.csv(file =estimated_file_path)
density.color = "#AFCCA0"
vline.color = "#AAB7D9"
P1 = ggplot(estimatedSolution, aes(x=x1)) +
geom_density(aes(x = x1, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x1")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(r[on]), y = "Posteriori") +

xlim(0.6,15)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P1
P2 = ggplot(estimatedSolution) +
geom_density(aes(x = x2, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x2")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(k[on]), y =  NULL) +
xlim(0.6,15)+
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
P2
P3 = ggplot(estimatedSolution) +
geom_density(aes(x = x3, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x3")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(r[off]), y = NULL) +

xlim(6,15)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.title.y = element_text(colour = 'black', size = 1),

        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )

P4 = ggplot(estimatedSolution) +
geom_density(aes(x = x4, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x4")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(k[off]), y =NULL) +

xlim(2.5,7.5)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.title.y = element_text(colour = 'black', size = 1),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )

P5 = ggplot(estimatedSolution) +
geom_density(aes(x = x5, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x5")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(mu), y =NULL) +
xlim(50,80)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.title.y = element_text(colour = 'black', size = 1),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P5

picD_sub_one <- ggarrange(P1,P2,P3,P4,P5, ncol = 5, nrow = 1, widths=c(10), heights=c(2), align = "v")
picD_sub_one

ggsave(str_glue(figure_path,'figureS11.pdf'), width = 6.3, height = 1, useDingbats = FALSE)
ggsave(str_glue(figure_path,'figureS11.jpg'), width = 7, height = 1)

