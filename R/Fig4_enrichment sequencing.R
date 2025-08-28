#Environmental air monitoring in international airports
#Fig 4 - Enrichment sequencing
##########################################################################

#PRIOR TO RUNNING - DEFINE PATHS TO DATA FILES
#All data available in supplementary data
seq_data_path <- ("../Supplementary_Data_Air_Paper/enrichment_sequencing_data.csv")

library(tidyverse)
library(cowplot)

###############################################################
#Prepare Data

#import data
seq_data <- read.csv(seq_data_path)

#structure date
seq_data$Collection.Date <- as.Date(seq_data$Collection.Date, format="%Y-%m-%d")

#structure airports
summary(as.factor(seq_data$Airport))

##############################################################################
#plot RPKM heatmap by sample


#for plotting, make zero > na
seq_data$RPKM[which(seq_data$RPKM==0)] <- NA

# structure Virus.Name
seq_data$Virus.Name <- as.factor(seq_data$Virus.Name)

#remove samples without collection date
seq_data <- seq_data[-which(is.na(seq_data$Collection.Date)==TRUE),]

#change dates to ordered factor
seq_data$date.fac <- factor(seq_data$Collection.Date, ordered=TRUE)
summary(seq_data$date.fac)
levels(seq_data$date.fac)

#for plotting axis labels function
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

#plot
rpkm_plot <- ggplot(seq_data, aes(x=date.fac,
                                 y=ordered(Virus.Name, levels=rev(levels(Virus.Name))),
                                 fill= RPKM)) +
  geom_tile()+
  scale_fill_distiller(palette = "OrRd", direction=1,
                       na.value="white",
                       name="RPKM") +
  scale_x_discrete(breaks = every_nth(n = 3))+
  facet_wrap(~Airport, nrow=3) +
  #annotate(geom="point", shape=17, x= as.Date("2023-12-18"), y=0.1, color="red", size=5)+
  geom_vline(xintercept="2023-12-02", color="red", size=2)+
  coord_cartesian(clip = "off")+
  theme_cowplot()+
  theme(panel.background = element_rect(fill = 'gray95', colour = 'gray90'),
        axis.text.x = element_text(angle=45, hjust=1, size=26),
        axis.text.y = element_text(size=20),
        axis.title = element_blank(),
        strip.text = element_text(size=32),
        strip.background = element_blank(),
        legend.text = element_text(size=26),
        legend.title = element_text(size=32),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        legend.key.size = unit(1.5, 'cm'))
rpkm_plot


save_plot(filename="Fig4_sequencing.png", rpkm_plot, base_height = 24, base_width = 16)
