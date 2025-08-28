#Environmental air monitoring in international airports
#Fig 2 - weekly positivity by sample type
##########################################################################

#PRIOR TO RUNNING - DEFINE PATHS TO DATA FILES
#All data available in supplementary data
air_data_path <- ("../Supplementary_Data_Air_Paper/air_PCR_data.csv")
nasal_data_path <- ("../Supplementary_Data_Air_Paper/nasal_PCR_data.csv")
ww_data_path <- ("../Supplementary_Data_Air_Paper/wastewater_PCR_data.csv")
NREVSS_data_path <- ("../Supplementary_Data_Air_Paper/clinical_NREVSS_data.csv")

#load libraries
library(tidyverse)
library(cowplot)

##############################################################################
#Prepare Air PCR Data

#import air PCR data
air_data <- read.csv(air_data_path)

#change Airport column name
colnames(air_data)[which(colnames(air_data)=="Sampler.Name")] <- "Airport"

#format collection date
air_data$collection.date <- as.Date(air_data$Start.Time, format="%Y-%m-%d")

#add collection week
air_data <- air_data %>%
  mutate(collection.week = collection.date - wday(collection.date) + 1)

#rename pathogens
air_data <- air_data %>%
  mutate(Assay.Target = case_when(
    Assay.Target == "SARS-CoV-2" ~ "SARS-CoV-2",
    Assay.Target == "RNaseP" ~ "RNaseP",
    Assay.Target == "Influenza A virus M gene" ~ "Influenza A",
    TRUE ~ Assay.Target
  ))

#subset to SFO and IAD (sites with overlapping data from nasal and WW)
air_data$Airport[which(air_data$Airport=="SFO1")] <- "SFO"
air_data$Airport[which(air_data$Airport=="SFO2")] <- "SFO"
air_data <- air_data[which(air_data$Airport%in%c("SFO", "IAD")),]

#summarize percent positivity by week
summary_air_data <- air_data %>%
  group_by(collection.week, Assay.Target) %>%
  summarize(
    total_samples = n(),
    positive_samples = sum(Qualitative.Result == "Positive", na.rm = TRUE),
    percent_positive = (positive_samples / total_samples) * 100,
    .groups = "drop"
  )

#add sample type
summary_air_data$sample_type <- rep("air", nrow(summary_air_data))

#subset air data to flu and sars-cov-2
air_flu <- summary_air_data[which(summary_air_data$Assay.Target=="Influenza A"),]
air_SC2 <- summary_air_data[which(summary_air_data$Assay.Target=="SARS-CoV-2"),]


#####################################################################################
#Prepare Wastewater PCR Data

#import WW pcr data
WW_data <- read.csv(ww_data_path)
head(WW_data)

#format date
WW_data$SAMPLE_COLLECTED_UTC <- as.Date(WW_data$SAMPLE_COLLECTED_UTC, format="%Y-%m-%d")

#rename pathogens
WW_data <- WW_data %>%
  mutate(Assay.Target = case_when(
    PATHOGEN_NAME == "SARS-CoV-2" ~ "SARS-CoV-2",
    PATHOGEN_NAME == "FLUAV" ~ "Influenza A",
    TRUE ~ PATHOGEN_NAME
  ))

#add collection week
WW_data <- WW_data %>%
  mutate(collection.week = SAMPLE_COLLECTED_UTC - wday(SAMPLE_COLLECTED_UTC) + 1)

#rename result column
WW_data <- WW_data %>%
  mutate(Qualitative.Result = case_when(
    FINAL_PCR_RESULT == "Positive" ~ "Positive",
    TRUE ~ FINAL_PCR_RESULT
  ))


#summarize percent positivity by week
summary_ww_data <- WW_data %>%
  group_by(collection.week, Assay.Target) %>%
  summarize(
    total_samples = n(),
    positive_samples = sum(Qualitative.Result == "Positive", na.rm = TRUE),
    percent_positive = (positive_samples / total_samples) * 100,
    .groups = "drop"
  )

#add sample type
summary_ww_data$sample_type <- rep("ww", nrow(summary_ww_data))

#subset ww data by pathogen
ww_flu <- summary_ww_data[which(summary_ww_data$Assay.Target=="Influenza A"),]
ww_SC2 <- summary_ww_data[which(summary_ww_data$Assay.Target=="SARS-CoV-2"),]

#####################################################################################
#Prepare Nasal PCR Data

#import individual nasal pcr data
nasal_data <- read.csv(nasal_data_path)

#rename pathogens
nasal_data <- nasal_data %>%
  mutate(Assay.Target = case_when(
    PCR_PATHOGEN == "Severe acute respiratory syndrome coronavirus 2" ~ "SARS-CoV-2",
    PCR_PATHOGEN == "Influenza A virus" ~ "Influenza A",
    TRUE ~ PCR_PATHOGEN
  ))
head(nasal_data)

#format collection date
nasal_data$collection.date <- as.Date(nasal_data$COLLECTION_DATE, format="%Y-%m-%d")

#add collection week
nasal_data <- nasal_data %>%
  mutate(collection.week = collection.date - wday(collection.date) + 1)

#rename result column
nasal_data <- nasal_data %>%
  mutate(Qualitative.Result = case_when(
    PCR_STATUS == "positive" ~ "Positive",
    TRUE ~ PCR_STATUS
  ))
head(nasal_data)

#summarize percent positivity by week
summary_nasal_data <- nasal_data %>%
  group_by(collection.week, Assay.Target) %>%
  summarize(
    total_samples = n(),
    positive_samples = sum(Qualitative.Result == "Positive", na.rm = TRUE),
    percent_positive = (positive_samples / total_samples) * 100,
    .groups = "drop"
  )

#add sample type
summary_nasal_data$sample_type <- rep("nasal", nrow(summary_nasal_data))

#subset nasal data by pathogen
nasal_flu <- summary_nasal_data[which(summary_nasal_data$Assay.Target=="Influenza A"),]
nasal_SC2 <- summary_nasal_data[which(summary_nasal_data$Assay.Target=="SARS-CoV-2"),]

#concatenate data
flu_data <- rbind(nasal_flu, air_flu, ww_flu)

#####################################################################################
#Prepare NREVSS data

#import clinical pcr data
clinical_data <- read.csv(NREVSS_data_path)

#format date
clinical_data$collection.week <- as.Date(clinical_data$collection.week)

#subset to flu and SC2
clinical_flu <- clinical_data[which(clinical_data$Assay.Target=="Influenza A"),]
clinical_SC2 <- clinical_data[which(clinical_data$Assay.Target=="SARS-CoV-2"),]

#concatenate data for plotting
flu_data <- rbind(nasal_flu, air_flu, ww_flu, clinical_flu)
SC2_data <- rbind(nasal_SC2, air_SC2, ww_SC2, clinical_SC2)

#################################################################
#Figure 2A: Flu and SC2 Weekly Percent Positive

#define color palette
colors <- c("#f28d2b", "#e15758", "#75b8b2", "#58a14e")


#Figure 2A: Flu and SC2 Weekly Percent Positive
all_data <- rbind(SC2_data, flu_data)
flu_SC2_plot <- ggplot(data=all_data, 
                          aes(x=collection.week, y=percent_positive,
                              color=sample_type)) +
  facet_wrap(~Assay.Target)+
  geom_line(size=1.2)+
  scale_color_manual(values=colors, name="Sample Type",
                     labels=c("airport air (TGS)", "clinical (NREVSS)", "traveler nasal (TGS)", "aviation wastewater (TGS)")
  )+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  scale_y_continuous(name="Weekly Percent Positive", limits = c(0,100))+
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_blank(),
        plot.title = element_text(size=25, hjust=0.5),
        strip.text = element_text(size=20), 
        strip.background = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
flu_SC2_plot

###########################################################################
#Compare flu % positivity across modalities

#pivot_wider
flu_comp <- pivot_wider(data=flu_data, id_cols = collection.week,
                    names_from = sample_type, values_from = percent_positive)

#check correction air and ww
cor(flu_comp$ww, flu_comp$air, method = "pearson", use = "complete.obs")
cor.test(flu_comp$ww, flu_comp$air, method = "pearson")

#check correction air and clinical
cor(flu_comp$clinical, flu_comp$air, method = "pearson", use = "complete.obs")
cor.test(flu_comp$clinical, flu_comp$air, method = "pearson")

#check correction air and nasal
cor(flu_comp$nasal, flu_comp$air, method = "pearson", use = "complete.obs")
cor.test(flu_comp$nasal, flu_comp$air, method = "pearson")

#Fig 2E: flu air and ww
flu_corr <- ggplot(data=flu_comp, aes(x=air, y=ww)) +
  geom_point(size=2, shape=1, stroke=1.5) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  scale_x_continuous(name="Airport Air % Positivity")+
  scale_y_continuous(name="Aviation WW % Positivity")+
  ggtitle("Influenza A")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        plot.title = element_text(size=17, hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
flu_corr

#Fig 2F: flu air and nasal
flu_corr2 <- ggplot(data=flu_comp, aes(x=air, y=clinical)) +
  geom_point(size=2, shape=1, stroke=1.5) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  scale_x_continuous(name="Airport Air % Positivity")+
  scale_y_continuous(name="Clinical % Positivity (NREVSS)")+
  ggtitle("Influenza A")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        plot.title = element_text(size=17, hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
flu_corr2

#Fig 2G: flu air and clinical
flu_corr3 <- ggplot(data=flu_comp, aes(x=air, y=nasal)) +
  geom_point(size=2, shape=1, stroke=1.5) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  scale_x_continuous(name="Airport Air % Positivity")+
  scale_y_continuous(name="Traveler Nasal % Positivity")+
  ggtitle("Influenza A")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        plot.title = element_text(size=17, hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
flu_corr3

#create flu plot grid
flu_correlations <- plot_grid(flu_corr, flu_corr3, flu_corr2, nrow=1, labels = c("E", "F", "G"),
                              label_size=17)
flu_correlations


###########################################################################
#Compare SC2 values across modalities

#Prepare data frame with PCR values for each sample type

#subset SC2 for air data
air_SC2_2 <- air_data[which(air_data$Assay.Target=="SARS-CoV-2"),]
#exclude negative samples
air_SC2_2 <- air_SC2_2[which(air_SC2_2$Qualitative.Result=="Positive"),]

#summarize air data mean PCR value by week
summary2_air_data <- air_SC2_2 %>%
  group_by(collection.week, Assay.Target) %>%
  summarize(
    mean_PCR = mean(Average.Result),
    .groups = "drop"
  )
head(summary2_air_data)

colnames(summary2_air_data)[which(colnames(summary2_air_data)=="mean_PCR")] <- "Air.PCR.Value"

#subset SC2 for ww data
WW_SC2_2 <- WW_data[which(WW_data$Assay.Target=="SARS-CoV-2"),]
#subset WW to positive samples
WW_SC2_2 <- WW_SC2_2[which(WW_SC2_2$Qualitative.Result=="Positive"),]

#summarize WW data mean PCR value by week
summary2_WW_data <- WW_SC2_2 %>%
  group_by(collection.week, Assay.Target) %>%
  summarize(
    mean_PCR = mean(PCR_TARGET_CONC),
    .groups = "drop"
  )

colnames(summary2_WW_data)[which(colnames(summary2_WW_data)=="mean_PCR")] <- "WW.PCR.Value"

#clincal
clin <- clinical_SC2
colnames(clin)[which(colnames(clin)=="percent_positive")] <- "Clinical.Percent.Pos"

#nasal
nas_ind <- summary_nasal_data[which(summary_nasal_data$Assay.Target=="SARS-CoV-2"),]
colnames(nas_ind)[which(colnames(nas_ind)=="percent_positive")] <- "Nasal.Percent.Pos"
nas_ind <- nas_ind[,c(1,5)]

#merge
airport_data <- merge(summary2_air_data, summary2_WW_data,  by="collection.week")
clin$collection.week <- as.Date(clin$collection.week)
nas_ind$collection.week <- as.Date(nas_ind$collection.week)
airport_data <- merge(airport_data, clin, by="collection.week")
airport_data <- merge(airport_data, nas_ind, by="collection.week")

#check correction air and ww
cor(airport_data$Air.PCR.Value, log10(airport_data$WW.PCR.Value),  method = "pearson", use = "complete.obs")
cor.test(airport_data$WW.PCR.Value, log10(airport_data$WW.PCR.Value), method = "pearson")

#check correction air and clinical
cor(airport_data$Air.PCR.Value, airport_data$Clinical.Percent.Pos, method = "pearson", use = "complete.obs")
cor.test(airport_data$Air.PCR.Value, airport_data$Clinical.Percent.Pos, method = "pearson")

#check correction air and nasal
cor(airport_data$Air.PCR.Value, airport_data$Nasal.Percent.Pos, method = "pearson", use = "complete.obs")
cor.test(airport_data$Air.PCR.Value, airport_data$Nasal.Percent.Pos, method = "pearson")

#Fig 2B: SC2 air and ww
SC2_corr <- ggplot(data=airport_data, aes(x=Air.PCR.Value, y=log10(WW.PCR.Value))) +
  geom_point(size=2, shape=1, stroke=1.5) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  scale_x_reverse(name="Airport Air PCR (Ct)", limits=c(35.2, 29.8))+
  scale_y_continuous(name="Aviation WW PCR log(cps/mL)")+
  ggtitle("SARS-CoV-2")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        plot.title = element_text(size=17, hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
SC2_corr

#Fig 2D: SC2 air and clinical
SC2_corr2 <- ggplot(data=airport_data, aes(x=Air.PCR.Value, y=Clinical.Percent.Pos)) +
  geom_point(size=2, shape=1, stroke=1.5) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  scale_x_reverse(name="Airport Air PCR (Ct)", limits=c(35.2, 29.8))+
  scale_y_continuous(name="Clinical % Positivity (NREVSS)")+
  ggtitle("SARS-CoV-2")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        plot.title = element_text(size=17, hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
SC2_corr2

#Fig 2C: SC2 air and nasal
SC2_corr3 <- ggplot(data=airport_data, aes(x=Air.PCR.Value, y=Nasal.Percent.Pos)) +
  geom_point(size=2, shape=1, stroke=1.5) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  scale_x_reverse(name="Airport Air PCR (Ct)", limits=c(35.2, 29.8))+
  scale_y_continuous(name="Traveler Nasal % Positivity")+
  ggtitle("SARS-CoV-2")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        plot.title = element_text(size=17, hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
SC2_corr3

#plot grid for SC2
SC2_correlations <- plot_grid(SC2_corr, SC2_corr3, SC2_corr2, nrow=1, labels = c("B", "C", "D"),
                              label_size=20)

######################################################################
#Create and export final Figure 2 plot

Fig2<- plot_grid(flu_SC2_plot,  SC2_correlations, flu_correlations,
                 nrow=3, labels = c("A", "", ""),
                              label_size=20)
Fig2

#export plot
save_plot(Fig2, filename="Fig2_Flu_SC2_positivity.png", dpi=300,
          base_height=15, base_width=12)
