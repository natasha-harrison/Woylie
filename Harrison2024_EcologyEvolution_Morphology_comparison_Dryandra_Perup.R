
##############################################################################################

# Compare morphology of havened/non-havened individuals at Dryandra and Perup

##############################################################################################

library(dplyr)
library(lme4)
library(ciTools)
library(DHARMa)
library(effsize)
library(performance)
library(lmerTest)
library(ggplot2)
library(ggpubr)


#############################
#       data prep
#############################

#   UPPER WARREN
ps <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Trapping_data\\PerupSanctuary_2022.csv")# MORPHOLOGY
boy <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Trapping_data\\Boyicup_2022.csv")# MORPHOLOGY
moo <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Trapping_data\\Moopinup_2022.csv")# MORPHOLOGY

#join them
warren <- bind_rows(ps, boy, moo)

warren_subset <- warren %>%
  filter(AGE == "A") %>%
  dplyr::group_by(ANIMAL_NO) %>%
  dplyr::summarise(Sex = first(SEX),
            Weight = first(WEIGHT),
            Pes = first(PES_LONG),
            Head = first(HEAD_LENGTH),
            Population = first(SSI_LABEL)) %>%
  mutate(Haven = ifelse(Population %in% c("51BOY/01", "51POS/01"), "Non-havened", "Havened"))

warren_subset$Population <- gsub("51PSF11O20", "PerupSanctuary" , warren_subset$Population)
warren_subset$Population <- gsub("51BOY/01", "Boyicup" , warren_subset$Population)
warren_subset$Population <- gsub("51POS/01", "Moopinup" , warren_subset$Population)


#   DRYANDRA
dryandra <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Dryandra Experimental Release\\Translocated Individuals.csv")

dryandra_subset <- dryandra %>%
  filter(Source %in% c("DNWS", "Dryandra_main")) %>%
  dplyr::select(Individual_ID, Sex, Weight, Pes, Head, Source) %>%
  mutate(Haven = ifelse(Source == "DNWS", "Havened", "Non-havened"))


#join them together
names(warren_subset) <- c("Individual_ID", "Sex", "Weight" ,"Pes", "Head", "Population", "Haven")
names(dryandra_subset) <- c("Individual_ID", "Sex", "Weight" ,"Pes", "Head", "Population", "Haven")

all_woylies <- rbind(warren_subset, dryandra_subset)


all_woylies <- all_woylies %>%
  filter(Sex != "U", !is.na(Weight)) %>%
  mutate(Region = ifelse(Population %in% c("Dryandra_main", "DNWS"), "Dryandra", "UpperWarren"))

#write.csv(all_woylies, "C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\all_woylie.csv")
#all_woylies <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\all_woylie.csv")
all_woylies <- read.csv("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\all_woylie.csv")

all_woylies <- all_woylies %>%
  mutate(Woylie_density = ifelse(Population == "Dryandra_main", 1.39,
                        ifelse(Population == "DNWS", 1.185,
                               ifelse(Population == "Boyicup", 0.494,
                                      ifelse(Population == "Moopinup", 0.346,
                                             ifelse(Population== "PerupSanctuary", 0.830, NA))))))


###mean weight and pes 
all_woylies %>%
  group_by(Region, Haven) %>%
  summarise(mean_weight = mean(Weight, na.rm=T),
            mean_pes = mean(Pes, na.rm=T))



##############################

#     model weight

##############################
all_woylies$Population <- as.factor(all_woylies$Population)
all_woylies$Haven <- as.factor(all_woylies$Haven)
all_woylies$Region <- as.factor(all_woylies$Region)

#weight_all <- lmer(log(Weight) ~ Sex + Haven*Region + (1|Region/Population), data=all_woylies)

weight_all <- lm(log(Weight) ~ Sex + Haven*Region, data=all_woylies)
summary(weight_all)



##post hoc tukey pairwise test
library(emmeans)

tuk <- emmeans(weight_all, ~ Haven*Region)

#fit2.emm.a
pairs(tuk, adjust="tukey")


#variance inflation factor to check for collinearity
check_collinearity(weight_all)

#check model fit
simulationOutput <- simulateResiduals(fittedModel = weight_all)
testResiduals(simulationOutput) #beautiful


#loglikelihood tests
weight_all_sex <- lm(log(Weight) ~ Haven*Region , data=all_woylies)
weight_all_haven <- lm(log(Weight) ~ Sex  + Region, data=all_woylies)
weight_all_region <- lm(log(Weight) ~ Sex  + Haven, data=all_woylies)
weight_all_int <- lm(log(Weight) ~ Sex  + Haven + Region, data=all_woylies)
#weight_all_density <- lm(log(Weight) ~ Sex  + Haven*Region, data=all_woylies)

#anova(weight_all_density, weight_all)
anova(weight_all_sex, weight_all)
anova(weight_all_haven, weight_all)
anova(weight_all_region, weight_all)
anova(weight_all_int, weight_all)


#get cohen's d 

dry_weight <- all_woylies %>%
  filter(Region == "Dryandra") %>%
  select(Weight, Haven)

cohen.d(Weight ~ Haven, dry_weight)


per_weight <- all_woylies %>%
  filter(Region == "UpperWarren") %>%
  select(Weight, Haven)

cohen.d(Weight ~ Haven, per_weight)





#plot model predictions
new_data_weight <- expand.grid(Haven= c("Havened", "Non-havened"),
                               Region= c("UpperWarren", "Dryandra"),
                               Sex = "M",
                               Population = c("Boyicup", "Moopinup", "PerupSanctuary", "DNWS", "Dryandra_main"))

#this ciTools function gives predictions and confidence intervals
predframe_weight <- add_ci(new_data_weight, weight_all, alpha = 0.05, names = c("lwr", "upr"), includeRanef=TRUE)

#back transform
predframe_weight <- predframe_weight %>%
  mutate(pred_exp = exp(pred),
         lwr_exp = exp(lwr),
         upr_exp = exp(upr))

#get the ones that I want
predframe_weight_subset <- predframe_weight %>% 
  mutate(ID = seq.int(nrow(predframe_weight))) %>%
  filter(ID %in% c(2, 6, 9, 15, 20))

predframe_weight_subset$Region <- relevel(predframe_weight_subset$Region, ref = "Dryandra")

#plot predictions
weight_region <- ggplot() +
  geom_errorbar(data=predframe_weight_subset, mapping=aes(x=Region, ymin=lwr_exp, ymax = upr_exp, colour=Haven, group=Population), position=position_dodge(width=0.5), linewidth=1, width=0.4, alpha=0.5) +
  geom_point(data=predframe_weight_subset, mapping=aes(x=Region, y=pred_exp, colour=Haven, group=Population), position=position_dodge(width=0.5), size=3) +
  scale_colour_manual(values= c("darkorange3", "springgreen4")) +
  xlab("Region") +
  ylab("Weight (g)") + 
  #ylim(1000,1350) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
                text = element_text(size=12),
                legend.position= "none")




##############################

#     model pes

##############################

#get more data from PS
#capture <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\DBCA_woylie\\trapping_woylie.csv")
capture <- read.csv("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\Data\\DBCA_woylie\\trapping_woylie.csv")


capture <- capture %>%
  mutate(haven = ifelse(SSI_LABEL %in% 
                          c("51PSD23F26", "51PSK17M19", "51PSK6M8", "51PWS/01"), "Haven", "Outside"))

capture$TRP_DATE <- as.Date(capture$TRP_DATE, format= "%Y-%m-%d")

capture_pes <- capture %>%
  filter(haven == "Haven", !is.na(PES_LONG), !is.na(HEAD_LENGTH), TRP_DATE > "2017-01-01", AGE=="A") %>%
  group_by(ANIMAL_NO, year(TRP_DATE)) %>%
  summarise(Sex = first(SEX), 
            Weight = first(WEIGHT), 
            Pes = first(PES_LONG),
            Head = first(HEAD_LENGTH))


capture_pes <- capture_pes %>%
  mutate(Population="PerupSanctuary",
         Haven="Havened",
         Region = "UpperWarren") %>%
  dplyr::select(ANIMAL_NO, Sex, Weight, Pes, Head, Population, Haven, Region)

##all_woylies <- all_woylies[,2:9] ##have to run this is read in all_woylies from csv

names(capture_pes) <- names(all_woylies)

all_woylies_pes <- rbind(all_woylies, capture_pes)

all_woylies_pes <- all_woylies_pes %>%
  filter(Sex %in% c("M", "F")) %>%
  mutate(log_head = log(Head),
         Woylie_density = ifelse(Population == "Dryandra_main", 1.39,
                                 ifelse(Population == "DNWS", 1.185,
                                        ifelse(Population == "Boyicup", 0.494,
                                               ifelse(Population == "Moopinup", 0.346,
                                                      ifelse(Population== "PerupSanctuary", 0.830, NA))))))

all_woylies_pes$Sex <- as.factor(all_woylies_pes$Sex)
all_woylies_pes$Region <- as.factor(all_woylies_pes$Region)
all_woylies_pes$Population <- as.factor(all_woylies_pes$Population)
all_woylies_pes$log_head <- as.numeric(as.character(all_woylies_pes$log_head))


pes_all <- lm(log(Pes) ~ Sex + log_head + Haven*Region, data=all_woylies_pes)
summary(pes_all)


##post hoc tukey pairwise test
tuk <- emmeans(pes_all, ~ Haven*Region)

#fit2.emm.a
pairs(tuk, adjust="tukey")




#check model fit
simulationOutput <- simulateResiduals(fittedModel = pes_all)
testResiduals(simulationOutput) #good


#logliklihoods
data_pes <- all_woylies_pes %>%
  filter(!is.na(log_head))

pes_all_sex <- lm(log(Pes) ~ log_head + Haven*Region, data=all_woylies_pes)
pes_all_haven <- lm(log(Pes) ~ Sex + log_head + Region, data=all_woylies_pes)
pes_all_region <- lm(log(Pes) ~ Sex + log_head + Haven, data=all_woylies_pes)
pes_all_loghead <- lm(log(Pes) ~ Sex + Haven + Region, data=data_pes)
pes_all_int <- lm(log(Pes) ~ Sex + log_head + Haven + Region, data=all_woylies_pes)

anova(pes_all_sex, pes_all)
anova(pes_all_haven, pes_all)
anova(pes_all_region, pes_all)
anova(pes_all_loghead, pes_all) #have to run original model on smaller dataset
anova(pes_all_int, pes_all)





#get cohen's d 

dry_pes <- all_woylies_pes %>%
  filter(Region == "Dryandra") %>%
  select(Pes, Haven)

cohen.d(Pes ~ Haven, dry_pes)


per_pes <- all_woylies_pes %>%
  filter(Region == "UpperWarren") %>%
  select(Pes, Haven)

cohen.d(Pes ~ Haven, per_pes)





#plot model predictions
new_data_pes <- expand.grid(Haven= unique(all_woylies_pes$Haven),
                               Region= unique(all_woylies_pes$Region),
                               Sex = "F",
                               Population = unique(all_woylies_pes$Population),
                              log_head = mean(all_woylies_pes$log_head, na.rm=T))

#this ciTools function gives predictions and confidence intervals
predframe_pes <- add_ci(new_data_pes, pes_all, alpha = 0.05, names = c("lwr", "upr"), includeRanef=F)

#back transform
predframe_pes <- predframe_pes %>%
  mutate(pred_exp = exp(pred),
         lwr_exp = exp(lwr),
         upr_exp = exp(upr))

#get the ones that I want
predframe_pes_subset <- predframe_pes %>% 
  mutate(ID = seq.int(nrow(predframe_pes))) %>%
  filter(ID %in% c(1, 6, 16, 19))


#factor relevel to look nice
predframe_pes_subset$Region <- fct_relevel(predframe_pes_subset$Region, "UpperWarren", "Dryandra")
predframe_pes_subset$Haven <- fct_relevel(predframe_pes_subset$Haven, "Havened", "Non-havened")

#plot predictions for ICTC 
pes_region <- ggplot() +
  geom_errorbar(data=predframe_pes_subset, mapping=aes(x=Region, ymin=lwr_exp, ymax = upr_exp, colour=Haven, group=Haven), position=position_dodge(width=0.5), linewidth=1, width=0.4, alpha=0.5) +
  geom_point(data=predframe_pes_subset, mapping=aes(x=Region, y=pred_exp, colour=Haven, group=Haven), position=position_dodge(width=0.5), size=3) +
  scale_colour_manual(values= c("darkorange3", "springgreen4")) +
  xlab("Region") +
  ylab("Pes (mm)") + 
  #ylim(1000,1350) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
                text = element_text(size=12),
                legend.position= "none")



weight_pes_region <- ggarrange(weight_region, pes_region, nrow =1, ncol = 2, labels =c("a", "b"))

ggsave("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\R_outputs_woylie_comparison\\weight_pes_region.pdf",
       plot = weight_pes_region,
       scale = 1,
       width = 18,
       height = 7,
       units = "cm")



##############################

#     Rain plots!

##############################
#https://labs.ala.org.au/posts/2023-08-28_alternatives-to-box-plots/post.html

library(ggdist)
library(gghalves)
library(pilot)

all_woylies_weight_plot <- all_woylies %>%
  filter(!is.na(Weight)) %>%
  mutate(Pop_merged = paste(Region,Haven, sep="_"))

##WEIGHT

weight_rain <- ggplot(data = all_woylies_weight_plot, 
       aes(x = Pop_merged %>% stringr::str_wrap(10) %>% reorder(Weight), 
           y = Weight, 
           colour = Haven, 
           fill = Haven)) +
  scale_colour_manual(values= c("darkorange3", "springgreen4")) +
  scale_fill_manual(values= c("darkorange3", "springgreen4")) +
  ggdist::stat_halfeye(adjust = .4, # smoothness of distribution
                       width = .87, # height of distribution
                       colour = NA,
                       alpha=0.6) +
  gghalves::geom_half_point(side = "l", # choose right or left side
                            range_scale = .3, # spread of points
                            alpha = .6,
                            size = 2.2) +
  #scale_colour_manual(values= c("springgreen4", "darkorange3")) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(colour, .7))),
    width = .12,        # adjust box width
    fill = NA,
    size = 1.1,         # size of box line
    outlier.shape = NA  # remove outlier points
  ) +
  scale_x_discrete(labels=c("Havened, \nDryandra",   "Non-havened, \nDryandra", "Havened, \nUpper Warren", "Non-havened, \nUpper Warren")) +
  coord_flip() +
  labs(x = "Population",
       y = "Weight (g)") +
  #scale_y_continuous(breaks = c(0, 100, 200, 300, 400),
  #                   labels = c(0, 100, 200, 300, 400),
  #                   limits = c(0, 400),
  #                   expand = c(0,0)) +
  #pilot::scale_color_pilot() +
  #pilot::scale_fill_pilot() +
  #pilot::theme_pilot(grid = "",
  #                   axes = "b") + 
  theme(legend.position = "none",
        axis.title.x = ggtext::element_markdown(),
        axis.text.y = element_text(face = "italic"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect("white"))
     


## PES

all_woylies_pes_plot <- all_woylies_pes %>% filter(!is.na(Pes)) %>%
  mutate(Pop_merged = paste(Region,Haven, sep="_"))

pes_rain <- ggplot(data = all_woylies_pes_plot, 
       aes(x = Pop_merged %>% stringr::str_wrap(10) %>% reorder(Pes), 
           y = Pes, 
           colour = Haven, 
           fill = Haven)) +
  scale_colour_manual(values= c("darkorange3", "springgreen4")) +
  scale_fill_manual(values= c("darkorange3", "springgreen4")) +
  ggdist::stat_halfeye(adjust = .4, # smoothness of distribution
                       width = .87, # height of distribution
                       colour = NA,
                       alpha=0.6) +
  gghalves::geom_half_point(side = "l", # choose right or left side
                            range_scale = .3, # spread of points
                            alpha = .6,
                            size = 2.2) +
  #scale_colour_manual(values= c("springgreen4", "darkorange3")) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(colour, .7))),
    width = .12,        # adjust box width
    fill = NA,
    size = 1.1,         # size of box line
    outlier.shape = NA  # remove outlier points
  ) +
  scale_x_discrete(labels=c("Havened, \nDryandra",   "Non-havened, \nDryandra", "Havened, \nUpper Warren", "Non-havened, \nUpper Warren")) +
  coord_flip() +
  labs(x = "Population",
       y = "Pes (mm)") +
  #scale_y_continuous(breaks = c(0, 100, 200, 300, 400),
  #                   labels = c(0, 100, 200, 300, 400),
  #                   limits = c(0, 400),
  #                   expand = c(0,0)) +
  #pilot::scale_color_pilot() +
  #pilot::scale_fill_pilot() +
  #pilot::theme_pilot(grid = "",
  #                   axes = "b") + 
  theme(legend.position = "none",
        axis.title.x = ggtext::element_markdown(),
        axis.text.y = element_text(face = "italic"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect("white"))



#combine them
rain_weight_pes <- ggarrange(weight_rain +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(r = 1) ), 
          pes_rain + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.margin = margin(l = 1) ), 
          ncol=2,nrow=1, widths = c(1, 0.80), labels = c("(a)", "(b)"))


ggsave("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\R_outputs_woylie_comparison\\weight_pes_rainplot.pdf",
       plot = rain_weight_pes,
       scale = 1,
       width = 22,
       height = 10,
       units = "cm")
