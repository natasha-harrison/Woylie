############################################################


#     Woylie stations models Perup and Dryandra


############################################################

library(dplyr)
library(ggplot2)
library(lubridate)
library(glmmTMB)
library(forcats)
library(wesanderson)
library(ggpubr)
library(lmerTest)
library(png)
library(ggpattern)
library(tidyr)
library(effsize)







######################################################
#   Look ar ratio of dine in to take away
######################################################
#woylies_behave_md <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\woylies_behave_md_20240109.csv")
woylies_behave_md <- read.csv("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\woylies_behave_md_20240109.csv")


peanuts_eaten <- woylies_behave_md %>%
  filter(Behavior %in% c("Eat (dine in)", "Eat (take away)"),
         !is.na(Night),
         Station_Ammend != "") %>%
  group_by(Station_Ammend, Night) %>%
  summarise(Site = first(Site),
            Treatment = first(Treatment),
            Dine_in = sum(Behavior == "Eat (dine in)"),
            Take_away = sum(Behavior == "Eat (take away)"),
            Total = n())

#write.csv(peanuts_eaten, "C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\peanuts_eaten_20240109.csv")

peanuts_eaten_7 <- peanuts_eaten %>%
  filter(Total <8)

peanuts_to_plot <- gather(peanuts_eaten_7, key = DI_TA, value = "Peanuts_eaten", -Site, -Station_Ammend, -Night, -Treatment)

peanuts_to_plot <- peanuts_to_plot %>%
  mutate(Region = ifelse(Site %in% c("DNWS", "Dryandra_main"), "Dryandra", "Upper Warren"),
         Haven = ifelse(Site %in% c("DNWS", "Perup_Sanctuary"), "Havened", "Non-havened"))

#DI_TA <- ggplot(data = filter(peanuts_to_plot, DI_TA != "Total")) + 
#  geom_col(aes(x = Treatment, y = Peanuts_eaten, fill = DI_TA, pattern=DI_TA, colour=Haven), position = "fill") +
#  facet_grid(~ Region + Haven, switch = "x") +
#  ylab("Peanuts consumed") +
#  xlab("Treatment") +
#  scale_colour_manual(breaks = c("Havened", "Non-havened"),
#                      values= c("darkorange3", "springgreen4")) +
#  theme(        panel.border = element_blank(),
#                panel.grid.major = element_blank(),
#                panel.grid.minor = element_blank(),
#                axis.line = element_line(colour = "black"),
#                panel.background = element_rect("white"),
#                # axis.title.x = element_text(vjust = -3.5),
#                axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
#                text = element_text(size=12),
#                legend.position = "none")


#this takes ages to plot
DI_TA <-  ggplot(data = filter(peanuts_to_plot, DI_TA != "Total", !is.na(Treatment))) + 
  geom_col_pattern(aes(x = Treatment, y = Peanuts_eaten, fill = Haven, pattern = DI_TA), alpha=0.6, position = "fill", pattern_density = 0.1) +
  facet_grid(~ Region + Haven, switch = "x") +
  ylab("Mode of peanut consumption") +
  xlab("Treatment") +
  scale_fill_manual(values = c("Havened" = "darkorange3", "Non-havened" = "springgreen4")) +
  scale_pattern_manual(values = c("Dine_in" = "stripe", "Take_away" = "none")) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_rect("white"),
    axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
    text = element_text(size = 12),
    legend.position = "none")


##model proportion of dine in? - only model one as they are complimentary, use binomial for proportions
peanuts_eaten_model <- peanuts_eaten %>%
  mutate(Prop_DI = Dine_in/Total,
         Region = ifelse(Site %in% c("DNWS", "Dryandra_main"), "Dryandra", "UpperWarren"),
         Haven = ifelse(Site %in% c("DNWS", "Perup_Sanctuary"), "Havened", "Non-havened"),
         Woylie_density = ifelse(Site == "Dryandra_main", 1.39,
                                 ifelse(Site == "DNWS", 1.185,
                                        ifelse(Site == "Boyicup", 0.494,
                                               ifelse(Site == "Moopinup", 0.346,
                                                      ifelse(Site== "Perup_Sanctuary", 0.830, NA))))))


peanut_dita_1 <- glmmTMB(Prop_DI ~ Night + Woylie_density + Region*Haven + Treatment  + (1|Station_Ammend), 
                    family="binomial", data=peanuts_eaten_model)

summary(peanut_dita_1)

#check model fit
simulationOutput <- simulateResiduals(fittedModel = peanut_dita_1)
testResiduals(simulationOutput) #good

##liklihood ratio
peanut_dita_int <- glmmTMB(Prop_DI ~ Night + Woylie_density + Region + Haven + Treatment  + (1|Station_Ammend), 
                         family="binomial", data=peanuts_eaten_model)

peanut_dita_region <- glmmTMB(Prop_DI ~ Night + Woylie_density + Haven + Treatment  + (1|Station_Ammend), 
                         family="binomial", data=peanuts_eaten_model)

peanut_dita_haven <- glmmTMB(Prop_DI ~ Night + Woylie_density + Region + Treatment  + (1|Station_Ammend), 
                         family="binomial", data=peanuts_eaten_model)

peanut_dita_night <- glmmTMB(Prop_DI ~  Woylie_density + Region*Haven + Treatment  + (1|Station_Ammend), 
                         family="binomial", data=peanuts_eaten_model)

peanut_dita_density <- glmmTMB(Prop_DI ~ Night + Region*Haven + Treatment  + (1|Station_Ammend), 
                         family="binomial", data=filter(peanuts_eaten_model, !is.na(Woylie_density), !is.na(Treatment)))

peanut_dita_treatment <- glmmTMB(Prop_DI ~ Night + Woylie_density + Region*Haven  + (1|Station_Ammend), 
                         family="binomial", data=filter(peanuts_eaten_model, !is.na(Woylie_density), !is.na(Treatment)))


peanut_dita_LL <- glmmTMB(Prop_DI ~ Night + Woylie_density + Region*Haven + Treatment  + (1|Station_Ammend), 
                         family="binomial", data=filter(peanuts_eaten_model, !is.na(Woylie_density), !is.na(Treatment)))

anova(peanut_dita_1, peanut_dita_int)
anova(peanut_dita_region, peanut_dita_int)
anova(peanut_dita_haven, peanut_dita_int)
anova(peanut_dita_LL, peanut_dita_density)
anova(peanut_dita_1, peanut_dita_night)
anova(peanut_dita_LL, peanut_dita_treatment)




######################################################
#   model behaviour proportions at pop level
######################################################

##maybe change the data input to not exclude competition behaviours?

all_woylies_behave_ID <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\woylies_behave_ID_20240109.csv")

#condense to only first hour after frst visit per station per night
all_woylies_behave_ID$Date.modified <- as.POSIXct(all_woylies_behave_ID$Date.modified)

all_woylies_behave_first_hour <- all_woylies_behave_ID %>%
  group_by(Station_Ammend, Night, Site) %>%
  filter(Date.modified <= first(Date.modified) + 3600)



##group by category within 'detection'
category_by_detection <- all_woylies_behave_first_hour %>%
  dplyr::group_by(group_id, Station_Ammend, Behavioral.category) %>%
  dplyr::summarise(behave_sum = sum(Duration..s.),
                   group_sum = first(group_sum),
                   Station = first(Station_Ammend),
                   Site = first(Site),
                   Observer = first(Observer),
                   RFID_DateTime = first(RFID_DateTime),
                   Treatment = first(Treatment),
                   Night = first(Night))

##create a new data frame with all behaviour possibilities and join in what we have (assume the rest is 0)
category_by_detection_All_Possibilities <- expand.grid(Station =unique(category_by_detection$Station),
                                                       group_id =unique(category_by_detection$group_id),
                                                       Behavioral.category =unique(category_by_detection$Behavioral.category))

#first inner join with group info
category_by_detection_groups <- category_by_detection %>%
  dplyr::group_by(Station, group_id) %>%
  dplyr::summarise(group_sum = first(group_sum),
                   Station = first(Station),
                   Site = first(Site),
                   Observer = first(Observer),
                   RFID_DateTime = first(RFID_DateTime),
                   Treatment = first(Treatment),
                   Night = first(Night))


category_by_detection_All_Possibilities <- inner_join(category_by_detection_All_Possibilities, category_by_detection_groups,
                                                      by=c("Station" = "Station",
                                                           "group_id" = "group_id"))

##then join on the actual values and make the rest 0s
category_by_detection_behave <- category_by_detection %>%
  select(Station, group_id, Behavioral.category, behave_sum)


category_by_detection_All_Possibilities <- left_join(category_by_detection_All_Possibilities, category_by_detection_behave,
                                                     by= c("Station" = "Station",
                                                           "group_id" = "group_id",
                                                           "Behavioral.category" = "Behavioral.category"))

category_by_detection_All_Possibilities$behave_sum <- ifelse(!is.na(category_by_detection_All_Possibilities$behave_sum), category_by_detection_All_Possibilities$behave_sum, 0)

category_by_detection_All_Possibilities <- category_by_detection_All_Possibilities %>%
  mutate(behave_prop = behave_sum/group_sum,
         Haven = ifelse(Site %in% c("DNWS", "Perup_Sanctuary"), "Havened", "Non-havened"),
         Region = ifelse(Site %in% c("DNWS", "Dryandra_main"), "Dryandra", "Upper_Warren"),
         Woylie_density = ifelse(Site == "Dryandra_main", 1.39,
                                        ifelse(Site == "DNWS", 1.185,
                                               ifelse(Site == "Boyicup", 0.494,
                                                      ifelse(Site == "Moopinup", 0.346,
                                                             ifelse(Site== "Perup_Sanctuary", 0.830, NA))))))







write.csv(category_by_detection_All_Possibilities, "C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\category_by_detection_All_Possibilities_20240109.csv")



##LETS MODEL
str(category_by_detection_All_Possibilities)


category_by_detection_All_Possibilities %>%
  summarise(sum(behave_sum))


##########################################

###   EXPLORATION

##########################################

cat_explore <- category_by_detection_All_Possibilities %>%
  filter(Behavioral.category == "Exploration", !is.na(Station))

explore1 <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_explore)
summary(explore1)

#check model fit
simulationOutput <- simulateResiduals(fittedModel = explore1)
testResiduals(simulationOutput) #good


      #liklihood ratio
explore_int <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Haven + Region + (1|Observer) + (1|Station), family='binomial', data= cat_explore)
explore_region <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Haven + (1|Observer) + (1|Station), family='binomial', data= cat_explore)
explore_haven <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Region + (1|Observer) + (1|Station), family='binomial', data= cat_explore)
explore_night <- glmmTMB(behave_prop ~  Treatment + Woylie_density + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_explore)
explore_density <- glmmTMB(behave_prop ~ Night + Treatment  + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_explore)
explore_treatment <- glmmTMB(behave_prop ~ Night + Woylie_density + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= filter(cat_explore, !is.na(Treatment)))

anova(explore1, explore_int)
anova(explore_region, explore_int)
anova(explore_haven, explore_int)
anova(explore1, explore_night)
anova(explore1, explore_density)
anova(explore1, explore_treatment)


dry_e_h <- cat_explore %>%
  filter(Region == "Dryandra", Haven == "Havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, dry_e_h)

dry_e_nh <- cat_explore %>%
  filter(Region == "Dryandra", Haven == "Non-havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, dry_e_nh)


per_e_h <- cat_explore %>%
  filter(Region == "Upper_Warren", Haven == "Havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, per_e_h)

per_e_nh <- cat_explore %>%
  filter(Region == "Upper_Warren", Haven == "Non-havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, per_e_nh)



##########################################

###       FORAGING

##########################################

cat_forage <- category_by_detection_All_Possibilities %>%
  filter(Behavioral.category == "Foraging")

forage1 <- glmmTMB(behave_prop ~ scale(Night) + Treatment + scale(Woylie_density) + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= filter(cat_forage, !is.na(Treatment)))
summary(forage1)

#check model fit
simulationOutput <- simulateResiduals(fittedModel = forage1)
testResiduals(simulationOutput) #good


#liklihood ratio
f_int <- glmmTMB(behave_prop ~ scale(Night) + Treatment + scale(Woylie_density) + Haven + Region + (1|Observer) + (1|Station), family='binomial', data= cat_forage)
f_region <- glmmTMB(behave_prop ~ scale(Night) + Treatment + scale(Woylie_density) + Haven + (1|Observer) + (1|Station), family='binomial', data= cat_forage)
f_haven <- glmmTMB(behave_prop ~ scale(Night) + Treatment + scale(Woylie_density) + Region + (1|Observer) + (1|Station), family='binomial', data= cat_forage)
f_night <- glmmTMB(behave_prop ~  Treatment + scale(Woylie_density) + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_forage)
f_density <- glmmTMB(behave_prop ~ scale(Night) + Treatment  + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_forage)
f_treatment <- glmmTMB(behave_prop ~ scale(Night) + scale(Woylie_density) + Haven + Region + (1|Observer) + (1|Station), family='binomial', data= filter(cat_forage, !is.na(Treatment)))

anova(forage1, f_int)
anova(f_region, f_int)
anova(f_haven, f_int)
anova(forage1, f_night)
anova(forage1, f_density)
anova(forage1, f_treatment)



dry_f_h <- cat_forage %>%
  filter(Region == "Dryandra", Haven == "Havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, dry_f_h)

dry_f_nh <- cat_forage %>%
  filter(Region == "Dryandra", Haven == "Non-havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, dry_f_nh)


per_f_h <- cat_forage %>%
  filter(Region == "Upper_Warren", Haven == "Havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, per_f_h)

per_f_nh <- cat_forage %>%
  filter(Region == "Upper_Warren", Haven == "Non-havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, per_f_nh)

##########################################

###       VIGILANCE

##########################################

cat_vig <- category_by_detection_All_Possibilities %>%
  filter(Behavioral.category == "Vigilance")

vig1 <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= filter(cat_vig, !is.na(Treatment)))
summary(vig1)

#check model fit
simulationOutput <- simulateResiduals(fittedModel = vig1)
testResiduals(simulationOutput) #good

#liklihood ratio
v_int <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Haven + Region + (1|Observer) + (1|Station), family='binomial', data= cat_vig)
v_region <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Haven + (1|Observer) + (1|Station), family='binomial', data= cat_vig)
v_haven <- glmmTMB(behave_prop ~ Night + Treatment + Woylie_density + Region + (1|Observer) + (1|Station), family='binomial', data= cat_vig)
v_night <- glmmTMB(behave_prop ~  Treatment + Woylie_density + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_vig)
v_density <- glmmTMB(behave_prop ~ Night + Treatment  + Haven*Region + (1|Observer) + (1|Station), family='binomial', data= cat_vig)
v_treatment <- glmmTMB(behave_prop ~ Night + Woylie_density + Haven + Region + (1|Observer) + (1|Station), family='binomial', data= filter(cat_vig, !is.na(Treatment)))

anova(vig1, v_int)
anova(v_region, v_int)
anova(v_haven, v_int)
anova(vig1, v_night)
anova(vig1, v_density)
anova(vig1, v_treatment)



dry_v_h <- cat_vig %>%
  filter(Region == "Dryandra", Haven == "Havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, dry_v_h)

dry_v_nh <- cat_vig %>%
  filter(Region == "Dryandra", Haven == "Non-havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, dry_v_nh)


per_v_h <- cat_vig %>%
  filter(Region == "Upper_Warren", Haven == "Havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, per_v_h)

per_v_nh <- cat_vig %>%
  filter(Region == "Upper_Warren", Haven == "Non-havened") %>%
  select(behave_prop, Treatment)

cohen.d(behave_prop ~ Treatment, per_v_nh)
