
##########################################################################

##      Compare peanuts remaining in Fox vs Control stations

##########################################################################

library(dplyr)
library(ggplot2)
library(glmmTMB)

#Dryandra <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Dryandra Experimental Release\\DryandraTranslocation_RFID_Stations_backUp.csv")
#Perup <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\RFID_Stations\\RFID_deployment.csv")

Dryandra <- read.csv("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\Data\\Dryandra Experimental Release\\DryandraTranslocation_RFID_Stations_backUp.csv")
Perup <- read.csv("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\Data\\RFID_Stations\\RFID_deployment.csv")



GUD_all <- rbind(Dryandra, Perup[,1:16])

GUD_all <- GUD_all %>%
  mutate(Region = ifelse(Site %in% c("DNWS", "Dryandra_main"), "Dryandra", "UpperWarren"),
         Haven = ifelse(Site %in% c("DNWS", "Perup_Sanctuary"), "Havened", "Non-havened"),
         Woylie_density = ifelse(Site == "Dryandra_main", 1.39,
                                 ifelse(Site == "DNWS", 1.185,
                                        ifelse(Site == "Boyicup", 0.494,
                                               ifelse(Site == "Moopinup", 0.346,
                                                      ifelse(Site== "Perup_Sanctuary", 0.830, NA))))))
  

###################################################

#               GUD peanut models                 #

###################################################

GUD_all$Peanuts_left <- as.integer(GUD_all$Peanuts_left)
GUD_all$Night <- as.numeric(GUD_all$Night)
GUD_all$Site <- as.factor(GUD_all$Site)
GUD_all$Treatment <- as.factor(GUD_all$Treatment)
GUD_all$Trap <- as.factor(GUD_all$Trap)
GUD_all$Region <- as.factor(GUD_all$Region)
GUD_all$Woylie_density <- as.numeric(GUD_all$Woylie_density)

peanut_1 <- glmmTMB(Peanuts_left ~ Night + Woylie_density + Region*Haven + Treatment  + (1|Trap), 
                    family="poisson", data=GUD_all)
summary(peanut_1)

#loglikilihoods
peanut_Night <- glmmTMB(Peanuts_left ~  Region*Haven + Treatment + Woylie_density + (1|Trap), 
        family="poisson", data=GUD_all)
peanut_region <- glmmTMB(Peanuts_left ~ Night + Haven + Treatment + Woylie_density + (1|Trap), 
                        family="poisson", data=GUD_all)
peanut_haven <- glmmTMB(Peanuts_left ~ Night + Region + Treatment + Woylie_density + (1|Trap), 
                        family="poisson", data=GUD_all)
peanut_treatment <- glmmTMB(Peanuts_left ~ Night + Region*Haven + Woylie_density + (1|Trap), 
                        family="poisson", data=GUD_all)
peanut_int <- glmmTMB(Peanuts_left ~ Night + Region + Haven + Treatment + Woylie_density + (1|Trap), 
                        family="poisson", data=GUD_all)
peanut_density <- glmmTMB(Peanuts_left ~ Night + Region*Haven + Treatment  + (1|Trap), 
                    family="poisson", data=GUD_all)



anova(peanut_1, peanut_Night)
anova(peanut_1, peanut_density)
anova(peanut_1, peanut_region)
anova(peanut_1, peanut_haven)
anova(peanut_1, peanut_treatment)
anova(peanut_1, peanut_int)






##model predictions
new_dataGUD <- expand.grid(Night = 1,
                           Site = unique(GUD_all$Site),
                           Region = unique(GUD_all$Region),
                           Woylie_density = mean(GUD_all$Woylie_density),
                           Treatment = c("F", "C"),
                           Haven = c("Havened", "Non-havened"),
                           Trap = "NA")


new_dataGUD2 <- data.frame(Site = c("DNWS", "DNWS", "Dryandra_main", "Dryandra_main", "Boyicup", "Boyicup", "Perup_Sanctuary", "Perup_Sanctuary"),
                           Night=c(1,1,1,1,1,1,1,1),
                           Region=c("Dryandra", "Dryandra","Dryandra","Dryandra", "UpperWarren", "UpperWarren", "UpperWarren", "UpperWarren"),
                           Woylie_density=c(1.185, 1.185, 1.39, 1.39, 0.494, 0.494, 0.83, 0.83),
                           Treatment=c("F", "C", "F", "C", "F", "C", "F", "C"),
                           Haven=c("Havened", "Havened", "Non-havened", "Non-havened", "Non-havened", "Non-havened","Havened", "Havened"),
                           Trap=c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"))


#predict values
predGUD <- predict(peanut_1, newdata = new_dataGUD2,
                   random.only=FALSE, 
                   se.fit =T, # <- this is from a different predict function not predict.merMOD
                   type = "response", 
                   allow.new.levels = TRUE, 
                   na.action = na.pass)

#combine dataframes
predframeGUD <- data.frame(new_dataGUD2, preds = predGUD$fit, se = predGUD$se.fit)

#get the ones that I want if I use newdataGUD
#predframeGUD_subset <- predframeGUD %>% 
#  mutate(ID = seq.int(nrow(predframeGUD))) %>%
#  filter(ID %in% c(1, 8, 11, 18, 22, 32, 29, 30, 39, 40))


##plot it
ggplot() +
  geom_point(data=predframeGUD_subset, mapping=aes(x=Region:Haven, y=preds, shape=Treatment, group=Haven:Treatment, colour=Haven), size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(data=predframeGUD_subset, mapping=aes(x=Region:Haven, ymin=preds-se, ymax=preds+se, group=Haven:Treatment, colour=Haven), linewidth=1, width=0.4, alpha=0.5, position=position_dodge(width=0.5)) +
  ylab("Peanuts left") +
  xlab("Region") +
  scale_colour_manual(breaks = c("Havened", "Non-havened"),
                      values= c("darkorange3", "springgreen4")) +
  scale_x_discrete(labels= c(expression(paste("Havened \nDryandra")), 
                                        expression(paste("Non-havened \nDryandra")), 
                                                   expression(paste("Havened \nUpper Warren")), 
                                                              expression(paste("Non-havened \nUpper Warren")))) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
               # axis.title.x = element_text(vjust = -3.5),
                axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
                text = element_text(size=12),
                legend.position = "none")

#can think about customising labels here:
#https://stackoverflow.com/questions/13223846/ggplot2-two-line-label-with-expression



##try new plot

ggplot() +
  geom_point(data=predframeGUD, mapping=aes(x=Treatment, y=preds, shape=Treatment, colour=Haven), size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(data=predframeGUD, mapping=aes(x=Treatment, ymin=preds-se, ymax=preds+se, colour=Haven), linewidth=1, width=0.4, alpha=0.5, position=position_dodge(width=0.5)) +
  facet_grid(~ Site, switch = "x") +
  ylab("Peanuts left") +
  xlab("Site") +
  scale_colour_manual(breaks = c("Havened", "Non-havened"),
                      values= c("darkorange3", "springgreen4")) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
                # axis.title.x = element_text(vjust = -3.5),
                axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
                text = element_text(size=12),
                legend.position = "none")






#get cohen's d 

dry_gud_fox <- GUD_all %>%
  filter(Region == "Dryandra", Treatment == "F") %>%
  select(Peanuts_left, Haven)

cohen.d(Peanuts_left ~ Haven, dry_gud_fox)

dry_gud_con <- GUD_all %>%
  filter(Region == "Dryandra", Treatment == "C") %>%
  select(Peanuts_left, Haven)

cohen.d(Peanuts_left ~ Haven, dry_gud_con)


per_gud_fox <- GUD_all %>%
  filter(Region == "UpperWarren", Treatment == "F") %>%
  select(Peanuts_left, Haven)

cohen.d(Peanuts_left ~ Haven, per_gud_fox)

per_gud_con <- GUD_all %>%
  filter(Region == "UpperWarren", Treatment == "C") %>%
  select(Peanuts_left, Haven)

cohen.d(Peanuts_left ~ Haven, per_gud_con)






#get cohen's d 

dry_gud_h <- GUD_all %>%
  filter(Region == "Dryandra", Haven == "Havened") %>%
  select(Peanuts_left, Treatment)

cohen.d(Peanuts_left ~ Treatment, dry_gud_h)

dry_gud_nh <- GUD_all %>%
  filter(Region == "Dryandra", Haven == "Non-havened") %>%
  select(Peanuts_left, Treatment)

cohen.d(Peanuts_left ~ Treatment, dry_gud_nh)


per_gud_h <- GUD_all %>%
  filter(Region == "UpperWarren", Haven == "Havened") %>%
  select(Peanuts_left, Treatment)

cohen.d(Peanuts_left ~ Treatment, per_gud_h)

per_gud_nh <- GUD_all %>%
  filter(Region == "UpperWarren", Haven == "Non-havened") %>%
  select(Peanuts_left, Treatment)

cohen.d(Peanuts_left ~ Treatment, per_gud_nh)


###############################################################

#               GUD peanut models as binomial                 #
#
###############################################################
GUD_all <- GUD_all %>%
  mutate(All_gone = ifelse(Peanuts_left== 0, 1, 0))

peanut_b1 <- glmmTMB(All_gone ~ Night + Region*Haven + Treatment + Woylie_density + (1|Trap), 
                    family="binomial", data=GUD_all)
summary(peanut_b1)


#loglikilihoods
peanut_Nightb <- glmmTMB(All_gone ~  Region*Haven + Treatment + Woylie_density + (1|Trap), 
                        family="binomial", data=GUD_all)
peanut_regionb <- glmmTMB(All_gone ~ Night + Haven + Treatment + Woylie_density + (1|Trap), 
                         family="binomial", data=GUD_all)
peanut_havenb <- glmmTMB(All_gone ~ Night + Region + Treatment + Woylie_density  + (1|Trap), 
                        family="binomial", data=GUD_all)
peanut_treatmentb <- glmmTMB(All_gone ~ Night + Region*Haven + Woylie_density + (1|Trap), 
                            family="binomial", data=GUD_all)
peanut_intb <- glmmTMB(All_gone ~ Night + Region + Haven + Treatment + Woylie_density + (1|Trap), 
                      family="binomial", data=GUD_all)
peanut_densityb <- glmmTMB(All_gone ~ Night + Region*Haven + Treatment + (1|Trap), 
                     family="binomial", data=GUD_all)



anova(peanut_b1, peanut_Nightb)
anova(peanut_b1, peanut_regionb)
anova(peanut_b1, peanut_havenb)
anova(peanut_b1, peanut_treatmentb)
anova(peanut_b1, peanut_intb)
anova(peanut_b1, peanut_densityb)




##model predictions
new_dataGUD_b <- expand.grid(Night = 1,
                           Site = unique(GUD_all$Site),
                           Region = unique(GUD_all$Region),
                           Treatment = c("F", "C"),
                           Woylie_density = mean(GUD_all$Woylie_density),
                           Haven = c("Havened", "Non-havened"),
                           Trap = "NA")

new_dataGUD2b <- data.frame(Site = c("DNWS", "DNWS", "Dryandra_main", "Dryandra_main", "Boyicup", "Boyicup", "Perup_Sanctuary", "Perup_Sanctuary"),
                           Night=c(1,1,1,1,1,1,1,1),
                           Region=c("Dryandra", "Dryandra","Dryandra","Dryandra", "UpperWarren", "UpperWarren", "UpperWarren", "UpperWarren"),
                           Woylie_density=c(1.185, 1.185, 1.39, 1.39, 0.494, 0.494, 0.83, 0.83),
                           Treatment=c("F", "C", "F", "C", "F", "C", "F", "C"),
                           Haven=c("Havened", "Havened", "Non-havened", "Non-havened", "Non-havened", "Non-havened","Havened", "Havened"),
                           Trap=c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"))



#predict values
predGUD_b <- predict(peanut_b1, newdata = new_dataGUD_b,
                   #random.only=FALSE, 
                   se.fit =T, # <- this is from a different predict function not predict.merMOD
                   type = "response", 
                   allow.new.levels = TRUE, 
                   na.action = na.pass)

#combine dataframes
predframeGUD_b <- data.frame(new_dataGUD_b, preds = predGUD_b$fit, se = predGUD_b$se.fit)

#get the ones that I want
#predframeGUD_subset_b <- predframeGUD_b %>% 
#  mutate(ID = seq.int(nrow(predframeGUD_b))) %>%
#  filter(ID %in% c(1, 8, 11, 18, 22, 32, 29, 30, 39, 40))


##plot it
ggplot() +
  geom_point(data=predframeGUD_subset_b, mapping=aes(x=Region:Haven, y=preds, shape=Treatment, group=Haven:Treatment, colour=Haven), size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(data=predframeGUD_subset_b, mapping=aes(x=Region:Haven, ymin=preds-se, ymax=preds+se, group=Haven:Treatment, colour=Haven), linewidth=1, width=0.4, alpha=0.5, position=position_dodge(width=0.5)) +
  ylab("Probability of all peanuts being taken") +
  xlab("Region") +
  scale_colour_manual(breaks = c("Havened", "Non-havened"),
                      values= c("darkorange3", "springgreen4")) +
  scale_x_discrete(labels= c(expression(paste("Havened \nDryandra")), 
                             expression(paste("Non-havened \nDryandra")), 
                             expression(paste("Havened \nUpper Warren")), 
                             expression(paste("Non-havened \nUpper Warren")))) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
                # axis.title.x = element_text(vjust = -3.5),
                axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
                text = element_text(size=12),
                legend.position = "none")





##try new plot

GUD_all_plot <- GUD_all %>%
  filter(!is.na(Haven), !is.na(Site), !is.na(Treatment))

GUD_all_plot$Site <- as.factor(GUD_all_plot$Site)
GUD_all_plot$Haven <- as.factor(GUD_all_plot$Haven)
GUD_all_plot$Treatment <- as.factor(GUD_all_plot$Treatment)

GUD <- ggplot() +
  geom_violin(data=GUD_all_plot, mapping=aes(x=Treatment, y=All_gone, group=Region:Haven:Treatment, fill=Haven), position=position_dodge(width=0.5), alpha=0.2) +
  geom_point(data=predframeGUD_b, mapping=aes(x=Treatment, y=preds, shape=Treatment, group=Region, colour=Haven), size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(data=predframeGUD_b, mapping=aes(x=Treatment, ymin=preds-se, ymax=preds+se, group=Region, colour=Haven), linewidth=1, width=0.4, alpha=0.5, position=position_dodge(width=0.5)) +
  facet_grid(~ Region + Haven, switch = "x") +
  ylab("Probability of all peanuts being taken") +
  xlab("Treatment") +
  scale_colour_manual(breaks = c("Havened", "Non-havened"),
                      values= c("darkorange3", "springgreen4")) +

  scale_fill_manual(breaks = c("Havened", "Non-havened"),
                      values= c("darkorange3", "springgreen4")) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
                # axis.title.x = element_text(vjust = -3.5),
                axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
                text = element_text(size=12),
                legend.position = "none")




GUD_plot <- ggarrange(GUD, DI_TA, nrow =1, ncol = 2, labels =c("a", "b"))


#ggsave("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\R_outputs_woylie_comparison\\GUD_DITA_draft.pdf",
#       plot = GUD_plot,
#       scale = 1.2,
#       width = 18,
#       height = 9,
#       units = "cm")


ggsave("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\PhD Chapters\\Woylie comparison\\R_outputs_woylie_comparison\\GUD_DITA_violin.pdf",
       plot = GUD_plot,
       scale = 1.2,
       width = 18,
       height = 9,
       units = "cm")
