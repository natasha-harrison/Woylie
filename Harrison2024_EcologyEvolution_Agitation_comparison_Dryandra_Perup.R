##############################################################################################

# Compare agitation scores of havened/non-havened individuals at Dryandra and Perup

##############################################################################################

library(dplyr)
library(lme4)
library(ciTools)
library(DHARMa)
library(lmerTest)
library(performance)
library(data.table)


#############################
#       data prep
#############################


##################   UPPER WARREN   #####################
ps<- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Trapping_data\\PerupSanctuary_2022.csv")# MORPHOLOGY
boy <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Trapping_data\\Boyicup_2022.csv")# MORPHOLOGY
moo <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Trapping_data\\Moopinup_2022.csv")# MORPHOLOGY

#join them
warren <- bind_rows(ps, boy, moo)

#subset of relevant columns

warren_subset <- warren %>%
  filter(AGE != "J") %>%
  select(TRP_DATE, SPT_LABEL, SSI_LABEL, ANIMAL_NO, IMPLANT_NO, LEFT_ID, RIGHT_ID, SEX, REPRO_STATUS, PY_CR_LEN, AGE)

warren_subset$TRP_DATE <- as.Date(warren_subset$TRP_DATE, format = "%d/%m/%Y")

#join to behaviour by trap
warren_ctb <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Cage trap behaviour\\PerupWoylieCageBehaviour_2022.csv")

warren_ctb$Date <- as.Date(warren_ctb$Date, format = "%d/%m/%Y")

warren_ctb <- warren_ctb %>%
  mutate(Trap_ammended = ifelse(Population == "Boyicup", paste("BOY", Trap, sep=""),
                                ifelse(Population== "Moopinup", paste("POS", Trap, sep=""), Trap)))

##join by trap and date
perup_agitation <- inner_join(warren_subset, warren_ctb, by= c("SPT_LABEL" = "Trap_ammended",
                                                               "TRP_DATE" = "Date"))

perup_agitation_subset <- perup_agitation %>%
  filter(!is.na(Approach),
         !is.na(Bag_on),
         !is.na(Door_open),
         !is.na(Bag_before),
         !is.na(Bag_during)) %>%
  select(TRP_DATE, SPT_LABEL, SSI_LABEL, ANIMAL_NO, IMPLANT_NO, LEFT_ID, RIGHT_ID, SEX, REPRO_STATUS, PY_CR_LEN, AGE, 
         Time, Handler, Approach, Bag_on, Door_open, Bag_before, Bag_during, Vocalise, Heavy_breathing, Trap_damage, Joey_eject)

names(perup_agitation_subset) <- c("Date", "Trap", "Site", "Animal_No", "Microchip", "Left_ear", "Right_ear", "Sex", "Repro_status", "PY_CR", "Age",
  "Time", "Handler", "Approach", "Bag_on", "Door_open", "Bag_before", "Bag_during", "Vocalise", "Heavy_breathing", "Trap_damage", "Joey_eject") 

#add in number of captures before this session
capture_perup_count <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Dryandra Experimental Release\\capture_perup_count.csv")

perup_agitation_subset <- left_join(perup_agitation_subset, capture_perup_count[,2:3], by=c("Animal_No" = "ANIMAL_NO"))


#make the NAs 0
perup_agitation_subset$No_Caps[is.na(perup_agitation_subset$No_Caps)] <- 0

###################    DRYANDRA   #####################

Dryandra_trap_data <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Dryandra Experimental Release\\Harvest_data_DNWS_DM.csv")
#get the date right
Dryandra_trap_data$Date <- as.Date(Dryandra_trap_data$Date, format = "%d/%m/%Y")

d_ctb <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Dryandra Experimental Release\\DryandraTranslocation_CageTrapBehaviour2.csv")
#get the date right
d_ctb$Date <- as.Date(d_ctb$Date, format = "%Y-%m-%d")


##join by trap and date

dryandra_agitation <- inner_join(Dryandra_trap_data, d_ctb, by=c("Date" = "Date",
                                                                 "Trap" = "Trap"))

dryandra_agitation_subset <- dryandra_agitation %>%
  filter(!is.na(Approach),
         !is.na(Bag_on),
         !is.na(Door_open),
         !is.na(Bag_before),
         !is.na(Bag_during)) %>%
  select(Date, Trap, Site, Animal_No, Microchip.x, Left_ear.x, Right_ear.x, Sex, Repro_status, PY_CR, Age, 
         Time, Handler, Approach, Bag_on, Door_open, Bag_before, Bag_during, Vocalise, Heavy_breathing, Trap_damage, Joey_eject)

names(dryandra_agitation_subset) <- c("Date", "Trap", "Site", "Animal_No", "Microchip", "Left_ear", "Right_ear", "Sex", "Repro_status", "PY_CR", "Age",
                                   "Time", "Handler", "Approach", "Bag_on", "Door_open", "Bag_before", "Bag_during", "Vocalise", "Heavy_breathing", "Trap_damage", "Joey_eject") 

#add in number of captures before this session
Dryandra_all_captures_counts <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\Data\\Dryandra Experimental Release\\Previous_capture_counts.csv")

Dryandra_all_captures_counts <- Dryandra_all_captures_counts %>%
  filter(nchar(ANIMAL_NO) >0) %>%
  select(ANIMAL_NO, No_Caps)


dryandra_agitation_subset <- left_join(dryandra_agitation_subset, Dryandra_all_captures_counts, by= c("Animal_No" = "ANIMAL_NO"), na.matches = "never")

#make the NAs 0
dryandra_agitation_subset$No_Caps[is.na(dryandra_agitation_subset$No_Caps)] <- 0




#JOIN THEM BOTH
agitation_all <- rbind(perup_agitation_subset, dryandra_agitation_subset)


#NO CAPTURSE WITHIN THIS SESSION
##need to sum the number in the session, remove those from total, then add the code to the no_caps
is.data.table(agitation_all) #this needs to be true

agitation_all <- data.table(agitation_all) #make it a data table first

#then run this (count previous entries/captures), takes 30 secs to run
agitation_all[,.ct:=mapply(function(x,y) agitation_all[(Date<=x & Animal_No ==y),.N],agitation_all$Date,agitation_all$Animal_No)]
##this line sometimes errors because the dataframe has been "copied by R", can fix by changing the order of the commands with capture <- 

names(agitation_all)[names(agitation_all) == '.ct'] <- 'previous_captures'

agitation_all <- as.data.frame(agitation_all)


##ADD CAPTURES WITHIN THIS SESSION TO PREVIOUS CAPTURES 

agitation_all <- agitation_all %>%
  mutate(Total_prev_caps = No_Caps + previous_captures)





##  First work with the cumulative scores ##

#need to convert L, M, H to number scores

## Create replace functions
fn1 <- function(x) gsub("L", "1", x)
fn2 <- function(x) gsub("M", "2", x)
fn3 <- function(x) gsub("H", "3", x)

## Turn it into a column-wise function; only apply to cols "a" and "b"
fncol1 <- plyr::colwise(fn1, .cols=c("Approach", "Bag_on", "Door_open", "Bag_before", "Bag_during"))
fncol2 <- plyr::colwise(fn2, .cols=c("Approach", "Bag_on", "Door_open", "Bag_before", "Bag_during"))
fncol3 <- plyr::colwise(fn3, .cols=c("Approach", "Bag_on", "Door_open", "Bag_before", "Bag_during"))

## Apply the function
cumulative <- fncol1(agitation_all)
cumulative <- fncol2(cumulative)
cumulative <- fncol3(cumulative)

#join back with other data
agitation_all[,14:18] <- cumulative

#make them numeric
agitation_all$Approach <- as.numeric(as.character(agitation_all$Approach))
agitation_all$Bag_on <- as.numeric(as.character(agitation_all$Bag_on))
agitation_all$Door_open <- as.numeric(as.character(agitation_all$Door_open))
agitation_all$Bag_before <- as.numeric(as.character(agitation_all$Bag_before))
agitation_all$Bag_during <- as.numeric(as.character(agitation_all$Bag_during))

##add in cumulative column 
agitation_all <- agitation_all %>% mutate(Cumulative_score = rowSums(agitation_all[,14:18]))


## change Y/N to 0 and 1
agitation_all$Vocalise <- if_else(agitation_all$Vocalise == "Y", gsub("Y", "1", agitation_all$Vocalise), gsub("N", "0", agitation_all$Vocalise))
agitation_all$Heavy_breathing <- if_else(agitation_all$Heavy_breathing == "Y", gsub("Y", "1", agitation_all$Heavy_breathing), gsub("N", "0", agitation_all$Heavy_breathing))
agitation_all$Trap_damage <- if_else(agitation_all$Trap_damage == "Y", gsub("Y", "1", agitation_all$Trap_damage), gsub("N", "0", agitation_all$Trap_damage))
agitation_all$Joey_eject <- if_else(agitation_all$Joey_eject == "Y", gsub("Y", "1", agitation_all$Joey_eject), gsub("N", "0", agitation_all$Joey_eject))


agitation_all$Vocalise <- as.numeric(as.character(agitation_all$Vocalise))
agitation_all$Heavy_breathing <- as.numeric(as.character(agitation_all$Heavy_breathing))
agitation_all$Trap_damage <- as.numeric(as.character(agitation_all$Trap_damage))
agitation_all$Joey_eject <- as.factor(as.character(agitation_all$Joey_eject))

#replace the sites with legible words
agitation_all$Site <- gsub("51PSF11O20", "PerupSanctuary", agitation_all$Site)
agitation_all$Site <- gsub("51BOY/01", "Boyicup", agitation_all$Site)
agitation_all$Site <- gsub("51POS/01", "Moopinup", agitation_all$Site)


agitation_all[agitation_all == ""] <- NA

agitation_all <- agitation_all %>%
  mutate(Haven = ifelse(Site %in% c("DNWS", "PerupSanctuary"), "Havened", "Non-havened"),
         Region = ifelse(Site %in% c("PerupSanctuary", "Moopinup", "Boyicup"), "UpperWarren", "Dryandra"),
         Animal_ID = ifelse(!is.na(Animal_No), Animal_No,
                            ifelse(!is.na(Microchip), Microchip, Left_ear))) %>%
  filter(!is.na(Approach),
         !is.na(Bag_during),
         Sex != "U")

#write.csv(agitation_all, "C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\agitation_all.csv")
#agitation_all <- read.csv("C:\\Users\\21116718\\OneDrive - The University of Western Australia\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\agitation_all.csv")
agitation_all <- read.csv("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\PhD Chapters\\Woylie comparison\\Data_woylie_comparison\\agitation_all.csv")

###get mean for EPI

agitation_all_EPI <- agitation_all %>%
  mutate(Site_epi = ifelse(Site %in% c("LOL", "VE", "VN", "VW", "VS"), "Dryandra", Site))

agitation_all_EPI %>%
  group_by(Site_epi) %>%
  summarise(mean = mean(Cumulative_score, na.rm=T), sd= sd(Cumulative_score, na.rm=T), n=n(), se=sd/sqrt(n))



#####################################
#       let the modelling begin
#####################################

agitation1 <- lmer(Cumulative_score ~ Sex + Haven*Region + Total_prev_caps + (1|Handler) + (1|Animal_ID), data=agitation_all, na.action="na.omit")

summary(agitation1)

agitation_all %>%
  group_by(Haven, Region) %>%
  summarise(mean(Cumulative_score, na.rm=T))


##post hoc tukey pairwise test
tuk <- emmeans(agitation1, ~ Haven*Region)

#fit2.emm.a
pairs(tuk, adjust="tukey")








#variance inflation factor to check for collinearity
check_collinearity(agitation1)

#test dispersion
testDispersion(agitation1) #fine

#Check residuals 
simulationOutput <- simulateResiduals(fittedModel = agitation1, plot = F)
plot(simulationOutput) #looks good

plotResiduals(agitation1)


##loglikilihood for variables significance
agitation_sex <- lmer(Cumulative_score ~ Haven*Region + Total_prev_caps + (1|Handler) + (1|Animal_ID), data=agitation_all)

agitation_pc <- lmer(Cumulative_score ~ Sex + Haven*Region +  (1|Handler) + (1|Animal_ID), data=agitation_all)

agitation_haven <- lmer(Cumulative_score ~ Sex + Region + Total_prev_caps + (1|Handler) + (1|Animal_ID), data=agitation_all)

agitation_reg <- lmer(Cumulative_score ~ Sex + Haven + Total_prev_caps + (1|Handler) + (1|Animal_ID), data=agitation_all)

agitation_int <- lmer(Cumulative_score ~ Sex + Haven + Region + Total_prev_caps + (1|Handler) + (1|Animal_ID), data=agitation_all)



anova(agitation1, agitation_sex)
anova(agitation1, agitation_pc)
anova(agitation1, agitation_haven)
anova(agitation1, agitation_reg)
anova(agitation1, agitation_int)









#model predictions

newdata_cumscore <- expand.grid(Region = unique(agitation_all$Region),
                                 Sex= c("F", "M"),
                                 Haven = c("Havened", "Non-havened"),
                                 Age ="A",
                                Total_prev_caps = mean(agitation_all$Total_prev_caps, na.rm=T),
                                 Handler = "Tash Harrison",
                                 Microchip="956000014811258",
                               # Animal_ID="19296ANEAA"
                                Animal_ID="19365ANEAA"
)

#this ciTools function gives predictions and confidence intervals
predframe_cumscore <- add_ci(newdata_cumscore, agitation1, alpha = 0.05, names = c("lwr", "upr"), includeRanef=F, allow.new.levels=T)


agitation_toplot <- agitation_all %>%
  filter(!is.na(Sex), !is.na(Haven), !is.na(Region)) 
  
agitation_toplot$Haven <- as.factor(agitation_toplot$Haven)
agitation_toplot$Sex <- as.factor(agitation_toplot$Sex)
agitation_toplot$Region <- as.factor(agitation_toplot$Region)

agitation_violin <- ggplot() +
  #geom_point(data=agitation_toplot, mapping=aes(x=Region, y=Cumulative_score, colour=Haven, group=Haven:Sex, shape=Sex), position=position_jitterdodge(dodge.width=0.5, jitter.width=0.2, jitter.height=0.2), size=1, alpha=0.4) +
  geom_violin(data=agitation_toplot, mapping=aes(x=Region, y=Cumulative_score, fill=Haven, group=Region:Haven:Sex), position=position_dodge(width=0.95), alpha=0.2)+
  geom_errorbar(data=predframe_cumscore, mapping=aes(x=Region, ymin=lwr, ymax = upr, colour=Haven, group=Haven:Sex), position=position_dodge(width=0.95), linewidth=1, width=0.4, alpha=0.7) +
  geom_point(data=predframe_cumscore, mapping=aes(x=Region, y=pred, colour=Haven, group=Haven:Sex, shape=Sex), position=position_dodge(width=0.95), size=3) +
  scale_colour_manual(values= c("darkorange3", "springgreen4")) +
  scale_fill_manual(values= c("darkorange3", "springgreen4")) +
  scale_shape_manual(values= c(15, 16)) +
  xlab("Region") +
  ylab("Cumulative agitation score") + 
  #ylim(1000,1350) +
  theme(        panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect("white"),
                text = element_text(size=12),
                legend.position= "top")




ggsave("C:\\Users\\TashHarrison\\OneDrive - Department of Biodiversity, Conservation and Attractions\\UWA\\PhD Chapters\\Woylie comparison\\R_outputs_woylie_comparison\\Agitation.pdf",
       plot = agitation_violin,
       scale = 1,
       width = 18,
       height = 10,
       units = "cm")




#get cohen's d 

dry_ag <- agitation_all %>%
  filter(Region == "Dryandra") %>%
  select(Cumulative_score, Haven)

cohen.d(Cumulative_score ~ Haven, dry_ag)


per_ag <- agitation_all %>%
  filter(Region == "UpperWarren") %>%
  select(Cumulative_score, Haven)

cohen.d(Cumulative_score ~ Haven, per_ag)








#####################################
#       joey ejections     NOT CONVERGING
#####################################


joey_data <- agitation_all %>%
  filter(PY_CR > 0)

joey_data$Joey_eject <- as.factor(joey_data$Joey_eject)

#ideal model
joey1 <- glmmTMB(Joey_eject ~ PY_CR + Total_prev_caps + Haven*Region + (1|Animal_No) + (1|Handler),
                 family="binomial", data=joey_data)

summary(joey1)


##get proportions of each
table(joey_data$Haven, joey_data$Region, joey_data$Joey_eject)

#ideal model minus handler
joey1 <- glmmTMB(as.numeric(Joey_eject) ~ PY_CR + Total_prev_caps + Haven*Region + (1|Animal_No),
                 family="nbinom1", data=joey_data)

summary(joey1)


##remove each variable
joey_PY <- glmmTMB(as.numeric(Joey_eject) ~ Total_prev_caps + Haven*Region + (1|Animal_No),
                 family="nbinom1", data=joey_data)

joey_pc <- glmmTMB(as.numeric(Joey_eject) ~ PY_CR + Haven*Region + (1|Animal_No),
                   family="nbinom1", data=joey_data)

joey_haven <- glmmTMB(as.numeric(Joey_eject) ~ PY_CR + Total_prev_caps +Region + (1|Animal_No),
                   family="nbinom1", data=joey_data)


joey_region <- glmmTMB(as.numeric(Joey_eject) ~ PY_CR + Total_prev_caps + Haven + (1|Animal_No),
                   family="nbinom1", data=joey_data)

joey_int <- glmmTMB(as.numeric(Joey_eject) ~ PY_CR + Total_prev_caps + Haven +Region + (1|Animal_No),
                   family="nbinom1", data=joey_data)



anova(joey1, joey_PY)
anova(joey1, joey_pc)
anova(joey1, joey_haven)
anova(joey1, joey_region)
anova(joey1, joey_int)





#get cohen's d 

ej_d <- joey_data %>%
  filter(Region == "Dryandra") %>%
  select(Joey_eject, Haven)

cohen.d(as.numeric(Joey_eject) ~ Haven, ej_d)


ej_uw <- joey_data %>%
  filter(Region == "UpperWarren") %>%
  select(Joey_eject, Haven)

cohen.d(as.numeric(Joey_eject) ~ Haven, ej_uw)



