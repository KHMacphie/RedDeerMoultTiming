## Moult Timing as a predictor of Fecundity, Birth Weight, Calf Survival, and Hind Survival ##

rm(list = ls())

library(MCMCglmm)
library(ggplot2)
library(readr)
library(dplyr)
library(colourpicker)
library(gridExtra)


#Import dataset
CoatChange_asPredictor <- read_csv("Data_MoultAsPredictor.csv")


#Organize dataset

#Check to see if there is any missing data in any of the other columns
sum(is.na(CoatChange_asPredictor$JulyCoatStatus_y1)) #2506 missing
sum(is.na(CoatChange_asPredictor$NovCoatStatus_y1)) #411 missing
sum(is.na(CoatChange_asPredictor$Year1))
sum(is.na(CoatChange_asPredictor$ReproStatus_y2)) #102 missing
sum(is.na(CoatChange_asPredictor$HindAge_y1))


#Delete missing data
CoatChange_asPredictor_new <- CoatChange_asPredictor[!is.na(CoatChange_asPredictor$JulyCoatStatus_y1), ] 

CoatChange_asPredictor_new <- CoatChange_asPredictor_new[!is.na(CoatChange_asPredictor_new$NovCoatStatus_y1), ] 

CoatChange_asPredictor_new <- CoatChange_asPredictor_new[!is.na(CoatChange_asPredictor_new$ReproStatus_y2), ] 

#Center Year1 and HindAge_y1
CoatChange_asPredictor_new$HindAge_y1_cd <- CoatChange_asPredictor_new$HindAge_y1 - 7
CoatChange_asPredictor_new$Year1_cd <- CoatChange_asPredictor_new$Year1 - 2012


#Convert First Calf category in reproductive status to immature in both winter and summer 
CoatChange_asPredictor_new$ReproStat_y2a <- ifelse(as.character(CoatChange_asPredictor_new$ReproStatus_y2) == "FirstCalf",
                                                   "Immature", as.character(CoatChange_asPredictor_new$ReproStatus_y2))

#Change the reference level for ReproStat_y2a to be Milk instead of Immature 
CoatChange_asPredictor_new$ReproStat_y2a <- as.factor(CoatChange_asPredictor_new$ReproStat_y2a) #Convert it to an unordered factor so that the relevel function works
CoatChange_asPredictor_new$ReproStat_y2a <- relevel(CoatChange_asPredictor_new$ReproStat_y2a, ref = "Milk")

#Change reference level for July and November coat status
CoatChange_asPredictor_new$JulyCoatStatus_y1 <- as.factor(CoatChange_asPredictor_new$JulyCoatStatus_y1) #Convert it to an unordered factor so that the relevel function works
CoatChange_asPredictor_new$JulyCoatStatus_y1 <- relevel(CoatChange_asPredictor_new$JulyCoatStatus_y1, ref = "summer")

CoatChange_asPredictor_new$NovCoatStatus_y1 <- as.factor(CoatChange_asPredictor_new$NovCoatStatus_y1) #Convert it to an unordered factor so that the relevel function works
CoatChange_asPredictor_new$NovCoatStatus_y1 <- relevel(CoatChange_asPredictor_new$NovCoatStatus_y1, ref = "winter")


###########################################################################################################################################################################################
##FECUNDITY##


#Delete all missing data for the fecundity model
CoatChange_asPredictor_new1 <- CoatChange_asPredictor_new[!is.na(CoatChange_asPredictor_new$Fecundity_y2), ] 

#Make fecundity is treated as a factor instead of numeric 
CoatChange_asPredictor_new1$Fecundity_y2 <- factor(CoatChange_asPredictor_new1$Fecundity_y2)
str(CoatChange_asPredictor_new1$Fecundity_y2)


k<-10000
prior_f<-list(R=list(V=1,fix=1),
              G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

MCMCFe <- MCMCglmm(
  fixed = Fecundity_y2 ~ JulyCoatStatus_y1 + NovCoatStatus_y1 + Year1_cd + ReproStat_y2a + HindAge_y1_cd + I(HindAge_y1_cd^2),
  random = ~Year1 + HindCode,
  family = "threshold", 
  prior = prior_f,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(CoatChange_asPredictor_new1)
)

summary(MCMCFe)
plot(MCMCFe)


#Fecundity figures

pred_dat <- expand.grid(JulyCoatStatus_y1 = levels(CoatChange_asPredictor_new1$JulyCoatStatus_y1),
                        NovCoatStatus_y1 = levels(CoatChange_asPredictor_new1$NovCoatStatus_y1),
                        ReproStat_y2a = levels(CoatChange_asPredictor_new1$ReproStat_y2a),
                        HindAge_y1_cd=c(-5:15),
                        Year1_cd=c(-11:11)) 
pred_dat$Fecundity_y2<- c(rep(0,nrow(pred_dat)-1),1) 
pred_dat$Year1 <- as.factor(1)  
pred_dat$HindCode <- as.factor(1)   

pred_new <- predict.MCMCglmm(object=MCMCFe, newdata=pred_dat, interval="confidence") 
pred_new <- cbind(pred_dat,pred_new) 
View(pred_new)



m <- subset(pred_new, HindAge_y1_cd == 0 & Year1_cd == 0 & ReproStat_y2a == "Milk" & NovCoatStatus_y1 == "winter")
m <- distinct(m)
m$coatstatus <- m$JulyCoatStatus_y1
m$Census <- "July"

m2 <- subset(pred_new, HindAge_y1_cd == 0 & Year1_cd == 0 & ReproStat_y2a == "Milk" & JulyCoatStatus_y1 == "summer")
m2 <- distinct(m2)
m2$coatstatus <- m2$NovCoatStatus_y1
m2$Census <- "November"

F_both <- rbind(m, m2)

(FPlot <- ggplot(F_both, aes(x = coatstatus, y = fit, colour = Census)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(title = "A",
         x = "Moult Status in Year 1",
         y = "Predicted Fecundity (Year 2)") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)

(FPlot2 <- ggplot(F_both, aes(x = Census, y = fit, colour = Census, shape = coatstatus)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(title = "C",
         x = "Moult status in year 1",
         y = "Fecundity in year 2",
         shape = "Coat Status") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    scale_shape_manual(values = c("summer" = 16, "half" = 17, "winter" = 18)) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)

###############################################################################################################################
##BIRTH WEIGHT##


#Check for NAs in birth weight
sum(is.na(CoatChange_asPredictor_new$CalfBirthWt_y2)) #1700 missing
CoatChange_asPredictor_new2 <- CoatChange_asPredictor_new[!is.na(CoatChange_asPredictor_new$CalfBirthWt_y2), ] 

#Check for NAs in birth date
sum(is.na(CoatChange_asPredictor_new2$CalfBirthDate_y2))#2 missing
CoatChange_asPredictor_new2 <- CoatChange_asPredictor_new2[!is.na(CoatChange_asPredictor_new2$CalfBirthDate_y2), ] 

#Center birth date
mean(CoatChange_asPredictor_new2$CalfBirthDate_y2)
CoatChange_asPredictor_new2$CalfBirthDate_y2cd <- CoatChange_asPredictor_new2$CalfBirthDate_y2 - 155

#Make sure there isn't any NAs and that it is treated as a factor instead of numeric 
sum(is.na(CoatChange_asPredictor_new2$CalfSex_y2))
CoatChange_asPredictor_new2$CalfSex_y2 <- factor(CoatChange_asPredictor_new2$CalfSex_y2)
str(CoatChange_asPredictor_new2$CalfSex_y2)

#Define the prior
k<-10000
priorBW<-list(R=list(V=1,fix=1),
              G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


MCMCBW <- MCMCglmm(
  fixed = CalfBirthWt_y2 ~ JulyCoatStatus_y1 + NovCoatStatus_y1 + Year1_cd + ReproStat_y2a + HindAge_y1_cd + I(HindAge_y1_cd^2) +
    CalfBirthDate_y2cd + CalfSex_y2,
  random = ~Year1 + HindCode,
  family = "gaussian", 
  prior = priorBW,
  nitt=100000, thin=75, burnin=7000,
  data = CoatChange_asPredictor_new2
)

summary(MCMCBW)
plot(MCMCBW)


#Summer birth weight table

intercept <- MCMCBW$Sol[,1]
halfcc <- rowSums(MCMCBW$Sol[,1:2])      # summer + half
wintercc <- rowSums(MCMCBW$Sol[,c(1,3)]) # summer + winter

# Prepare data frame
BW <- data.frame(
  coatstatus = c("summer", "half", "winter"),
  mean = c(mean(intercept), mean(halfcc), mean(wintercc)),
  lowerCI = c(HPDinterval(mcmc(intercept))[1],
              HPDinterval(mcmc(halfcc))[1],
              HPDinterval(mcmc(wintercc))[1]),
  higherCI = c(HPDinterval(mcmc(intercept))[2],
               HPDinterval(mcmc(halfcc))[2],
               HPDinterval(mcmc(wintercc))[2])
)

BW$Census <- "July"

View(BW)

# Converting the categories to factor allowed me to re-order the levels to be
#consistent with the other graphs 
BW$Coatstatus <- factor(BW$coatstatus, levels = c("summer", "half", "winter"))
str(BW$coatstatus)

#Winter birth weight table

intercept <- MCMCBW$Sol[,1]
halfcc <- rowSums(MCMCBW$Sol[,c(1,4)])      # winter + half
summercc <- rowSums(MCMCBW$Sol[,c(1,5)]) # winter + summer

# Prepare data frame
WBW <- data.frame(
  coatstatus = c("winter", "half", "summer"),
  mean = c(mean(intercept), mean(halfcc), mean(summercc)),
  lowerCI = c(HPDinterval(mcmc(intercept))[1],
              HPDinterval(mcmc(halfcc))[1],
              HPDinterval(mcmc(summercc))[1]),
  higherCI = c(HPDinterval(mcmc(intercept))[2],
               HPDinterval(mcmc(halfcc))[2],
               HPDinterval(mcmc(summercc))[2])
)

WBW$Census <- "November"

View(WBW)

# Converting the categories to factor allowed me to re-order the levels to be
#consistent with the other graphs 
WBW$Coatstatus <- factor(WBW$coatstatus, levels = c("winter", "half", "summer"))
str(WBW$coatstatus)

#Merge tables for the figure

BW_both <- rbind(BW, WBW)
View(BW_both)

(BWPlot <- ggplot(BW_both, aes(x = coatstatus, y = mean, colour = Census)) +
    geom_errorbar(aes(ymin = lowerCI, ymax = higherCI), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(5.5,7.0) +
    labs(title = "B",
         x = "Moult Status in Year 1",
         y = "Predicted Birth Weight in Year 2 (kg)") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)



(BWPlot2 <- ggplot(BW_both, aes(x = Census, y = mean, colour = Census, shape = coatstatus)) +
    geom_errorbar(aes(ymin = lowerCI, ymax = higherCI), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(5.5,7.0) +
    labs(title = "D",
         x = "Moult status in year 1",
         y = "Birth weight in year 2 (kg)",
         shape= "Coat Status") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    scale_shape_manual(values = c("summer" = 16, "half" = 17, "winter" = 18)) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)


#############################################################################################################################
##CALF SURVIVAL##

#Delete missing Calf survival data
sum(is.na(CoatChange_asPredictor_new$CalfSurvive_y1)) #1526 missing
CoatChange_asPredictor_new3 <- CoatChange_asPredictor_new[!is.na(CoatChange_asPredictor_new$CalfSurvive_y1), ] 

#Ensure calf survival binomial is a two level factor
CoatChange_asPredictor_new3$CalfSurvive_y1 <- as.factor(CoatChange_asPredictor_new3$CalfSurvive_y1)
str(CoatChange_asPredictor_new3$CalfSurvive_y1)

#Here I’ve just set all the FYS=1 to be N for deathtype so I can keep all rows with deathtype=N (N means natural)
CoatChange_asPredictor_new3$CalfDeath_SurvN <- ifelse(CoatChange_asPredictor_new3$CalfSurvive_y1==1,"N",CoatChange_asPredictor_new3$CalfDeathType_y1) 
CoatChange_asPredictor_new3 <- CoatChange_asPredictor_new3[-which(CoatChange_asPredictor_new3$CalfDeath_SurvN!="N"|is.na(CoatChange_asPredictor_new3$CalfDeath_SurvN)==T),]

#Delete missing values for calf birth date and birth weight
sum(is.na(CoatChange_asPredictor_new3$CalfBirthDate_y1)) #25 missing
sum(is.na(CoatChange_asPredictor_new3$CalfBirthWt_y1)) #168 missing

#Delete missing values 
CoatChange_asPredictor_new3 <- CoatChange_asPredictor_new3[!is.na(CoatChange_asPredictor_new3$CalfBirthDate_y1), ] 
CoatChange_asPredictor_new3 <- CoatChange_asPredictor_new3[!is.na(CoatChange_asPredictor_new3$CalfBirthWt_y1), ] 

#Center birth date and birth weight year 1
CoatChange_asPredictor_new3$CalfBirthDate_y1cd <- CoatChange_asPredictor_new3$CalfBirthDate_y1 - 155
CoatChange_asPredictor_new3$CalfBirthWt_y1cd <- CoatChange_asPredictor_new3$CalfBirthWt_y1 - 6.503008


#Define the prior
k<-10000
priorCS<-list(R=list(V=1,fix=1),
              G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

MCMCCAS <- MCMCglmm(
  fixed = CalfSurvive_y1 ~ JulyCoatStatus_y1 + NovCoatStatus_y1 + HindAge_y1_cd + I(HindAge_y1_cd^2) + 
    Year1_cd + CalfBirthDate_y1cd + CalfBirthWt_y1cd,
  random = ~Year1 + HindCode,
  family = "threshold", 
  prior = priorCS,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(CoatChange_asPredictor_new3)
)

summary(MCMCCAS)
par(mar = c(4, 4, 2, 2))
plot(MCMCCAS)


#Summer Calf Survival Graph

pred_dat2 <- expand.grid(JulyCoatStatus_y1 = levels(CoatChange_asPredictor_new3$JulyCoatStatus_y1),
                         NovCoatStatus_y1 = levels(CoatChange_asPredictor_new3$NovCoatStatus_y1),
                         HindAge_y1_cd=c(-1:1),
                         Year1_cd=c(-1:1),
                         CalfBirthDate_y1cd = c(-1:1),
                         CalfBirthWt_y1cd = c(-1:1)) 
pred_dat2$CalfSurvive_y1 <- c(rep(0,nrow(pred_dat2)-1),1)
pred_dat2$Year1 <- as.factor(1) 
pred_dat2$HindCode <- as.factor(1)   

pred_new2 <- predict.MCMCglmm(object=MCMCCAS, newdata=pred_dat2, interval="confidence") 
pred_new2 <- cbind(pred_dat2,pred_new2)
View(pred_new2)



m3 <- subset(pred_new2, HindAge_y1_cd == 0 & Year1_cd == 0 & CalfBirthDate_y1cd == 0 & CalfBirthWt_y1cd == 0 & NovCoatStatus_y1 == "winter")
m3 <- distinct(m3)
m3$coatstatus <- m3$JulyCoatStatus_y1
m3$Census <- "July"
View(m3)

m4 <- subset(pred_new2, HindAge_y1_cd == 0 & Year1_cd == 0 & CalfBirthDate_y1cd == 0 & CalfBirthWt_y1cd == 0 & JulyCoatStatus_y1 == "summer")
m4 <- distinct(m4)
m4$coatstatus <- m4$NovCoatStatus_y1
m4$Census <- "November"
View(m4)

CS_both <- rbind(m3, m4)

(CSPlot <- ggplot(CS_both, aes(x = coatstatus, y = fit, colour = Census)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(title = "C",
         x = "Moult Status in Year 1",
         y = "Probability of Calf Survival in Year 1") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)

(CSPlot2 <- ggplot(CS_both, aes(x = Census, y = fit, colour = Census, shape = coatstatus)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(title = "B",
         x = "Moult status in year 1",
         y = "Probability of calf survival in year 1",
         shape = "Coat Status") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    scale_shape_manual(values = c("summer" = 16, "half" = 17, "winter" = 18)) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)


############################################################################################################################################################################
##HIND SURVIVAL##


#Delete missing hind survival data
sum(is.na(CoatChange_asPredictor_new$HindSurvive_y1)) # 1 missing
CoatChange_asPredictor_new4 <- CoatChange_asPredictor_new[!is.na(CoatChange_asPredictor_new$HindSurvive_y1), ] 

#Ensure hind survival binomial is a two level factor
CoatChange_asPredictor_new4$HindSurvive_y1 <- as.factor(CoatChange_asPredictor_new4$HindSurvive_y1)
str(CoatChange_asPredictor_new4$HindSurvive_y1)

#Change the reference level for ReproStat_y2a to be True yeld instead of Milk 
CoatChange_asPredictor_new4$ReproStat_y2a <- as.factor(CoatChange_asPredictor_new4$ReproStat_y2a) #Convert it to an unordered factor so that the relevel function works
CoatChange_asPredictor_new4$ReproStat_y2a <- relevel(CoatChange_asPredictor_new4$ReproStat_y2a, ref = "Winter yeld")

#Define the prior
k<-10000
priorHS<-list(R=list(V=1,fix=1),
              G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

MCMCHS <- MCMCglmm(
  fixed = HindSurvive_y1 ~ JulyCoatStatus_y1 + NovCoatStatus_y1 + ReproStat_y2a + HindAge_y1_cd + I(HindAge_y1_cd^2) + 
    Year1_cd, 
  random = ~Year1,
  family = "threshold", 
  prior = priorHS,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(CoatChange_asPredictor_new4)
)

summary(MCMCHS)
par(mar = c(4, 4, 2, 2))
plot(MCMCHS)

save(MCMCHS, file = "Data/MCMCHS.RData")
load("Data/MCMCHS.RData")


#Adult Survival Graphs

pred_dat3 <- expand.grid(JulyCoatStatus_y1 = levels(CoatChange_asPredictor_new4$JulyCoatStatus_y1),
                         NovCoatStatus_y1 = levels(CoatChange_asPredictor_new4$NovCoatStatus_y1),
                         ReproStat_y2a = levels(CoatChange_asPredictor_new4$ReproStat_y2a),
                         HindAge_y1_cd=c(-1:1),
                         Year1_cd=c(-1:1)) 
pred_dat3$HindSurvive_y1 <- c(rep(0,nrow(pred_dat3)-1),1) 
pred_dat3$Year1 <- as.factor(1)  

pred_new3 <- predict.MCMCglmm(object=MCMCHS, newdata=pred_dat3, interval="confidence") 
pred_new3 <- cbind(pred_dat3,pred_new3) 
View(pred_new3)

m5 <- subset(pred_new3, HindAge_y1_cd == 0 & Year1_cd == 0 & ReproStat_y2a == "Winter yeld" & NovCoatStatus_y1 == "winter")
m5 <- distinct(m5)
m5$coatstatus <- m5$JulyCoatStatus_y1
m5$Census <- "July"

m6 <- subset(pred_new3, HindAge_y1_cd == 0 & Year1_cd == 0 & ReproStat_y2a == "Winter yeld" & JulyCoatStatus_y1 == "summer")
m6 <- distinct(m6)
m6$coatstatus <- m6$NovCoatStatus_y1
m6$Census <- "November"

HS_both <- rbind(m5, m6)


(HSPlot <- ggplot(HS_both, aes(x = coatstatus, y = fit, colour = Census)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(title = "D",
         x = "Moult Status in Year 1",
         y = "Probability of Female Survival in Year 1") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)

(HSPlot2 <- ggplot(HS_both, aes(x = Census, y = fit, colour = Census, shape = coatstatus)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 5) +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(title = "A",
         x = "Moult status in year 1",
         y = "Probability of female survival in year 1",
         shape = "Coat Status") +
    scale_color_manual(values = c("July" = "#81C784", "November" = "steelblue3")) +
    scale_shape_manual(values = c("summer" = 16, "half" = 17, "winter" = 18)) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)
    )
)



############################################################################################################################################################################

#Places all of the graphs into one figure 


#Version 1: X-axis says each moult status (summer, half, winter) with each census color coded side by side
grid.arrange(FPlot, BWPlot, CSPlot, HSPlot, ncol = 2, nrow = 2)


#Version 2: X-axis says each census (July and November) with each moult status side by side organised by shape. Census still color coded 
grid.arrange(HSPlot2, CSPlot2, FPlot2, BWPlot2, ncol = 2, nrow = 2)


##END##