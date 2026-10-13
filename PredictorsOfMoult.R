## Predictors of moult timing ##

rm(list = ls())

#Packages installed

library(lme4)
library(MCMCglmm)
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)

##SUMMER MODELS##

#Import dataset

Summer_CoatChange2 <- read_csv("Data_SummerMoultPrediction.csv")

#Convert summer coat status from categorical to binomial
sum(is.na(Summer_CoatChange2$CoatStatus))
Summer_CoatChange2$CoatStatus_binomial <- ifelse(Summer_CoatChange2$CoatStatus %in% c("summer", "half"), 1, 0)

#Ensure CoatStatus_binomial is considered a factor
Summer_CoatChange2$CoatStatus_binomial <- factor(Summer_CoatChange2$CoatStatus_binomial)
str(Summer_CoatChange2$CoatStatus_binomial)

#Center CensusYear and HindAge
Summer_CoatChange2$HindAge_cd <- Summer_CoatChange2$HindAge - 7
Summer_CoatChange2$CensusYear_cd <- Summer_CoatChange2$CensusYear - 2012

#Check to see if any columns have any missing values
sum(is.na(Summer_CoatChange2$CensusYear))  
sum(is.na(Summer_CoatChange2$ReproStatus_SameYear)) #5 missing
sum(is.na(Summer_CoatChange2$HindAge)) 
sum(is.na(Summer_CoatChange2$ReproStatus_PrevYear)) #7 missing

#Delete missing values 
Summer_CoatChange2_new <- Summer_CoatChange2[!is.na(Summer_CoatChange2$ReproStatus_SameYear), ] 
Summer_CoatChange2_new <- Summer_CoatChange2_new[!is.na(Summer_CoatChange2_new$ReproStatus_PrevYear), ] 

#Change the reference level for ReproStatus_SameYear and ReproStatus_PrevYear to be Milk instead of Immature 
Summer_CoatChange2_new$ReproStatus_SameYear <- as.factor(Summer_CoatChange2_new$ReproStatus_SameYear) #Convert it to an unordered factor so that the relevel function works
Summer_CoatChange2_new$ReproStatus_SameYear <- relevel(Summer_CoatChange2_new$ReproStatus_SameYear, ref = "Milk")

Summer_CoatChange2_new$ReproStatus_PrevYear <- as.factor(Summer_CoatChange2_new$ReproStatus_PrevYear) #Convert it to an unordered factor so that the relevel function works
Summer_CoatChange2_new$ReproStatus_PrevYear <- relevel(Summer_CoatChange2_new$ReproStatus_PrevYear, ref = "Milk")



#Summer Baysian Binomial GLMM 

#Fix character format by changing NaÃ¯ve to say Naive (uses dplyr package) 
Summer_CoatChange2_new <- Summer_CoatChange2_new %>%
  mutate(ReproStatus_PrevYear = recode(ReproStatus_PrevYear, "NaÃ¯ve" = "Naïve")) 

#Anytime it says Naive, change it to Immature 
Summer_CoatChange2_new$PrevYrRepro2 <- ifelse(as.character(Summer_CoatChange2_new$ReproStatus_PrevYear) == "Naïve",
                                              "Immature", as.character(Summer_CoatChange2_new$ReproStatus_PrevYear))


Summer_CoatChange2_new$Fec_SameYr <- ifelse(Summer_CoatChange2_new$ReproStatus_SameYear=="Milk",1,0)
Summer_CoatChange2_new$Fec_SameYr <- as.factor(Summer_CoatChange2_new$Fec_SameYr)
str(Summer_CoatChange2_new$Fec_SameYr)

#Re-level PrevYrRepro2
Summer_CoatChange2_new$PrevYrRepro2 <- as.factor(Summer_CoatChange2_new$PrevYrRepro2) #Convert it to an unordered factor so that the relevel function works
Summer_CoatChange2_new$PrevYrRepro2 <- relevel(Summer_CoatChange2_new$PrevYrRepro2, ref = "Milk")
str(Summer_CoatChange2_new$PrevYrRepro2)

#Define the prior
k<-10000
summerprior<-list(R=list(V=1,fix=1),
                  G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                         G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#Model predictors of moult
MCMCSummer_prev2 <- MCMCglmm(
  fixed = CoatStatus_binomial ~ CensusYear_cd + PrevYrRepro2*Fec_SameYr + HindAge_cd + I(HindAge_cd^2),
  random = ~CensusYear + Code,
  family = "threshold", 
  prior = summerprior,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(Summer_CoatChange2_new)
)

summary(MCMCSummer_prev2)
plot(MCMCSummer_prev2)



#######################################################################################################################
##WINTER MODEL##

Winter_CoatChange2 <- read_csv("Data_WinterMoultPrediction.csv")

#Convert winter coat status from categorical to binomial
sum(is.na(Winter_CoatChange2$CoatStatus)) #2 missing
Winter_CoatChange2_new <- Winter_CoatChange2[!is.na(Winter_CoatChange2$CoatStatus), ] #Delete NAs
Winter_CoatChange2_new$CoatStatus_binomial <- ifelse(Winter_CoatChange2_new$CoatStatus %in% c("winter", "half"), 1, 0)

#Ensure it is treated as a factor and not numerical
Winter_CoatChange2_new$CoatStatus_binomial <- as.factor(Winter_CoatChange2_new$CoatStatus_binomial)
str(Winter_CoatChange2_new$CoatStatus_binomial)

#Center CensusYear and HindAge
Winter_CoatChange2_new$HindAge_cd <- Winter_CoatChange2_new$HindAge - 7
Winter_CoatChange2_new$CensusYear_cd <- Winter_CoatChange2_new$CensusYear - 2005

#Check to see if any columns have any missing values
sum(is.na(Winter_CoatChange2_new$CensusYear))  
sum(is.na(Winter_CoatChange2_new$ReproStatus_SameYear))  
sum(is.na(Winter_CoatChange2_new$HindAge))  

#Delete missing values 
Winter_CoatChange2_new <- Winter_CoatChange2_new[!is.na(Winter_CoatChange2_new$ReproStatus_SameYear), ]

# Set the reference level to "Milk"
Winter_CoatChange2_new$ReproStatus_SameYear <- as.factor(Winter_CoatChange2_new$ReproStatus_SameYear) #Convert it to an unordered factor so that the relevel function works
Winter_CoatChange2_new$ReproStatus_SameYear <- relevel(Winter_CoatChange2_new$ReproStatus_SameYear, ref = "Milk")


#Winter Baysian Binomial GLMM 
#Define the prior
k<-10000
winterprior<-list(R=list(V=1,fix=1),
                  G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                         G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


MCMCWinter <- MCMCglmm(
  fixed = CoatStatus_binomial ~ CensusYear_cd + ReproStatus_SameYear + HindAge_cd + I(HindAge_cd^2),
  random = ~CensusYear + Code,
  family = "threshold", 
  prior = winterprior,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(Winter_CoatChange2_new)
)

summary(MCMCWinter)
plot(MCMCWinter)

#######################################################################################################################
##SUMMER AND WINTER PARASITE LOAD MODELS##

SummerParasite_new <- Summer_CoatChange2_new[!is.na(Summer_CoatChange2_new$StrongylesSpring), ] 

WinterParasite_new <- Winter_CoatChange2_new[!is.na(Winter_CoatChange2_new$StrongylesSummer), ]
WinterParasite_new <- WinterParasite_new[!is.na(WinterParasite_new$StrongylesAutumn), ]

#Center CensusYear, HindAge, and Strongyles
SummerParasite_new$CensusYear_cd <- SummerParasite_new$CensusYear - 2020
SummerParasite_new$StrongylesSpring_cd <- SummerParasite_new$StrongylesSpring - 15

WinterParasite_new$CensusYear_cd <- WinterParasite_new$CensusYear - 2019
WinterParasite_new$StrongylesSummer_cd <- WinterParasite_new$StrongylesSummer - 6
WinterParasite_new$StrongylesAutumn_cd <- WinterParasite_new$StrongylesAutumn - 1

str(SummerParasite_new)
str(WinterParasite_new)


#Summer Parasite model

k<-10000
summerprior<-list(R=list(V=1,fix=1),
                  G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                         G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


MCMCSummerPS <- MCMCglmm(
  fixed = CoatStatus_binomial ~ CensusYear_cd + PrevYrRepro2*Fec_SameYr + HindAge_cd + I(HindAge_cd^2) + StrongylesSpring_cd,
  random = ~CensusYear + Code,
  family = "threshold", 
  prior = summerprior,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(SummerParasite_new)
)

summary(MCMCSummerPS)
plot(MCMCSummerPS)



#Winter Parasite model
k<-10000
winterprior<-list(R=list(V=1,fix=1),
                  G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                         G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


MCMCWinterPS <- MCMCglmm(
  fixed = CoatStatus_binomial ~ CensusYear_cd + ReproStatus_SameYear + HindAge_cd + I(HindAge_cd^2) + 
    StrongylesSummer_cd + StrongylesAutumn_cd,
  random = ~CensusYear + Code,
  family = "threshold", 
  prior = winterprior,
  nitt=100000, thin=75, burnin=7000,
  data = as.data.frame(WinterParasite_new)
)

summary(MCMCWinterPS)
plot(MCMCWinterPS)


#######################################################################################################################

## Figures##


# MCMC values for summer census year variance
SDYear <- data.frame(var1 = MCMCSummer_prev2$VCV[,1])
#View(SDYear)

ggplot(SDYear, aes(x = var1)) +
  geom_density(fill = "lightgreen", alpha = 0.6) +
  labs(
    title = "A",
    x = "Variance",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)


# MCMC values for summer individual hind variance
SDHind <- data.frame(var1 = MCMCSummer_prev2$VCV[,2])
#View(SDHind)

ggplot(SDHind, aes(x = var1)) +
  geom_density(fill = "lightgreen", alpha = 0.6) +
  labs(
    title = "B",
    x = "Variance",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)


# MCMC values for summer census year variance
WDYear <- data.frame(var1 = MCMCWinter$VCV[,1])
#View(WDYear)

ggplot(WDYear, aes(x = var1)) +
  geom_density(fill = "steelblue3", alpha = 0.6) +
  labs(
    title = "C",
    x = "Variance",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)


# MCMC values for summer census year variance
WDHind <- data.frame(var1 = MCMCWinter$VCV[,2])
#View(WDHind)

ggplot(WDHind, aes(x = var1)) +
  geom_density(fill = "steelblue3", alpha = 0.6) +
  labs(
    title = "D",
    x = "Variance",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)



#Now develop a graph that represents all random effects from each model
#all in one graph 

# Create labeled data frames for each set of values
summer1 <- data.frame(
  value = (MCMCSummer_prev2$VCV[,1]),
  model = "Summer",
  component = "Census Year"
)

summer2 <- data.frame(
  value = (MCMCSummer_prev2$VCV[,2]),
  model = "Summer",
  component = "Female ID"
)

winter1 <- data.frame(
  value = (MCMCWinter$VCV[,1]),
  model = "Winter",
  component = "Census Year"
)

winter2 <- data.frame(
  value = (MCMCWinter$VCV[,2]),
  model = "Winter",
  component = "Female ID"
)

# Combine all into one data frame
AllDensity <- bind_rows(summer1, summer2, winter1, winter2)


#Combined group label
AllDensity <- AllDensity %>%
  mutate(group = paste(model, component, sep = " "))

#View(AllDensity)



ggplot(AllDensity, aes(x = var1, fill = group)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "coral") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") + 
  scale_fill_manual(values = c("Summer Census Year" = "lightgreen",
                               "Summer Female ID" = "lightgreen",
                               "Winter Census Year" = "steelblue3",
                               "Winter Female ID" = "steelblue3")) +
  labs(
    x = "Variance",
    y = "Probability Density"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    strip.text = element_text(size = 12, face = "bold"), 
    legend.position = "none",  
    title = element_text(size = 16),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 13),  
    axis.text.y = element_text(size = 13)
  )
# saved as 7x7.5" portrait




##########################################################################################################

#SUMMER PLOTS

pred_dat <- expand.grid(PrevYrRepro2 = levels(Summer_CoatChange2_new$PrevYrRepro2),
                        Fec_SameYr = levels(Summer_CoatChange2_new$Fec_SameYr), 
                        HindAge_cd=c(-5:17),
                        CensusYear_cd=c(-11:12)) 
pred_dat$CoatStatus_binomial<- c(rep(0,nrow(pred_dat)-1),1) 
pred_dat$CensusYear <- as.factor(1)  
pred_dat$Code <- as.factor(1)   

pred_new <- predict.MCMCglmm(object=MCMCSummer_prev2, newdata=pred_dat, interval="confidence") 
pred_new <- cbind(pred_dat,pred_new) # combine the predictions with the dataframe they are predicted from
#View(pred_new)


#July Census Year Graph
m1 <- subset(pred_new, HindAge_cd == 0 & PrevYrRepro2 == "Milk"&Fec_SameYr==1)
m1 <- distinct(m1)
m1$CensusYear2 <- m1$CensusYear_cd + 2012

#The graph requires CoatStatus_binomial to be numeric for the geom_jitter layer, so I have to re-run
#the function that converted CoatStatus to CoatStatus_binomial, so that it reverts back to numeric instead of factor.
Summer_CoatChange2_new$CoatStatus_binomial <- ifelse(Summer_CoatChange2_new$CoatStatus %in% c("summer", "half"), 1, 0)
str(Summer_CoatChange2_new$CoatStatus_binomial)


(SCensus <- ggplot(m1, aes(x = CensusYear2, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#81C784", alpha = 0.1) +
    geom_line(color = "darkgreen") +
    geom_jitter(data = Summer_CoatChange2_new, aes(x = CensusYear, y = CoatStatus_binomial),
                color = "darkgreen", alpha = 0.2, size = 1.5, width = 1.0, height = 0.05,) +
    labs(title = "E", x = "Year", y = "") +
    theme_minimal() +
    ylim(0.0,1.0) +
    xlim(1987, 2024) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)  
    ))


##July Census Age Graph

m3 <- subset(pred_new, CensusYear_cd == 0 & PrevYrRepro2 == "Milk"&Fec_SameYr==1)
m3 <- distinct(m3)
m3$HindAge <- m3$HindAge_cd + 7
#View(m3)


(SAge <- ggplot(m3, aes(x = HindAge, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#81C784", alpha = 0.1) +
    geom_line(color = "darkgreen") +
    geom_jitter(data = Summer_CoatChange2_new, aes(x = HindAge, y = CoatStatus_binomial),
                color = "darkgreen", alpha = 0.2, size = 1.5, width = 1.0, height = 0.05,) +
    labs(title = "A", x = "Age (years)", y = "") +
    theme_minimal() +
    ylim(0.0,1.0) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)   
    )
)


#Previous year reproductive status plus same year fecundity summer graph

#Change CoatStatus_binomial back to a factor since the jitter is no longer needed 
Summer_CoatChange2$CoatStatus_binomial <- factor(Summer_CoatChange2$CoatStatus_binomial)
str(Summer_CoatChange2$CoatStatus_binomial)



pred_dat3 <- expand.grid(PrevYrRepro2 = levels(Summer_CoatChange2_new$PrevYrRepro2),
                         HindAge_cd=c(-5:17),
                         CensusYear_cd=c(-11:12),
                         Fec_SameYr= levels(Summer_CoatChange2_new$Fec_SameYr)) # a grid of combinations of the fixed effects, so include the range of each variable in your model
pred_dat3$CoatStatus_binomial<- c(rep(0,nrow(pred_dat3)-1),1) # it needs to include 0 and 1 for the Bernoulli models but doesn’t seem to matter how many
pred_dat3$CensusYear <- as.factor(1)  # Random term
pred_dat3$Code <- as.factor(1)   # Random term

pred_new3 <- predict.MCMCglmm(object=MCMCSummer_prev2, newdata=pred_dat3, interval="confidence") # using pred_dat as newdata
pred_new3 <- cbind(pred_dat3,pred_new3) # combine the predictions with the dataframe they are predicted from
#View(pred_new3)


m5 <- subset(pred_new3, HindAge_cd == 0 & CensusYear_cd ==0) 
m5 <- distinct(m5)
#View(m5)



(SPrevF <- ggplot(m5, aes(x = PrevYrRepro2, y = fit, color = Fec_SameYr)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), 
                  position = position_dodge(width = 0.5),
                  width = 0.3, size = 1) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    scale_color_manual(values = c("0" = "darkgreen", "1" = "#81C784")) +
    theme_minimal() +
    labs(
      title = "C",
      x = "Reproductive status",
      y = "Probability of moulting by census",
      color = "Fecundity") +
    ylim(0.0,1.0) +
    theme(
      title = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      legend.position = c(0.2,0.2)
    ))

#####################################################################################################################
#WINTER PLOTS


pred_dat4 <- expand.grid(ReproStatus_SameYear = levels(Winter_CoatChange2_new$ReproStatus_SameYear),
                         HindAge_cd=c(-5:16),
                         CensusYear_cd=c(-18:18)) 
pred_dat4$CoatStatus_binomial<- c(rep(0,nrow(pred_dat4)-1),1) 
pred_dat4$CensusYear <- as.factor(1)  
pred_dat4$Code <- as.factor(1)   

pred_new4 <- predict.MCMCglmm(object=MCMCWinter, newdata=pred_dat4, interval="confidence") 
pred_new4 <- cbind(pred_dat4,pred_new4) 
#View(pred_new4)



#November Census Year Graph
m6 <- subset(pred_new4, HindAge_cd == 0 & ReproStatus_SameYear == "Milk")
m6 <- distinct(m6)
m6$CensusYear2 <- m6$CensusYear_cd + 2005



#The graph requires CoatStatus_binomial to be numeric for the geom_jitter layer, so I have to re-run
#the function that converted CoatStatus to CoatStatus_binomial, so that it reverts back to numeric instead of factor.
Winter_CoatChange2_new$CoatStatus_binomial <- ifelse(Winter_CoatChange2_new$CoatStatus %in% c("winter", "half"), 1, 0)
str(Winter_CoatChange2_new$CoatStatus_binomial)

(WCensus <- ggplot(m6, aes(x = CensusYear2, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "steelblue", alpha = 0.1) +
    geom_line(color = "steelblue4") +
    geom_jitter(data = Winter_CoatChange2_new, aes(x = CensusYear, y = CoatStatus_binomial),
                color = "steelblue3", alpha = 0.2, size = 1.5, width = 1.0, height = 0.07,) +
    labs(title = "F", x = "Year", y = "") +
    theme_minimal() +
    ylim(0.0,1.0) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)   
    ))



#November Census Same Year Reproductive Status Graph
m7 <- subset(pred_new4, HindAge_cd == 0 & CensusYear_cd == 0)
m7 <- distinct(m7)

levels(m7$ReproStatus_SameYear) <- c( "Milk","Immature","SummerYeld","True Yeld" )

(WRepro <- ggplot(m7, aes(x = ReproStatus_SameYear, y = fit)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15, color = "steelblue3", size =1) +
    geom_point(size = 3, color = "steelblue3") +
    theme_minimal() +
    ylim(0.0,1.0) +
    labs(
      title = "D",
      x = "Reproductive status",
      y = "") +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14,angle = 45, hjust = 1),  
      axis.text.y = element_text(size = 14)  
    ))



##November Census Age Graph

m8 <- subset(pred_new4, CensusYear_cd == 0 & ReproStatus_SameYear == "Milk")
m8 <- distinct(m8)
m8$HindAge <- m8$HindAge_cd + 7
#View(m8)


(WAge <- ggplot(m8, aes(x = HindAge, y = fit)) +
    geom_line(color = "steelblue4") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "steelblue", alpha = 0.1) +
    geom_jitter(data = Winter_CoatChange2_new, aes(x = HindAge, y = CoatStatus_binomial),
                color = "steelblue3", alpha = 0.1, size = 1.5, width = 1.0, height = 0.04,) +
    labs(title = "B", x = "Age (years)", y = "") +
    theme_minimal() +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)   
    ))


###############################################################################################################################
##SUMMER AND WINTER PARASITE LOAD GRAPHS##


#Summer Parasite Graph

pred_dat5 <- expand.grid(PrevYrRepro2 = levels(SummerParasite_new$PrevYrRepro2),
                         Fec_SameYr = levels(SummerParasite_new$Fec_SameYr), 
                         HindAge_cd=c(-1:1),
                         CensusYear_cd=c(-1:1),
                         StrongylesSpring_cd=c(-15:219)) 
pred_dat5$CoatStatus_binomial<- c(rep(0,nrow(pred_dat5)-1),1) 
pred_dat5$CensusYear <- as.factor(1)  
pred_dat5$Code <- as.factor(1)   

pred_new5 <- predict.MCMCglmm(object=MCMCSummerPS, newdata=pred_dat5, interval="confidence") 
pred_new5 <- cbind(pred_dat5,pred_new5) 
#View(pred_new5)


m9 <- subset(pred_new5, HindAge_cd == 0 & PrevYrRepro2 == "Milk" & Fec_SameYr==1 & CensusYear_cd == 0)
m9 <- distinct(m9)
m9$StrongylesSpring2 <- m9$StrongylesSpring_cd + 15



SummerParasite_new$CoatStatus_binomial <- ifelse(SummerParasite_new$CoatStatus %in% c("summer", "half"), 1, 0)
str(SummerParasite_new$CoatStatus_binomial)


(SParasite <- ggplot(m9, aes(x = StrongylesSpring2, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#81C784", alpha = 0.2) +
    geom_line(color = "darkgreen") +
    geom_jitter(data = SummerParasite_new, aes(x = StrongylesSpring, y = CoatStatus_binomial),
                color = "darkgreen", alpha = 0.2, size = 1.5, width = 1.0, height = 0.05,) +
    labs(title = "A", x = "April FEC (eggs/g)", y = "Probability of moulting by July census") +
    theme_minimal() +
    ylim(0.0,1.0) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)  
    ))




#Winter Parasite Graph

pred_dat6 <- expand.grid(ReproStatus_SameYear = levels(WinterParasite_new$ReproStatus_SameYear),
                         HindAge_cd=c(-1:1),
                         CensusYear_cd=c(-1:1),
                         StrongylesSummer_cd=seq(-6, 84, 5),
                         StrongylesAutumn_cd=seq(-1, 95, 5)) 
pred_dat6$CoatStatus_binomial<- c(rep(0,nrow(pred_dat6)-1),1) 
pred_dat6$CensusYear <- as.factor(1)  
pred_dat6$Code <- as.factor(1)   

pred_new6 <- predict.MCMCglmm(object=MCMCWinterPS, newdata=pred_dat6, interval="confidence")
pred_new6 <- cbind(pred_dat6,pred_new6) 
#View(pred_new6)




m10 <- subset(pred_new6, HindAge_cd == 0 & ReproStatus_SameYear == "Milk" & CensusYear_cd == 0 & StrongylesAutumn_cd == 4)
m10 <- distinct(m10)
m10$Strongyles <- m10$StrongylesSummer_cd + 6
m10$Strongyles_S <- m10$StrongylesSummer_cd + 6
m10$Season <- "Summer"


m11 <- subset(pred_new6, HindAge_cd == 0 & ReproStatus_SameYear == "Milk" & CensusYear_cd == 0 & StrongylesSummer_cd == 4)
m11 <- distinct(m11)
m11$Strongyles <- m11$StrongylesAutumn_cd + 1
m11$Strongyles_A <- m11$StrongylesAutumn_cd + 1
m11$Season <- "Autumn"

#The graph requires CoatStatus_binomial to be numeric for the geom_jitter layer, so I have to re-run
#the function that converted CoatStatus to CoatStatus_binomial, so that it reverts back to numeric instead of factor.
WinterParasite_new$CoatStatus_binomial <- ifelse(WinterParasite_new$CoatStatus %in% c("winter", "half"), 1, 0)
str(WinterParasite_new$CoatStatus_binomial)


(WParasite_S <- ggplot(m10, aes(x = Strongyles_S, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "steelblue", alpha = 0.2) +
    geom_line(color = "steelblue4") +
    geom_jitter(data = WinterParasite_new, aes(x = StrongylesSummer, y = CoatStatus_binomial),
                color = "steelblue3", alpha = 0.2, size = 1.5, width = 1.0, height = 0.05,) +
    labs(title = "B", x = "August FEC (eggs/g)", y = "Probability of moulting by November census") +
    theme_minimal() +
    ylim(0.0,1.0) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)  
    ))


(WParasite_A <- ggplot(m11, aes(x = Strongyles_A, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "steelblue", alpha = 0.2) +
    geom_line(color = "steelblue4") +
    geom_jitter(data = WinterParasite_new, aes(x = StrongylesAutumn, y = CoatStatus_binomial),
                color = "steelblue3", alpha = 0.2, size = 1.5, width = 1.0, height = 0.05,) +
    labs(title = "C", x = "November FEC (eggs/g)", y = "Probability of moulting by November census") +
    theme_minimal() +
    ylim(0.0,1.0) +
    theme(
      title = element_text(size = 16),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14)  
    ))


###############################################################################################################################


#Place graphs into the same plot panel

library(gridExtra)

#Temporal, reproductive status and age
space <- ggplot()+theme_void()
grid.arrange(SAge,WAge,space,space,SPrevF,WRepro,space,space,SCensus,WCensus, ncol = 2, nrow = 5,heights=c(1,0.05,1.2,0.05,1))  #10*14"


#Summer and Winter Parasite Graphs
grid.arrange(SParasite, WParasite_S, WParasite_A, ncol = 3) #7x12"
