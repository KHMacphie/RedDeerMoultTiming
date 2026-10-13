#### Organise data for multivariate model ####


rm(list = ls())

library(MCMCglmm)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)



#Import dataframes

Summer_CoatChange2 <- read_csv("Data_SummerMoultPrediction.csv")

Winter_CoatChange2 <- read_csv("Data_WinterMoultPrediction.csv")

Parturition <- read_csv("Data_Parturition.csv") # NEW DATAFRAME WITH ALL PARTURITION DATES
Parturition <- subset(Parturition,OrdinalBirthDate>120&OrdinalBirthDate<213) #use May-Jul births

#Convert CoatStatus to binomial
sum(is.na(Summer_CoatChange2$CoatStatus)) #Just to make sure I am not converting any NAs to 0
Summer_CoatChange2$CoatStatus_binomial <- ifelse(Summer_CoatChange2$CoatStatus %in% c("summer", "half"), 1, 0)

sum(is.na(Winter_CoatChange2$CoatStatus)) #Just to make sure I am not converting any NAs to 0
Winter_CoatChange2_new <- Winter_CoatChange2[!is.na(Winter_CoatChange2$CoatStatus), ] 
Winter_CoatChange2_new$CoatStatus_binomial <- ifelse(Winter_CoatChange2_new$CoatStatus %in% c("winter", "half"), 1, 0)

#Check for missing values
sum(is.na(Summer_CoatChange2$CensusYear))  
sum(is.na(Summer_CoatChange2$HindAge)) 
sum(is.na(Summer_CoatChange2$ReproStatus_PrevYear)) #7 missing
sum(is.na(Summer_CoatChange2$ReproStatus_SameYear)) #5 missing

sum(is.na(Winter_CoatChange2_new$CensusYear))  
sum(is.na(Winter_CoatChange2_new$HindAge)) 
sum(is.na(Winter_CoatChange2_new$ReproStatus_SameYear)) #7 missing

#Delete rows with NA for reproductive status  
Summer_CoatChange2_new <- Summer_CoatChange2[!is.na(Summer_CoatChange2$ReproStatus_PrevYear), ] 
Summer_CoatChange2_new <- Summer_CoatChange2_new[!is.na(Summer_CoatChange2_new$ReproStatus_SameYear), ] 
Winter_CoatChange2_new <- Winter_CoatChange2_new[!is.na(Winter_CoatChange2_new$ReproStatus_SameYear), ] 

#Set CoatStatus_binomial as a 2 level factor for the threshold family to work
Summer_CoatChange2_new$CoatStatus_binomial <- as.factor(Summer_CoatChange2_new$CoatStatus_binomial)
str(Summer_CoatChange2_new$CoatStatus_binomial)

Winter_CoatChange2_new$CoatStatus_binomial <- as.factor(Winter_CoatChange2_new$CoatStatus_binomial)
str(Winter_CoatChange2_new$CoatStatus_binomial)

#Fix character format by changing NaÃ¯ve to say Naive (uses dplyr package) 
Summer_CoatChange2_new <- Summer_CoatChange2_new %>%
  mutate(ReproStatus_PrevYear = recode(ReproStatus_PrevYear, "NaÃ¯ve" = "Naive"))

#Add Fecundity and the updated version of previous reproductive status where Naive and Immature are combined
Summer_CoatChange2_new$Fec_SameYr <- ifelse(Summer_CoatChange2_new$ReproStatus_SameYear=="Milk",1,0)
Summer_CoatChange2_new$PrevYrRepro2 <- ifelse(as.character(Summer_CoatChange2_new$ReproStatus_PrevYear) == "Naive",
                                              "Immature", as.character(Summer_CoatChange2_new$ReproStatus_PrevYear))

### Bringing the different responses into one dataframe ###

# We want to maximise the data we keep across the responses

# Years are limited by summer moult census
table(Summer_CoatChange2_new$CensusYear) #2001-2024
table(Winter_CoatChange2_new$CensusYear) #1987-2023
table(Parturition$BirthYear) #1968-2024

Summer <- Summer_CoatChange2_new[,c("Code","CensusYear","HindAge","CoatStatus_binomial","Fec_SameYr","PrevYrRepro2")]

Winter <- Winter_CoatChange2_new[,c("Code","CensusYear","HindAge","CoatStatus_binomial","ReproStatus_SameYear")]
Winter <- subset(Winter,CensusYear>=2000) # inc. 2000 because previous winter for 2001 summer/parturition

Parturition <- subset(Parturition,BirthYear>=2001)

# Identify deer present in all three
hinds <- as.data.frame(table(c(unique(Summer$Code),
                               unique(Winter$Code),
                               unique(Parturition$MumCode)))) # count of how many of the dataframe each individual appears in
table(hinds$Freq)

hinds <- subset(hinds,Freq==3) # hinds present in all 3

Summer <- Summer[-which(Summer$Code %in% hinds$Var1==F),]
Winter <- Winter[-which(Winter$Code %in% hinds$Var1==F),]
Parturition <- Parturition[-which(Parturition$MumCode %in% hinds$Var1==F),]

# All relevant deer:year combinations for winter of same year and winter from previous year
tri_bothwinter <- data.frame(Code=c(Summer$Code,Winter$Code,Winter$Code,Parturition$MumCode),
                             Year=c(Summer$CensusYear,Winter$CensusYear,Winter$CensusYear+1,Parturition$BirthYear),
                             HindAge=c(Summer$HindAge,Winter$HindAge,Winter$HindAge+1,Parturition$MumAge))
tri_bothwinter <- distinct(tri_bothwinter) # Keep all distinct deer:year rows
tri_bothwinter$CodeYear <- paste(tri_bothwinter$Code,tri_bothwinter$Year)

# Rename columns and merge the three responses
Summer$CodeYear <- paste(Summer$Code,Summer$CensusYear)
Winter$CodeYear <- paste(Winter$Code,Winter$CensusYear)
Parturition$CodeYear <- paste(Parturition$MumCode,Parturition$BirthYear)

colnames(Summer) <- c("Code_S","SummerCensusYear","SummerHindAge","Summer_CoatCensus","Fec_SameYr","S_PrevYrRepro","CodeYear")
colnames(Winter) <- c("Code_W","WinterCensusYear","WinterHindAge","Winter_CoatCensus","SameYrRepro","CodeYear","CodePrevYear")
colnames(Parturition) <- c("CalfCode","Code_P","BirthYear","MumAge","P_PrevYrRepro","BirthDate","CodeYear")

tri_bothwinter <- merge(tri_bothwinter,Summer,by=c("CodeYear"), all=T)
tri_bothwinter <- merge(tri_bothwinter,Parturition,by=c("CodeYear"), all=T)
tri_bothwinter <- merge(tri_bothwinter,Winter,by=c("CodeYear"), all=T)
Winter$CodePrevYear <- paste(Winter$Code_W,Winter$WinterCensusYear) # actually code:year repeat but called prev to match other dataframe
tri_bothwinter$CodePrevYear <- paste(tri_bothwinter$Code,tri_bothwinter$Year-1)
tri_bothwinter <- merge(tri_bothwinter,Winter,by=c("CodePrevYear"), all=T)

#Center CensusYear and HindAge
tri_bothwinter$HindAge_cd <- tri_bothwinter$HindAge - 7
tri_bothwinter$Year_cd <- tri_bothwinter$Year - 2012

# make repro stat factor with milk as reference and missing as 'unnecessary' - 
# the unecessary wont show up in model output as all known in the relevant columns with rows of data for each response
length(which(is.na(tri_bothwinter$S_PrevYrRepro)))
length(which(is.na(tri_bothwinter$S_PrevYrRepro[which(is.na(tri_bothwinter$Summer_CoatCensus)==F)])))
tri_bothwinter$S_PrevYrRepro <- ifelse(is.na(tri_bothwinter$S_PrevYrRepro)==T,"Unnecessary",tri_bothwinter$S_PrevYrRepro)
tri_bothwinter$S_PrevYrRepro <- relevel(as.factor(tri_bothwinter$S_PrevYrRepro), ref = "Milk")

length(which(is.na(tri_bothwinter$P_PrevYrRepro)))
length(which(is.na(tri_bothwinter$P_PrevYrRepro[which(is.na(tri_bothwinter$BirthDate)==F)])))
tri_bothwinter$P_PrevYrRepro <- ifelse(is.na(tri_bothwinter$P_PrevYrRepro)==T,"Unnecessary",tri_bothwinter$P_PrevYrRepro)
tri_bothwinter$P_PrevYrRepro <- relevel(as.factor(tri_bothwinter$P_PrevYrRepro), ref = "Milk")

colnames(tri_bothwinter)[which(colnames(tri_bothwinter)=="SameYrRepro.x")] <- "SameYrRepro_sameW"
colnames(tri_bothwinter)[which(colnames(tri_bothwinter)=="Winter_CoatCensus.x")] <- "Winter_CoatCensus_sameW"
length(which(is.na(tri_bothwinter$SameYrRepro_sameW)))
length(which(is.na(tri_bothwinter$SameYrRepro_sameW[which(is.na(tri_bothwinter$Winter_CoatCensus_sameW)==F)])))
tri_bothwinter$SameYrRepro_sameW <- ifelse(is.na(tri_bothwinter$SameYrRepro_sameW)==T,"Unnecessary",tri_bothwinter$SameYrRepro_sameW)
tri_bothwinter$SameYrRepro_sameW <- relevel(as.factor(tri_bothwinter$SameYrRepro_sameW), ref = "Milk")

colnames(tri_bothwinter)[which(colnames(tri_bothwinter)=="SameYrRepro.y")] <- "SameYrRepro_prevW"
colnames(tri_bothwinter)[which(colnames(tri_bothwinter)=="Winter_CoatCensus.y")] <- "Winter_CoatCensus_prevW"
length(which(is.na(tri_bothwinter$SameYrRepro_prevW)))
length(which(is.na(tri_bothwinter$SameYrRepro_prevW[which(is.na(tri_bothwinter$Winter_CoatCensus_prevW)==F)])))
tri_bothwinter$SameYrRepro_prevW <- ifelse(is.na(tri_bothwinter$SameYrRepro_prevW)==T,"Unnecessary",tri_bothwinter$SameYrRepro_prevW)
tri_bothwinter$SameYrRepro_prevW <- relevel(as.factor(tri_bothwinter$SameYrRepro_prevW), ref = "Milk")

length(which(is.na(tri_bothwinter$Fec_SameYr)))
length(which(is.na(tri_bothwinter$Fec_SameYr[which(is.na(tri_bothwinter$Summer_CoatCensus)==F)])))
tri_bothwinter$Fec_SameYr <- ifelse(is.na(tri_bothwinter$Fec_SameYr)==T,"Unnecessary",tri_bothwinter$Fec_SameYr)
tri_bothwinter$Fec_SameYr <- as.factor(tri_bothwinter$Fec_SameYr)

# Check data types
str(tri_bothwinter)

### Looking at sample size

# Sample size for each response
length(which(is.na(tri_bothwinter$BirthDate)==F)) #2033 
length(which(is.na(tri_bothwinter$Summer_CoatCensus)==F)) #3319 
length(which(is.na(tri_bothwinter$Winter_CoatCensus_sameW)==F)) #3422 
length(which(is.na(tri_bothwinter$Winter_CoatCensus_prevW)==F)) #3422 



# Priors

k<-10000
prior2<-list(R=list(V=diag(4), nu=0, fix=2),
             G=list(G1=list(V=diag(4),nu=4,alpha.mu=c(0,0,0,0),alpha.V=diag(k,4,4)),
                    G1=list(V=diag(4),nu=4,alpha.mu=c(0,0,0,0),alpha.V=diag(k,4,4))))


# Trivariate model using same year winter coat change

MCMCMultiVar <- MCMCglmm(
  fixed = cbind(BirthDate, Summer_CoatCensus, Winter_CoatCensus_sameW, Winter_CoatCensus_prevW) ~ trait - 1 + 
    trait:Year_cd + trait:HindAge_cd + trait:I(HindAge_cd^2) + at.level(trait,1):P_PrevYrRepro +
    at.level(trait,2):S_PrevYrRepro*at.level(trait,2):Fec_SameYr + at.level(trait,3):SameYrRepro_sameW + at.level(trait,4):SameYrRepro_prevW,
  random = ~ us(trait):Year + us(trait):Code,
  rcov = ~ cors(trait):units,
  family = c("gaussian", "threshold", "threshold", "threshold"),
  prior = prior2,
  pr = TRUE,
  nitt = 450000, thin = 250, burnin = 50000, 
  data = tri_bothwinter
)

summary(MCMCMultiVar)
plot(MCMCMultiVar)


#This is supposed to save the model to the environment I think? We'll see if it works with Github
saveRDS(MCMCMultiVar, file = "MCMCMultiVar.rds")
MCMCTriM <- readRDS("MCMCMultiVar.rds")



#####################################################################################################################################

##MCMCTriM Figure##

MCMCTriM<- MCMCMultiVar

#Correlation Calculations:

#Code

SW_code_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitSummer_CoatCensus.Code"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.Code"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.Code"]))

SpW_code_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitSummer_CoatCensus.Code"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.Code"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.Code"]))

SC_code_cor <- MCMCTriM$VCV[,"traitBirthDate:traitSummer_CoatCensus.Code"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.Code"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.Code"]))

WC_code_cor <- MCMCTriM$VCV[,"traitBirthDate:traitWinter_CoatCensus_sameW.Code"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.Code"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.Code"]))

pWC_code_cor <- MCMCTriM$VCV[,"traitBirthDate:traitWinter_CoatCensus_prevW.Code"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.Code"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.Code"]))


pWW_code_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_prevW.Code"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.Code"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.Code"]))

#Year

SW_year_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitSummer_CoatCensus.Year"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.Year"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.Year"]))

SpW_year_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitSummer_CoatCensus.Year"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.Year"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.Year"]))

SC_year_cor <- MCMCTriM$VCV[,"traitBirthDate:traitSummer_CoatCensus.Year"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.Year"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.Year"]))

WC_year_cor <- MCMCTriM$VCV[,"traitBirthDate:traitWinter_CoatCensus_sameW.Year"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.Year"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.Year"]))

pWC_year_cor <- MCMCTriM$VCV[,"traitBirthDate:traitWinter_CoatCensus_prevW.Year"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.Year"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.Year"]))

pWW_year_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_prevW.Year"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.Year"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.Year"]))


#Residual (units)

SW_units_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitSummer_CoatCensus.units"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.units"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.units"]))

SpW_units_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitSummer_CoatCensus.units"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.units"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.units"]))

SC_units_cor <- MCMCTriM$VCV[,"traitBirthDate:traitSummer_CoatCensus.units"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.units"])*
     sqrt(MCMCTriM$VCV[,"traitSummer_CoatCensus:traitSummer_CoatCensus.units"]))

WC_units_cor <- MCMCTriM$VCV[,"traitBirthDate:traitWinter_CoatCensus_sameW.units"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.units"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.units"]))

pWC_units_cor <- MCMCTriM$VCV[,"traitBirthDate:traitWinter_CoatCensus_prevW.units"]/
  (sqrt(MCMCTriM$VCV[,"traitBirthDate:traitBirthDate.units"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.units"]))

pWW_units_cor <- MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_prevW.units"]/
  (sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_sameW:traitWinter_CoatCensus_sameW.units"])*
     sqrt(MCMCTriM$VCV[,"traitWinter_CoatCensus_prevW:traitWinter_CoatCensus_prevW.units"]))


#Table:

Table1 <- data.frame(
  Traits = c(
    "Winter Moult/Next Calving Birth Dates",
    "Winter Moult/Next Calving Birth Dates",
    "Winter Moult/Next Calving Birth Dates",
    "Winter Moult/Next Summer Moult",
    "Winter Moult/Next Summer Moult",
    "Winter Moult/Next Summer Moult",
    "Winter Moult/Next Winter Moult",
    "Winter Moult/Next Winter Moult",
    "Winter Moult/Next Winter Moult",
    "Calving Birth Date/Summer Moult",
    "Calving Birth Date/Summer Moult",
    "Calving Birth Date/Summer Moult",
    "Calving Birth Date/Next Winter Moult",
    "Calving Birth Date/Next Winter Moult",
    "Calving Birth Date/Next Winter Moult",
    "Summer Moult/Next Winter Moult",
    "Summer Moult/Next Winter Moult",
    "Summer Moult/Next Winter Moult"
  ),
  Subcategory = c(
    "Individual",
    "Year",
    "Residual",
    "Individual",
    "Year",
    "Residual",
    "Individual",
    "Year",
    "Residual",    
    "Individual",
    "Year",
    "Residual",
    "Individual",
    "Year",
    "Residual",
    "Individual",
    "Year",
    "Residual"
  ),
  Estimate = c(
    mean(pWC_code_cor),
    mean(pWC_year_cor),
    mean(pWC_units_cor),
    mean(SpW_code_cor),
    mean(SpW_year_cor),
    mean(SpW_units_cor),
    mean(pWW_code_cor),
    mean(pWW_year_cor),
    mean(pWW_units_cor),
    mean(SC_code_cor),
    mean(SC_year_cor),
    mean(SC_units_cor),
    mean(WC_code_cor),
    mean(WC_year_cor),
    mean(WC_units_cor),
    mean(SW_code_cor),
    mean(SW_year_cor),
    mean(SW_units_cor)
  ),
  Lower = c(
    HPDinterval(pWC_code_cor)[, "lower"],
    HPDinterval(pWC_year_cor)[, "lower"],
    HPDinterval(pWC_units_cor)[, "lower"],
    HPDinterval(SpW_code_cor)[, "lower"],
    HPDinterval(SpW_year_cor)[, "lower"],
    HPDinterval(SpW_units_cor)[, "lower"],
    HPDinterval(pWW_code_cor)[, "lower"],
    HPDinterval(pWW_year_cor)[, "lower"],
    HPDinterval(pWW_units_cor)[, "lower"],
    HPDinterval(SC_code_cor)[, "lower"],
    HPDinterval(SC_year_cor)[, "lower"],
    HPDinterval(SC_units_cor)[, "lower"],
    HPDinterval(WC_code_cor)[, "lower"],
    HPDinterval(WC_year_cor)[, "lower"],
    HPDinterval(WC_units_cor)[, "lower"],
    HPDinterval(SW_code_cor)[, "lower"],
    HPDinterval(SW_year_cor)[, "lower"],
    HPDinterval(SW_units_cor)[, "lower"]
  ),
  Upper = c(
    HPDinterval(pWC_code_cor)[, "upper"],
    HPDinterval(pWC_year_cor)[, "upper"],
    HPDinterval(pWC_units_cor)[, "upper"],
    HPDinterval(SpW_code_cor)[, "upper"],
    HPDinterval(SpW_year_cor)[, "upper"],
    HPDinterval(SpW_units_cor)[, "upper"],
    HPDinterval(pWW_code_cor)[, "upper"],
    HPDinterval(pWW_year_cor)[, "upper"],
    HPDinterval(pWW_units_cor)[, "upper"],
    HPDinterval(SC_code_cor)[, "upper"],
    HPDinterval(SC_year_cor)[, "upper"],
    HPDinterval(SC_units_cor)[, "upper"],
    HPDinterval(WC_code_cor)[, "upper"],
    HPDinterval(WC_year_cor)[, "upper"],
    HPDinterval(WC_units_cor)[, "upper"],
    HPDinterval(SW_code_cor)[, "upper"],
    HPDinterval(SW_year_cor)[, "upper"],
    HPDinterval(SW_units_cor)[, "upper"])
)


#Provides the precise correlation values 
View(Table1)


#Changes the order of the data points
Table1$Subcategory <- factor(
  Table1$Subcategory,
  levels = c( "Residual", "Individual","Year")
)

Table1$Traits <- factor(Table1$Traits,levels =   c("Winter Moult/Next Calving Birth Dates",
                                                   "Winter Moult/Next Summer Moult",
                                                   "Winter Moult/Next Winter Moult",
                                                   "Calving Birth Date/Summer Moult",
                                                   "Calving Birth Date/Next Winter Moult",
                                                   "Summer Moult/Next Winter Moult"))
#Graph:

(MultiVarFig <- ggplot(Table1, aes(x = Traits, y = Estimate, color = Subcategory)) +
    geom_pointrange(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.5), size=0.7, linewidth=1) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.3) +
    scale_color_manual(values = c("Individual" = "goldenrod", "Year" = "darkorange", "Residual" = "lightgoldenrod"),guide = guide_legend(reverse = TRUE)) +
    scale_x_discrete(limits=rev)+
    labs(
      
      x = "Trait Combinations",
      y = "Correlation (Estimate ± 95% CIs)",
      color = "Random Term"
    ) +
    coord_cartesian(ylim = c(-1, 1)) +
    coord_flip() +
    theme_classic() +
    theme(
      title = element_text(size = 16),
      axis.title.x = element_text(size = 16), 
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 11), 
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    )) #7*7



