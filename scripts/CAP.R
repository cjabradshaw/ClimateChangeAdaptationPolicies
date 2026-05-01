## Climate change adaptation policy (CAP) dataset analysis
## Corey Bradshaw & Maddy King
## Apr 2026
## Github repository: https://github.com/cjabradshaw/ClimateChangeAdaptationPolicies

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(lubridate)
library(ozmaps)
library(purrr)
library(rnaturalearth)
library(scales)
library(sf)
library(terra)
library(tidyr)
library(viridis)
library(units)

# functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# import policy data
setwd("~/Documents/GitHub/ClimateChangeAdaptationPolicies/data/")

## read comma-delimited text file
data <- read.table("CAPdatV2.csv", header=TRUE, sep=",", dec=".", strip.white=TRUE, quote="\"")
head(data)

## recode legbody to Australian State/Territory abbreviations
data$legbodyCODE <- recode(data$legbody, "New South Wales"="NSW",
                       "Victoria"="VIC",
                       "Queensland"="QLD",
                       "South Australia"="SA",
                       "Western Australia"="WA",
                       "Tasmania"="TAS",
                       "Northern Territory"="NT",
                       "Australian Capital Territory"="ACT")
table(data$legbodyCODE)


##################################
## summary by political ideology

## which columns start with 'Political_Party'
politics_colnames <- colnames(data)[grepl("^Party", colnames(data))]
politics_colnames

## create left-right categories
data$pol_LEFT <- ifelse(data$Party_Labor == 1 | data$Party_Labor_Greens == 1, 1, 0)
data$pol_RIGHT <- ifelse(data$Party_Liberal == 1 | data$Party_Nationals == 1 |
                           data$Party_Liberal_National == 1, 1, 0)

## plot proportion of records with left or right political affiliation by state
sum_pol <- sum(data$pol_LEFT, na.rm=TRUE) + sum(data$pol_RIGHT, na.rm=TRUE)
state_num <- data %>%
  group_by(legbodyCODE) %>%
  summarise(n_LEFT = sum(pol_LEFT, na.rm=TRUE),
            n_RIGHT = sum(pol_RIGHT, na.rm=TRUE),
            p_LEFT = n_LEFT / (n_LEFT + n_RIGHT),
            p_RIGHT = n_RIGHT / (n_LEFT + n_RIGHT))

## remove untitled category
state_num <- subset(state_num, legbodyCODE != "")

## pivot proportions
df_prop_long <- state_num %>%
  pivot_longer(
    cols = c(p_LEFT, p_RIGHT),
    names_to = "ideology",
    values_to = "p"
  ) %>%
  mutate(
    ideology = recode(ideology,
                  p_LEFT = "LEFT",
                  p_RIGHT = "RIGHT")
  )

## plot proportion LEFT/RIGHT by state in ggplot
ggplot(df_prop_long, aes(x = legbodyCODE, y = p, fill = ideology)) +
  geom_col() +
  scale_fill_manual(values = c("LEFT" = "red", "RIGHT" = "blue")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "State / Territory",
    y = "proportion",
    fill = "ideology"
  ) +
  theme_minimal()

## summarise across all states/territories
total_pol <- data %>%
  summarise(n_LEFT = sum(pol_LEFT, na.rm=TRUE),
            n_RIGHT = sum(pol_RIGHT, na.rm=TRUE),
            p_LEFT = n_LEFT / (n_LEFT + n_RIGHT),
            p_RIGHT = n_RIGHT / (n_LEFT + n_RIGHT))

## pivot proportions
df_total_long <- total_pol %>%
  pivot_longer(
    cols = c(p_LEFT, p_RIGHT),
    names_to = "ideology",
    values_to = "p"
  ) %>%
  mutate(
    ideology = recode(ideology,
                      p_LEFT = "LEFT",
                      p_RIGHT = "RIGHT")
  )

## plot total proportion LEFT/RIGHT in ggplot
ggplot(df_total_long, aes(x = ideology, y = p, fill = ideology)) +
  geom_col() +
  scale_fill_manual(values = c("LEFT" = "red", "RIGHT" = "blue")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.7)) +
  labs(
    x = "ideology",
    y = "proportion",
    fill = "ideology"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


## summary by IPCC category
data$IPCC_BehavCult <- data$IPCCA_Behavioural_and_Cultural
data$IPCC_Knowl <- data$IPCCA_Knowledge
data$IPCC_InfrTech <- data$IPCCA_Infrastructural_and_Technological
data$IPCC_NatBased <- data$IPCCA_Nature_Based
data$IPCC_Inst <- data$IPCCA_Institutional
head(data)

## plot proportion of records with each IPCC category by state
ippc_sum <- data %>%
  group_by(legbodyCODE) %>%
  summarise(n_BehavCult = sum(IPCC_BehavCult, na.rm=TRUE),
            n_Knowl = sum(IPCC_Knowl, na.rm=TRUE),
            n_InfrTech = sum(IPCC_InfrTech, na.rm=TRUE),
            n_NatBased = sum(IPCC_NatBased, na.rm=TRUE),
            n_Inst = sum(IPCC_Inst, na.rm=TRUE),
            p_BehavCult = n_BehavCult / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_Knowl = n_Knowl / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_InfrTech = n_InfrTech / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_NatBased = n_NatBased / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_Inst = n_Inst / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst))
ippc_sum

## remove untitled category
ippc_sum <- subset(ippc_sum, legbodyCODE != "")
ippc_sum

## pivot counts
df_ippc_long <- ippc_sum %>%
  pivot_longer(
    cols = c(p_BehavCult, p_Knowl, p_InfrTech, p_NatBased, p_Inst),
    names_to = "IPCC_category",
    values_to = "p"
  ) %>%
  mutate(
    IPCC_category = recode(IPCC_category,
                          p_BehavCult = "Behavioural and Cultural",
                          p_Knowl = "Knowledge",
                          p_InfrTech = "Infrastructural and Technological",
                          p_NatBased = "Nature-Based",
                          p_Inst = "Institutional")
  )
df_ippc_long

## plot
ggplot(df_ippc_long, aes(x = legbodyCODE, y = p, fill = IPCC_category)) +
  geom_col() +
  scale_fill_manual(values = c("Behavioural and Cultural" = "darkblue",
                               "Knowledge" = "green",
                               "Infrastructural and Technological" = "lightpink",
                               "Nature-Based" = "darkgrey",
                               "Institutional" = "red")) +
  labs(
    x = "State / Territory",
    y = "proportion",
    fill = "IPCC category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "top")

## summarise across all states/territories
total_ippc <- data %>%
  summarise(n_BehavCult = sum(IPCC_BehavCult, na.rm=TRUE),
            n_Knowl = sum(IPCC_Knowl, na.rm=TRUE),
            n_InfrTech = sum(IPCC_InfrTech, na.rm=TRUE),
            n_NatBased = sum(IPCC_NatBased, na.rm=TRUE),
            n_Inst = sum(IPCC_Inst, na.rm=TRUE),
            p_BehavCult = n_BehavCult / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_Knowl = n_Knowl / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_InfrTech = n_InfrTech / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_NatBased = n_NatBased / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst),
            p_Inst = n_Inst / sum(n_BehavCult, n_Knowl, n_InfrTech, n_NatBased, n_Inst))
total_ippc

## order from highest to lowest proportion
total_ippc_long <- total_ippc %>%
  pivot_longer(
    cols = c(p_BehavCult, p_Knowl, p_InfrTech, p_NatBased, p_Inst),
    names_to = "IPCC_category",
    values_to = "p"
  ) %>%
  mutate(
    IPCC_category = recode(IPCC_category,
                          p_BehavCult = "Behavioural and Cultural",
                          p_Knowl = "Knowledge",
                          p_InfrTech = "Infrastructural and Technological",
                          p_NatBased = "Nature-Based",
                          p_Inst = "Institutional")
  ) %>%
  arrange(desc(p))

## plot
ggplot(total_ippc_long, aes(x = reorder(IPCC_category, -p), y = p, fill = IPCC_category)) +
  geom_col() +
  scale_fill_manual(values = c("Behavioural and Cultural" = "darkblue",
                               "Knowledge" = "green",
                               "Infrastructural and Technological" = "lightpink",
                               "Nature-Based" = "darkgrey",
                               "Institutional" = "red")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "IPCC category",
    y = "proportion",
    fill = "IPCC category"
  ) +
  theme_minimal() +
  theme(legend.position = "none")



##################################
## summary by hazard

## columns about hazard
hazard_colnames <- colnames(data)[grepl("^Hazard_", colnames(data))]
hazard_colnames

data$haz_fire <- data$Hazard_Fire
data$haz_flood <- ifelse(data$Hazard_Coastal_erosion_flooding_change == 1 | data$Hazard_Riverine_Flooding == 1, 1, 0)
data$haz_tmpprcphum <- ifelse(data$Hazard_Temperature_Extremes == 1 | data$Hazard_Changes_Precipitation == 1 |
                              data$Hazard_Extreme_wind == 1, 1, 0)
data$haz_wind <- data$Hazard_Extreme_wind
data$haz_drought <- data$Hazard_Drought
data$haz_gen <- data$Hazard_General_Climate_Change
head(data)

## plot proportions by each hazard by state
haz_sum <- data %>%
  group_by(legbodyCODE) %>%
  summarise(n_fire = sum(haz_fire, na.rm=TRUE),
            n_flood = sum(haz_flood, na.rm=TRUE),
            n_tmpprcphum = sum(haz_tmpprcphum, na.rm=TRUE),
            n_wind = sum(haz_wind, na.rm=TRUE),
            n_drought = sum(haz_drought, na.rm=TRUE),
            n_gen = sum(haz_gen, na.rm=TRUE),
            p_fire = n_fire / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_flood = n_flood / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_tmpprcphum = n_tmpprcphum / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_wind = n_wind / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_drought = n_drought / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_gen = n_gen / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen))
haz_sum

## remove untitled category
haz_sum <- subset(haz_sum, legbodyCODE != "")
haz_sum

## pivot counts
df_haz_long <- haz_sum %>%
  pivot_longer(
    cols = c(p_fire, p_flood, p_tmpprcphum, p_wind, p_drought, p_gen),
    names_to = "hazard",
    values_to = "p"
  ) %>%
  mutate(
    hazard = recode(hazard,
                    p_fire = "fire",
                    p_flood = "flood",
                    p_tmpprcphum = "temperature/precipitation/humidity",
                    p_wind = "wind",
                    p_drought = "drought",
                    p_gen = "general")
  )
df_haz_long

## plot
ggplot(df_haz_long, aes(x = legbodyCODE, y = p, fill = hazard)) +
  geom_col() +
  scale_fill_manual(values = c("fire" = "red",
                               "flood" = "green",    
                               "temperature/precipitation/humidity" = "blue",
                               "wind" = "purple",
                               "drought" = "brown",    
                               "general" = "darkgrey")) +
  labs(
    x = "State / Territory",
    y = "proportion",
    fill = "hazard"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "top")


## summarise across all states/territories
total_haz <- data %>%
  summarise(n_fire = sum(haz_fire, na.rm=TRUE),
            n_flood = sum(haz_flood, na.rm=TRUE),
            n_tmpprcphum = sum(haz_tmpprcphum, na.rm=TRUE),
            n_wind = sum(haz_wind, na.rm=TRUE),
            n_drought = sum(haz_drought, na.rm=TRUE),
            n_gen = sum(haz_gen, na.rm=TRUE),
            p_fire = n_fire / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_flood = n_flood / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_tmpprcphum = n_tmpprcphum / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_wind = n_wind / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_drought = n_drought / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen),
            p_gen = n_gen / sum(n_fire, n_flood, n_tmpprcphum, n_wind, n_drought, n_gen))
total_haz

## order from highest to lowest proportion
total_haz_long <- total_haz %>%
  pivot_longer(
    cols = c(p_fire, p_flood, p_tmpprcphum, p_wind, p_drought, p_gen),
    names_to = "hazard",
    values_to = "p"
  ) %>%
  mutate(
    hazard = recode(hazard,
                    p_fire = "fire",
                    p_flood = "flood",
                    p_tmpprcphum = "temperature/precipitation/humidity",
                    p_wind = "wind",
                    p_drought = "drought",
                    p_gen = "general")
  ) %>%
  arrange(desc(p))
total_haz_long

## plot
ggplot(total_haz_long, aes(x = reorder(hazard, -p), y = p, fill = hazard)) +
  geom_col() +
  scale_fill_manual(values = c("fire" = "red",
                               "flood" = "green",    
                               "temperature/precipitation/humidity" = "blue",
                               "wind" = "purple",
                               "drought" = "brown",    
                               "general" = "darkgrey")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "hazard",
    y = "proportion",
    fill = "hazard"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


###############################################
## summary by climate risk assessment category
data$risk_infrastr <- data$Climate_Risk_Assessment_Infrastructure_Built_Environment
data$risk_health <- data$Climate_Risk_Assessment_Health_Social_Support
data$risk_nature <- data$Climate_Risk_Assessment_Natural_Environment
data$risk_primprod <- data$Climate_Risk_Assessment_Primary_Industries_Food_System
data$risk_regions <- data$Climate_Risk_Assessment_Regional_Remote_Communities
data$risk_economy <- data$Climate_Risk_Assessment_Economy_Trade_Financial_System
data$risk_defence <- data$Climate_Risk_Assessment_Defence_National_Security
data$risk_multi <- data$Climate_Risk_Assessment_Multi_System
colnames(data)

## plot proportions by each risk category by state, if a category is missing in 1 state, ignore warning and set to zero
risk_sum <- data %>%
  group_by(legbodyCODE) %>%
  summarise(n_infrastr = sum(risk_infrastr, na.rm=TRUE),
            n_health = sum(risk_health, na.rm=TRUE),
            n_nature = sum(risk_nature, na.rm=TRUE),
            n_primprod = sum(risk_primprod, na.rm=TRUE),
            n_regions = sum(risk_regions, na.rm=TRUE),
            n_economy = sum(risk_economy, na.rm=TRUE),
            n_defence = sum(risk_defence, na.rm=TRUE),
            n_multi = sum(risk_multi, na.rm=TRUE),
            p_infrastr = n_infrastr / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_health = n_health / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_nature = n_nature / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_primprod = n_primprod / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_regions = n_regions / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_economy = n_economy / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_defence = ifelse(sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi) == 0, 0,
                               ifelse(is.na(n_defence), 0 ,n_defence / sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi))),
            p_multi = ifelse(sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi) == 0, 0,
                             ifelse(is.na(n_multi), 0, n_multi / sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi))))
risk_sum

## pivot counts
df_risk_long <- risk_sum %>%
  pivot_longer(
    cols = c(p_infrastr, p_health, p_nature, p_primprod, p_regions, p_economy, p_defence, p_multi),
    names_to = "risk_category",
    values_to = "p"
  ) %>%
  mutate(
    risk_category = recode(risk_category,
                           p_infrastr = "infrastructure/built environment",
                           p_health = "health/social support",
                           p_nature = "natural environment",
                           p_primprod = "primary industries/food system",
                           p_regions = "regional/remote communities",
                           p_economy = "economy/trade/financial system",
                           p_defence = "defence/national security",
                           p_multi = "multi-system")
  )
df_risk_long
                           
## plot
ggplot(df_risk_long, aes(x = legbodyCODE, y = p, fill = risk_category)) +
  geom_col() +
  scale_fill_manual(values = c("infrastructure/built environment" = "darkgrey",
                               "health/social support" = "red",                         
                               "natural environment" = "darkgreen",
                               "primary industries/food system" = "lightgreen",
                               "regional/remote communities" = "purple",
                               "economy/trade/financial system" = "brown",
                               "defence/national security" = "blue",
                               "multi-system" = "black")) +
  labs(
    x = "State / Territory",
    y = "proportion",
    fill = "climate risk assessment category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "top")

## summarise across all states/territories
total_risk <- data %>%
  summarise(n_infrastr = sum(risk_infrastr, na.rm=TRUE),
            n_health = sum(risk_health, na.rm=TRUE),
            n_nature = sum(risk_nature, na.rm=TRUE),
            n_primprod = sum(risk_primprod, na.rm=TRUE),
            n_regions = sum(risk_regions, na.rm=TRUE),
            n_economy = sum(risk_economy, na.rm=TRUE),
            n_defence = sum(risk_defence, na.rm=TRUE),
            n_multi = sum(risk_multi, na.rm=TRUE),
            p_infrastr = n_infrastr / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_health = n_health / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_nature = n_nature / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_primprod = n_primprod / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_regions = n_regions / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_economy = n_economy / sum(n_infrastr, n_health, n_nature, n_primprod, n_regions, n_economy, n_defence, n_multi),
            p_defence = ifelse(sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi) == 0, 0,
                               ifelse(is.na(n_defence), 0 ,n_defence / sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi))),
            p_multi = ifelse(sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi) == 0, 0,
                             ifelse(is.na(n_multi), 0, n_multi / sum(n_infrastr,n_health,n_nature,n_primprod,n_regions,n_economy,n_defence,n_multi))))
total_risk
                                          
## plot
df_total_risk_long <- total_risk %>%
  pivot_longer(
    cols = c(p_infrastr, p_health, p_nature, p_primprod, p_regions, p_economy, p_defence, p_multi),
    names_to = "risk_category",
    values_to = "p"
  ) %>%
  mutate(
    risk_category = recode(risk_category,
                           p_infrastr = "infrastructure/built environment",
                           p_health = "health/social support",
                           p_nature = "natural environment",
                           p_primprod = "primary industries/food system",
                           p_regions = "regional/remote communities",
                           p_economy = "economy/trade/financial system",
                           p_defence = "defence/national security",
                           p_multi = "multi-system")
  )
df_total_risk_long
                               
ggplot(df_total_risk_long, aes(x = risk_category, y = p, fill = risk_category)) +
  geom_col() +
  scale_fill_manual(values = c("infrastructure/built environment" = "darkgrey",
                               "health/social support" = "red",                         
                               "natural environment" = "darkgreen",
                               "primary industries/food system" = "lightgreen",
                               "regional/remote communities" = "purple",
                               "economy/trade/financial system" = "brown",
                               "defence/national security" = "blue",
                               "multi-system" = "black")) +
  labs(
    x = "",
    y = "proportion",
    fill = "climate risk assessment category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")                               
                               
## order from highest to lowest proportion
df_total_risk_long$risk_category <- factor(df_total_risk_long$risk_category, levels = df_total_risk_long$risk_category[order(df_total_risk_long$p, decreasing = TRUE)])

# replot
ggplot(df_total_risk_long, aes(x = risk_category, y = p, fill = risk_category)) +
  geom_col() +
  scale_fill_manual(values = c("infrastructure/built environment" = "darkgrey",
                               "health/social support" = "red",                         
                               "natural environment" = "darkgreen",
                               "primary industries/food system" = "lightgreen",
                               "regional/remote communities" = "purple",
                               "economy/trade/financial system" = "brown",
                               "defence/national security" = "blue",
                               "multi-system" = "black")) +
  labs(
    x = "",
    y = "proportion",
    fill = "climate risk assessment category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
                  

## summarise by purpose, amalgamating categories into 4 broad groups:
## 1. 'governance' = governance instruments + institutional change
## 2. 'research' = information gathering + research + knowledge building
## 3. 'coordination' = stakeholder engagement + coordination
## 4. 'intervention' = intervention
data$governance <- ifelse(data$Purpose_Governance_Instruments == 1 | 
                          data$Purpose_Institutional_Change == 1, 1, 0)
data$research <- ifelse(data$Purpose_Information_Gathering == 1 | 
                        data$Purpose_Research == 1 | 
                        data$Purpose_Knowledge_Building == 1, 1, 0)
data$coordination <- ifelse(data$Purpose_Stakeholder_Engagement == 1 | 
                            data$Purpose_Coordination == 1, 1, 0)
data$intervention <- data$Purpose_Intervention
                            
## plot proportions by each purpose category by state, if a category is missing in 1 state, ignore warning and set to zero
purpose_sum <- data %>%
  group_by(legbodyCODE) %>%
  summarise(n_governance = sum(governance, na.rm=TRUE),
            n_research = sum(research, na.rm=TRUE),
            n_coordination = sum(coordination, na.rm=TRUE),
            n_intervention = sum(intervention, na.rm=TRUE),
            p_governance = n_governance / sum(n_governance, n_research, n_coordination, n_intervention),
            p_research = n_research / sum(n_governance, n_research, n_coordination, n_intervention),
            p_coordination = n_coordination / sum(n_governance, n_research, n_coordination, n_intervention),
            p_intervention = ifelse(sum(n_governance,n_research,n_coordination,n_intervention) == 0, 0,
                                    ifelse(is.na(n_intervention), 0 ,n_intervention / sum(n_governance,n_research,n_coordination,n_intervention))))
purpose_sum

## pivot counts
df_purpose_long <- purpose_sum %>%
  pivot_longer(
    cols = c(p_governance, p_research, p_coordination, p_intervention),
    names_to = "purpose_category",
    values_to = "p"
  ) %>%
  mutate(
    purpose_category = recode(purpose_category,
                           p_governance = "governance",
                            p_research = "research",
                            p_coordination = "coordination",
                            p_intervention = "intervention")
  )
df_purpose_long

# plot
ggplot(df_purpose_long, aes(x = legbodyCODE, y = p, fill = purpose_category)) +
  geom_col() +
  scale_fill_manual(values = c("governance" = "darkgrey",
                               "research" = "red",                         
                               "coordination" = "blue",
                               "intervention" = "green")) +
  labs(
    x = "state/territory",
    y = "proportion",
    fill = "purpose category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "top")

## summarise across all states/territories
total_purpose <- data %>%
  summarise(n_governance = sum(governance, na.rm=TRUE),
            n_research = sum(research, na.rm=TRUE),
            n_coordination = sum(coordination, na.rm=TRUE),
            n_intervention = sum(intervention, na.rm=TRUE),
            p_governance = n_governance / sum(n_governance, n_research, n_coordination, n_intervention),
            p_research = n_research / sum(n_governance, n_research, n_coordination, n_intervention),
            p_coordination = n_coordination / sum(n_governance, n_research, n_coordination, n_intervention),
            p_intervention = ifelse(sum(n_governance,n_research,n_coordination,n_intervention) == 0, 0,
                                    ifelse(is.na(n_intervention), 0 ,n_intervention / sum(n_governance,n_research,n_coordination,n_intervention))))
total_purpose

# plot
df_total_purpose_long <- total_purpose %>%
  pivot_longer(
    cols = c(p_governance, p_research, p_coordination, p_intervention),
    names_to = "purpose_category",
    values_to = "p"
  ) %>%
  mutate(
    purpose_category = recode(purpose_category,
                              p_governance = "governance",
                              p_research = "research",
                              p_coordination = "coordination",
                              p_intervention = "intervention")
  )
df_total_purpose_long

ggplot(df_total_purpose_long, aes(x = purpose_category, y = p, fill = purpose_category)) +
  geom_col() +
  scale_fill_manual(values = c("governance" = "darkgrey",
                               "research" = "red",                         
                               "coordination" = "blue",
                               "intervention" = "green")) +
  labs(
    x = "",
    y = "proportion",
    fill = "purpose category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

## reorder from highest to lowest proportion
df_total_purpose_long$purpose_category <- factor(df_total_purpose_long$purpose_category,
                                                levels = df_total_purpose_long$purpose_category[order(df_total_purpose_long$p, decreasing = TRUE)])
# replot
ggplot(df_total_purpose_long, aes(x = purpose_category, y = p, fill = purpose_category)) +
  geom_col() +
  scale_fill_manual(values = c("governance" = "darkgrey",
                               "research" = "red",
                               "coordination" = "blue",
                               "intervention" = "green")) +
  labs(
    x = "",
    y = "proportion",
    fill = "purpose category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")



#######################
## occurrence by year
## first, create dataset with expanding range between Start_Year and End_Year  by year
data_years <- data %>%
  filter(!is.na(Start_Year)) %>%
  mutate(End_Year = ifelse(is.na(End_Year), Start_Year, End_Year)) %>%
  rowwise() %>%
  mutate(year = list(seq(Start_Year, End_Year))) %>%
  unnest(cols = c(year))
head(data_years)
colnames(data_years)
data_years$year

## count number of records per year
year_sum <- data_years %>%
  group_by(year) %>%
  summarise(n = n())
tail(year_sum)

## identify peak year
peak_year <- year_sum$year[which.max(year_sum$n)]

## plot number of records per year, add vertical dashed line at peak and current (2006) years
recyr_plot <- ggplot(year_sum, aes(x = year, y = n)) +
  geom_line() +
  geom_vline(xintercept = peak_year, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 2026, linetype = "dashed", color = "red") +
  labs(x = "year", y = "number of records") +
  theme_minimal()

## plot cumulative records by year
cumrecyr_plot <- ggplot(year_sum, aes(x = year, y = cumsum(n))) +
  geom_line() +
  #scale_y_log10() +
  geom_vline(xintercept = peak_year, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 2026, linetype = "dashed", color = "red") +
  labs(x = "year", y = "cumulative number of records") +
  theme_minimal()

## plot together
ggarrange(recyr_plot, cumrecyr_plot, ncol = 1, nrow = 2)


##############################
## number of records by state
state_sum <- data %>%
  group_by(legbodyCODE) %>%
  summarise(n = n())
state_sum

## average population size by state
state_pop <- data %>%
  group_by(legbodyCODE) %>%
  summarise(pop = mean(Population_Abundance, na.rm=TRUE))

## number of records / population size by state (millions of people)
state_sum_pop <- merge(state_sum, state_pop, by = "legbodyCODE")
state_sum_pop$nrecXpop <- state_sum_pop$n / state_sum_pop$pop * 1e6
state_sum_pop

## import median income data
## https://www.abs.gov.au/statistics/labour/earnings-and-working-conditions/personal-income-australia/2021-22#data-downloads
income_data <- read.table("earningMedSA4.csv", header=TRUE, sep=",", dec=".", strip.white=TRUE, quote="\"")
head(income_data)
state_income <- income_data[which(is.na(income_data$SA4CODE)==T),]
state_income <- na.omit(state_income[,c("region", "medEarn22")])
state_income                            

## merge with state_sum_pop
state_sum_pop_income <- merge(state_sum_pop, state_income, by.x = "legbodyCODE", by.y = "region")
state_sum_pop_income

## plot nrecXpop by medEarn22
ggplot(state_sum_pop_income, aes(x = medEarn22, y = nrecXpop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "median income 2021-2022", y = "number of records / million people") +
  theme_minimal() +
  geom_text_repel(aes(label = legbodyCODE), size = 4, box.padding = 0.5, point.padding = 0.5)

## plot nrec by medEarn22
ggplot(state_sum_pop_income, aes(x = medEarn22, y = n)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "median income 2021-2022", y = "number of records") +
  theme_minimal() +
  geom_text_repel(aes(label = legbodyCODE), size = 4, box.padding = 0.5, point.padding = 0.5)

## SA4 income
sa4_income <- income_data[which(is.na(income_data$SA4CODE)==F),]

## SA4 population estimates
## https://dataexplorer.abs.gov.au/vis?tm=ABS_ANNUAL_ERP_ASGS2021&pg=0&snb=1&df%5Bds%5D=PEOPLE_TOPICS&df%5Bid%5D=ABS_ANNUAL_ERP_ASGS2021&df%5Bag%5D=ABS&df%5Bvs%5D=1.2.0&dq=.GCCSA..A&pd=2015,&to%5BTIME_PERIOD%5D=false
sa4pop <- read.csv("popSA4.csv", header=TRUE, sep=",", dec=".", strip.white=TRUE, quote="\"")
head(sa4pop)
sa4pop$pop_mean <- apply(sa4pop[,2:ncol(sa4pop)], 1, mean, na.rm=TRUE)
head(sa4pop)

## merge sa4pop with sa4income
sa4_income_pop <- merge(sa4_income, sa4pop[,c("SA4CODE", "pop_mean")], by.x = "SA4CODE", by.y = "SA4CODE")
head(sa4_income_pop)

## tabulate number of records by SA4 region in data
data$SA4CODE_2021_2022

## for records with > 1 SA4 region, split into multiple rows (one per SA4 region)
data_sa4 <- data %>%
  filter(!is.na(SA4CODE_2021_2022)) %>%
  rowwise() %>%
  mutate(SA4CODE_2021_2022 = strsplit(as.character(SA4CODE_2021_2022), ";")) %>%
  unnest(SA4CODE_2021_2022)
head(data_sa4)
                                                   
colnames(data_sa4)
data_sa4$SA4CODE_2021_2022

## create new numeric SA4CODE variable from the numbers in the character SA4CODE_2021_2022 variable
data_sa4$SA4CODE <- as.numeric(data_sa4$SA4CODE_2021_2022)
data_sa4$SA4CODE

## tabulate number of records by SA4CODE
sa4_sum <- data_sa4 %>%
  group_by(SA4CODE) %>%
  summarise(n = n())
sa4_sum

## merge sa4_income_pop with sa4_sum by SA4CODE & retain only relevant columns
sa4_income_pop_sum <- merge(sa4_income_pop, sa4_sum, by.x = "SA4CODE", by.y = "SA4CODE", all.x = TRUE)
sa4_income_pop_sum

## create nrecords by population (per 1 million people) variable)
sa4_income_pop_sum$nrecXpop <- sa4_income_pop_sum$n / (sa4_income_pop_sum$pop_mean / 1000000)
head(sa4_income_pop_sum)

## plot nrecXpop by medEarn22
ggplot(sa4_income_pop_sum, aes(x = medEarn22, y = nrecXpop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "median income 2021-2022", y = "number of records / million people") +
  theme_minimal() # +
  #geom_text_repel(aes(label = SA4CODE), size = 2.5, box.padding = 0.5, point.padding = 0.5)

## plot nrec by medEarn22
ggplot(sa4_income_pop_sum, aes(x = medEarn22, y = n)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "median income 2021-2022", y = "number of records") +
  theme_minimal() #+
  #geom_text_repel(aes(label = SA4CODE), size = 2.5, box.padding = 0.5, point.padding = 0.5)




####################################################
## state/territory temperature anomaly time series
## SILO open data https://www.longpaddock.qld.gov.au/silo/
setwd("~/Documents/GitHub/ClimateChangeAdaptationPolicies/")

## download (in Terminal) SILO data directly to disk (creates relevant subdirectories)
## /silo/... included in .gitignore to avoid syncing huge data files with Github
## warning: each year's .nc file is 419.3 MB, so full download (max & min temp) = ~ 113.2 GB
## code also available as 'SILOncDownload.txt' in the Github repository (https://github.com/cjabradshaw/ClimateChangeAdaptationPolicies)
  ## for y in $(seq 1891 2025); do
     ## aws s3 cp \
     ## s3://silo-open-data/Official/annual/max_temp/${y}.max_temp.nc \
     ## data/silo/max_temp/ \
     ## --no-sign-request
  ## done

  ## and

  ## for y in $(seq 1891 2025); do
     ## aws s3 cp \
     ## s3://silo-open-data/Official/annual/min_temp/${y}.min_temp.nc \
     ## data/silo/min_temp/ \
     ## --no-sign-request
  ## done

## retrieve data from downloaded files
tmax <- rast(list.files("data/silo/max_temp", full.names = TRUE))
tmin <- rast(list.files("data/silo/min_temp", full.names = TRUE))

## calculate daily mean
tmean_daily <- (tmax + tmin) / 2

## calculate annual means
dates <- time(tmean_daily)

tmean_annual <- tapp(
  tmean_daily,
  index = lubridate::year(dates),
  fun   = mean,
  na.rm = TRUE
)
names(tmean_annual)

## define crs
albers_wkt <- "PROJCRS[\"GDA94 / Australian Albers\",
  BASEGEOGCRS[\"GDA94\",
    DATUM[\"Geocentric Datum of Australia 1994\",
      ELLIPSOID[\"GRS 1980\",6378137,298.257222101]],
    PRIMEM[\"Greenwich\",0],
    CS[ellipsoidal,2],
    AXIS[\"latitude\",north],
    AXIS[\"longitude\",east],
    UNIT[\"degree\",0.0174532925199433]],
  CONVERSION[\"Australian Albers\",
    METHOD[\"Albers Equal Area\"],
    PARAMETER[\"Latitude of false origin\",0],
    PARAMETER[\"Longitude of false origin\",132],
    PARAMETER[\"Latitude of 1st standard parallel\",-18],
    PARAMETER[\"Latitude of 2nd standard parallel\",-36],
    PARAMETER[\"Easting at false origin\",0],
    PARAMETER[\"Northing at false origin\",0]],
  CS[Cartesian,2],
  AXIS[\"easting\",east],
  AXIS[\"northing\",north],
  UNIT[\"metre\",1]]"

Sys.setenv(PROJ_LIB = "/usr/share/proj")

## assign crs
crs(tmean_annual) <- albers_wkt

## state/territory boundaries
names(ozmaps::abs_ste)
states_sf <- ozmaps::abs_ste

## terra format
states_v <- states_sf %>%
  select(
    name = all_of("NAME"),
    geometry
  ) %>%
  st_transform(3577) %>%     # Australian Albers
  vect() %>%
  makeValid()

## check
class(states_v)
nrow(states_v)
table(states_v$name)
unique(states_v$name)

# ensure CRS match
crs(states_v)
crs(tmean_annual) <- "GEOGCRS[\"WGS 84\",
  DATUM[\"World Geodetic System 1984\",
    ELLIPSOID[\"WGS 84\",6378137,298.257223563]],
  CS[ellipsoidal,2],
  AXIS[\"latitude\",north],
  AXIS[\"longitude\",east],
  UNIT[\"degree\",0.0174532925199433]]"
tmean_annual_aea <- terra::project(tmean_annual, albers_wkt)

# repair geometries
states_v <- terra::makeValid(states_v)

## state mean temps
suppressWarnings({
  state_means <- terra::extract(
    tmean_annual_aea,
    states_v,
    fun     = mean,
    weights = TRUE,
    na.rm   = TRUE,
    ID      = TRUE
  )
})
names(state_means)

## long format
state_df <- state_means |>
  dplyr::rename(poly_id = ID) |>
  dplyr::left_join(
    data.frame(
      poly_id = seq_len(nrow(states_v)),
      state   = states_v$name
    ),
    by = "poly_id"
  ) |>
  tidyr::pivot_longer(
    cols = matches("^X\\d{4}$"),
    names_to  = "year",
    values_to = "tmean"
  ) |>
  dplyr::mutate(
    year = as.integer(sub("^X", "", year))
  ) |>
  dplyr::select(state, year, tmean)

state_df

state_df$state <- recode(state_df$state, "New South Wales"="NSW",
                           "Victoria"="VIC",
                           "Queensland"="QLD",
                           "South Australia"="SA",
                           "Western Australia"="WA",
                           "Tasmania"="TAS",
                           "Northern Territory"="NT",
                           "Australian Capital Territory"="ACT",
                           "Other Territories"="other")
table(state_df$state)
state_df

## plot time series of annual mean temperature by state
ggplot(state_df, aes(x = year, y = tmean, color = state)) +
  geom_line() +
  labs(x = "year", y = "annual mean temperature (°C)", color = "state/territory") +
  theme_minimal() +
  theme(legend.position = "top")

## baseline for temperature anomalies
baseline_state <- state_df %>%
  filter(year >= 1900, year <= 1950) %>%
  group_by(state) %>%
  summarise(
    baseline = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

## anomalies
anomalies_state <- state_df %>%
  left_join(baseline_state, by = "state") %>%
  mutate(anomaly = tmean - baseline)

## remove 'other' category from 'anomalies_state'
anomalies_state <- subset(anomalies_state, state != "other")

## plot time series of temperature anomalies by state
ggplot(anomalies_state, aes(x = year, y = anomaly, color = state)) +
  geom_line() +
  labs(x = "year", y = "temperature anomaly (°C)", color = "state/territory") +
  theme_minimal() +
  theme(legend.position = "top")

## create map of temperature anomalies by state for most recent year (2025)
# filter anomalies for 2025
state_2025 <- anomalies_state %>%
  filter(year == 2025) %>%
  select(state, anomaly)

# join to state geometries
# NOTE: adjust the join key if your state name column differs
names(states_sf)
unique(state_2025$state)
state_name_col <- grep("NAME", names(states_sf), value = TRUE)[1]
unique(states_sf[[state_name_col]])

## create lookup table for join
state_lookup <- tibble::tibble(
  state = c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"),
  name  = c(
    "New South Wales",
    "Victoria",
    "Queensland",
    "South Australia",
    "Western Australia",
    "Tasmania",
    "Northern Territory",
    "Australian Capital Territory"
  )
)

state_2025_fixed <- state_2025 %>%
  left_join(state_lookup, by = "state")

map_sf <- states_sf %>%
  left_join(
    state_2025_fixed,
    by = setNames("name", state_name_col)
  )

summary(map_sf$anomaly)
str(map_sf$anomaly)
map_sf$anomaly <- as.numeric(map_sf$anomaly)

ggplot(map_sf) +
  geom_sf(aes(fill = anomaly), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal()

map_sf_ae <- st_transform(map_sf, 3577)
state_2025_limits <- c(min(map_sf_ae$anomaly, na.rm = TRUE),
                          max(map_sf_ae$anomaly, na.rm = TRUE))

ggplot(map_sf_ae) +
  geom_sf(aes(fill = anomaly), colour = "grey30", linewidth = 0.2) +
    scale_fill_gradient2(
    name     = "2025 anomaly (°C)",
    low      = "#2C7BB6",
    mid      = "white",
    high     = "#D7191C",
    midpoint = 0,
    limits   = state_2025_limits,
    oob      = squish
  ) +
  coord_sf(crs = 3577) +
  labs(
    caption  = "source: SILO (BoM-derived)"
  ) +
  theme_minimal()


## calculate anomalies by SA4 region
## download shapefile from ABS: https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/SA4_2021_AUST_SHP_GDA2020.zip

## load SA4 boundaries
setwd("~/Documents/GitHub/ClimateChangeAdaptationPolicies/data/")

## sf object
sa4_sf <- st_read(
  "shapefiles/SA4_2021_AUST_SHP_GDA2020/SA4_2021_AUST_GDA2020.shp",
  quiet = TRUE
)
st_crs(sa4_sf)

## clean
sa4_sf_clean <- sa4_sf |>
  select(
    code = SA4_CODE21,
    name = SA4_NAME21,
    geometry
  ) |>
  st_make_valid() |>
  st_zm(drop = TRUE)

## Albers projection
sa4_sf_albers <- st_transform(sa4_sf_clean, albers_wkt)

sa4_sf_albers_fixed <- st_sf(
  data.frame(
    code = sa4_sf_albers$code,
    name = sa4_sf_albers$name
  ),
  geometry = st_geometry(sa4_sf_albers)
)
stopifnot(
  nrow(sa4_sf_albers_fixed) == length(st_geometry(sa4_sf_albers_fixed))
)

## dissolve
sa4_sf_albers_dissolved <- sa4_sf_albers |>
  group_by(code, name) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

## sanity check
stopifnot(
  nrow(sa4_sf_albers_dissolved) == length(st_geometry(sa4_sf_albers_dissolved))
)

## check
st_is_valid(sa4_sf_albers_dissolved)
# should be all TRUE

exists("sa4_sf_albers_dissolved")
# should return TRUE

## remove empty geometries if any
which(sf::st_is_empty(sa4_sf_albers_dissolved))
sa4_sf_albers_clean <- sa4_sf_albers_dissolved |>
  filter(!st_is_empty(geometry))

## convert to sp
sa4_sp <- as(sa4_sf_albers_clean, "Spatial")

## convert for terra
sa4_v  <- vect(sa4_sp)
sa4_v  <- makeValid(sa4_v)

## sanity check
nrow(sa4_v)        # ~90 SA4s
length(sa4_v$name)

## area-weighted extraction by SA4
sa4_means <- suppressWarnings(
  terra::extract(
    tmean_annual_aea,
    sa4_v,
    fun     = mean,
    weights = TRUE,
    na.rm   = TRUE,
    ID      = TRUE
  )
)
summary(sa4_means[ , -1])
names(sa4_means)[1:10]

# long format
names(sa4_v)
sa4_v$ID

sa4_df <- sa4_means %>%
  rename(poly_id = ID) %>%
  left_join(
    data.frame(
      poly_id = seq_len(nrow(sa4_v)),
      sa4_code = sa4_v$code,
      sa4_name = sa4_v$name
    ),
    by = "poly_id"
  ) %>%
  tidyr::pivot_longer(
    cols = tidyselect::matches("^X\\d{4}$"),
    names_to  = "year",
    values_to = "tmean"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year))
  ) %>%
  select(sa4_code, sa4_name, year, tmean)


## baseline
baseline_sa4 <- sa4_df %>%
  filter(year >= 1900, year <= 1950) %>%
  group_by(sa4_code) %>%
  summarise(
    baseline = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

## anomalies
anomalies_sa4 <- sa4_df %>%
  left_join(baseline_sa4, by = "sa4_code") %>%
  mutate(anomaly = tmean - baseline)

summary(anomalies_sa4$anomaly)

## plot an example time series
## find Adelaide Hills in sa4_name
anomalies_sa4
grep("Adelaide", anomalies_sa4$sa4_name, value = TRUE)

one_sa4 <- anomalies_sa4 %>%
  filter(sa4_name == "Adelaide - Central and Hills")

ggplot(one_sa4, aes(x = year, y = anomaly)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.9, colour = "#2C7BB6") +
  labs(
    title = "temperature anomalies — Adelaide (Central & Hills)",
    subtitle = "1900–1950 baseline",
    x = "year",
    y = "anomaly (°C)"
  ) +
  theme_minimal(base_size = 12)

## plot time series for all South Australia SA4s
SA_SA4_anomalies <- anomalies_sa4[which(anomalies_sa4$sa4_code > 400 & anomalies_sa4$sa4_code < 500),]

ggplot(SA_SA4_anomalies, aes(x = year, y = anomaly)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(colour = "#D7191C", linewidth = 0.5) +
  facet_wrap(~ sa4_name, scales = "free_y") +
  labs(
    title = "SA4 temperature anomalies — South Australia",
    subtitle = "1900-1950 baseline",
    x = "year",
    y = "anomaly (°C)"
  ) +
  theme_minimal(base_size = 10)

## create map of 2025 anomalies by SA4
sa4_2025 <- anomalies_sa4 %>%
  filter(year == 2025) %>%
  select(sa4_code, sa4_name, anomaly)
summary(sa4_2025$anomaly)
names(sa4_sf_albers_dissolved)

map_sa4 <- sa4_sf_albers_dissolved %>%
  left_join(
    sa4_2025,
    by = c("code" = "sa4_code")
  )
sum(!is.na(map_sa4$anomaly))
summary(sa4_2025$anomaly)
sa4_2025_limits <- as.numeric(c(summary(sa4_2025$anomaly)[1], summary(sa4_2025$anomaly)[6]))

ggplot(map_sa4) +
  geom_sf(aes(fill = anomaly), colour = NA) +
  scale_fill_gradient2(
    name     = "2025 anomaly (°C)",
    low      = "#2C7BB6",
    mid      = "white",
    high     = "#D7191C",
    midpoint = 0,
    limits   = sa4_2025_limits,
    oob      = squish
  ) +
  coord_sf(crs = 3577) +
  labs(
    caption  = "source: SILO (BoM-derived)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    legend.position  = "right"
  )


################################################################
## split record time series (data_years) by state (legbodyCODE)
colnames(data_years)
table(data_years$legbodyCODE)
str(data_years)

data_yearsXstate <- data_years %>%
  group_split(legbodyCODE, .keep = TRUE) %>%
  setNames(sort(unique(data_years$legbodyCODE)))


## count number of records per year per state
state_year_sum <- data_years %>%
                    count(legbodyCODE, year, name = "n")

ggplot(state_year_sum,
       aes(x = year, y = n, colour = legbodyCODE)) +
  geom_line() +
  labs(
    x = "year",
    y = "number of records"
  ) +
  theme_minimal()

## plot number of records time series with state-specific temperature anomalies
## SA
head(anomalies_state)
SA_anom <- 
  anomalies_state %>%
  filter(state == "SA") %>%
  select(year, anomaly)
head(state_year_sum)
SA_rec <- subset(state_year_sum, legbodyCODE == "SA") %>%
  select(year, n)

SA_plot_df <- SA_rec %>%
  left_join(SA_anom, by = "year")
summary(SA_plot_df)

scale_factor <- max(cumsum(SA_plot_df$n), na.rm = TRUE) /
  max(abs(SA_plot_df$anomaly), na.rm = TRUE)

ggplot(SA_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(SA_plot_df)
ggplot(SA_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
SA.linreg <- linreg.ER(x=SA_plot_df$anomaly, y=cumsum(SA_plot_df$n))
SA.linreg

## QLD
head(anomalies_state)
QLD_anom <- 
  anomalies_state %>%
  filter(state == "QLD") %>%
  select(year, anomaly)
head(state_year_sum)
QLD_rec <- subset(state_year_sum, legbodyCODE == "QLD") %>%
  select(year, n)

QLD_plot_df <- QLD_rec %>%
  left_join(QLD_anom, by = "year")
summary(QLD_plot_df)

scale_factor <- max(cumsum(QLD_plot_df$n), na.rm = TRUE) /
  max(abs(QLD_plot_df$anomaly), na.rm = TRUE)

ggplot(QLD_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(QLD_plot_df)
ggplot(QLD_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
QLD.linreg <- linreg.ER(x=QLD_plot_df$anomaly, y=cumsum(QLD_plot_df$n))
QLD.linreg

## NSW
head(anomalies_state)
NSW_anom <- 
  anomalies_state %>%
  filter(state == "NSW") %>%
  select(year, anomaly)
head(state_year_sum)
NSW_rec <- subset(state_year_sum, legbodyCODE == "NSW") %>%
  select(year, n)

NSW_plot_df <- NSW_rec %>%
  left_join(NSW_anom, by = "year")
summary(NSW_plot_df)

scale_factor <- max(cumsum(NSW_plot_df$n), na.rm = TRUE) /
  max(abs(NSW_plot_df$anomaly), na.rm = TRUE)

ggplot(NSW_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(NSW_plot_df)
ggplot(NSW_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
NSW.linreg <- linreg.ER(x=NSW_plot_df$anomaly, y=cumsum(NSW_plot_df$n))
NSW.linreg

## VIC
head(anomalies_state)
VIC_anom <- 
  anomalies_state %>%
  filter(state == "VIC") %>%
  select(year, anomaly)
head(state_year_sum)
VIC_rec <- subset(state_year_sum, legbodyCODE == "VIC") %>%
  select(year, n)

VIC_plot_df <- VIC_rec %>%
  left_join(VIC_anom, by = "year")
summary(VIC_plot_df)

scale_factor <- max(cumsum(VIC_plot_df$n), na.rm = TRUE) /
  max(abs(VIC_plot_df$anomaly), na.rm = TRUE)

ggplot(VIC_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(VIC_plot_df)
ggplot(VIC_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
VIC.linreg <- linreg.ER(x=VIC_plot_df$anomaly, y=cumsum(VIC_plot_df$n))
VIC.linreg

## TAS
head(anomalies_state)
TAS_anom <- 
  anomalies_state %>%
  filter(state == "TAS") %>%
  select(year, anomaly)
head(state_year_sum)
TAS_rec <- subset(state_year_sum, legbodyCODE == "TAS") %>%
  select(year, n)

TAS_plot_df <- TAS_rec %>%
  left_join(TAS_anom, by = "year")
summary(TAS_plot_df)

scale_factor <- max(cumsum(TAS_plot_df$n), na.rm = TRUE) /
  max(abs(TAS_plot_df$anomaly), na.rm = TRUE)

ggplot(TAS_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(TAS_plot_df)
ggplot(TAS_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
TAS.linreg <- linreg.ER(x=TAS_plot_df$anomaly, y=cumsum(TAS_plot_df$n))
TAS.linreg

# NT
head(anomalies_state)
NT_anom <- 
  anomalies_state %>%
  filter(state == "NT") %>%
  select(year, anomaly)
head(state_year_sum)
NT_rec <- subset(state_year_sum, legbodyCODE == "NT") %>%
  select(year, n)

NT_plot_df <- NT_rec %>%
  left_join(NT_anom, by = "year")
summary(NT_plot_df)

scale_factor <- max(cumsum(NT_plot_df$n), na.rm = TRUE) /
  max(abs(NT_plot_df$anomaly), na.rm = TRUE)

ggplot(NT_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(NT_plot_df)
ggplot(NT_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
NT.linreg <- linreg.ER(x=NT_plot_df$anomaly, y=cumsum(NT_plot_df$n))
NT.linreg

## WA
head(anomalies_state)
WA_anom <- 
  anomalies_state %>%
  filter(state == "WA") %>%
  select(year, anomaly)
head(state_year_sum)
WA_rec <- subset(state_year_sum, legbodyCODE == "WA") %>%
  select(year, n)

WA_plot_df <- WA_rec %>%
  left_join(WA_anom, by = "year")
summary(WA_plot_df)

scale_factor <- max(cumsum(WA_plot_df$n), na.rm = TRUE) /
  max(abs(WA_plot_df$anomaly), na.rm = TRUE)

ggplot(WA_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(WA_plot_df)
ggplot(WA_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
WA.linreg <- linreg.ER(x=WA_plot_df$anomaly, y=cumsum(WA_plot_df$n))
WA.linreg


## ACT
head(anomalies_state)
ACT_anom <- 
  anomalies_state %>%
  filter(state == "ACT") %>%
  select(year, anomaly)
head(state_year_sum)
ACT_rec <- subset(state_year_sum, legbodyCODE == "ACT") %>%
  select(year, n)

ACT_plot_df <- ACT_rec %>%
  left_join(ACT_anom, by = "year")
summary(ACT_plot_df)

scale_factor <- max(cumsum(ACT_plot_df$n), na.rm = TRUE) /
  max(abs(ACT_plot_df$anomaly), na.rm = TRUE)

ggplot(ACT_plot_df, aes(x = year)) +
  # records (left axis)
  geom_line(
    aes(y = cumsum(n)),
    colour = "#2C7BB6",
    linewidth = 1
  ) +
  
  # temperature anomaly (right axis, rescaled)
  geom_line(
    aes(y = anomaly * scale_factor),
    colour = "#D7191C",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "cumulative number of records",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "anomaly (°C)"
    )
  ) +
  
  labs(
    x        = "year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#2C7BB6"),
    axis.text.y.left   = element_text(color = "#2C7BB6"),
    axis.title.y.right = element_text(color = "#D7191C"),
    axis.text.y.right  = element_text(color = "#D7191C")
  )

## bivariate plot of cumulative records and temperature anomalies by year
head(ACT_plot_df)
ggplot(ACT_plot_df, aes(x = anomaly, y = cumsum(n))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(x = "temperature anomaly (°C)", y = "cumulative number of records") +
  theme_minimal()
ACT.linreg <- linreg.ER(x=ACT_plot_df$anomaly, y=cumsum(ACT_plot_df$n))
ACT.linreg




state_2025_anomalies <- na.omit(data.frame(state=map_sf_ae$state, anomaly=map_sf_ae$anomaly))
state_2025_anomalies
linreg_R2 <- c(NSW.linreg[2], VIC.linreg[2], QLD.linreg[2], SA.linreg[2],
                     WA.linreg[2], TAS.linreg[2], NT.linreg[2], ACT.linreg[2])
linreg_R2
state_linreg_R2 <- data.frame(
  state = c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"),
  R2 = linreg_R2)
state_linreg_R2

## merge
state_2025_anomalies_linreg <- state_2025_anomalies %>%
  left_join(state_linreg_R2, by = c("state" = "state"))

## plot relationship
ggplot(state_2025_anomalies_linreg, aes(x = anomaly, y = R2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = state), size = 3) +
  labs(x = "2025 temperature anomaly (°C)", y = "R² of records-anomaly relationship") +
  theme_minimal()
linreg.ER(x=state_2025_anomalies_linreg$anomaly, y=state_2025_anomalies_linreg$R2)


## year of first record by state
first_record_year <- data_years %>%
  group_by(legbodyCODE) %>%
  summarise(first_year = min(year, na.rm = TRUE), .groups = "drop")

## average anomaly from year of initial record to 2025 by state
head(anomalies_state)

state_anom_joined <- anomalies_state %>%
  left_join(
    first_record_year,
    by = c("state" = "legbodyCODE")
  )

mean_anomaly_by_state <- state_anom_joined %>%
  filter(
    year >= first_year,
    year <= 2025
  ) %>%
  group_by(state) %>%
  summarise(
    first_year = first(first_year),
    mean_anomaly = mean(anomaly, na.rm = TRUE),
    sd_anomaly = sd(anomaly, na.rm = TRUE),
    se_anomaly = sd_anomaly / sqrt(n()),
    n_years = n(),
    .groups = "drop"
  )

## merge with state_linreg_R2
state_anomalies_linreg <- mean_anomaly_by_state %>%
  left_join(state_linreg_R2, by = c("state" = "state"))

## plot relationship
ggplot(state_anomalies_linreg, aes(x = mean_anomaly, y = R2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = state), size = 3) +
  labs(x = "mean temperature anomaly (°C)", y = "R² of records-anomaly relationship") +
  theme_minimal()
linreg.ER(x=state_anomalies_linreg$mean_anomaly, y=state_anomalies_linreg$R2)

## scale records time series by state for better comparison to temperature anomalies
## SA
SA_plot_df_scaled <- SA_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(SA_plot_df_scaled$anomaly, cumsum(SA_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
SA.lm <- lm(cumsum(SA_plot_df_scaled$n_scaled) ~ SA_plot_df_scaled$anomaly)
SA.lm.summ <- summary(SA.lm)
SA.lm.summ
SA.slope <- as.numeric(coefficients(SA.lm)[2])
SA.slope
SA.slope.se <- as.numeric(coef(summary(SA.lm))[2, 2])
SA.slope.se

## NSW
NSW_plot_df_scaled <- NSW_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(NSW_plot_df_scaled$anomaly, cumsum(NSW_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
NSW.lm <- lm(cumsum(NSW_plot_df_scaled$n_scaled) ~ NSW_plot_df_scaled$anomaly)
NSW.lm.summ <- summary(NSW.lm)
NSW.lm.summ
NSW.slope <- as.numeric(coefficients(NSW.lm)[2])
NSW.slope
NSW.slope.se <- as.numeric(coef(summary(NSW.lm))[2, 2])
NSW.slope.se

## QLD
QLD_plot_df_scaled <- QLD_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(QLD_plot_df_scaled$anomaly, cumsum(QLD_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
QLD.lm <- lm(cumsum(QLD_plot_df_scaled$n_scaled) ~ QLD_plot_df_scaled$anomaly)
QLD.lm.summ <- summary(QLD.lm)
QLD.lm.summ
QLD.slope <- as.numeric(coefficients(QLD.lm)[2])
QLD.slope
QLD.slope.se <- as.numeric(coef(summary(QLD.lm))[2, 2])
QLD.slope.se

## VIC
VIC_plot_df_scaled <- VIC_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(VIC_plot_df_scaled$anomaly, cumsum(VIC_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
VIC.lm <- lm(cumsum(VIC_plot_df_scaled$n_scaled) ~ VIC_plot_df_scaled$anomaly)
VIC.lm.summ <- summary(VIC.lm)
VIC.lm.summ
VIC.slope <- as.numeric(coefficients(VIC.lm)[2])
VIC.slope
VIC.slope.se <- as.numeric(coef(summary(VIC.lm))[2, 2])
VIC.slope.se

## TAS
TAS_plot_df_scaled <- TAS_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(TAS_plot_df_scaled$anomaly, cumsum(TAS_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
TAS.lm <- lm(cumsum(TAS_plot_df_scaled$n_scaled) ~ TAS_plot_df_scaled$anomaly)
TAS.lm.summ <- summary(TAS.lm)
TAS.lm.summ
TAS.slope <- as.numeric(coefficients(TAS.lm)[2])
TAS.slope
TAS.slope.se <- as.numeric(coef(summary(TAS.lm))[2, 2])
TAS.slope.se

## NT
NT_plot_df_scaled <- NT_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(NT_plot_df_scaled$anomaly, cumsum(NT_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
NT.lm <- lm(cumsum(NT_plot_df_scaled$n_scaled) ~ NT_plot_df_scaled$anomaly)
NT.lm.summ <- summary(NT.lm)
NT.lm.summ
NT.slope <- as.numeric(coefficients(NT.lm)[2])
NT.slope
NT.slope.se <- as.numeric(coef(summary(NT.lm))[2, 2])
NT.slope.se

## ACT
ACT_plot_df_scaled <- ACT_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(ACT_plot_df_scaled$anomaly, cumsum(ACT_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
ACT.lm <- lm(cumsum(ACT_plot_df_scaled$n_scaled) ~ ACT_plot_df_scaled$anomaly)
ACT.lm.summ <- summary(ACT.lm)
ACT.lm.summ
ACT.slope <- as.numeric(coefficients(ACT.lm)[2])
ACT.slope
ACT.slope.se <- as.numeric(coef(summary(ACT.lm))[2, 2])
ACT.slope.se

## WA
WA_plot_df_scaled <- WA_plot_df %>%
  mutate(n_scaled = n / max(cumsum(n), na.rm = TRUE))
plot(WA_plot_df_scaled$anomaly, cumsum(WA_plot_df_scaled$n_scaled), pch=19, xlab="temperature anomaly (°C)", ylab="scaled cumulative number of records")
WA.lm <- lm(cumsum(WA_plot_df_scaled$n_scaled) ~ WA_plot_df_scaled$anomaly)
WA.lm.summ <- summary(WA.lm)
WA.lm.summ
WA.slope <- as.numeric(coefficients(WA.lm)[2])
WA.slope
WA.slope.se <- as.numeric(coef(summary(WA.lm))[2, 2])
WA.slope.se

state.lm.slopes <- c(NSW.slope, VIC.slope, QLD.slope, SA.slope, WA.slope, TAS.slope, NT.slope, ACT.slope)
state.lm.slope.ses <- c(NSW.slope.se, VIC.slope.se, QLD.slope.se, SA.slope.se, WA.slope.se, TAS.slope.se, NT.slope.se, ACT.slope.se)  
state.lm.coeff <- data.frame(state=c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"),
                             slope = state.lm.slopes, slope_se = state.lm.slope.ses)
state.lm.coeff

## merge with state_anomalies_linreg
state_2025_anomalies_linreg_lm <- state_anomalies_linreg %>%
  left_join(
    state.lm.coeff,
    by = c("state" = "state")
  )
state_2025_anomalies_linreg_lm

## plot relationship between anomalies and slope of records-anomaly relationship
ggplot(state_2025_anomalies_linreg_lm, aes(x = mean_anomaly, y = slope)) +
  geom_point() +
  geom_errorbar(aes(ymin = slope - slope_se,ymax = slope + slope_se), width = 0.01,
                linewidth=0.2, color="grey7") +
  geom_errorbar(aes(xmin = mean_anomaly - se_anomaly, xmax = mean_anomaly + se_anomaly),
                width = 0.01, linewidth=0.2, color="grey7") +
  geom_smooth(method = "lm", fullrange=T, linetype="dashed") +
  geom_text_repel(aes(label = state), size = 4, color="grey7") +
  labs(x = "mean temperature anomaly (°C)", y = "slope of scaled cumulative records-anomaly relationship") +
  theme_minimal()
linreg.ER(x=state_2025_anomalies_linreg_lm$mean_anomaly, y=state_2025_anomalies_linreg_lm$slope)
lin.mod <- lm(slope ~ mean_anomaly, data = state_2025_anomalies_linreg_lm)
summary(lin.mod)


