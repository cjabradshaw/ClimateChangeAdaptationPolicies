## Climate change adaptation policy (CAP) dataset analysis
## Corey Bradshaw & Maddy King
## Apr 2026

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggrepel)

# import data
setwd("~/Documents/Students/King")

## read comma-delimited text file
data <- read.table("CAPdat.csv", header=TRUE, sep=",", dec=".", strip.white=TRUE, quote="\"")
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
politics_colnames <- colnames(data)[grepl("^Political_Party", colnames(data))]

## create left-right categories
data$pol_LEFT <- data$Political_Party_Labor
data$pol_RIGHT <- ifelse(data$Political_Party_Coalition | data$Political_Party_Liberal, 1, 0)

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
                               

#######################
## occurrence by year

## first, create dataset with expanding range between Start_Year and End_Year
## by year; if no 'End_Year', just use Start_Year; remove NAs first
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
