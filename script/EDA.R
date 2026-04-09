library(Certara.RsNLME) # data
library(dplyr)

# DATA ----
## Dataset ----
df <- pkData
head(df)


## Summary Statistics ----
# mean and sd of numerical variables
df %>%
  distinct(Subject, Age, BodyWeight) %>%
  summarise(mean_age = mean(Age), mean_weight = mean(BodyWeight),
            sd_age = sd(Age), sd_weight = sd(BodyWeight))

# frequency table of categorical variables
df %>% 
  distinct(Subject, Gender) %>%
  count(Gender)

# concentration by nominal time
df %>% 
  select(c(Nom_Time, Conc)) %>%
  group_by(Nom_Time) %>%
  summarise(avg = mean(Conc), sd = sd(Conc)) 


## EDA ----




