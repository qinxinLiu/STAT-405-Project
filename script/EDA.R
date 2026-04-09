library(Certara.RsNLME) # data
library(dplyr)
library(ggplot2)

# DATA ----
## Dataset ----
df <- pkData
head(df)

df$Gender <- as.factor(df$Gender)
df$Subject <- as.factor(df$Subject)


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
# concentration by actual time
ggplot(data = df, aes(y = Conc, x = Act_Time, colour = Subject, group = Subject))+
  geom_point() +
  geom_line() +
  labs(x = "Actual Time (h)", y = "Concentration") +
  theme_light()

# log-concentration by nominal time
ggplot(data = df, aes(y = log(Conc), x = Nom_Time, colour = Subject, group = Subject))+
  geom_point() +
  geom_line() +
  theme_light()

# log-concentration by actual time and body weight
ggplot(data = df, aes(y = log(Conc), x = Act_Time, colour = BodyWeight, group = Subject))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "Concentration") +
  scale_colour_viridis_c() +
  theme_light()

# log-concentration by actual time, age and gender
ggplot(data = df, aes(y = log(Conc), x = Act_Time, colour = Age, shape = Gender, group = Subject))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "Concentration") +
  scale_colour_viridis_c() +
  theme_light() 
