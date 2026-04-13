library(dplyr)
library(ggplot2)

# DATA ----
## Dataset ----
load(file = "data/warfarin.rda")
df <- warfarin
head(df)

df <- df %>%
  mutate(id = as.factor(id))


## Summary Statistics ----
# mean and sd of numerical variables
df %>%
  distinct(id, age, wt) %>%
  summarise(mean_age = mean(age), mean_weight = mean(wt),
            sd_age = sd(age), sd_weight = sd(wt))

# frequency table of categorical variables
df %>% 
  distinct(id, sex) %>%
  count(sex)

# concentration by nominal time
df %>% 
  select(c(time, dv)) %>%
  group_by(time) %>%
  summarise(avg = mean(dv), sd = sd(dv)) 


## EDA ----
# concentration by actual time
ggplot(data = df, aes(y = dv, x = time, colour = id, group = id))+
  geom_point() +
  geom_line() +
  labs(x = "Actual Time (h)", y = "dventration") +
  theme_light()

# log-dventration by nominal time
ggplot(data = df, aes(y = log(dv), x = time, colour = id, group = id))+
  geom_point() +
  geom_line() +
  theme_light()

# log-dventration by actual time and body weight
ggplot(data = df, aes(y = log(dv), x = time, colour = wt, group = id))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "dventration") +
  scale_colour_viridis_c() +
  theme_light()

# log-dventration by actual time, age and sex
ggplot(data = df, aes(y = log(dv), x = time, colour = age, shape = sex, group = id))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "dventration") +
  scale_colour_viridis_c() +
  theme_light() 
