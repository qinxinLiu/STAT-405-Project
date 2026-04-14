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

df %>%
  distinct(id, time) %>%
  count(time)

# concentration by time
df %>% 
  select(c(time, dv, dvid)) %>%
  filter(dvid == "cp") %>%
  group_by(time) %>%
  summarise(avg = mean(dv), sd = sd(dv)) 


## EDA ----
# concentration by time
ggplot(data = subset(df, dvid == "cp"), aes(y = dv, x = time, colour = id, group = id))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "Concentration (mg/L)", colour = "Patient ID") +
  theme_light()

# concentration by time and body weight
ggplot(data = subset(df, dvid == "cp"), aes(y = dv, x = time, colour = wt, group = id))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "Concentration (mg/L)", colour = "Weight (kg)") +
  scale_colour_viridis_c() +
  theme_light()

# concentration by time, age and sex
ggplot(data = subset(df, dvid == "cp"), aes(y = dv, x = time, colour = age, shape = sex, linetype = sex, group = id))+
  geom_point() +
  geom_line() +
  labs(x = "Time (h)", y = "Concentrationv(mg/L)", colour = "Age (years)", linetype = "Sex", shape = "Sex") +
  scale_linetype_manual(values = c("female" = "dashed", "male" = "solid")) +
  scale_shape_manual(values = c("female" = 17, "male" = 19)) +
  scale_colour_viridis_c() +
  theme_light() 
