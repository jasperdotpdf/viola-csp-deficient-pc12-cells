library(tidyr)
library(dplyr)
library(afex)
library(emmeans)

citation("tidyr")
citation("dplyr")


#### ATP ####

df_ATP <- read.csv("Directory", h=T)

df_ATP$Time        <- factor(df_ATP$Time)
df_ATP$Cell.Line   <- factor(df_ATP$Cell.Line)
df_ATP$Substrate   <- factor(df_ATP$Substrate)
df_ATP$Antimycin.A <- factor(df_ATP$Antimycin.A)
df_ATP$Replicate   <- factor(df_ATP$Replicate)



### Full Model: 0, 15 and 60 min ###

## Four-Way ANOVA ## 

ATP_rm_model <- aov_ez(
  id      = "Replicate",
  dv      = "ATP",
  data    = df_ATP,
  within  = c("Cell.Line", "Substrate", "Time", "Antimycin.A"),
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov=TRUE
)
print(ATP_rm_model)



### Restricted Model: 15 and 60 min

str(df_ATP)
df_ATP_0<- df_ATP %>%
  filter(Time != "0") %>%
  droplevels()  

ATP_rm_model0 <- aov_ez(
  id      = "Replicate",
  dv      = "ATP",
  data    = df_ATP_0,
  within  = c("Cell.Line", "Substrate", "Time", "Antimycin.A"),
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov=TRUE
)
print(ATP_rm_model0)

## Cell Line post-hoc ##

emmeans(ATP_rm_model0,~Cell.Line | Substrate * Time * Antimycin.A,
        model="univariate")%>%
  contrast(method="pairwise")%>%
  summary(infer=TRUE,adjust="tukey")

## Antimycin A post-hoc ##

emmeans(ATP_rm_model0,~Antimycin.A | Substrate * Time * Cell.Line,
        model="univariate")%>%
  contrast(method="pairwise")%>%
  summary(infer=TRUE,adjust="tukey")


### Endpoint Substrate Analysis: Percentage of Initial ###

df_perc <- read.csv("Directory", h=T) 

df_perc$Time        <- factor(df_perc$Time)
df_perc$Cell.Line   <- factor(df_perc$Cell.Line)
df_perc$Substrate   <- factor(df_perc$Substrate)
df_perc$Antimycin.A <- factor(df_perc$Antimycin.A)
df_perc$Replicate   <- factor(df_perc$Replicate)
df_perc$ATP <- as.numeric(df_perc$ATP)

# subset data for 60 min endpoint

df_perc_60 <- droplevels(subset(df_perc, Time == "60"))


# subset antimycin A

df_noAntiA<-droplevels(subset(df_ATP_0,Antimycin.A=="0"))

rm_noAntiA <- aov_ez(
  id      = "Replicate",
  dv      = "ATP",
  data    = df_noAntiA,
  within  = c("Substrate","Cell.Line","Time"),  
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov = TRUE)

print(rm_noAntiA)

emmeans(rm_noAntiA, ~ Substrate | Cell.Line*Time,
        model="univariate") %>%
  contrast(method = "pairwise") %>%
  summary(infer = TRUE, adjust = "tukey")

df_AntiA<-droplevels(subset(df_ATP_0,Antimycin.A=="1"))

rm_AntiA<- aov_ez(
  id      = "Replicate",
  dv      = "ATP",
  data    = df_AntiA,
  within  = c("Substrate","Cell.Line","Time"),  
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov = TRUE)


print(rm_AntiA)

emmeans(rm_AntiA, ~ Substrate | Cell.Line*Time,
        model="univariate") %>%
  contrast(method = "pairwise") %>%
  summary(infer = TRUE, adjust = "tukey")




#### LACTATE ####

lac_df <- read.csv("Directory", h = T)

lac_df$Time        <- factor(lac_df$Time)
lac_df$Cell.Line   <- factor(lac_df$Cell.Line)
lac_df$Substrate   <- factor(lac_df$Substrate)
lac_df$Antimycin.A <- factor(lac_df$Antimycin.A)
lac_df$Replicate   <- factor(lac_df$Replicate)
str(lac_df)
View(lac_df)


str(df_ATP)
df_lac_0<- lac_df %>%
  filter(Time != "0") %>%
  droplevels()  

lac_rm_model <- aov_ez(
  id      = "Replicate",
  dv      = "Lactate",
  data    = df_lac_0,
  within  = c("Cell.Line", "Substrate", "Time", "Antimycin.A"),
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov=TRUE
)
print(lac_rm_model)

## 60 min Endpoint ##

lac_df_60<-droplevels(subset(lac_df,Time=="60"))


lac_60_rm_model <- aov_ez(
  id      = "Replicate",
  dv      = "Lactate",
  data    = lac_df_60,
  within  = c("Cell.Line", "Substrate", "Antimycin.A"),
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov="TRUE"
)
print(lac_60_rm_model)

# Cell line post-hoc # 

emmeans(lac_60_rm_model,~Cell.Line | Substrate * Antimycin.A,
        model="univariate")%>%
  contrast(method="pairwise")%>%
  summary(infer=TRUE,adjust="tukey")

# Antimycin A post-hoc #

emmeans(lac_60_rm_model,~Antimycin.A | Substrate * Cell.Line,
        model="univariate")%>%
  contrast(method="pairwise")%>%
  summary(infer=TRUE,adjust="tukey")
