library(tidyr)
library(dplyr)
library(afex)
library(emmeans)

lac_data<-read.csv("Directory", h=T)


lac_data$Time        <- factor(lac_data$Time)
lac_data$Cell.Line   <- factor(lac_data$Cell.Line)
lac_data$Glucose   <- factor(lac_data$Glucose)
lac_data$Inhibitor <- factor(lac_data$Inhibitor)
lac_data$Replicate   <- factor(lac_data$Replicate)


### Exclude t = 0 (Set at 0) ###

lac_data_noT0 <- lac_data %>%
  filter(Time != "0") %>%
  droplevels()  

lac_mod_noT0 <- aov_ez(
  id = "Replicate",
  dv = "Lactate",
  data = lac_data_noT0,
  within = c("Cell.Line","Glucose","Time","Inhibitor"),
  type = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov = TRUE
)

print(lac_mod_noT0)


# 2 h Endpoint Subset #

lac_data_2h<-subset(lac_data,Time=="2")
str(lac_data_2h)

lac_2h_model <- aov_ez(
  id      = "Replicate",
  dv      = "Lactate",
  data    = lac_data_2h,
  within  = c("Cell.Line", "Glucose", "Inhibitor"),
  type    = 3,
  anova_table = list(correction = "none", es = "pes"),
  include_aov = TRUE
)
print(lac_2h_model)

# Cell Line post-hoc #

emmeans(lac_2h_model,~Cell.Line | Glucose * Inhibitor,
                         model = "univariate")%>%
  contrast(method="pairwise")%>%
  summary(infer=TRUE,adjust="tukey")

# Antimycin A post-hoc #

emmeans(lac_2h_model,~Inhibitor | Glucose * Cell.Line,
                      model="univariate")%>%
  contrast(method="pairwise")%>%
  summary(infer=TRUE,adjust="tukey")
