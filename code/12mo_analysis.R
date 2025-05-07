
## required packages ##
library(pcnetmeta)
library(dplyr)

## set wd ##
getwd()
setwd("path/to/dir")
df <- read.csv("data.csv", header=TRUE) # dataset is available in the data dir

## data preparation ##
df <- df %>%
  mutate(Treatment = case_when(
    Treatment == "3_Small_Lesion" ~ "3smalllesions",
    Treatment == "3_Large_Lesion" ~ "3largelesions",
    Treatment == ">3_Large_Lesion" ~ "Morethan3largelesions",
    Treatment == "OTHER" ~ "Other",
    TRUE ~ Treatment # Keep other values as they are
  ))
df
df <- subset(df, followup == 12) # subset 12-month timepoint data
df <- df %>%
  filter(!is.na(X50oroverpainreduction_n))
result <- df %>%
  group_by(studies) %>%
  summarise(TreatmentCount = n_distinct(Treatment)) %>%
  count(TreatmentCount, name = "NumberOfStudies")
(result)
unique_studies <- unique(df$studies)
(unique_studies)

## network plot ##
nma.networkplot(studies, Treatment, data = df, title = "12-Month",
                #trtname = c("NC", "SH", "IC", "GC"),
                alphabetic = FALSE,
                weight.edge = TRUE, adjust.thick = 3, weight.node = TRUE,
                adjust.node.size = 10, node.col = "orange", edge.col = "lightblue",
                text.cex = 1.1, adjust.figsizex = 1.2, adjust.figsizey = 1.2)

## model building ##
set.seed(12345)
model.out <- nma.ab.bin(studies, Treatment, X50oroverpainreduction_n , X50oroverpainreduction_N, data = df, 
                        trtname = c("Morethan3largelesions", "3largelesions","3smalllesions",  "Other"
                        ), 
                        param = c("AR", "OR", "RR", "LOR", "LRR", "RD",
                                  "rank.prob"), model = "het_cor", higher.better = TRUE, digits = 3,
                        n.adapt = 20000, n.iter = 200000, n.thin = 2, conv.diag = TRUE,
                        dic = TRUE, trace = "LOR", postdens = TRUE)

model.out$OddsRatio$Median_CI # check for 1 to see which trt is more effective or differ significantly
model.out$TrtRankProb
model.out$DIC # prefer the model w lower dic

## model results ##
# treatment-specific effect 95% CI risk plot
absolute.plot(model.out, width = 5, height = 1.5, 
              network.name="12mo", alphabetic = FALSE)

# treatment-specific relative effect, pairwise 95% CI contrast plot
contrast.plot(model.out, reference = "Other", width = 5, height = 1.5, 
              network.name="12mo",
              digit=2) # 

# rank probability 
rank.prob(model.out, cex.axis = 1, cex.lab = 2)
title("12-Month Rank Probability")

# patients' response rate data
final_df <- df %>%
  group_by(Treatment) %>%
  summarise(
    Total_50oroverpainreduction_n = sum(X50oroverpainreduction_n, na.rm = TRUE),
    Total_50oroverpainreduction_N = sum(X50oroverpainreduction_N, na.rm = TRUE)
  ) %>%
  mutate(
    Percent_n = (Total_50oroverpainreduction_n / Total_50oroverpainreduction_N) * 100
  ) %>%
  ungroup()
# save the resulting dataframe to a CSV file
write.csv(final_df, "12mo_final_data.csv", row.names = FALSE)





