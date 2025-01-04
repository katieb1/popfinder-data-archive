# Create function for plotting confusion matrix based on predicted/actual
# Assuming columns = sampleID, predicted, actual

library(tidyverse)

# AssignPOP ----
dataPath <- "F:/PopFinder/Results/assignPOP/mod1_test_output/assignment_predictions.csv"
outputPath <- "F:/PopFinder/Results/assignPOP/mod1_test_output/LS_assignPOP.png"
df <- read_csv(dataPath)
empirical = FALSE

# TODO: uncomment when running empirical
datasetInput <- "lesp_all_54"
datasetOutput <- datasetInput
empirical = TRUE
assignPOPResultsFolder <- file.path("F:/PopFinder/Results/assignPOP/", datasetOutput)
dataPath <- file.path(assignPOPResultsFolder, "test_results.csv")
outputPath <- file.path(assignPOPResultsFolder, paste0(datasetInput, "_assignPOP.png"))
df <- read_csv(dataPath)

highColour <- "orangered3"
lowColour <- '#FFFFFF'

# Simulations
if (!empirical){
  pop_xwalk <- data.frame(p1 = c(2),
                          p2 = c(3),
                          p3 = c(4),
                          p4 = c(5),
                          p5 = c(6))
} else if (datasetInput == "lesp_all_54"){
  pop_xwalk <- data.frame(popname = c("Hernyken", "Baccalieu", "Corrosol", 
                                      "Kent", "Green"),
                          popID = c(1,2,3,4,5))
} else if (datasetInput == "tbmu"){
  pop_xwalk <- data.frame(popname = c("Akpatok", "Digges Island", "Gannet Island",
                                      "Hornoya", "Latrabjarg", "Minarets",
                                      "Prince Leopold Island"),
                          popID = c(2,6,8,9,10,11,12))
} else if (datasetInput == "nofu4_all"){
  pop_xwalk <- data.frame(popname = c("Cape Vera", "Faeroe Islands", "Qaquluit", 
                                      "Prince Leopold Island"),
                          popID = c(1,2,3,4))
}





if (empirical){
  numPops <- nrow(pop_xwalk)
  df <- df %>% dplyr::rename(Actual = origin, 
                             Predicted = assignment,
                             Frequency = Freq)
  df$Actual <- pop_xwalk$popname[match(unlist(df$Actual), pop_xwalk$popID)]
  df$Predicted <- pop_xwalk$popname[match(unlist(df$Predicted), pop_xwalk$popID)]
} else {
  numPops <- length(pop_xwalk)
  pop_xwalk_mod <- pop_xwalk %>% gather(key = "popname", value = "popID")
  df$Actual <- pop_xwalk_mod$popname[match(unlist(df$actual), pop_xwalk_mod$popID)]
  df$Predicted <- pop_xwalk_mod$popname[match(unlist(df$predicted), pop_xwalk_mod$popID)]
  
  cm <- table(df$Actual, df$Predicted)
  
  # normalize
  for(i in 1:numPops){
    cm[i,] <- round(cm[i,]/sum(cm[i,]), digits=2)
  }
  
  cm <- as.data.frame(cm)
  
  cm <- cm %>% mutate(Freq = case_when(is.nan(Freq) ~ 0,
                                       TRUE ~ Freq))
  cm <- cm %>% dplyr::rename("Frequency" = "Freq",
                             "Actual" = "Var1",
                             "Predicted" = "Var2")
}



# Prediction = y axis = Var2
# Actual = x axis = Var 1
# TODO: data = df when empirical, data = cm when simulations
if (empirical){
  plotData <- df
} else {
  plotData <- cm
}

p <- ggplot(data = plotData,
       mapping = aes(x = Actual,
                     y = Predicted)) +
  geom_tile(aes(fill = Frequency)) +
  geom_text(aes(label = sprintf("%.2f", Frequency)), vjust = 1) +
  scale_fill_gradientn(limits = c(0, 1),
                      breaks=c(0, 0.25, 0.5, 0.75, 1), 
                      colours = c(lowColour, highColour)) +
  xlab("Origin Population") +
  ylab("Predicted Population") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(outputPath, p, width=1712, height=1296, units="px")

# For test results LS sim
output_folder = file.path(outputFolderPath, "mod1_test_output")
assign_mat <- assign.matrix_mod(dir=paste0(output_folder, "/"), k.fold=5)
assign_dat <- data.frame(assign_mat[[1]])
assign_dat <- assign_dat %>%
  rename("Frequency" = Freq)

p <- ggplot(data = assign_dat,
            mapping = aes(x = origin,
                          y = assignment)) +
  geom_tile(aes(fill = Frequency)) +
  geom_text(aes(label = sprintf("%.2f", Frequency)), vjust = 1) +
  scale_fill_gradientn(limits = c(0, 1),
                       breaks=c(0, 0.25, 0.5, 0.75, 1), 
                       colours = c('#FFFFFF', '#008080')) +
  xlab("Origin Population") +
  ylab("Predicted Population") +
  theme_classic()

outputPath <- "F:/PopFinder/Results/assignPOP/mod1_test_output/LS_sim_test_cm.png"
ggsave(outputPath, p)

# Popfinder ----
# Modify to generate confusion matrices for the various datasets
# Options = LS_outputs, MS_outputs, HS_outputs, lesp, tbmu, nofu
datasetInput <- "oct27_lesp" 
datasetOutput <- "oct27_lesp"
empirical <- T
popmapName <- "popmap_54.txt"

# Don't change below
if (empirical) {
  popfinderInputsFolder <- "F:/PopFinder/Data/empirical data/"
} else {
  popfinderInputsFolder <- "F:/PopFinder/Data/simulated data/"
}
popfinderResultsFolder <- file.path("F:/PopFinder/Results/popfinder/", datasetOutput)
if (empirical){
  dataPath <- file.path(popfinderResultsFolder, "classifier_test_results.csv")
} else {
  dataPath <- file.path(popfinderResultsFolder, "classifier_assignment_results.csv")
}
outputPath <- file.path(popfinderResultsFolder, paste0(datasetInput, "_popfinder.png"))

highColour <- "darkslateblue"
lowColour <- '#FFFFFF'

df <- read_csv(dataPath)

if (empirical){
  # popmap <- read_tsv(file.path(popfinderInputsFolder, datasetInput, popmapName))
  df <- df %>%
    # merge(popmap, by = c("sampleID")) %>%
    dplyr::rename(predicted = "pred_pop", 
           actual = "true_pop") %>%
    dplyr::select(predicted, actual)
  
  if (datasetInput == "lesp_all_54"){
    pop_xwalk <- data.frame(popname = c("Hernyken", "Baccalieu", "Corrosol", 
                                        "Kent", "Green"),
                            popID = c("Hernyken", "Baccalieu", "Corrosol", 
                                      "Kent", "Green"))
  } else if (datasetInput == "tbmu"){
    pop_xwalk <- data.frame(popname = c("Akpatok", "Digges Island", "Gannet Island",
                                        "Hornoya", "Latrabjarg", "Minarets",
                                        "Prince Leopold Island"),
                            popID = c("akp", "dig", "gan", "hor", "lat", "min", "pli"))
  } else if (datasetInput == "nofu4_all"){
    pop_xwalk <- data.frame(popname = c("Cape Vera", "Faeroe Islands", "Qaquluit", 
                                        "Prince Leopold Island"),
                            popID = c("CV","FI","QU", "PLI"))
  }
  
  # df$actual <- pop_xwalk$popname[match(unlist(df$actual), pop_xwalk$popID)]
  # df$predicted <- pop_xwalk$popname[match(unlist(df$predicted), pop_xwalk$popID)]
 # df <- df %>% filter(!actual %in% c("cog", "fun", "coa"))
} else {
  df <- df %>%
    rename(predicted = "most_assigned_pop_across_models", 
           actual = "real_pop") %>%
    select(sampleID, predicted, actual)
}

pops <- unique(df$actual)
numPops <- length(pops)

df$actual <- factor(df$actual, levels=levels(factor(pops)))
df$predicted <- factor(df$predicted, levels=levels(factor(pops)))
cm <- table(df$predicted, df$actual)

# normalize
for(i in 1:numPops){
  cm[,i] <- round(cm[,i]/sum(cm[,i]), digits=2)
}

cm <- as.data.frame(cm)
cm <- cm %>% mutate(Freq = case_when(is.nan(Freq) ~ 0,
                               TRUE ~ Freq))

cm <- cm %>% dplyr::rename("Frequency" = "Freq",
                    "Actual" = "Var2",
                    "Predicted" = "Var1")

# Prediction = y axis = Var2
# Actual = x axis = Var 1
p <- ggplot(data = cm,
            mapping = aes(x = Actual,
                          y = Predicted)) +
  geom_tile(aes(fill = Frequency)) +
  geom_text(aes(label = sprintf("%.2f", Frequency)), vjust = 1) +
  scale_fill_gradientn(limits = c(0, 1),
                       breaks=c(0, 0.25, 0.5, 0.75, 1), 
                       colours = c(lowColour, highColour)) +
  xlab("Origin Population") +
  ylab("Predicted Population") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(outputPath, p, width=1712, height=1296, units="px")

# Rubias ----
datasetInput <- "lesp" 
datasetOutput <- "lesp"
empirical <- T
highColour <- "seagreen"
lowColour <- '#FFFFFF'

# Don't change below

rubiasResultsFolder <- file.path("F:/PopFinder/Results/rubias/", datasetOutput)
if (empirical){
  dataPath <- file.path(rubiasResultsFolder, "rubias_self_assignments.csv")
} else {
  dataPath <- file.path(rubiasResultsFolder, "assignment_predictions.csv")
}
outputPath <- file.path(rubiasResultsFolder, paste0(datasetInput, "_rubias.png"))

df <- read_csv(dataPath)


# Simulations
if (!empirical){
  pop_xwalk <- data.frame(p1 = c(1),
                          p2 = c(2),
                          p3 = c(3),
                          p4 = c(4),
                          p5 = c(5))
} else   if (datasetInput == "lesp"){
  pop_xwalk <- data.frame(popname = c("Hernyken", "Baccalieu", "Corrosol", 
                                      "Kent", "Green"),
                          popID = c(1, 2, 3, 4, 5))
} else if (datasetInput == "tbmu"){
  pop_xwalk <- data.frame(popname = c("Akpatok", "Digges Island", "Gannet Island",
                                      "Hornoya", "Latrabjarg", "Minarets",
                                      "Prince Leopold Island"),
                          popID = c(1,4,6,7,8,9,10))
} else if (datasetInput == "nofu4"){
  pop_xwalk <- data.frame(popname = c("Cape Vera", "Faeroe Islands", "Qaquluit", 
                                      "Prince Leopold Island"),
                          popID = c(1, 2, 3, 4))
}

# assign_mat <- createConfusionMatrix(df, pop_xwalk)
# colnames(assign_mat) <- colnames(pop_xwalk)
# rownames(assign_mat) <- colnames(pop_xwalk)
# cm <- assign_mat

# Don't change below


if (empirical){
  numPops <- nrow(pop_xwalk)
  pops <- pop_xwalk$popname
  df$Actual <- pop_xwalk$popname[match(unlist(df$real.pop), pop_xwalk$popID)]
  df$Predicted <- pop_xwalk$popname[match(unlist(df$pred.pop), pop_xwalk$popID)]
} else {
  numPops <- ncol(pop_xwalk)
  pops <- colnames(pop_xwalk)
  pop_xwalk_mod <- pop_xwalk %>% gather(key = "popname", value = "popID")
  df$Actual <- pop_xwalk_mod$popname[match(unlist(df$actual), pop_xwalk_mod$popID)]
  df$Predicted <- pop_xwalk_mod$popname[match(unlist(df$predicted), pop_xwalk_mod$popID)]
}


# Change counts of assignment to rates of assignment
df$Actual <- factor(df$Actual, levels=levels(factor(pops)))
df$Predicted <- factor(df$Predicted, levels=levels(factor(pops)))

cm <- table(df$Actual, df$Predicted)

# normalize
for(i in 1:numPops){
  cm[i,] <- round(cm[i,]/sum(cm[i,]), digits=2)
}

cm <- as.data.frame(cm)

# if (!empirical){
#   cm <- pop_xwalk %>%
#     gather() %>%
#     merge(cm, by.x=c("value"), by.y=c("Var1")) %>%
#     dplyr::rename("Actual" = "key") %>%
#     dplyr::select(-value)
#   cm <- pop_xwalk %>%
#     gather() %>%
#     merge(cm, by.x=c("value"), by.y=c("Var2")) %>%
#     dplyr::rename("Predicted" = "key") %>%
#     dplyr::select(-value)
#   cm <- cm %>% dplyr::rename("Frequency" = "Freq")
# } else {
cm <- cm %>% mutate(Freq = case_when(is.nan(Freq) ~ 0,
                                     TRUE ~ Freq))
cm <- cm %>% dplyr::rename("Frequency" = "Freq",
                    "Actual" = "Var1",
                    "Predicted" = "Var2")
# }

# Prediction = y axis = Var2
# Actual = x axis = Var 1
p <- ggplot(data = cm,
            mapping = aes(x = Actual,
                          y = Predicted)) +
  geom_tile(aes(fill = Frequency)) +
  geom_text(aes(label = sprintf("%.2f", Frequency)), vjust = 1) +
  scale_fill_gradientn(limits = c(0, 1),
                       breaks=c(0, 0.25, 0.5, 0.75, 1), 
                       colours = c(lowColour, highColour)) +
  xlab("Origin Population") +
  ylab("Predicted Population") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(outputPath, p, width=1712, height=1296, units="px")
