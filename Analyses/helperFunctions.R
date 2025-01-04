library(tidyverse)
library(rubias)
library(sf)
library(ggspatial)
library(tibble)
library(dplyr)
library(adegenet)
library(mltools)
library(caret)

# Function to load and format data for rubias
loadDataForRubiasAll <- function(structureFolder, strucPrefix, numInds, numLoci, 
                                 popmap, modOutputFolder){
  
  # Load all data
  dat <- read.structure(
    file.path(structureFolder, paste0(strucPrefix, ".str")), 
    onerowperind=FALSE, n.ind=numInds, n.loc=numLoci, col.lab=1, 
    col.pop=2, ask=F, NA.char = "-9")
  datPops <- as.numeric(dat@pop)
  datSamples <- indNames(dat)
  
  # Convert genind to dataframe
  dat <- as.data.frame(dat)
  
  # Use popmap to determine whether an individual is known or unknown
  knownIndiv <-  popmap$sampleID[!is.na(popmap$pop)]
  unknownIndiv <-  popmap$sampleID[is.na(popmap$pop)]
  
  # Need to add relevant columns for rubias assignment
  # (1) sample_type: "reference" if known, "mixture" if unknown
  # (2) repunit: population if known, NA if unknown - must be character
  # (3) collection: population if known, name of sample if unknown?
  # (4) indiv: sample ID
  datMeta <- data.frame(
    sample_type = rep("reference", nrow(dat)), 
    repunit = datPops,
    collection = datPops, 
    indiv = datSamples)
  
  datMeta <- datMeta %>%
    mutate(sample_type = case_when(indiv %in% unknownIndiv ~ "mixture",
                                   TRUE ~ as.character(sample_type)),
           repunit = case_when(indiv %in% unknownIndiv ~ NA,
                               TRUE ~ as.character(repunit)),
           collection = case_when(indiv %in% unknownIndiv ~ "unknowns",
                                  TRUE ~ as.character(collection)))
  
  # Append metadata to genetic data, then combine knowns and unknowns
  dat <- cbind(datMeta, dat)
  
  # Convert repunit and collection columns to character
  dat$repunit <- as.character(dat$repunit)
  dat$collection <- as.character(dat$collection)
  
  # Save to output folder
  write.csv(dat, file=file.path(modOutputFolder, "rubias_reference.csv")) 
  
  return(dat)
}

# Function to load and format data for rubias
loadDataForRubias <- function(structureFolder, strucPrefix, numIndsKnown, 
                              numIndsUnknown, numLoci, modOutputFolder,
                              unknownSuffix = "_unknown"){
  
  # Load known data
  dat <- read.structure(
    file.path(structureFolder, paste0(strucPrefix, ".str")), 
    onerowperind=FALSE, n.ind=numIndsKnown, n.loc=numLoci, col.lab=1, 
    col.pop=2, ask=F, NA.char = "-9")
  datPops <- as.numeric(dat@pop)
  datSamples <- indNames(dat)
  
  # Load unknowns
  datUnknown <- read.structure(
    file.path(structureFolder, paste0(strucPrefix, unknownSuffix, ".str")), 
    onerowperind=FALSE, n.ind=numIndsUnknown, n.loc=numLoci, col.lab=1, 
    col.pop=2, ask=F, NA.char = "-9")
  datUnknownSamples <- indNames(datUnknown)
  
  # Convert genind to dataframe
  dat <- as.data.frame(dat)
  datUnknown <- as.data.frame(datUnknown)
  
  # Need to add relevant columns for rubias assignment
  # (1) sample_type: "reference" if known, "mixture" if unknown
  # (2) repunit: population if known, NA if unknown - must be character
  # (3) collection: population if known, name of sample if unknown?
  # (4) indiv: sample ID
  datMeta <- data.frame(
    sample_type = rep("reference", nrow(dat)), 
    repunit = datPops,
    collection = datPops, 
    indiv = datSamples)
  
  datUnknownMeta <- data.frame(
    sample_type = rep("mixture", nrow(datUnknown)), 
    repunit = rep(NA, nrow(datUnknown)),
    collection = rep("unknowns", nrow(datUnknown)), 
    indiv = datUnknownSamples)
  
  # Append metadata to genetic data, then combine knowns and unknowns
  dat <- cbind(datMeta, dat)
  datUnknown <- cbind(datUnknownMeta, datUnknown)
  
  # Drop any loci that differ between the two datasets
  lociToDropUnknown <- setdiff(names(datUnknown), names(dat))
  datUnknown <- datUnknown %>% dplyr::select(-all_of(lociToDropUnknown))
  lociToDropKnown <- setdiff(names(dat), names(datUnknown))
  dat <- dat %>% dplyr::select(-all_of(lociToDropKnown))
  
  # Convert repunit and collection columns to character
  dat$repunit <- as.character(dat$repunit)
  dat$collection <- as.character(dat$collection)
  datUnknown$repunit <- as.character(datUnknown$repunit)
  datUnknown$collection <- as.character(datUnknown$collection)
  
  # Bind knowns and unknowns
  datAll <- bind_rows(dat, datUnknown)
  
  # Save to output folder
  write.csv(dat, file=file.path(modOutputFolder, "rubias_reference.csv")) 
  
  return(datAll)
}

# Function to calculate performance metrics
reportPerformanceMetrics <- function(mat, pop_xwalk, unknown_dat=NA){
  
  precision <- diag(mat) / rowSums(mat)
  recall <- diag(mat) / colSums(mat)
  f1 <- (2 * precision * recall) / (precision + recall)
  cm <- matrix(mat, nrow=length(pop_xwalk))
  
  macroPrecision <- mean(precision)
  macroRecall <- mean(recall)
  macroF1 <- mean(f1)
  macroMCC <- mcc(confusionM=cm)
  
  if (!is.data.frame(unknown_dat)) {
    pop_accuracies <- data.frame(diag(mat))
    pop_accuracies  <- cbind(popID = rownames(pop_accuracies ), pop_accuracies)
    pop_accuracies <- pop_xwalk %>% 
      gather() %>% 
      merge(pop_accuracies, by.x = "value", by.y = "popID") %>%
      dplyr::select(key, diag.mat.) %>%
      dplyr::rename(value = diag.mat.)
  } else {
    unknown_dat <- unknown_dat %>% 
      mutate(success = case_when(pred.pop == real.pop ~ 1,
                                 pred.pop != real.pop ~ 0))
    pop_accuracies <- unknown_dat %>%
      group_by(real.pop) %>%
      summarize(accuracy = mean(success))
    pop_accuracies <- pop_xwalk %>% 
      gather() %>%
      merge(pop_accuracies, by.x = "value", by.y = "real.pop") %>%
      dplyr::select(key, accuracy) %>%
      rename(value = accuracy)
  }
  
  macroAccuracy <- mean(pop_accuracies$value)
  
  results <- data.frame(accuracy = c(macroAccuracy),
                        precision = c(macroPrecision),
                        recall = c(macroRecall),
                        f1 = c(macroF1),
                        mcc = c(macroMCC))
  
  results <- rbind(results %>% gather(), pop_accuracies)
  results$value <- round(results$value, 3)
  
  return(results)
}

# Function to calculate performance metrics
reportPerformanceMetricsEmpirical <- function(mat, pop_xwalk, unknown_dat=NA, popID=TRUE){
  
  precision <- (diag(mat) / colSums(mat)) %>% replace_na(0)
  recall <- (diag(mat) / rowSums(mat)) %>% replace_na(0)
  f1 <- ((2 * precision * recall) / (precision + recall)) %>% replace_na(0)
  numPops <- dim(mat)[1]
  cm <- matrix(mat, nrow=numPops)
  
  macroPrecision <- mean(precision)
  macroRecall <- mean(recall)
  macroF1 <- mean(f1)
  macroMCC <- mcc(confusionM=cm)
  
  if (!is.data.frame(unknown_dat)) {
    pop_accuracies <- data.frame(diag(mat))
    pop_accuracies  <- cbind(popID = rownames(pop_accuracies ), pop_accuracies)
    
    if (popID){
      pop_accuracies <- pop_xwalk %>% 
        gather() %>% 
        merge(pop_accuracies, by.x = "key", by.y = "popID") %>%
        dplyr::select(key, diag.mat.) %>%
        dplyr::rename(value = diag.mat.)
    } else {
      pop_accuracies <- pop_xwalk %>% 
        gather() %>% 
        merge(pop_accuracies, by.x = "value", by.y = "popID") %>%
        dplyr::select(key, diag.mat.) %>%
        dplyr::rename(value = diag.mat.)
    }

  } else {
    unknown_dat <- unknown_dat %>% 
      mutate(success = case_when(real.pop == pred.pop ~ 1,
                                 real.pop != pred.pop ~ 0))
    pop_accuracies <- unknown_dat %>%
      group_by(real.pop) %>%
      dplyr::summarize(accuracy = mean(success))
    pop_accuracies <- pop_xwalk %>% 
      gather() %>%
      merge(pop_accuracies, by.x = "key", by.y = "real.pop") %>%
      dplyr::select(key, accuracy) %>%
      dplyr::rename(value = accuracy)
  }
  
  macroAccuracy <- mean(pop_accuracies$value)
  
  results <- data.frame(accuracy = c(macroAccuracy),
                        precision = c(macroPrecision),
                        recall = c(macroRecall),
                        f1 = c(macroF1),
                        mcc = c(macroMCC))
  
  results <- rbind(results %>% gather(), pop_accuracies)
  results$value <- round(results$value, 3)
  
  return(results)
}

# Function for turning assignment results into a confusion matrix
createConfusionMatrix <- function(df, pop_xwalk, popID = TRUE){
 
  numPops <- ncol(pop_xwalk)
  
  if (popID){
    pops <- pop_xwalk %>%
      gather() %>%
      pull(value)
  } else {
    pops <- pop_xwalk %>%
      gather() %>%
      pull(key)
  }
  
  # Change counts of assignment to rates of assignment
  df$real.pop <- factor(df$real.pop, levels=levels(factor(pops)))
  df$pred.pop <- factor(df$pred.pop, levels=levels(factor(pops)))
  assign_tab <- table(df$pred.pop, df$real.pop)
  assign_mat <- matrix(assign_tab, ncol = length(pop_xwalk)) 
  colnames(assign_mat) <- colnames(assign_tab)
  rownames(assign_mat) <- rownames(assign_tab)
  
  for(i in 1:numPops){
    assign_mat[,i] <- round(assign_mat[,i]/sum(assign_mat[,i]), digits=2)
  }
  
  assign_mat[is.nan(assign_mat)] <- 0
  
  return(assign_mat)
}
