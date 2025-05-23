#!/usr/bin/env Rscript
require(tidyverse)
require(data.table)
require(GenomicRanges)

#############
# ARGUMENTS #
#############
# [1] = project directory
#setwd("/projects/b1059/projects/Ryan/telomere_variation/NemaScan")
#args <- list("Analysis_Results-20220714")
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
baseDir <- paste0(paste(strsplit(getwd(), split = "/")[[1]][1:(length(strsplit(eval(getwd()), split = "/")[[1]])-2)],collapse = "/"),"/")
today <- format(Sys.time(), '%Y%m%d')

pop_designation <- args[2]

# Simulated QTLs and Effects
effects <- list.files(pattern = "sims.par",recursive = T)
print("Gathered Simulated QTL")

if(pop_designation == TRUE){
   iterations <- purrr::map(effects, .f = function(x){
      paste(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "_")[[1]][1:7], collapse = "_") # QUEST
})
} else{
   iterations <- purrr::map(effects, .f = function(x){
      paste(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "_")[[1]][1:6], collapse = "_") # QUEST
})
}
# Assessing Mapping Performance
print("Measuring Performance")
simulation.metrics <- function(x, pop_designation = FALSE){
   
   #Pull out params from simulation file name
   # I decided to add the population designations to the file name instead of having 
   # to do a more complicated join later on in analysis. Since old populations don't have
   # this syntax I made this optional
   if(pop_designation == TRUE){
      nQTL <- strsplit(x,split = "_")[[1]][1]
      rep <- strsplit(x,split = "_")[[1]][2]
      h2 <- strsplit(x,split = "_")[[1]][3]
      MAF <- strsplit(x,split = "_")[[1]][4]
      effect.range <- strsplit(x,split = "_")[[1]][5]
      sample.designation <- strsplit(x,split = "_")[[1]][6]
      sample.population <- strsplit(x,split = "_")[[1]][7]

      geno_file <- paste(sample.designation, sample.population,MAF,"Genotype_Matrix.tsv",sep = "_")
   }else{
      nQTL <- strsplit(x,split = "_")[[1]][1]
      rep <- strsplit(x,split = "_")[[1]][2]
      h2 <- strsplit(x,split = "_")[[1]][3]
      MAF <- strsplit(x,split = "_")[[1]][4]
      effect.range <- strsplit(x,split = "_")[[1]][5]
      sample.population <- strsplit(x,split = "_")[[1]][6]

      geno_file <- paste(sample.population,MAF,"Genotype_Matrix.tsv",sep = "_")
   }
   print(x)
   
   # Effects
   effects <- data.table::fread(paste(paste(effect.range,
                                     nQTL,
                                     "Phenotypes", sep = "/"),
                                     paste(x,"sims.par",sep = "_"), sep = "/"),
                                header = T)
   effects <- effects %>% 
      tidyr::separate(QTL, c("CHROM","POS"), sep = ":", remove = F)
   
   # Phenotypes
   phenos <- data.table::fread(paste(paste(effect.range,
                                            nQTL,
                                            "Phenotypes", sep = "/"),
                                      paste(x,"sims.phen",sep = "_"), sep = "/"),
                                header = F)  %>%
      `colnames<-`(c("strain","strain2","trait.value")) %>%
      dplyr::select(-strain2)
   
   # Genotype Matrix
   
   complete.effects <- data.table::fread(paste(baseDir, 
                                               args[1],
                                               "/Genotype_Matrix/",
                                               geno_file, sep = ""),
                                         header = T) %>%
      tidyr::unite("QTL",c(CHROM, POS), sep = ":", remove = F) %>%
      dplyr::filter(QTL %in% effects$QTL) %>%
      dplyr::select(-CHROM, -POS) %>%
      tidyr::pivot_longer(cols = !c(QTL, REF, ALT), names_to = "strain", values_to = "allele") %>%
      dplyr::full_join(.,effects) %>%
      dplyr::full_join(.,phenos) %>%
      dplyr::group_by(QTL) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n != 1,
                    !is.na(QTL))
   
   
   if(nrow(complete.effects) == nrow(effects)){
      genos.effects <- data.table::fread(paste(baseDir,
                                               args[1],
                                               "/Genotype_Matrix/",
                                                geno_file, sep = ""),
                                         header = T) %>%
         tidyr::unite("QTL",c(CHROM, POS), sep = ":", remove = F) %>%
         dplyr::filter(QTL %in% effects$QTL) %>%
         dplyr::select(-CHROM, -POS) %>%
         tidyr::pivot_longer(cols = !c(QTL, REF, ALT), names_to = "strain", values_to = "allele") %>%
         dplyr::full_join(.,effects) %>%
         dplyr::full_join(.,phenos) %>%
         dplyr::group_by(QTL) %>%
         tidyr::nest()
      
   } else {
      
      effects.2 <- effects %>%
         dplyr::filter(QTL %in% complete.effects$QTL) %>%
         droplevels()
      
      genos.effects <- data.table::fread(paste(baseDir,
                                               args[1],
                                               "/Genotype_Matrix/",
                                               geno_file, sep = ""),
                                         header = T) %>%
         tidyr::unite("QTL",c(CHROM, POS), sep = ":", remove = F) %>%
         dplyr::filter(QTL %in% effects.2$QTL) %>%
         dplyr::select(-CHROM, -POS) %>%
         tidyr::pivot_longer(cols = !c(QTL, REF, ALT), names_to = "strain", values_to = "allele") %>%
         dplyr::left_join(.,effects.2) %>%
         dplyr::left_join(.,phenos) %>%
         dplyr::group_by(QTL) %>%
         tidyr::nest()
   }
   # ("-1" = "REF", "1" = "ALT")
   
   
   # Simulated Variance Explained
   var.exp <- function(data,QTL){
      aov.out <- summary(aov(data = data, trait.value ~ allele + 1))
      SST <- sum(aov.out[[1]][[2]])
      SSallele <- aov.out[[1]][[2]][[1]]
      simulated.variance.exp <- SSallele/SST
      data.frame(QTL, simulated.variance.exp)
   }
   simQTL.variance.explained <- purrr::map2(genos.effects$data, 
                                            genos.effects$QTL, 
                                            var.exp) %>%
      Reduce(rbind,.)
   
   effects <- effects %>%
      dplyr::full_join(.,simQTL.variance.explained) %>%
      dplyr::rename(Simulated.QTL.VarExp = simulated.variance.exp)
      
   # safe.aggregate <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
   #                                                                               paste(x,"processed_aggregate_mapping.tsv",sep = "_"), sep = "/"),header = T)}, 
   #                                      otherwise = "No successful mapping matching simulation parameters :( ")
   
   safe.lmm.exact.inbred <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                   paste(x,"processed_LMM-EXACT-INBRED_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                          otherwise = "No successful mapping matching simulation parameters :( ")
   safe.lmm.exact.inbred.pca <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                   paste(x,"processed_LMM-EXACT-INBRED_PCA_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                          otherwise = "No successful mapping matching simulation parameters :( ")   
   safe.lmm.exact.loco <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                    paste(x,"processed_LMM-EXACT-LOCO_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                          otherwise = "No successful mapping matching simulation parameters :( ")
   safe.lmm.exact.loco.pca <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                    paste(x,"processed_LMM-EXACT-LOCO_PCA_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                          otherwise = "No successful mapping matching simulation parameters :( ")
   
   mapping.lmm.exact.inbred <- safe.lmm.exact.inbred(x)
   mapping.lmm.exact.inbred.pca <- safe.lmm.exact.inbred.pca(x)
   mapping.lmm.exact.loco <- safe.lmm.exact.loco(x)
   mapping.lmm.exact.loco.pca <- safe.lmm.exact.loco.pca(x)
   
   mappings <- list("LMM-EXACT-INBRED" = mapping.lmm.exact.inbred$result, "LMM-EXACT-INBRED-PCA" = mapping.lmm.exact.inbred.pca$result, 
                     "LMM-EXACT-LOCO" = mapping.lmm.exact.loco$result, "LMM-EXACT-LOCO-PCA" = mapping.lmm.exact.loco.pca$result) 
   #purrr::set_names(c("LMM-EXACT-INBRED", "LMM-EXACT-INBRED-PCA", "LMM-EXACT-LOCO", "LMM-EXACT-LOCO-PCA"))
   #print(mappings)
   #print(mappings)
   #Function to process all processed mapping.tsv 
   generate_sim_data <- function(map_obj){
  ###
  # map obj - a mapping data table read in safely with purr, names of the list are the Algorithm ID
  # 
  # 
  ###


  #Check if map obj was sucessfully read-in 
   if(is.character(mapping.lmm.exact.inbred$result)){ 
      lmm.exact.inbred <- c("No mapping for simulation parameters")} 
      else {

         peak.info <- map_obj %>%
            dplyr::filter(!is.na(peak_id)) %>%
            dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, var.exp) %>%
            dplyr::filter(!duplicated(.)) %>%
            dplyr::mutate(detected.peak = marker)
         
         simulated.mapping.results.scores <- map_obj %>%
            dplyr::rename(QTL = marker) %>%
            dplyr::filter(QTL %in% effects$QTL) %>%
            dplyr::select(QTL, log10p, aboveBF)
         
         effects.scores <- effects  %>%
            dplyr::full_join(., simulated.mapping.results.scores) %>%
            dplyr::filter(!duplicated(QTL),
                        !is.na(log10p))
         
         peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
                                          ranges = IRanges::IRanges(start = peak.info$startPOS,
                                                                  end = peak.info$endPOS),
                                          peakPOS = peak.info$peakPOS,
                                          detected.peak = peak.info$detected.peak)
         real.effects <- GenomicRanges::GRanges(seqnames = effects.scores$CHROM,
                                                ranges = IRanges::IRanges(start = as.numeric(effects.scores$POS),
                                                                           end = as.numeric(effects.scores$POS)),
                                                QTL = effects.scores$QTL)
         
         overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
            as.data.frame() %>%
            dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
            `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
            dplyr::right_join(., peak.info) %>%
            dplyr::mutate(QTL = if_else(is.na(QTL),
                                       true = detected.peak,
                                       false = QTL)) %>%
            dplyr::rename(interval.log10p = log10p,
                        interval.var.exp = var.exp,
                        interval.Frequency = AF1) %>%
            dplyr::select(-c(CHROM, marker, POS))
         
         all.QTL <- data.frame(c(effects.scores$QTL, overlap$QTL)) %>%
            `colnames<-`(c("QTL")) %>%
            dplyr::filter(!duplicated(QTL)) %>%
            dplyr::mutate(QTL = as.character(QTL),
                        Simulated = (QTL %in% effects.scores$QTL),
                        Detected = (QTL %in% overlap$QTL)) %>%
            dplyr::full_join(.,effects.scores, by = "QTL") %>%
            dplyr::full_join(.,overlap, by = "QTL") %>%
            dplyr::mutate(
                        #algorithm = algorithm_id ,
                        top.hit = QTL == detected.peak,
                        sim = x)
         
         all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
         all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
         
         return(all.QTL)

         }
         
      }  
   all_algorithims <- plyr::ldply(mappings, generate_sim_data, .id = "algorithm")
   return(all_algorithims)

   #### AGGREGATE ####
   # mapping.aggregate <- safe.aggregate(x)
   
   # if(is.character(mapping.aggregate$result)){
   #    all.QTL <- c("No mapping for simulation parameters")
   # } else {
   #    # QTL Interval Information
   #    peak.info <- mapping.aggregate$result %>%
   #       dplyr::filter(!is.na(peak_id)) %>%
   #       dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, var.exp, algorithm) %>%
   #       dplyr::filter(!duplicated(.)) %>%
   #       dplyr::mutate(detected.peak = marker)
   #    
   #    # Simulated Variant Mapping Results
   #    simulated.mapping.results.scores <- mapping.aggregate$result %>%
   #       dplyr::rename(QTL = marker) %>%
   #       dplyr::filter(QTL %in% effects$QTL) %>%
   #       dplyr::select(QTL, log10p, aboveBF)
   #    
   #    # Join Interval Metrics with Simulated Variant Mapping Results
   #    effects.scores <- effects %>%
   #       dplyr::full_join(., simulated.mapping.results.scores) %>%
   #       dplyr::filter(!duplicated(QTL),
   #                     !is.na(log10p))
   #    
   #    # Genomic Ranges: Interval Information
   #    peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
   #                                    ranges = IRanges::IRanges(start = peak.info$startPOS, 
   #                                                              end = peak.info$endPOS),
   #                                    peakPOS = peak.info$peakPOS,
   #                                    detected.peak = peak.info$detected.peak)
   #    
   #    # Genomic Ranges: Simulated Variants
   #    real.effects <- GenomicRanges::GRanges(seqnames = effects.scores$CHROM,
   #                                           ranges = IRanges::IRanges(start = as.numeric(effects.scores$POS), 
   #                                                                     end = as.numeric(effects.scores$POS)),
   #                                           QTL = effects.scores$QTL)
   #    
   #    overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
   #       as.data.frame() %>%
   #       dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
   #       `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
   #       dplyr::right_join(., peak.info) %>%
   #       dplyr::mutate(QTL = if_else(is.na(QTL), 
   #                                   true = detected.peak, 
   #                                   false = QTL)) %>%
   #       dplyr::rename(interval.log10p = log10p,
   #                     interval.var.exp = var.exp,
   #                     interval.Frequency = AF1) %>%
   #       dplyr::select(-c(CHROM, marker, POS))
   #    
   #    all.QTL <- data.frame(c(effects.scores$QTL, overlap$QTL)) %>%
   #       `colnames<-`(c("QTL")) %>%
   #       dplyr::filter(!duplicated(QTL)) %>%
   #       dplyr::mutate(QTL = as.character(QTL),
   #                     Simulated = (QTL %in% effects.scores$QTL),
   #                     Detected = (QTL %in% overlap$QTL)) %>%
   #       dplyr::full_join(.,effects.scores, by = "QTL") %>%
   #       dplyr::full_join(.,overlap, by = "QTL") %>%
   #       dplyr::mutate(algorithm = "MIXED",
   #                     top.hit = QTL == detected.peak,
   #                     sim = x)
   #    
   #    all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
   #    all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
   #    result.list[[1]] <- all.QTL
   # }
   #####
   
   #### GCTA: LMM-EXACT w/ INBRED GRM ####


 }


if(pop_designation == TRUE){
   simulation.metrics.list <- purrr::map(iterations, simulation.metrics,  pop_designation = TRUE)
   check_col <- function(x){ncol(x)<= 2}
   filtered.simulation.metrics.list <- discard(simulation.metrics.list, check_col)
   simulation.metrics.df <- Reduce(rbind, filtered.simulation.metrics.list)
   save(simulation.metrics.df, file = paste("NemaScan_Performance",args[1],today,"RData", sep = "."))
}
if(pop_designation == FALSE){
   
   simulation.metrics.list <- purrr::map(iterations, simulation.metrics, pop_designation = FALSE)
   check_col <- function(x){ncol(x)<= 2}
   filtered.simulation.metrics.list <- discard(simulation.metrics.list, check_col)
   simulation.metrics.df <- Reduce(rbind, filtered.simulation.metrics.list)
   save(simulation.metrics.df, file = paste("NemaScan_Performance",args[1],today,"RData", sep = "."))
}

















### GRAVEYARD ###
# safe.lmm.exact <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
#                                                                                 paste(x,"processed_LMM_EXACT_mapping.tsv",sep = "_"), sep = "/"),header = T)}, 
#                                        otherwise = "No successful mapping matching simulation parameters :( ")
# 
# safe.lmm.exact.loco.inbred <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
#                                                                                 paste(x,"processed_LMM_EXACT_INBRED_LOCO_mapping.tsv",sep = "_"), sep = "/"),header = T)}, 
#                                        otherwise = "No successful mapping matching simulation parameters :( ")
# safe.EMMA <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
#                                                                     paste(x,"processed_mapping.tsv",sep = "_"), sep = "/"),header = T) %>%
#       dplyr::mutate(marker = gsub(marker, pattern = "_", replacement = ":"))}, 
#       otherwise = "No successful mapping matching simulation parameters :( ")

#### GCTA: LMM-EXACT ####
# mapping.lmm.exact <- safe.lmm.exact(x)
# if(is.character(mapping.lmm.exact$result)){
#    lmm.exact <- c("No mapping for simulation parameters")
# } else {
#    peak.info <- mapping.lmm.exact$result %>%
#       dplyr::mutate(causal.variant = as.factor(marker %in% effects$QTL)) %>%
#       dplyr::filter(!is.na(peak_id)) %>%
#       dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, causal.variant, var.exp) %>%
#       dplyr::filter(!duplicated(.)) %>%
#       dplyr::mutate(detected.peak = marker)
#    
#    simulated.mapping.results.scores <- mapping.lmm.exact$result %>%
#       dplyr::rename(QTL = marker) %>%
#       dplyr::filter(QTL %in% effects$QTL) %>%
#       dplyr::select(QTL, log10p)
#    effects.scores <- effects %>%
#       dplyr::full_join(., simulated.mapping.results.scores) %>%
#       dplyr::filter(!duplicated(QTL))
#    
#    peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
#                                    ranges = IRanges::IRanges(start = peak.info$startPOS, 
#                                                              end = peak.info$endPOS),
#                                    peakPOS = peak.info$peakPOS,
#                                    detected.peak = peak.info$detected.peak)
#    real.effects <- GenomicRanges::GRanges(seqnames = effects.scores$CHROM,
#                                           ranges = IRanges::IRanges(start = as.numeric(effects.scores$POS), 
#                                                                     end = as.numeric(effects.scores$POS)),
#                                           QTL = effects.scores$QTL)
#    
#    overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
#       as.data.frame() %>%
#       dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
#       `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
#       dplyr::right_join(., peak.info) %>%
#       dplyr::mutate(QTL = if_else(is.na(QTL), true = detected.peak, false = QTL))
#    
#    all.QTL <- data.frame(c(effects.scores$QTL, overlap$QTL)) %>%
#       `colnames<-`(c("QTL")) %>%
#       dplyr::filter(!duplicated(QTL)) %>%
#       dplyr::mutate(QTL = as.character(QTL),
#                     Simulated = (QTL %in% effects.scores$QTL),
#                     Detected = (QTL %in% overlap$QTL)) %>%
#       dplyr::full_join(.,effects.scores, by = "QTL") %>%
#       dplyr::full_join(.,overlap, by = "QTL") %>%
#       dplyr::mutate(sim = x,
#                     var.exp = as.character(var.exp),
#                     BETA = as.character(BETA),
#                     log10p = if_else(condition = is.na(log10p.y), 
#                                      true = log10p.x, 
#                                      false = log10p.y),
#                     algorithm = "LMM-EXACT") %>%
#       dplyr::select(-CHROM.y, -marker, -POS.y, -AF1, -log10p.x, -log10p.y)
#    
#    all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
#    all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
#    all.QTL.lmm.exact <- all.QTL
#    result.list[[2]] <- all.QTL.lmm.exact
# }
#### rrBLUP: EMMAx or EMMA ####
# mapping.EMMA <- safe.EMMA(x)
# if(is.character(mapping.EMMA$result)){
#    EMMA <- c("No mapping for simulation parameters")
# } else {
#    peak.info <- mapping.EMMA$result %>%
#       dplyr::mutate(causal.variant = as.factor(marker %in% effects$QTL)) %>%
#       dplyr::filter(!is.na(peak_id)) %>%
#       dplyr::select(CHROM, marker, POS, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, causal.variant, var.exp) %>%
#       dplyr::filter(!duplicated(.)) %>%
#       dplyr::mutate(detected.peak = marker)
#    
#    simulated.mapping.results.scores <- mapping.EMMA$result %>%
#       dplyr::rename(QTL = marker) %>%
#       dplyr::filter(QTL %in% effects$QTL) %>%
#       dplyr::select(QTL, log10p)
#    effects.scores <- effects %>%
#       dplyr::full_join(., simulated.mapping.results.scores) %>%
#       dplyr::filter(!duplicated(QTL))
#    
#    peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
#                                    ranges = IRanges::IRanges(start = peak.info$startPOS, 
#                                                              end = peak.info$endPOS),
#                                    peakPOS = peak.info$peakPOS,
#                                    detected.peak = peak.info$detected.peak)
#    real.effects <- GenomicRanges::GRanges(seqnames = effects.scores$CHROM,
#                                           ranges = IRanges::IRanges(start = as.numeric(effects.scores$POS), 
#                                                                     end = as.numeric(effects.scores$POS)),
#                                           QTL = effects.scores$QTL)
#    overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
#       as.data.frame() %>%
#       dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
#       `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
#       dplyr::right_join(., peak.info) %>%
#       dplyr::mutate(QTL = if_else(is.na(QTL), true = detected.peak, false = QTL))
#    
#    all.QTL <- data.frame(c(effects.scores$QTL, overlap$QTL)) %>%
#       `colnames<-`(c("QTL")) %>%
#       dplyr::filter(!duplicated(QTL)) %>%
#       dplyr::mutate(QTL = as.character(QTL),
#                     Simulated = (QTL %in% effects.scores$QTL),
#                     Detected = (QTL %in% overlap$QTL),
#                     BETA = NA) %>%
#       dplyr::full_join(.,effects.scores, by = "QTL") %>%
#       dplyr::full_join(.,overlap, by = "QTL") %>%
#       dplyr::mutate(sim = x,
#                     var.exp = as.character(var.exp),
#                     BETA = as.character(BETA),
#                     log10p = if_else(condition = is.na(log10p.y), 
#                                      true = log10p.x, 
#                                      false = log10p.y),
#                     algorithm = "EMMA") %>%
#       dplyr::select(-CHROM.y, -marker, -POS.y, -log10p.x, -log10p.y)
#    
#    all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
#    all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
#    all.QTL.EMMA <- all.QTL
#    result.list[[5]] <- all.QTL.EMMA
# }
#### GCTA: LMM-EXACT w/ INBRED GRM, LOCO ####
# mapping.lmm.exact.inbred.loco <- safe.lmm.exact.loco.inbred(x)
# 
# if(is.character(mapping.lmm.exact.inbred.loco$result)){
#    lmm.exact.inbred <- c("No mapping for simulation parameters")
# } else {
#    peak.info <- mapping.lmm.exact.inbred.loco$result %>%
#       dplyr::mutate(causal.variant = as.factor(marker %in% effects$QTL)) %>%
#       dplyr::filter(!is.na(peak_id)) %>%
#       dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, causal.variant, var.exp) %>%
#       dplyr::filter(!duplicated(.)) %>%
#       dplyr::mutate(detected.peak = marker)
#    
#    simulated.mapping.results.scores <- mapping.lmm.exact.inbred.loco$result %>%
#       dplyr::rename(QTL = marker) %>%
#       dplyr::filter(QTL %in% effects$QTL) %>%
#       dplyr::select(QTL, log10p)
#    
#    effects.scores <- effects %>%
#       dplyr::full_join(., simulated.mapping.results.scores) %>%
#       dplyr::filter(!duplicated(QTL))
#    
#    peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
#                                    ranges = IRanges::IRanges(start = peak.info$startPOS, 
#                                                              end = peak.info$endPOS),
#                                    peakPOS = peak.info$peakPOS,
#                                    detected.peak = peak.info$detected.peak)
#    real.effects <- GenomicRanges::GRanges(seqnames = effects.scores$CHROM,
#                                           ranges = IRanges::IRanges(start = as.numeric(effects.scores$POS), 
#                                                                     end = as.numeric(effects.scores$POS)),
#                                           QTL = effects.scores$QTL)
#    
#    overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
#       as.data.frame() %>%
#       dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
#       `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
#       dplyr::right_join(., peak.info) %>%
#       dplyr::mutate(QTL = if_else(is.na(QTL), true = detected.peak, false = QTL))
#    
#    all.QTL <- data.frame(c(effects.scores$QTL, overlap$QTL)) %>%
#       `colnames<-`(c("QTL")) %>%
#       dplyr::filter(!duplicated(QTL)) %>%
#       dplyr::mutate(QTL = as.character(QTL),
#                     Simulated = (QTL %in% effects.scores$QTL),
#                     Detected = (QTL %in% overlap$QTL)) %>%
#       dplyr::full_join(.,effects.scores, by = "QTL") %>%
#       dplyr::full_join(.,overlap, by = "QTL") %>%
#       dplyr::mutate(sim = x,
#                     var.exp = as.character(var.exp),
#                     BETA = as.character(BETA),
#                     log10p = if_else(condition = is.na(log10p.y), 
#                                      true = log10p.x, 
#                                      false = log10p.y),
#                     algorithm = "LMM-EXACT-LOCO-INBRED") %>%
#       dplyr::select(-CHROM.y, -marker, -POS.y, -AF1, -log10p.x, -log10p.y)
#    
#    all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
#    all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
#    all.QTL.lmm.exact.inbred.loco <- all.QTL
#    result.list[[3]] <- all.QTL.lmm.exact.inbred.loco
# }
