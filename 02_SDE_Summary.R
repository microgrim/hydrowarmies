# 02_SDE_Summary.R

#this script uses the significant results from DESeq
# to compare the number of significantly differentially expressed genes
# in three comparison types:
# across temperatures within a vent
# between vents at a specific temperature
# within vents at specific temperatures between years
# preliminarily, we also explored effects of using as a baseline for expression either:
# the 13C-heavy fraction only, or using the 13C-heavy + 12C-light fractions 
#this script ultimately summarizes the results to be plotted


# Resource allocation time ------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_resource_allocation.R"))){
  cat("Preparing resource allocations.")
  suppressWarnings(
    source(paste0(getwd(), "/", "00_resource_allocation.R"), local = FALSE,
           verbose = getOption("verbose"), prompt.echo = getOption("prompt"),
           echo = FALSE),
    classes = "warnings")
  
} else {
  cat("You have to manually load these packages and set paths.")
  # #packages needed for this script:
  
  library(tidyverse)
  library(data.table)
  library(viridis)
  library(BiocManager)
  library(BiocParallel)
  library(DESeq2)
  library(cli)
  library(furrr)
  library(progressr)
  library(UpSetR)
  library(RVenn)
  library(ggVennDiagram)
  
  # 
  # 
}
try(f_projectpath())

# Do you needto set up a server or remote path? ---------------------------


try(f_remotepath())
try(f_serverpath())

if(!exists("projectpath", envir = .GlobalEnv)){
  projectpath <- getwd()
}




# Load custom functions for parsing ---------------------------------------

if(file.exists(paste0(getwd(), "/", "00_custom_functions.R"))){
  cat("Loading custom functions.")
  source(paste0(getwd(), "/", "00_custom_functions.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
}

# P_venn_upset <- function(vennlist, mdata_filter, prefix){
#   #vennlist: your list of simplified character vectors of genes/targets, named by contrast condition
#   #mdata_filter: optional selector of your data
#   #provide the filter a list of named vectors
#   # metadata_filter <- list(`Fraction` = c("12L", "13H"),
#   #                         `Temperature` = c("30; 55"),
#   #                         `baseline` = c("axial_temp_12L", "axial_temp_13H", "axial_temp_both"))
#   if(!missing(prefix)){
#     prefix <- rlang::parse_expr(prefix)
#   } else {
#     prefix <- "default"
#   }
#   vennlist <- vennlist
#   mdata_filter <- mdata_filter
#   
#   #params filter
#   if(exists("labels_vennlist", inherits = FALSE)){
#     rm(labels_vennlist)}
#   if(!is.null(mdata_filter)){
#     labels_idx <- imap(mdata_filter,
#                        ~vennlist[.x]) %>%
#       Filter(Negate(function(x) is.null(unlist(x))), .) %>%
#       purrr::list_flatten(., name_spec = "{inner}")
#     
#     for(i in seq(length(labels_idx))){
#       namevar <- names(labels_idx[i])
#       for(j in seq(pluck_depth(labels_idx[i]) -1)){
#         temp_list <- labels_idx[[namevar]][j] %>%
#           setNames(., paste0(namevar, "_", names(.))) %>%
#           Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
#           purrr::list_flatten(.) %>%
#           purrr::list_flatten(.)
#         if(exists("labels_vennlist")){
#           labels_vennlist <- c(labels_vennlist, temp_list)
#         } else {
#           labels_vennlist <- temp_list
#         }
#       }
#     }
#     keep <- mdata_filter %>%
#       map(., ~grepl(paste0(.x, collapse = "|"), names(labels_vennlist))) %>%
#       purrr::reduce(`&`)
#     labels_vennlist <- labels_vennlist[keep]
#   } else {
#     labels_vennlist <- map(vennlist,
#                            ~.x %>%
#                              purrr::list_flatten(., name_spec = "{outer}_{inner}"))
#   }
#   
#   params_list <- list("anemone" = grep(paste0(c("anemone"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#                       "marker113" = grep(paste0(c("marker113"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#                       "marker33" = grep(paste0(c("marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#                       "anemone \u2229 marker113" = grep(paste0(c("anemone", "marker113"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#                       "anemone \u2229 marker33" = grep(paste0(c("anemone", "marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#                       "marker113 \u2229 marker33" = grep(paste0(c("marker113", "marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#                       "anemone \u2229 marker113 \u2229 marker33" = grep(paste0(c("anemone", "marker113", "marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE))
#   
#   # params_list <- list("13H" = grep(paste0(c("temp_13H"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "12L" = grep(paste0(c("temp_12L"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "13H+12L" = grep(paste0(c("temp_both"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "13H \u2229 12L" = grep(paste0(c("temp_13H", "temp_12L"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "13H \u2229 13H+12L" = grep(paste0(c("temp_13H", "temp_both"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "12L \u2229 13H+12L" = grep(paste0(c("temp_12L", "temp_both"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "12L \u2229 13H \u2229 13H+12L" = grep(paste0(c("temp_13H", "temp_both", "temp_12L"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE))
#   # 
#   #13H: red
#   #12L: blue
#   #12L+13H: green
#   #13H vs 12L: purple
#   #13H vs 12L+13H: brown
#   #12L vs 12L+13H: teal
#   #all: white/black
#   
#   
#   g <- UpSetR::upset(UpSetR::fromList(labels_vennlist),
#                      main.bar.color = "black",
#                      mainbar.y.label = "Intersection size",
#                      sets.x.label = "Set size",
#                      sets = c(names(labels_vennlist)),
#                      empty.intersections = "keep",
#                      order.by = "degree", decreasing = F,
#                      queries = list(list(query = intersects, params = params_list[1], color = "red", active = T, query.name = paste(names(params_list[1]))),
#                                     list(query = intersects, params = params_list[2], color = "blue", active = T, query.name = paste(names(params_list[2]))),
#                                     list(query = intersects, params = params_list[3], color = "darkgreen", active = T, query.name = paste(names(params_list[3]))),
#                                     list(query = intersects, params = params_list[4], color = "purple", active = T, query.name = paste(names(params_list[4]))),
#                                     list(query = intersects, params = params_list[5], color = "brown", active = T, query.name = paste(names(params_list[5]))),
#                                     list(query = intersects, params = params_list[6], color = "turquoise", active = T, query.name = paste(names(params_list[6]))),
#                                     list(query = intersects, params = params_list[7], color = "black", active = T, query.name = paste(names(params_list[7])))),
#                      keep.order = F,
#                      query.legend = "bottom")
#   assign(paste0(prefix, "_upset_plot"), g, envir = .GlobalEnv)
# }
# 
# #which genes are picked up only when using either 12L or 13H+12L as baseline?
# F_venn_setdiffs <- function(vennlist, mdata_filter, prefix) {
#   if(!missing(prefix)){
#     prefix <- rlang::parse_expr(prefix)
#   } else {
#     prefix <- "default"
#   }
#   mdata_filter <- get("mdata_filter", inherits = TRUE)
#   title_temp <- mdata_filter[["Temperature"]]
#   vennlist <- get("vennlist", inherits = TRUE)
#   
#   #params filter
#   if(exists("labels_vennlist", inherits = FALSE)){
#     rm(labels_vennlist)}
#   if(exists("keep", inherits = FALSE)){
#     rm(keep)}
#   
#   if(exists("mdata_filter", inherits = TRUE)){
#     labels_idx <- imap(mdata_filter,
#                        ~vennlist[.x]) %>%
#       Filter(Negate(function(x) is.null(unlist(x))), .) %>%
#       purrr::list_flatten(., name_spec = "{inner}")
#     
#     for(i in seq(length(labels_idx))){
#       namevar <- names(labels_idx[i])
#       for(j in seq(pluck_depth(labels_idx[i]) -1)){
#         temp_list <- labels_idx[[namevar]][j] %>%
#           setNames(., paste0(namevar, "_", names(.))) %>%
#           Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
#           purrr::list_flatten(.) %>%
#           purrr::list_flatten(.)
#         if(exists("labels_vennlist")){
#           labels_vennlist <- c(labels_vennlist, temp_list)
#         } else {
#           labels_vennlist <- temp_list
#         }
#         rm(temp_list)
#       }
#     }
#     keep <- mdata_filter %>%
#       map(., ~grepl(paste0(.x, collapse = "|"), names(labels_vennlist))) %>%
#       purrr::reduce(`&`)
#     labels_vennlist <- labels_vennlist[keep]
#     rm(labels_idx)
#   } else {
#     labels_vennlist <- map(vennlist,
#                            ~.x %>%
#                              # purrr::list_flatten(., name_spec = "{inner}_{outer}") %>%
#                              purrr::list_flatten(., name_spec = "{outer}_{inner}"))
#   }
#   
#   if(any(grepl(paste0(vents, collapse = "|"), mdata_filter))){
#     idx <- map2(names(labels_vennlist), mdata_filter["Vents"], 
#                 ~grepl(paste0(.y, collapse = "|"), .x, ignore.case = TRUE, perl = TRUE) %>%
#                   purrr::reduce(`&`)) %>% 
#       unlist(., recursive = FALSE, use.names = FALSE)
#     labels_vennlist <- labels_vennlist[idx]
#   }
#   
#   # for(i in c("anemone")){
#   for(i in mdata_filter[["Vent"]]){
#     namevar <- i
#     idx <- grepl(i, names(labels_vennlist))
#     temp_list <- labels_vennlist[idx]
#     params_list <- list("13H \u2229 12L" = grep(paste0(c("temp_13H", "temp_12L"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE),
#                         "13H \u2229 13H+12L" = grep(paste0(c("temp_13H", "temp_both"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE),
#                         "12L \u2229 13H+12L" = grep(paste0(c("temp_12L", "temp_both"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE),
#                         "12L \u2229 13H \u2229 13H+12L" = grep(paste0(c("temp_13H", "temp_both", "temp_12L"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE))
#     params_list2 <- list("13H \u2284 12L" = grep(paste0(c("temp_13H", "temp_12L"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE),
#                          "13H \u2284 13H+12L" = grep(paste0(c("temp_13H", "temp_both"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE),
#                          "12L \u2284 13H+12L" = grep(paste0(c("temp_12L", "temp_both"), collapse = "|"), names(temp_list), value = TRUE, ignore.case = TRUE))
#     temp_venn_list1 <- params_list %>%
#       map(., ~temp_list[.x] %>%
#             purrr::reduce(intersect))
#     temp_venn_list2 <- params_list2 %>%
#       map(., ~temp_list[.x] %>%
#             purrr::reduce(setdiff))
#     temp_venn_list <- list(temp_venn_list1, temp_venn_list2) %>%
#       purrr::list_flatten(.)
#     
#     assign(paste(prefix, title_temp, namevar, "setdiffs", sep = "_"), temp_venn_list, inherits = TRUE, envir = .GlobalEnv)
#   }
# }

# Global labellers/renamers -----------------------------------------------


if(file.exists(paste0(getwd(), "/", "00_global_labellers.R"))){
  cat("Loading project-wide labellers and lookups.")
  source(paste0(getwd(), "/", "00_global_labellers.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 



# Import sample metadata --------------------------------------------------

if(!exists("sample_metadata", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/processed/", "sample_metadata.rds"))){
    cli::cli_alert_info("Reading in sample metadata object.")
    sample_metadata <- readr::read_rds(paste0(projectpath, "/data/processed/", "sample_metadata.rds"))
  } else {
    cli::cli_abort("Sample metadata has not been processed in step 1.")  
  }
}


# Import processed DESeq results ------------------------------------------

processed_deseq_res_idx <- c("coa_axial_temp_13H",
                             "coa_axial_temp_bothfrac", #supplemental: vent-specific, thermal comparison using 13C-heavy and 12C-light as baseline
                             "coa_axial_eachtemp_13H", #temp-specific, intervent comparison, consolidates coa_axial_*_13H_res_list
                             "coa_axial_eachtemp_bothfrac", #supplemental: consolidates coa_axial_*_bothfrac_res_list
                             "coa_axial_eachtemp_13H_vent_year") #temp- and vent-specific, temporal comparison, consolidates coa_axial_*_13H_vent_year_res_list

for(i in processed_deseq_res_idx){
  namevar <- i
  if(!exists(paste0(namevar, "_res_list"), envir = .GlobalEnv)){
      if(file.exists(paste0(projectpath, "/intermediate/deseq/", namevar, "_res_list", ".rds"))){
      cli::cli_alert_info("Reading in DESeq results list: {namevar}.")
      temp_list <- readr::read_rds(file = paste0(projectpath, "/intermediate/deseq/", namevar, "_res_list", ".rds"))
      assign(paste0(namevar, "_res_list"), temp_list, envir = .GlobalEnv, inherits = TRUE)
      rm(temp_list)
    }
    else {
      cli::cli_alert_info("DESeq results list does not exist with this name: {namevar}. Please process through DESeq.")
    }
  }
  rm(namevar)
}


# Plot number of SDE genes ------------------------------------------------

#working on just 13H as baseline now:
coa_axial_SDE.df <- grep("bothfrac", processed_deseq_res_idx, invert = TRUE, value = TRUE) %>%
  paste0(., "_res_list") %>%
  lapply(., get0) %>%
  setNames(., grep("bothfrac", processed_deseq_res_idx, invert = TRUE, value = TRUE)) %>%
  purrr::list_flatten(., name_spec =  "{inner}") %>%
  dplyr::bind_rows(., .id = "model")

#the vent-specific, thermal comparisons
coa_axial_temp_SDE_summary <- coa_axial_SDE.df %>%
  dplyr::filter(!grepl("2013|2014", contrast)) %>%
  dplyr::filter(!grepl("eachtemp", model)) %>%
  dplyr::select(c(model, contrast, gene_symbol, padj)) %>%
  dplyr::mutate(contrast = factor(contrast)) %>%
  tidyr::pivot_longer(., cols = !c("model", "contrast", "gene_symbol"),
                      names_to = "metric",
                      values_to = "padj") %>%
  tidyr::drop_na(.) %>%
  dplyr::filter(padj < set.alpha) %>%
  dplyr::summarise(Num_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
  tidyr::complete(contrast, model, fill = list(Num_SDE = 0)) %>%
  # dplyr::arrange(desc(Num_SDE)) %>%
  #adding in the "almost-significant" results
  dplyr::full_join(., 
                   (coa_axial_SDE.df %>%
                      dplyr::filter(!grepl("2013|2014", contrast)) %>%
                      dplyr::filter(!grepl("eachtemp", model)) %>%
                      dplyr::select(c(model, contrast, gene_symbol, pvalue)) %>%
                      dplyr::mutate(contrast = factor(contrast)) %>%
                      tidyr::pivot_longer(., cols = !c(model, "contrast", "gene_symbol"),
                                          names_to = "metric",
                                          values_to = "pvalue") %>%
                      tidyr::drop_na(.) %>%
                      dplyr::filter(pvalue < set.alpha) %>%
                      dplyr::summarise(Num_almost_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
                      tidyr::complete(contrast, model, fill = list(Num_almost_SDE = 0)) %>%
                      droplevels), by = join_by(contrast, model)) %>%
  droplevels

#temperature-specific, intervent comparisons
coa_axial_vent_SDE_summary <- coa_axial_SDE.df %>%
  dplyr::filter(!grepl("2013|2014", contrast)) %>%
  dplyr::filter(grepl("eachtemp", model)) %>%
  dplyr::select(c(model, contrast, gene_symbol, padj)) %>%
  dplyr::mutate(contrast = factor(contrast)) %>%
  tidyr::pivot_longer(., cols = !c("model", "contrast", "gene_symbol"),
                      names_to = "metric",
                      values_to = "padj") %>%
  tidyr::drop_na(.) %>%
  dplyr::filter(padj < set.alpha) %>%
  dplyr::summarise(Num_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
    tidyr::complete(contrast, model, fill = list(Num_SDE = 0)) %>%
  dplyr::arrange(desc(Num_SDE)) %>%
  #adding in the "almost-significant" results
  dplyr::full_join(., 
                   (coa_axial_SDE.df %>%
                      dplyr::filter(!grepl("2013|2014", contrast)) %>%
                      dplyr::filter(grepl("eachtemp", model)) %>%
                      dplyr::select(c(model, contrast, gene_symbol, pvalue)) %>%
                      dplyr::mutate(contrast = factor(contrast)) %>%
                      tidyr::pivot_longer(., cols = !c(model, "contrast", "gene_symbol"),
                                          names_to = "metric",
                                          values_to = "pvalue") %>%
                      tidyr::drop_na(.) %>%
                      dplyr::filter(pvalue < set.alpha) %>%
                      dplyr::summarise(Num_almost_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
                      tidyr::complete(contrast, model, fill = list(Num_almost_SDE = 0)) %>%
                      droplevels), by = join_by(contrast, model)) %>%
  droplevels

#vent- and temperature-specific, temporal comparisons
coa_axial_year_SDE_summary <- coa_axial_SDE.df %>%
  dplyr::filter(grepl("2013|2014", contrast)) %>%
  dplyr::select(c(model, contrast, gene_symbol, padj)) %>%
  dplyr::mutate(contrast = factor(contrast)) %>%
  tidyr::pivot_longer(., cols = !c("model", "contrast", "gene_symbol"),
                      names_to = "metric",
                      values_to = "padj") %>%
  tidyr::drop_na(.) %>%
  dplyr::filter(padj < set.alpha) %>%
  dplyr::summarise(Num_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
  tidyr::complete(contrast, model, fill = list(Num_SDE = 0)) %>%
  dplyr::arrange(desc(Num_SDE)) %>%
  #adding in the "almost-significant" results
  dplyr::full_join(., 
            (coa_axial_SDE.df %>%
               dplyr::filter(grepl("2013|2014", contrast)) %>%
               dplyr::select(c(model, contrast, gene_symbol, pvalue)) %>%
               dplyr::mutate(contrast = factor(contrast)) %>%
               tidyr::pivot_longer(., cols = !c(model, "contrast", "gene_symbol"),
                                   names_to = "metric",
                                   values_to = "pvalue") %>%
               tidyr::drop_na(.) %>%
               dplyr::filter(pvalue < set.alpha) %>%
               dplyr::summarise(Num_almost_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
               tidyr::complete(contrast, model, fill = list(Num_almost_SDE = 0)) %>%
               droplevels), by = join_by(contrast, model)) %>%
  droplevels


#export the table results:

STab_coa_axial_SDE_summary <- dplyr::bind_rows((coa_axial_temp_SDE_summary %>%
                                                  dplyr::mutate(Comparison = dplyr::case_when(grepl("vent_year", model) ~ "Temporal",
                                                                                              grepl("eachtemp", model) ~ "Spatial",
                                                                                              .default = "Thermal")) %>%
                                                  tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
                                                  dplyr::mutate(Year = stringr::str_extract(contrast, paste0(year_lookup, collapse = "|"))) %>%
                                                  dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
                                                  dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup), labels = vent_lookup)) %>%
                                                  dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
                                                                  stringr::str_remove(., "_")) %>%
                                                  dplyr::arrange(Vent, Temperature) %>%
                                                  dplyr::mutate(Hold = Vent,
                                                                Variable = Temperature) %>%
                                                  dplyr::select(Comparison, Hold, Variable, Num_SDE) %>%
                                                  dplyr::arrange(Hold, Variable)),
                                               (coa_axial_vent_SDE_summary %>%
                                                  dplyr::mutate(Comparison = dplyr::case_when(grepl("vent_year", model) ~ "Temporal",
                                                                                              grepl("eachtemp", model) ~ "Spatial",
                                                                                              .default = "Thermal")) %>%
                                                  tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
                                                  dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
                                                  dplyr::mutate(Year = stringr::str_extract(contrast, paste0(year_lookup, collapse = "|"))) %>%
                                                  dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup), labels = vent_lookup)) %>%
                                                  dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
                                                                  stringr::str_remove(., "_")) %>%
                                                  dplyr::arrange(Vent, Temperature) %>%
                                                  dplyr::mutate(Hold = Temperature,
                                                                Variable = Vent) %>%
                                                  dplyr::select(Comparison, Hold, Variable, Num_SDE) %>%
                                                  dplyr::arrange(Hold, Variable)),
                                               (coa_axial_year_SDE_summary %>%
                                                  dplyr::mutate(Comparison = dplyr::case_when(grepl("vent_year", model) ~ "Temporal",
                                                                                              grepl("eachtemp", model) ~ "Spatial",
                                                                                              .default = "Thermal")) %>%
                                                  tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
                                                  dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
                                                  dplyr::mutate(Year = stringr::str_extract(contrast, paste0(year_lookup, collapse = "|"))) %>%
                                                  dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup), labels = vent_lookup)) %>%
                                                  dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
                                                                  stringr::str_remove(., "_")) %>%
                                                  dplyr::arrange(Vent, Year, Temperature) %>%
                                                  dplyr::summarise(Num_SDE = sum(Num_SDE), .by = c(Comparison, Vent, Year)) %>%
                                                  dplyr::mutate(Hold = Vent,
                                                                Variable = Year) %>%
                                                  dplyr::distinct(Comparison, Hold, Variable, Num_SDE)) %>%
                                                 dplyr::arrange(Hold, Variable))
  
if(!file.exists(paste0(projectpath, "/output/tables/", "stab_10_coa_axial_SDE_genes_summary.tsv"))){
  readr::write_delim(STab_coa_axial_SDE_summary, 
                     paste0(projectpath, "/output/tables/", "stab_10-coa_axial_SDE_genes_summary", ".tsv"),
                     delim = "\t", col_names = TRUE, na = "NA")  
  
}



# Supplementary: baseline of 13H+12L --------------------------------------


#bonus:
##comparing to the baseline using 13H+12L:

coa_axial_SDE_bothfrac.df <- grep("bothfrac", processed_deseq_res_idx, invert = FALSE, value = TRUE) %>%
  paste0(., "_res_list") %>%
  lapply(., get0) %>%
  setNames(., grep("bothfrac", processed_deseq_res_idx, invert = FALSE, value = TRUE)) %>%
  purrr::list_flatten(., name_spec =  "{inner}") %>%
  dplyr::bind_rows(., .id = "model")


coa_axial_vent_SDE_baseline_summary <- coa_axial_SDE_bothfrac.df %>%
  dplyr::filter(!grepl("2013|2014", contrast)) %>%
  dplyr::filter(grepl("eachtemp", model)) %>%
  dplyr::select(c(model, contrast, gene_symbol, padj)) %>%
  dplyr::mutate(contrast = factor(contrast)) %>%
  tidyr::pivot_longer(., cols = !c("model", "contrast", "gene_symbol"),
                      names_to = "metric",
                      values_to = "padj") %>%
  tidyr::drop_na(.) %>%
  dplyr::filter(padj < set.alpha) %>%
  dplyr::summarise(Num_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
  tidyr::complete(contrast, model, fill = list(Num_SDE = 0)) %>%
  tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
  dplyr::select(sample, Num_SDE) %>%
  dplyr::rename(coa_axial_eachtemp_bothfrac = "Num_SDE") %>%
  dplyr::full_join(., coa_axial_vent_SDE_summary %>%
                     dplyr::select(contrast, Num_SDE) %>%
                     tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
                     dplyr::select(sample, Num_SDE) %>%
                     dplyr::rename(coa_axial_eachtemp_13H = "Num_SDE") %>%
                     droplevels,
                   by = join_by(sample)) %>%
  droplevels

coa_axial_temp_SDE_baseline_summary <- coa_axial_SDE_bothfrac.df %>%
  dplyr::filter(!grepl("2013|2014", contrast)) %>%
  dplyr::filter(!grepl("eachtemp", model)) %>%
  dplyr::select(c(model, contrast, gene_symbol, padj)) %>%
  dplyr::mutate(contrast = factor(contrast)) %>%
  tidyr::pivot_longer(., cols = !c("model", "contrast", "gene_symbol"),
                      names_to = "metric",
                      values_to = "padj") %>%
  tidyr::drop_na(.) %>%
  dplyr::filter(padj < set.alpha) %>%
  dplyr::summarise(Num_SDE = length(gene_symbol), .by = c(model, contrast)) %>%
  tidyr::complete(contrast, model, fill = list(Num_SDE = 0)) %>%
  tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
  dplyr::select(sample, Num_SDE) %>%
  dplyr::rename(coa_axial_temp_bothfrac = "Num_SDE") %>%
  dplyr::full_join(., coa_axial_temp_SDE_summary %>%
                     dplyr::select(contrast, Num_SDE) %>%
                     tidyr::separate_wider_delim(contrast,delim = " - ", names = c("sample", "baseline"), cols_remove = FALSE) %>%
                     dplyr::select(sample, Num_SDE) %>%
                     dplyr::rename(coa_axial_temp_13H = "Num_SDE") %>%
                     droplevels,
                   by = join_by(sample)) %>%
  droplevels


# Making barcharts for the summary ----------------------------------------

  temp_df1 <- coa_axial_temp_SDE_summary %>%
    dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
    dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup))) %>%
    dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
                    stringr::str_remove(., "_")) %>%
    dplyr::mutate(baseline_temperature = stringr::str_extract(contrast, "axial_30|axial_55|axial_80|axial_temp") %>%
                    stringr::str_remove_all(., "axial_") %>%
                    stringr::str_replace_all("temp", "all")) %>%
    dplyr::arrange(Vent, Temperature) %>%
    # dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)([[:digit:]]+)(.*)$", "\\2\\3\\4", contrast)) %>%
    dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)(.*)$", "\\2\\3", contrast)) %>%
    dplyr::mutate(baseline = factor(baseline, levels = names(contrast_lookup), labels = contrast_lookup)) %>%
    # dplyr::mutate(Almost_SDE = Num_almost_SDE - Num_SDE) %>%
    dplyr::mutate(label_y = dplyr::case_when( (Num_SDE < 10) ~ 10, 
                                              .default = (Num_SDE+10)*1.1)) %>%
    # .default = (Num_SDE+5)*0.9)) %>%
    droplevels
  
  g1 <- print(
    ggplot(data = temp_df1 %>%
             dplyr::filter(baseline_temperature == "all") %>%
             droplevels, aes(y = Num_SDE, 
                             # x = baseline, group = baseline, 
                             x = Temperature, group = Vent, 
                             fill = Vent))
    + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
               position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + theme_bw() 
    + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
    + scale_fill_manual(values = vent_colors, labels = vent_lookup, breaks = names(vent_lookup))
    # + scale_fill_viridis(discrete = TRUE, option = "plasma")
    + scale_x_discrete(labels = contrast_lookup)
    + geom_text(aes(label = Num_SDE, x = Temperature, vjust = "outward", hjust = "center",
                    y = label_y),
                # y = (Num_SDE+5)*0.9),
                # check_overlap = TRUE,
                colour = "black", fontface = "bold")
    + facet_grid(.~Vent,
                 # + facet_wrap(baseline_temperature~Temperature,
                 # ncol = 3, dir = "v", strip.position = "top",
                 labeller = labeller(baseline_temperature = baseline_lookup,
                                     Vent = vent_lookup),
                 scales = "free_x", drop = TRUE, shrink = TRUE)
    + labs(x = "Incubation temperature (˚C)", y = "Number of SDE genes")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Vent",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             color = "none")
    + ggtitle("Within vents across temperatures"))
  
  
  temp_df2 <- coa_axial_vent_SDE_summary %>%
    dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
    dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup))) %>%
    dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
                    stringr::str_remove(., "_")) %>%
    dplyr::mutate(baseline_temperature = stringr::str_extract(contrast, "axial_30|axial_55|axial_80|axial_temp") %>%
                    stringr::str_remove_all(., "axial_") %>%
                    stringr::str_replace_all("temp", "all")) %>%
    # dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)([[:digit:]]+)(.*)$", "\\2\\3\\4", contrast)) %>%
    dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)(.*)$", "\\2\\3", contrast)) %>%
    dplyr::mutate(baseline = factor(baseline, levels = names(contrast_lookup), labels = contrast_lookup)) %>%
    # dplyr::mutate(Almost_SDE = Num_almost_SDE - Num_SDE) %>%
    dplyr::mutate(label_y = dplyr::case_when( (Num_SDE < 10) ~ 10, 
                                              .default = (Num_SDE+10)*1.1)) %>%
    # .default = (Num_SDE+5)*0.9)) %>%
    droplevels 
  
  g2 <- print(
    ggplot(data = temp_df2 %>%
             dplyr::filter(baseline_temperature != "all") %>%
             droplevels, aes(y = Num_SDE, 
                             x = Vent, group = Vent, 
                             fill = Vent))
    + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
               position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + theme_bw() 
    + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
    + scale_fill_manual(values = vent_colors, labels = vent_lookup, breaks = names(vent_lookup))
    # + scale_fill_viridis(discrete = TRUE, option = "plasma")
    + scale_x_discrete(labels = contrast_lookup)
    + geom_text(aes(label = Num_SDE, x = Vent, vjust = "outward", hjust = "center", y = label_y),
                check_overlap = TRUE, colour = "black", fontface = "bold")
    + facet_wrap(baseline_temperature~.,
                 ncol = 3, dir = "v", strip.position = "top",
                 labeller = labeller(Temperature = baseline_lookup,
                                     baseline_temperature = baseline_lookup),
                 scales = "free_x", drop = TRUE, shrink = TRUE)
    + labs(x = "Vent", 
           y = "Number of SDE genes")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Vent",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)))
    + ggtitle("Between vents at each temperature"))
  
  temp_df3 <- coa_axial_year_SDE_summary %>%
    dplyr::mutate(Year = stringr::str_extract(contrast, paste0(year_lookup, collapse = "|"))) %>%
    dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
    dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup))) %>%
    dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
                    stringr::str_remove(., "_")) %>%
    dplyr::mutate(Year = stringr::str_extract(contrast, "_2013|_2014") %>%
                    stringr::str_remove(., "_")) %>%
    dplyr::mutate(baseline_temperature = stringr::str_extract(contrast, "axial_30|axial_55|axial_80|axial_temp") %>%
                    stringr::str_remove_all(., "axial_") %>%
                    stringr::str_replace_all("temp", "all")) %>%
    # dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)([[:digit:]]+)(.*)$", "\\2\\3\\4", contrast)) %>%
    dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)(.*)$", "\\2\\3", contrast)) %>%
    dplyr::mutate(baseline = factor(baseline, levels = names(contrast_lookup), labels = contrast_lookup)) %>%
    # dplyr::mutate(Almost_SDE = Num_almost_SDE - Num_SDE) %>%
    dplyr::mutate(label_y = dplyr::case_when( (Num_SDE < 10) ~ 10, 
                                              .default = (Num_SDE+5)*1.1)) %>%
    dplyr::mutate(Label_SDE = as.character(Num_SDE)) %>%
    dplyr::mutate(Label_SDE = dplyr::case_when( is.na(label_y) ~ "NA",
                                                .default = Label_SDE)) %>%
    dplyr::mutate(label_y = dplyr::case_when( is.na(label_y) ~ 10, 
                                              .default = label_y)) %>%
    dplyr::arrange(Vent, Temperature, Year) %>%
    dplyr::mutate(Temp_year = paste0(Temperature, "_", Year)) %>%
    droplevels
  
  temp_df3 <- temp_df3 %>%
    full_join(., (temp_df3 %>%
                    tidyr::pivot_wider(., id_cols = c(Vent),
                                       names_from = "Temp_year", values_from = "Num_SDE") %>%
                    tidyr::pivot_longer(., cols = !c("Vent"),
                                        names_to = "Temp_year",
                                        values_to = "Num_SDE")  %>%
                    tidyr::separate_wider_delim(., Temp_year, names = c("Temperature", "Year"), cols_remove = FALSE, delim = "_") %>%
                    dplyr::mutate(baseline_temperature = Temperature) %>%
                    dplyr::mutate(baseline = "13H samples as baseline") %>%
                    droplevels %>%
                    distinct(Vent, Temperature, .keep_all = TRUE)) ) %>%
    dplyr::mutate(Label_SDE = as.character(Num_SDE)) %>%
    dplyr::mutate(Label_SDE = dplyr::case_when( is.na(label_y) ~ "NA",
                                                .default = Label_SDE)) %>%
    dplyr::mutate(label_y = dplyr::case_when( is.na(label_y) ~ 10, 
                                              .default = label_y))
  
  
  #Temp on x-axis, facet by vent and color by year
  g3 <- print(
    ggplot(data = temp_df3, aes(y = Num_SDE, 
                                x = Temperature, fill = Year,
                                # x = Year, fill = Vent
                                group = Vent))
    + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, 
               position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
    + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
               position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
    + theme_bw() 
    + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
    + scale_fill_manual(values = colors_year, labels = year_lookup, breaks = names(year_lookup))
    # + scale_fill_manual(values = pals::ocean.oxy(n = 2), breaks = c(2013, 2014))
    # + scale_fill_viridis(discrete = TRUE, option = "cividis")
    # + scale_x_discrete(labels = temp_year_lookup)
    + geom_text(data = temp_df3 %>% tidyr::drop_na(Num_SDE),
                aes(label = Label_SDE, 
                    x = Temperature,
                    # x = Year, 
                    y = label_y), position = position_dodge2(0.9),
                check_overlap = TRUE, colour = "black", fontface = "bold")
    + geom_text(data = temp_df3 %>% dplyr::slice(which(is.na(Num_SDE))) %>% droplevels,
                aes(label = Label_SDE, 
                    x = Temperature,
                    # x = Year, 
                    y = label_y), position = position_dodge2(0.9),
                check_overlap = TRUE, colour = "black", fontface = "italic")
    
    + facet_wrap(vars(Vent), 
                 labeller = labeller(Vent = vent_lookup),
                 scales = "fixed", shrink = TRUE, drop = FALSE)
    + labs(x = "Incubation temperature (˚C)",y = "Number of SDE genes")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Year",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)))
    + ggtitle("Between years in each vent"))
  
  #temp on x-axis, facet by years, fill by vent
  g4 <- print(
    ggplot(data = temp_df3, aes(y = Num_SDE, 
                                # x = Temp_year, fill = Year,
                                x = Vent, fill = Year,
                                group = Temperature))
    + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
               position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + theme_bw() 
    + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
    + scale_fill_manual(values = colors_year, labels = year_lookup, breaks = names(year_lookup))
    # + scale_fill_manual(values = pals::ocean.oxy(n = 2), breaks = c(2013, 2014))
    # + scale_fill_viridis(discrete = TRUE, option = "plasma")
    + scale_x_discrete(labels = contrast_lookup)
    + geom_text(aes(label = Num_SDE, 
                    x = Vent,
                    y = label_y), position = position_dodge2(0.9),
                check_overlap = TRUE, colour = "black", fontface = "bold")
    + geom_text(data = data.frame(Temperature = c("30", "30", "55"), label_y = c(10, 10, 10), Vent = c("marker113", "anemone", "anemone"), Num_SDE = c("NA", "NA", "NA")),
                position = position_dodge2(0.9),
                check_overlap = TRUE,
                aes(label = Num_SDE, x = Vent, y = label_y, group = Temperature), color = "black", fontface = "italic", inherit.aes = FALSE)
    + facet_wrap(.~Temperature,
                 labeller = labeller(Temperature = baseline_lookup),
                 scales = "fixed", shrink = TRUE, drop = TRUE)
    # + facet_grid(rows = vars(Temperature), cols = vars(Vent),
    #              scales = "fixed", space = "free", shrink = TRUE, drop = TRUE,
    #              labeller = labeller(Temperature = baseline_lookup))
    + labs(x = "Vent", y = "Number of SDE genes")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Year",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)))
    + ggtitle("Between years in each vent"))


design <- "
1122
3334
"
# gpatch <- (g1 | g2 | (g3 / patchwork::guide_area())) + patchwork::plot_layout(guides = "collect") +
gpatch <- (g1 + g2 + g3 + patchwork::guide_area()) + 
  patchwork::plot_layout(guides = "collect", design = design) + 
  patchwork::plot_annotation(title = "Abundance of SDE genes, by vent, temperature, and year",
                             tag_level = "A")
gpatch <- gpatch & theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
gpatch

namevar <- "SFig-coa_axial_SDE_barcharts"
if(!all(grepl(paste0(namevar, ".png"), list.files(paste0(projectpath, "/figures/"), recursive = TRUE, include.dirs = TRUE)))){
  ggsave(paste0(projectpath, "/figures/exploratory/", namevar, ".png"),
         gpatch,
         width = 14, height = 10, units = "in")
} else {
  cli::cli_alert_info("Already generated a figure.")
}

if(!all(grepl(paste0(namevar, ".svg"), list.files(paste0(projectpath, "/figures/"), recursive = TRUE, include.dirs = TRUE)))){
  ggsave(paste0(projectpath, "/figures/exploratory/", namevar,  ".svg"),
         gpatch,
         width = 14, height = 10, units = "in")
} else {
  cli::cli_alert_info("Already generated a figure.")
}




# # Venn diagram to compare intersections of SDE genes ----------------------
# 
# #what genes are SDE in each temp at each vent, based on which group of samples form the baseline?
# 
# headers <- c("model", "contrast", "gene_symbol", "baseMean", "padj", "log2FoldChange_MMSE")
# headers <- intersect(headers, names(coa_axial_SDE.df))
# 
# p_rlist <- coa_axial_SDE.df %>%
#   dplyr::select(c(model, contrast, gene_symbol, pvalue)) %>%
#   # dplyr::select(c(contrast, gene_symbol, padj)) %>%
#   tidyr::pivot_longer(., cols = !c(model, "contrast", "gene_symbol"),
#                       names_to = "metric",
#                       values_to = "pvalue") %>%
#   tidyr::drop_na(.) %>%
#   dplyr::filter(pvalue < set.alpha) %>%
#   # droplevels %>%
#   dplyr::select(-c(metric, pvalue))
# 
# coa_axial_SDE_venn.df <- coa_axial_SDE.df %>%
#   dplyr::select(all_of(headers)) %>%
#   dplyr::mutate(contrast = factor(contrast)) %>%
#   dplyr::distinct(., .keep_all = TRUE) %>%
#   dplyr::filter(gene_symbol %in% unique(p_rlist[["gene_symbol"]])) %>%
#   dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
#   dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
#                   stringr::str_remove(., "_")) %>%
#   dplyr::mutate(Year = stringr::str_extract(contrast, "_2013|_2014") %>%
#                   stringr::str_remove(., "_")) %>%
#   dplyr::mutate(baseline_temperature = stringr::str_extract(contrast, "axial_30|axial_55|axial_80|axial_temp") %>%
#                   stringr::str_remove_all(., "axial_") %>%
#                   stringr::str_replace_all("temp", "all")) %>%
#   dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)(.*)$", "\\2\\3", contrast)) %>%
#   dplyr::rename(gene_id = "gene_symbol") %>%
#   droplevels %>%
#   dplyr::right_join(., p_rlist,
#              by = c("gene_id" = "gene_symbol", "contrast" = "contrast", "model" = "model"),
#              relationship = "many-to-many") %>%
#   droplevels %>%
#   dplyr::group_by(Vent, gene_id, baseline, model) %>%
#   dplyr::mutate(obs = length(contrast)) %>%
#   droplevels %>%
#   dplyr::select(-obs) %>%
#   ungroup
# 
# #temperature-specific, intervent comparisons
# coa_axial_eachtemp_venn_list <- coa_axial_SDE_venn.df %>%
#   dplyr::filter(!grepl("2013|2014", contrast)) %>%
#   dplyr::filter(grepl("eachtemp", model)) %>%
#   dplyr::filter(!grepl("all", baseline_temperature)) %>%
#   dplyr::select(-Year) %>%
#   droplevels %>%
#   split(., f = .$baseline_temperature) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$Vent, drop = FALSE)) %>%
#   purrr::list_flatten(., name_spec = "{inner}_{outer}") %>%
#   bind_rows(., .id = "temp_vent") %>%
#   split(., f = .$Temperature) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$temp_vent) %>%
#         map(., ~.x %>%
#               droplevels %>%
#               dplyr::select(gene_id) %>%
#               droplevels %>%
#               simplify))
#               # simplify) %>% purrr::list_flatten(., name_spec = "{inner}"))
#               # simplify)) %>%
#   # map(., ~.x %>% purrr::list_flatten(., name_spec = "{inner}_{outer}"))
#   # purrr::list_flatten(., name_spec = "{inner}")
# 
# coa_axial_eachtemp_year_venn_list <- coa_axial_SDE_venn.df %>%
#   dplyr::filter(grepl("2013|2014", contrast)) %>%
#   dplyr::distinct(Vent, Temperature, Year, baseline_temperature, baseline, .keep_all = FALSE) %>%
#   droplevels %>%
#   tidyr::expand(nesting(Vent, Temperature, baseline_temperature, baseline), Year) %>%
#   dplyr::mutate(contrast = paste0(Vent, "_", Temperature, "_13H_", Year, " - ", baseline)) %>%
#   dplyr::full_join(., coa_axial_SDE_venn.df %>%
#                      dplyr::filter(grepl("2013|2014", contrast)) %>%
#                      droplevels,
#                    by = join_by(contrast, Vent, Temperature, Year, baseline_temperature, baseline)) %>%
#   droplevels %>%
#   split(., f = .$baseline_temperature) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$Vent, drop = FALSE)) %>%
#   purrr::list_flatten(., name_spec = "{inner}_{outer}") %>%
#   bind_rows(., .id = "temp_vent") %>%
#   split(., f = .$Temperature) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$temp_vent) %>%
#         map(., ~.x %>%
#               droplevels %>%
#               split(., f = .$Year, drop = FALSE) %>%
#               map(., ~.x %>%
#                     dplyr::select(gene_id) %>%
#                     droplevels %>%
#                     simplify))) 
# 
# #the vent-specific, thermal comparisons
# coa_axial_temp_venn_list <- coa_axial_SDE_venn.df %>%
#   dplyr::filter(!grepl("2013|2014", contrast)) %>%
#   dplyr::filter(!grepl("eachtemp", model)) %>%
#   dplyr::filter(grepl("all", baseline_temperature)) %>%
#   dplyr::select(-Year) %>%
#   droplevels %>%
#   split(., f = .$Vent) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$Temperature) %>%
#         map(., ~.x %>%
#               droplevels %>%
#               dplyr::select(gene_id) %>%
#               droplevels %>%
#               simplify))
#   #             simplify)) %>% purrr::list_flatten(., name_spec = "{inner}_{outer}")
# 
# {
# # df <- rlist %>%
# #   # dplyr::select(c(contrast, gene_symbol, pvalue)) %>%
# #   dplyr::select(c(contrast, gene_symbol, padj)) %>%
# #   dplyr::mutate(contrast = factor(contrast)) %>%
# #   tidyr::pivot_longer(., cols = !c("contrast", "gene_symbol"),
# #                       names_to = "metric",
# #                       values_to = "pvalue") %>%
# #   tidyr::drop_na(.) %>%
# #   dplyr::filter(pvalue < set.alpha) %>%
# #   droplevels
# 
# # #do the btw-vent and btw-temp comparisons first
# # coa_axial_temp_venn_list <- rlist %>%
# #   dplyr::select(all_of(headers)) %>%
# #   dplyr::mutate(contrast = factor(contrast)) %>%
# #   dplyr::distinct(., .keep_all = TRUE) %>%
# #   dplyr::filter(gene_symbol %in% unique(df[["gene_symbol"]])) %>%
# #   dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
# #   dplyr::mutate(Vent = stringr::str_extract(contrast, paste0(vents, collapse = "|"))) %>%
# #   dplyr::mutate(Temperature = stringr::str_extract(contrast, "_30|_55|_80") %>%
# #                   stringr::str_remove(., "_")) %>%
# #   dplyr::mutate(Year = stringr::str_extract(contrast, "_2013|_2014") %>%
# #                   stringr::str_remove(., "_")) %>%
# #   dplyr::mutate(baseline_temperature = stringr::str_extract(contrast, "axial_30|axial_55|axial_80|axial_temp") %>%
# #                   stringr::str_remove_all(., "axial_") %>%
# #                   stringr::str_replace_all("temp", "all")) %>%
# #   dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)(.*)$", "\\2\\3", contrast)) %>%
# #   dplyr::rename(gene_id = "gene_symbol") %>%
# #   droplevels %>%
# #   right_join(., (df %>%
# #                    dplyr::select(-c(metric, pvalue))),
# #              by = c("gene_id" = "gene_symbol", "contrast" = "contrast"),
# #              relationship = "many-to-many") %>%
# #   droplevels %>%
# #   group_by(Vent, gene_id, baseline) %>%
# #   dplyr::mutate(obs = length(contrast)) %>%
# #   droplevels %>%
# #   dplyr::select(-obs) %>%
# #   ungroup
# 
# # coa_axial_eachtemp_venn_list <- coa_axial_temp_venn_list %>%
# #   dplyr::filter(!grepl("2013|2014", contrast)) %>%
# #   dplyr::filter(!grepl("all", baseline_temperature)) %>%
# #   dplyr::select(-Year) %>%
# #   droplevels %>%
# #   split(., f = .$baseline_temperature) %>%
# #   map(., ~.x %>%
# #         droplevels %>%
# #         split(., f = .$Vent) %>%
# #         map(., ~.x %>%
# #               droplevels %>%
# #                     dplyr::select(gene_id) %>%
# #                     droplevels %>%
# #                     simplify))
# # 
# # coa_axial_eachtemp_year_venn_list <- coa_axial_temp_venn_list %>%
# #   dplyr::filter(grepl("2013|2014", contrast)) %>%
# #   dplyr::distinct(Vent, Temperature, Year, baseline_temperature, baseline, .keep_all = FALSE) %>%
# #   droplevels %>%
# #   tidyr::expand(nesting(Vent, Temperature, baseline_temperature, baseline), Year) %>%
# #   dplyr::mutate(contrast = paste0(Vent, "_", Temperature, "_13H_", Year, " - ", baseline)) %>%
# #   dplyr::full_join(., coa_axial_temp_venn_list %>%
# #                      dplyr::filter(grepl("2013|2014", contrast)) %>%
# #                      droplevels,
# #                    by = join_by(contrast, Vent, Temperature, Year, baseline_temperature, baseline)) %>%
# # #   droplevels
# # # coa_axial_eachtemp_year_venn_list <- coa_axial_temp_venn_list %>%
# # #   dplyr::filter(grepl("2013|2014", contrast)) %>%
# # #   # dplyr::filter(!grepl("all", baseline_temperature)) %>%
# #   droplevels %>%
# #   split(., f = .$baseline_temperature) %>%
# #   map(., ~.x %>%
# #         droplevels %>%
# #         split(., f = .$Vent) %>%
# #         map(., ~.x %>%
# #               droplevels %>%
# #               split(., f = .$Year) %>%
# #               map(., ~.x %>%
# #                     dplyr::select(gene_id) %>%
# #                     droplevels %>%
# #                     simplify)))
# # 
# # coa_axial_temp_venn_list <- coa_axial_temp_venn_list %>%
# #   dplyr::filter(grepl("all", baseline_temperature)) %>%
# #   dplyr::filter(!grepl("2013|2014", contrast)) %>%
# #   dplyr::select(-Year) %>%
# #   droplevels %>%
# #   split(., f = .$Vent) %>%
# #   map(., ~.x %>%
# #         droplevels %>%
# #         split(., f = .$Temperature) %>%
# #         map(., ~.x %>%
# #               droplevels %>%
# #                     dplyr::select(gene_id) %>%
# #                     droplevels %>%
# #                     simplify))
# 
# if(!file.exists(paste0(projectpath, "/intermediate/deseq/", "coa_axial_temp_venn_list.rds"))){
#   # readr::write_rds(axial_temp_venn_list, paste0(getwd(), "/", "coa_axial_temp_venn_list.rds"),
#   #                  compress = "gz")
# }
# }
# 
# 
# #simple venns:
# {
#   #all of them:
#   # coa_axial_temp_venn_upset <- coa_axial_temp_venn_list %>%
#   #   purrr::list_flatten(., name_spec = "{outer}_{inner}") %>%
#   #   purrr::list_flatten(., name_spec = "{outer}_{inner}")
#   # g <- UpSetR::upset(UpSetR::fromList(coa_axial_temp_venn_upset),
#   #                    nsets = 27,
#   #                    keep.order = TRUE)
#   # 
#   # g <- ggplotify::as.ggplot(g)
#   # # if(is.null(grep("coa_axial_temp_frac_shared_SDE.png", list.files(getwd(), recursive = TRUE, include.dirs = TRUE), value = TRUE))){
#   # #   ggsave(paste0(getwd(), "/", Sys.Date(), "-", "coa_axial_temp_frac_shared_SDE.png"),
#   # #          g,
#   # #          width = 12, height = 8, units = "in")
#   # # } else {
#   # #   cli::cli_alert_info("Already generated a figure.")
#   # # }
#   # 
#   # #each temp:
#   # 
#   # coa_axial_eachtemp_venn_upset <- coa_axial_eachtemp_venn_list %>%
#   #   # purrr::list_flatten(., name_spec = "{outer}_{inner}") %>%
#   #   purrr::list_flatten(., name_spec = "{outer}_{inner}")
#   # g2 <- UpSetR::upset(UpSetR::fromList(coa_axial_eachtemp_venn_upset),
#   #                    nsets = 27,
#   #                    keep.order = TRUE)
#   # 
#   # g2 <- ggplotify::as.ggplot(g)
# }
# 
# 
# # metadata_filter <- list(`Vent` = c("anemone", "marker113", "marker33"),
# #                         `Temperature` = c("30", "55", "80"),
# #                         `baseline` = c("axial_temp_12L","axial_temp_both",  "axial_temp_13H"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "coa_axial_all_temp")
# 
# {
#   # mdata_filter <- metadata_filter
#   # vennlist <- coa_axial_temp_venn_list
#   # prefix <- "coa_axial_all_temp"
#   # 
#   # params_idx <- map(mdata_filter,
#   #                   ~unlist(., recursive = FALSE)) %>%
#   #   purrr::reduce(interaction, sep = "_") %>%
#   #   levels %>%
#   #   as.character
#   # 
#   # params_sets <- bind_cols(as_tibble_col(params_idx, column_name = "label"),
#   #                          as_tibble(stringr::str_split_fixed(params_idx, "_", n = length(mdata_filter)))) %>%
#   #   dplyr::filter(!grepl("12L", label)) %>%
#   #   tidyr::nest(data = label, .by = c(V1, V2)) 
#   # params_list <- params_sets %>%
#   #   dplyr::select(data) %>%
#   #   droplevels %>%
#   #   tibble::deframe(.) %>%
#   #   setNames(., paste0(params_sets[[1]], "_", params_sets[[2]])) %>%
#   #   map(., ~.x %>%
#   #         unlist(., recursive = FALSE, use.names = FALSE))
#   # 
#   # 
#   # # params_list <- list("13H" = grep(paste0(c("temp_13H"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   # #                     "12L" = grep(paste0(c("temp_12L"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   # #                     "13H+12L" = grep(paste0(c("temp_both"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   # #                     "13H \u2229 12L" = grep(paste0(c("temp_13H", "temp_12L"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   # #                     "13H \u2229 13H+12L" = grep(paste0(c("temp_13H", "temp_both"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   # #                     "12L \u2229 13H+12L" = grep(paste0(c("temp_12L", "temp_both"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   # #                     "12L \u2229 13H \u2229 13H+12L" = grep(paste0(c("temp_13H", "temp_both", "temp_12L"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE))
#   # 
#   # g <- UpSetR::upset(UpSetR::fromList(labels_vennlist), 
#   #                    main.bar.color = "black",
#   #                    mainbar.y.label = "Intersection size",
#   #                    sets.x.label = "Set size", 
#   #                    sets = c(params_list),
#   #                    # order.by = "freq",
#   #                    # nsets = length(labels_vennlist),
#   #                    # queries = list(list(query = intersects, params = params_list[1], color = "red", active = T, query.name = paste(names(params_list[1]))),
#   #                    #                list(query = intersects, params = params_list[2], color = "purple", active = T, query.name = paste(names(params_list[2]))),
#   #                    #                list(query = intersects, params = params_list[3], color = "brown", active = T, query.name = paste(names(params_list[3]))),
#   #                    #                list(query = intersects, params = params_list[4], color = "blue", active = T, query.name = paste(names(params_list[4]))),
#   #                    #                list(query = intersects, params = params_list[5], color = "turquoise", active = T, query.name = paste(names(params_list[5]))),
#   #                    #                list(query = intersects, params = params_list[6], color = "navy", active = T, query.name = paste(names(params_list[6]))),
#   #                    #                list(query = intersects, params = params_list[7], color = "orange", active = T, query.name = paste(names(params_list[7]))),
#   #                    #                list(query = intersects, params = params_list[8], color = "black", active = T, query.name = paste(names(params_list[8]))),
#   #                    #                list(query = intersects, params = params_list[9], color = "tan", active = T, query.name = paste(names(params_list[9])))),
#   #                    # query.legend = "bottom",
#   #                    keep.order = TRUE)
#   # g
# }
# 
# # ggsave(paste0(getwd(), "/", Sys.Date(), "-", "axial_temp_frac_shared_SDE.png"),
# #        g,
# #        width = 12, height = 8, units = "in")
# 
# 
# 
# # Individual temperatures by vent -----------------------------------------
# 
# 
# #not working at this time...
# 
# #Temperature-specific, between vents:
# metadata_filter <- list(`Temperature` = c("80"),
#                         `Vent` = c("anemone", "marker113", "marker33"))
# P_venn_upset(vennlist = coa_axial_eachtemp_venn_list,
#              mdata_filter = metadata_filter,
#              prefix = "baseline_80")
# {
#   # metadata_filter <- list(`Temperature` = c("80"),
#   #                         `Vent` = c("anemone", "marker113", "marker33"))
#   # mdata_filter <- metadata_filter
#   # vennlist <- coa_axial_eachtemp_venn_list
#   # # vennlist <- get("vennlist", inherits = TRUE)
#   # # mdata_filter <- get("mdata_filter", inherits = TRUE)
#   # 
#   # #params filter
#   # if(exists("labels_vennlist", inherits = FALSE)){
#   #   rm(labels_vennlist)}
#   # 
#   # if(exists("mdata_filter", inherits = TRUE)){
#   #   labels_idx <- imap(mdata_filter,
#   #                      ~vennlist[.x]) %>%
#   #     Filter(Negate(function(x) is.null(unlist(x))), .) %>%
#   #     purrr::list_flatten(., name_spec = "{inner}")
#   # 
#   #   for(i in seq(length(labels_idx))){
#   #     namevar <- names(labels_idx[i])
#   #     for(j in seq(pluck_depth(labels_idx[i]) -1)){
#   #       temp_list <- labels_idx[[namevar]][j] %>%
#   #         setNames(., paste0(namevar, "_", names(.))) %>%
#   #         Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
#   #         # purrr::list_flatten(., name_spec = "{inner}_{outer}") %>%
#   #         # purrr::list_flatten(.) %>%
#   #         purrr::list_flatten(.)
#   #       if(exists("labels_vennlist")){
#   #         labels_vennlist <- c(labels_vennlist, temp_list)
#   #       } else {
#   #         labels_vennlist <- temp_list
#   #       }
#   #     }
#   #   }
#   #   keep <- mdata_filter %>%
#   #     map(., ~grepl(paste0(.x, collapse = "|"), names(labels_vennlist))) %>%
#   #     purrr::reduce(`&`)
#   #   labels_vennlist <- labels_vennlist[keep]
#   # 
#   # } else {
#   #   labels_vennlist <- map(vennlist,
#   #                          ~.x %>%
#   #                            # purrr::list_flatten(., name_spec = "{inner}_{outer}") %>%
#   #                            purrr::list_flatten(., name_spec = "{outer}_{inner}"))
#   # }
#   # 
#   # params_list <- list("anemone" = grep(paste0(c("anemone"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "marker113" = grep(paste0(c("marker113"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "marker33" = grep(paste0(c("marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "anemone \u2229 marker113" = grep(paste0(c("anemone", "marker113"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "anemone \u2229 marker33" = grep(paste0(c("anemone", "marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "marker113 \u2229 marker33" = grep(paste0(c("marker113", "marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE),
#   #                     "anemone \u2229 marker113 \u2229 marker33" = grep(paste0(c("anemone", "marker113", "marker33"), collapse = "|"), names(labels_vennlist), value = TRUE, ignore.case = TRUE))
#   # 
#   # #13H: red
#   # #12L: blue
#   # #12L+13H: green
#   # #13H vs 12L: purple
#   # #13H vs 12L+13H: brown
#   # #turquoise
#   # #all: white/black
#   # 
#   # g <- UpSetR::upset(UpSetR::fromList(labels_vennlist),
#   #                    main.bar.color = "black",
#   #                    mainbar.y.label = "Intersection size",
#   #                    sets.x.label = "Set size",
#   #                    sets = c(names(labels_vennlist)),
#   #                    # empty.intersections = "keep",
#   #                    order.by = "degree", decreasing = F,
#   #                    queries = list(list(query = intersects, params = params_list[1], color = "red", active = T, query.name = paste(names(params_list[1]))),
#   #                                   list(query = intersects, params = params_list[2], color = "blue", active = T, query.name = paste(names(params_list[2]))),
#   #                                   list(query = intersects, params = params_list[3], color = "darkgreen", active = T, query.name = paste(names(params_list[3]))),
#   #                                   list(query = intersects, params = params_list[4], color = "purple", active = T, query.name = paste(names(params_list[4]))),
#   #                                   list(query = intersects, params = params_list[5], color = "brown", active = T, query.name = paste(names(params_list[5]))),
#   #                                   list(query = intersects, params = params_list[6], color = "turquoise", active = T, query.name = paste(names(params_list[6]))),
#   #                                   list(query = intersects, params = params_list[7], color = "black", active = T, query.name = paste(names(params_list[7])))),
#   #                    keep.order = F,
#   #                    query.legend = "bottom")
#   # g
# }
# 
# #vent-specific, between temperatures
# metadata_filter <- list(`Vent` = c("anemone"),
#                         `Temperature` = c("30", "55", "80"))
# P_venn_upset(vennlist = coa_axial_temp_venn_list,
#              mdata_filter = metadata_filter,
#              prefix = "anemone") 
# 
# # metadata_filter <- list(`Temperature` = c("55"),
# #                         `Vent` = c("anemone", "marker113", "marker33"))
# # # `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_55") 
# 
# metadata_filter <- list(`Temperature` = c("30"),
#                         `Vent` = c("anemone", "marker113", "marker33"))
# # `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# P_venn_upset(vennlist = coa_axial_temp_venn_list,
#              mdata_filter = metadata_filter,
#              prefix = "baseline_30") 
# # g1 <- ggplotify::as.ggplot(baseline_anemone_80_upset_plot)
# # g2 <- ggplotify::as.ggplot(baseline_marker113_80_upset_plot)
# # g3 <- ggplotify::as.ggplot(baseline_marker33_80_upset_plot)
# # g4 <- ggplotify::as.ggplot(baseline_marker33_80_upset_plot[[2]])
# g5 <- ggplotify::as.ggplot(baseline_80_upset_plot[["legend"]])
# 
# # layout <- c(patchwork::area(1, 2, 2, 3),
# #             patchwork::area(3, 1, 3, 3))
# # plot(layout)
# # ggplotify::as.ggplot(baseline_marker33_80_upset_plot[[1]]) + ggplotify::as.ggplot(baseline_marker33_80_upset_plot[[2]]) + patchwork::plot_layout(design = layout)
# 
# design <- "
# #111
# 2222
# "
# # g1A <- ggplotify::as.ggplot(baseline_anemone_80_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_anemone_80_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # g2A <- ggplotify::as.ggplot(baseline_marker113_80_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_marker113_80_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # g3A <- ggplotify::as.ggplot(baseline_marker33_80_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_marker33_80_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # 
# # title_temp <- metadata_filter[["Temperature"]]
# # # gpatch <- ((g1 +ggtitle("Anemone")) / (g2 + ggtitle("Marker113")) / (g3 + ggtitle("Marker33")))
# # gpatch <- (g1A + ggtitle("Anemone")) / (g2A + ggtitle("Marker113")) / (g3A + ggtitle("Marker33"))
# # gpatch <- (gpatch / (patchwork::wrap_elements(g5, ignore_tag = TRUE))) + patchwork::plot_layout(ncol = 1, nrow = 4, heights = c(3,3,3, 1), byrow = FALSE, tag_level = "new")
# # gpatch <- gpatch + patchwork::plot_annotation(title = "Overlap in SDE gene sets",
# #                                               subtitle = paste0("At temperature ", title_temp, "\u00B0 C"),
# #                                               tag_level = "A")
# # gpatch
# # 
# # ggsave(paste0(getwd(), "/", Sys.Date(), "-","coa-", "axial_", title_temp, "_frac_shared_SDE.png"),
# #        gpatch,
# #        width = 6, height = 12, units = "in")
# g1A <- ggplotify::as.ggplot(baseline_80_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_80_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# g2A <- ggplotify::as.ggplot(baseline_55_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_55_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# g3A <- ggplotify::as.ggplot(baseline_30_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_30_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# 
# # gpatch <- ((g1 +ggtitle("Anemone")) / (g2 + ggtitle("Marker113")) / (g3 + ggtitle("Marker33")))
# gpatch <- (g1A + ggtitle("80")) / (g2A + ggtitle("55")) / (g3A + ggtitle("30"))
# gpatch <- (gpatch / (patchwork::wrap_elements(g5, ignore_tag = TRUE))) + patchwork::plot_layout(ncol = 1, nrow = 4, heights = c(3,3,3, 1), byrow = FALSE, tag_level = "new")
# gpatch <- gpatch + patchwork::plot_annotation(title = "Overlap in SDE gene sets between vents at each temperature",
#                                               subtitle = paste0("Using 13H as baseline expression"),
#                                               tag_level = "A")
# gpatch
# ggsave(paste0(projectpath, "/figures/exploratory/", Sys.Date(), "-","coa-", "axial_by_temp", "_frac_shared_SDE.png"),
#        gpatch,
#        width = 6, height = 12, units = "in")
# 
# # # 55 and 30 next ----------------------------------------------------------
# # 
# # 
# # metadata_filter <- list(`Vent` = c("anemone"),
# #                         `Temperature` = c("55"),
# #                         `baseline` = c("axial_temp_13H"))
# # # `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_anemone_55")
# # 
# # metadata_filter <- list(`Vent` = c("marker113"),
# #                         `Temperature` = c("55"),
# #                         `baseline` = c("axial_temp_13H"))
# # # `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_marker113_55") 
# # metadata_filter <- list(`Vent` = c("marker33"),
# #                         `Temperature` = c("55"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_marker33_55") 
# # 
# # 
# # design <- "
# # #111
# # 2222
# # "
# # g1A <- ggplotify::as.ggplot(baseline_anemone_55_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_anemone_55_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # g2A <- ggplotify::as.ggplot(baseline_marker113_55_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_marker113_55_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # g3A <- ggplotify::as.ggplot(baseline_marker33_55_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_marker33_55_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # 
# # title_temp <- metadata_filter[["Temperature"]]
# # gpatch <- (g1A + ggtitle("Anemone")) / (g2A + ggtitle("Marker113")) / (g3A + ggtitle("Marker33"))
# # gpatch <- (gpatch / (patchwork::wrap_elements(g5, ignore_tag = TRUE))) + patchwork::plot_layout(ncol = 1, nrow = 4, heights = c(3,3,3, 1), byrow = FALSE, tag_level = "new")
# # gpatch <- gpatch + patchwork::plot_annotation(title = "Overlap in SDE gene sets",
# #                                               subtitle = paste0("At temperature ", title_temp, "\u00B0 C"),
# #                                               tag_level = "A")
# # gpatch
# # 
# # ggsave(paste0(getwd(), "/", Sys.Date(), "-","coa-", "axial_", title_temp, "_frac_shared_SDE.png"),
# #        gpatch,
# #        width = 6, height = 12, units = "in")
# # 
# # 
# # 
# # # 30 now ------------------------------------------------------------------
# # 
# # metadata_filter <- list(`Vent` = c("anemone"),
# #                         `Temperature` = c("30"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_anemone_30")
# # 
# # metadata_filter <- list(`Vent` = c("marker113"),
# #                         `Temperature` = c("30"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_marker113_30") 
# # metadata_filter <- list(`Vent` = c("marker33"),
# #                         `Temperature` = c("30"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # P_venn_upset(vennlist = coa_axial_temp_venn_list,
# #              mdata_filter = metadata_filter,
# #              prefix = "baseline_marker33_30") 
# # 
# # design <- "
# # #111
# # 2222
# # "
# # g1A <- ggplotify::as.ggplot(baseline_anemone_30_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_anemone_30_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # g2A <- ggplotify::as.ggplot(baseline_marker113_30_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_marker113_30_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # g3A <- ggplotify::as.ggplot(baseline_marker33_30_upset_plot[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(baseline_marker33_30_upset_plot[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = design)
# # 
# # title_temp <- metadata_filter[["Temperature"]]
# # gpatch <- (g1A + ggtitle("Anemone")) / (g2A + ggtitle("Marker113")) / (g3A + ggtitle("Marker33"))
# # gpatch <- (gpatch / (patchwork::wrap_elements(g5, ignore_tag = TRUE))) + patchwork::plot_layout(ncol = 1, nrow = 4, heights = c(3,3,3, 1), byrow = FALSE, tag_level = "new")
# # gpatch <- gpatch + patchwork::plot_annotation(title = "Overlap in SDE gene sets",
# #                                               subtitle = paste0("At temperature ", title_temp, "\u00B0 C"),
# #                                               tag_level = "A")
# # gpatch
# # 
# # ggsave(paste0(getwd(), "/", Sys.Date(), "-","coa-", "axial_", title_temp, "_frac_shared_SDE.png"),
# #        gpatch,
# #        width = 6, height = 12, units = "in")
# 
# 
# # What are the diffs in SDE genes/ ----------------------------------------
# 
# #`13H ∩ 13H+12L`
# #In anemone 80: 26+108=134 genes are shared between contrasts using baseline expression of 13H and 13H+12L
# #marker113 80: 14+67=81 genes are shared
# #marker33: 63+171=234 genes are shared
# #30 genes in Anemone_80, 34 genes in Marker113_80, and 32 genes in Marker33_80 are SDE only when using only 12L as baseline (not seen elsewhere)
# 
# # metadata_filter <- list(`Vent` = c("anemone"),
# #                         `Temperature` = c("80"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # F_venn_setdiffs(vennlist = coa_axial_temp_venn_list,
# #                 mdata_filter = metadata_filter,
# #                 prefix = "baseline")
# # metadata_filter <- list(`Vent` = c("marker113"),
# #                         `Temperature` = c("80"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # F_venn_setdiffs(vennlist = coa_axial_temp_venn_list,
# #                 mdata_filter = metadata_filter,
# #                 prefix = "baseline")
# # metadata_filter <- list(`Vent` = c("marker33"),
# #                         `Temperature` = c("80"),
# #                         `baseline` = c("axial_temp_13H", "axial_temp_12L","axial_temp_both"))
# # F_venn_setdiffs(vennlist = coa_axial_temp_venn_list,
# #                 mdata_filter = metadata_filter,
# #                 prefix = "baseline")
# 
# # Filter SDES ---------------
# {
# # vent13H_axial12L_SDE_full <- rlist %>%
# #   dplyr::select(all_of(headers)) %>%
# #   dplyr::mutate(contrast = factor(contrast)) %>%
# #   dplyr::distinct(., .keep_all = TRUE) %>%
# #   dplyr::filter(gene_symbol %in% unique(df[["gene_symbol"]])) %>%
# #   dplyr::mutate(Vent = str_extract(contrast, paste0(vents, collapse = "|"))) %>%
# #   dplyr::rename(gene_id = "gene_symbol") %>%
# #   droplevels %>%
# #   ungroup %>%
# #   left_join(., df %>% 
# #               dplyr::select(-c(metric, padj)), 
# #             # relationship = "many-to-many",
# #             by = c("gene_id" = "gene_symbol", "contrast" = "contrast")) %>%
# #   dplyr::filter(grepl("temp_12L", contrast)) %>%
# #   droplevels %>%
# #   dplyr::mutate(originalcontrast = contrast) %>%
# #   dplyr::mutate(contrast = dplyr::case_when(grepl("30", originalcontrast) ~ "vent_13H_30",
# #                                             grepl("55", originalcontrast) ~ "vent_13H_55",
# #                                             grepl("80", originalcontrast) ~ "vent_13H_80")) %>%
# #   dplyr::mutate(contrast = paste0(contrast, " - axial_temp_12L") %>% factor(.)) %>%
# #   dplyr::mutate(Temperature = str_extract(contrast, "30|55|80")) %>%
# #   dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_temp_)(.*)$", "\\2\\3", contrast)) %>%
# #   droplevels 
# # vent13H_axialboth_SDE_full <- rlist %>%
# #   dplyr::select(all_of(headers)) %>%
# #   dplyr::mutate(contrast = factor(contrast)) %>%
# #   dplyr::distinct(., .keep_all = TRUE) %>%
# #   dplyr::filter(gene_symbol %in% unique(df[["gene_symbol"]])) %>%
# #   dplyr::mutate(Vent = str_extract(contrast, paste0(vents, collapse = "|"))) %>%
# #   dplyr::rename(gene_id = "gene_symbol") %>%
# #   droplevels %>%
# #   left_join(., df %>% 
# #               dplyr::select(-c(metric, padj)), 
# #             relationship = "many-to-many",
# #             by = c("gene_id" = "gene_symbol", "contrast" = "contrast")) %>%
# #   dplyr::filter(grepl("temp_both", contrast)) %>%
# #   dplyr::mutate(originalcontrast = contrast) %>%
# #   dplyr::mutate(contrast = dplyr::case_when(grepl("30", originalcontrast) ~ "vent_13H_30",
# #                                             grepl("55", originalcontrast) ~ "vent_13H_55",
# #                                             grepl("80", originalcontrast) ~ "vent_13H_80")) %>%
# #   dplyr::mutate(contrast = paste0(contrast, " - axial_temp_both") %>% factor(.)) %>%
# #   dplyr::mutate(Temperature = str_extract(contrast, "30|55|80")) %>%
# #   dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_temp_)(.*)$", "\\2\\3", contrast)) %>%
# #   # dplyr::mutate(padj = ifelse(padj < 0.05, padj, NA)) %>%
# #   droplevels 
# }
# 
# # vent13H_axial13H_SDE_full <- rlist %>%
# #   dplyr::select(all_of(headers)) %>%
# #   dplyr::mutate(contrast = factor(contrast)) %>%
# #   dplyr::distinct(., .keep_all = TRUE) %>%
# #   dplyr::filter(gene_symbol %in% unique(df[["gene_symbol"]])) %>%
# #   dplyr::mutate(Vent = str_extract(contrast, paste0(vents, collapse = "|"))) %>%
# #   dplyr::rename(gene_id = "gene_symbol") %>%
# #   droplevels %>%
# #   left_join(., df %>% 
# #               dplyr::select(-c(metric, pvalue)), 
# #             relationship = "many-to-many",
# #             by = c("gene_id" = "gene_symbol", "contrast" = "contrast")) %>%
# #   # dplyr::filter(grepl("temp_13H", contrast)) %>%
# #   droplevels %>%
# #   dplyr::mutate(originalcontrast = contrast) %>%
# #   dplyr::mutate(contrast = dplyr::case_when(grepl("30", originalcontrast) ~ "vent_13H_30",
# #                                             grepl("55", originalcontrast) ~ "vent_13H_55",
# #                                             grepl("80", originalcontrast) ~ "vent_13H_80")) %>%
# #   dplyr::mutate(contrast = paste0(contrast, " - axial_temp_13H") %>% factor(.)) %>%
# #   dplyr::mutate(Temperature = str_extract(contrast, "30|55|80")) %>%
# #   dplyr::mutate(baseline = gsub("([[:print:]]+)(axial_)(.*)$", "\\2\\3", contrast)) %>%
# #   droplevels 
# # 
# # 
# # # Plot the baseline expression patterns -----------------------------------
# # 
# # 
# # for(i in c("coa_axial_30_13H", "coa_axial_55_13H", "coa_axial_80_13H")){
# #   if(!exists(paste0(i, "_deseq_list"), envir = .GlobalEnv)){
# #     namevar <- i
# #     if(file.exists(paste0(getwd(), "/", namevar, "_deseq_list", ".rds"))){
# #       cli::cli_alert_info("Reading in DESeq results list: {namevar}.")
# #       temp_list <- readr::read_rds(file = paste0(getwd(), "/", namevar, "_deseq_list", ".rds"))
# #       assign(paste0(namevar, "_deseq_list"), temp_list, envir = .GlobalEnv, inherits = TRUE)
# #       rm(temp_list)
# #     }
# #     else {
# #       cli::cli_alert_info("DESeq results list does not exist with this name: {namevar}. Please process through DESeq.")
# #     }
# #   }
# #   rm(namevar)
# # }
# # 
# # # c("coa_axial_temp_13H_deseq_list", "coa_axial_temp_12L_deseq_list", "coa_axial_temp_bothfrac_deseq_list")
# # # if(!exists("coa_axial_temp_bothfrac_deseq_list", envir = .GlobalEnv)){
# # #   if(file.exists(paste0(getwd(), "/", "coa_axial_temp_bothfrac_deseq_list", ".rds"))){
# # #     cli::cli_alert_info("Reading in the coa-mtg DESeq formatted counts table...")
# # #     progressr::with_progress({
# # #       coa_axial_temp_bothfrac_deseq_list <- readr::read_rds(file = paste0(getwd(), "/", "coa_axial_temp_bothfrac_deseq_list", ".rds"))
# # #     })
# # #   } else {
# # #     NULL
# # #   }
# # # }
# # # if(!exists("coa_axial_temp_12L_deseq_list", envir = .GlobalEnv)){
# # #   if(file.exists(paste0(getwd(), "/", "coa_axial_temp_12L_deseq_list", ".rds"))){
# # #     cli::cli_alert_info("Reading in the coa-mtg DESeq formatted counts table...")
# # #     progressr::with_progress({
# # #       coa_axial_temp_12L_deseq_list <- readr::read_rds(file = paste0(getwd(), "/", "coa_axial_temp_12L_deseq_list", ".rds"))
# # #     })
# # #   } else {
# # #     NULL
# # #   }
# # # }
# # # if(!exists("coa_axial_temp_13H_deseq_list", envir = .GlobalEnv)){
# # #   if(file.exists(paste0(getwd(), "/", "coa_axial_temp_13H_deseq_list", ".rds"))){
# # #     cli::cli_alert_info("Reading in the coa-mtg DESeq formatted counts table...")
# # #     progressr::with_progress({
# # #       coa_axial_temp_13H_deseq_list <- readr::read_rds(file = paste0(getwd(), "/", "coa_axial_temp_13H_deseq_list", ".rds"))
# # #     })
# # #   } else {
# # #     NULL
# # #   }
# # # }
# # 
# # p_rlist
# # list2env(coa_axial_30_13H_deseq_list, envir = .GlobalEnv)
# # coa_axial_30_13H_SDE_counts <- counts(dds, normalized = TRUE)[match(unique(p_rlist[["gene_symbol"]]), rownames(counts(dds))),]
# # 
# # list2env(coa_axial_55_13H_deseq_list, envir = .GlobalEnv)
# # coa_axial_55_13H_SDE_counts <- counts(dds, normalized = TRUE)[match(unique(p_rlist[["gene_symbol"]]), rownames(counts(dds))),]
# # 
# # list2env(coa_axial_80_13H_deseq_list, envir = .GlobalEnv)
# # coa_axial_80_13H_SDE_counts <- counts(dds, normalized = TRUE)[match(unique(p_rlist[["gene_symbol"]]), rownames(counts(dds))),]
# # 
# # #use the metadata/mdata from axial_both
# # metadata <- sample_metadata
# # # mdata_filter <- list(`Fraction` = c("12L", "13H"))
# # # metadata <- map2(mdata_filter, names(mdata_filter),
# # #                  ~coa_axial_temp_bothfrac_deseq_list[["metadata"]] %>%
# # #                    dplyr::filter(
# # #                      grepl(paste0(.x, collapse = "|"), .data[[.y]], ignore.case = TRUE)
# # #                    )) %>%
# # #   purrr::reduce(inner_join) %>%
# # #   bind_rows %>%
# # #   distinct(sample_name, .keep_all = TRUE) %>%
# # #   droplevels %>%
# # #   dplyr::mutate(Fraction = relevel(Fraction, ref = "13H")) %>%
# # #   dplyr::arrange(Fraction, .by_group = TRUE) %>%
# # #   dplyr::mutate(base_12L_frac = ifelse(grepl("A", frac.3), 
# # #                                        as.character(temp_fraction),
# # #                                        as.character(vent_temp_fraction))) %>%
# # #   dplyr::mutate(across(contains("frac", ignore.case = TRUE), ~factor(.x, levels = unique(.x)))) %>%
# # #   droplevels
# # 
# # #norm counts from DESeq's dds:
# # coa_axial_temp_counts_tbl_df <- full_join((coa_axial_30_13H_SDE_counts %>%
# #                                          as.data.frame() %>%
# #                                          tibble::rownames_to_column(., var = "gene_id") %>%
# #                                          tidyr::pivot_longer(., cols = !c("gene_id"),
# #                                                              names_to = "sample_name",
# #                                                              values_to = "counts_norm") %>%
# #                                          dplyr::mutate(baseline = "axial_temp_13H") %>%
# #                                          # dplyr::mutate(baseline = "vent_13H - axial_temp_13H") %>%
# #                                          dplyr::mutate(sampleset = "only_13H")),
# #                                       (coa_axial_55_13H_SDE_counts %>%
# #                                          as.data.frame() %>%
# #                                          tibble::rownames_to_column(., var = "gene_id") %>%
# #                                          tidyr::pivot_longer(., cols = !c("gene_id"),
# #                                                              names_to = "sample_name",
# #                                                              values_to = "counts_norm") %>%
# #                                          dplyr::mutate(baseline = "axial_temp_13H") %>%
# #                                          # dplyr::mutate(baseline = "vent_13H - axial_temp_both") %>%
# #                                          dplyr::mutate(sampleset = "only_13H"))) %>%
# #   bind_rows(., (coa_axial_80_13H_SDE_counts %>%
# #                   as.data.frame() %>%
# #                   tibble::rownames_to_column(., var = "gene_id") %>%
# #                   tidyr::pivot_longer(., cols = !c("gene_id"),
# #                                       names_to = "sample_name",
# #                                       values_to = "counts_norm") %>%
# #                   dplyr::mutate(baseline = "axial_temp_13H") %>%
# #                   # dplyr::mutate(baseline = "vent_13H - axial_temp_12L") %>%
# #                   dplyr::mutate(sampleset = "only_13H"))) %>%
# #   dplyr::mutate(across(c(gene_id, sample_name, 
# #                          # contrast,
# #                          sampleset), ~factor(.x))) %>%
# #   left_join(., metadata %>%
# #               dplyr::select(sample_name, Vent, Temperature, Fraction) %>%
# #               droplevels) %>%
# #   left_join(., vent13H_axial13H_SDE_full,
# #   # left_join(., bind_rows(vent13H_axial13H_SDE_full, vent13H_axialboth_SDE_full, vent13H_axial12L_SDE_full),
# #             by = c("gene_id", "Vent", "Temperature", "baseline"),
# #             relationship = "many-to-many") %>%
# #   dplyr::relocate(c(Vent, Temperature, Fraction), .after = "originalcontrast") %>%
# #   # dplyr::mutate(SDE = ifelse(padj < set.alpha, "sig", NA)) %>%
# #   droplevels
# # 
# # {
# #   
# #   g3 <- ggplot(data = coa_axial_temp_counts_tbl_df %>%
# #                  dplyr::select(gene_id, baseMean, baseline, Temperature) %>%
# #                  distinct(., .keep_all = FALSE) %>%
# #                  tidyr::drop_na(baseMean) %>%
# #                  droplevels)
# #   g3 <- g3 + geom_point(aes(x = gene_id, y = baseMean, shape = Temperature), size = 3, color = "black", 
# #                         show.legend = FALSE,
# #                         alpha = 1, na.rm = TRUE)
# #   g3 <- g3 + scale_y_continuous(labels = scales::label_number_auto(),
# #                                 expand = expansion(mult = 0.2))
# #   g3 <- g3 + scale_shape_manual(values = c(13, 14, 15))
# #   g3 <- g3 + theme_bw()
# #   g3 <- g3 + facet_grid(rows = vars(Temperature),
# #                         shrink = FALSE,
# #                         scales = "free")
# #   # g3 <- g3 + facet_grid(rows = vars(baseline),
# #   #                       labeller = labeller(baseline = contrast_lookup),
# #                         # shrink = FALSE)
# #   g3 <- g3 + theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
# #                    strip.text.y = element_text(angle = 0),
# #                    legend.position = "right")
# #   g3 <- g3 + labs(x = "gene family", 
# #                   y = "normalized counts")
# #   g3 <- g3 + guides(shape = NULL)
# #   # guide_legend("Baseline",
# #   #                                        direction = "vertical",
# #   #                                        override.aes = list(shape = 13, size = 3, alpha = 1)))
# #   g3
# # }
# # ggsave(paste0(getwd(), "/", Sys.Date(), "-", "coa-", "axial_temp_baseline_exp.png"),
# #        g3,
# #        width = 10, height = 10, units = "in")
# # 
# 
# 
# 
# # Count how many genes have padj and pvalue < 0.05 ------------------------
# 
# #find how many SDE genes pass the following criteria:
# #adjusted p-value < 0.05
# #adjusted p-value < 0.001
# #adjusted p-value < 10e-6
# #unadjusted p-value < 10e-6
# 
# #also determine how many have abs(log2FC) > 2
# 
# # in this labeling scheme:
# #A Within vents: baselien all temps, each column is a vent (corresponds to panel A of the orriginal bar chart)
# #B Between vents: baseline each temp, each column is a temp (panel B)
# #C Between years: baseline each temp, each column is a vent/year pair (panel C)
# 
# temp_list <- rlist %>%
#   dplyr::mutate(pvalue_micron = dplyr::case_when((pvalue > 10^-6) ~ NA,
#                                                  .default = pvalue)) %>%
#   dplyr::mutate(padj_micron = dplyr::case_when((padj > 10^-6) ~ NA,
#                                                .default = padj)) %>%
#   dplyr::mutate(padj_milli = dplyr::case_when((padj > 10^-3) ~ NA,
#                                                .default = padj)) %>%
#   dplyr::filter(abs(log2FoldChange) >= 2) %>%
#   dplyr::filter((padj <= 0.05) | (pvalue_micron <= 10^-6)) %>%
#   droplevels %>%
#   dplyr::mutate(baseline = stringr::str_extract(contrast, "axial_.*13H"),
#                 Vent = stringr::str_extract(contrast, "(anemone|marker113|marker33)"),
#                 Temperature = stringr::str_extract(contrast, "(30|55|80)"),
#                 Group = stringr::str_remove(contrast, " - axial_.*13H")) %>%
#   dplyr::mutate(Comparison = dplyr::case_when(grepl("axial_temp_13H", contrast) ~ "A Within vents",
#                                               grepl("(2013|2014)", contrast) ~ "C Between years",
#                                               .default ="B Between vents")) %>%
#   split(., f = .$Comparison) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$Group) %>%
#   map(., ~.x %>%
#         dplyr::select(gene_symbol, contrast, contains("padj"), contains("micron")) %>%
#         tidyr::pivot_longer(., cols = !c(gene_symbol, contrast),
#                             names_to = "metric",
#                             values_to = "value") %>%
#         droplevels %>%
#         split(., f = .$metric) %>%
#         map(., ~.x %>%
#               droplevels %>%
#               tidyr::drop_na(value) %>%
#               dplyr::select(gene_symbol) %>%
#               droplevels %>%
#               simplify)) %>%
#   map2(., names(.),
#        ~.x %>%
#          setNames(., paste0(names(.))))
#   ) |>
#   # map_depth(1, ~.x %>% purrr::keep_at(~ grepl("anemone", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::list_flatten(.))
# 
# temp_vennlist <- temp_list |>
#   map_depth(2, ~RVenn::Venn(.x) %>% ggVennDiagram::process_set_data(.) %>% purrr::pluck("item"))
# 
# 
# 
# myBreaks <- purrr::flatten(temp_vennlist) %>%
#   purrr::flatten() %>%
#   sapply(., length) %>%
#   unlist %>%
#   max %>%
#   signif(., digits = 2)
# myBreaks <- c(seq(from = 0, to = myBreaks,
#                   by = (myBreaks / 100)))
# 
# patch_plots <- map(temp_vennlist %>% names,
#                  ~map2(., names(temp_vennlist[[.x]]),
#                        ~ggVennDiagram::ggVennDiagram(temp_vennlist[[.x]][[.y]], 
#                                                      set_size = NA, 
#                                                      label = "count",
#                                                      label_color = "black",
#                                                      label_geom = "label") +
#                          scale_fill_viridis(discrete = FALSE, option = "turbo",
#                                             limits = c(0, max(myBreaks)),
#                                             guide = "none",
#                                             name = "Gene count") +
#                          scale_x_continuous(expand = expansion(mult = 0.2)) +
#                          ggtitle(.y) +
#                          theme_void()) %>%
#                    purrr::reduce(`%+%`) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = paste0(.x))
# ) %>%
#   setNames(., names(temp_vennlist)) %>%
#   print(.)
# 
# 
# 
# for(i in 1:length(patch_plots)){
#   g <- patch_plots[[i]]
#   namevar <- LETTERS[i]
#   ggsave(paste0(projectpath, "/figures/exploratory/", Sys.Date(), "-", "venn_of_different_padj_", namevar, ".png"),
#          g,
#          width = 12, height = 8, units = "in")
# }
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("A Within", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("anemone", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   purrr::reduce(c) %>%
#   purrr::reduce(c) %>%
#   # length(.) #get the total occurrences of SDE genes
# unique %>% length(.) #get the number of unique genes
# 
# #equivalent to:
# # c(temp_vennlist[["A Within vents"]][["anemone_13H_30"]][["padj_milli"]],
# #           temp_vennlist[["A Within vents"]][["anemone_13H_55"]][["padj_milli"]],
# #           temp_vennlist[["A Within vents"]][["anemone_13H_80"]][["padj_milli"]])
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("A Within", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("113", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   purrr::reduce(c) %>%
#   purrr::reduce(c) %>%
#   # length(.)
#   unique %>% length(.)
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("A Within", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("33", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   purrr::reduce(c) %>%
#   purrr::reduce(c) %>%
#   length(.)
#   # unique %>% length(.)
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("B Between", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("30", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   purrr::reduce(c) %>%
#   purrr::reduce(c) %>%
#   length(.) #get the total occurrences of SDE genes
#   # unique %>% length(.) #get the number of unique genes
# 
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("B Between", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("80", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   purrr::reduce(c) %>%
#   purrr::reduce(c) %>%
#   # length(.) #get the total occurrences of SDE genes
# unique %>% length(.) #get the number of unique genes
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("C Between", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("marker33", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   map_depth(1, ~.x %>% purrr::list_flatten(., name_spec = "{outer}")) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   purrr::reduce(c) %>%
#   purrr::reduce(c) %>%
#   # print
#   # length(.) #get the total occurrences of SDE genes
#   unique %>% length(.) #get the number of unique genes
# 
# temp_vennlist |>
#   purrr::keep_at(~ grepl("C Between", .)) |>
#   map_depth(1, ~.x %>% purrr::keep_at(~ grepl("marker113_55", .))) |>  #if you want to just look at one vent
#   map_depth(2, ~.x %>% purrr::keep_at(~ grepl("padj_milli", .))) |>
#   map_depth(1, ~.x %>% purrr::list_flatten(., name_spec = "{outer}")) |>
#   c() %>%
#   purrr::list_flatten(.) %>%
#   # purrr::reduce(c) %>%
#   # purrr::reduce(c) %>%
#   print
#   # length(.) #get the total occurrences of SDE genes
# # unique %>% length(.) #get the number of unique genes
# 
# #what are theye?
# 
# full_functions_list <- readr::read_delim(paste0(projectpath, "/data/reference/", "full_functions_list", ".tsv"),
#                                          show_col_types = FALSE,
#                    delim = "\t", col_names = TRUE)
# setdiff(temp_vennlist[["A Within vents"]][["marker113_13H_30"]][["padj_milli"]], temp_vennlist[["A Within vents"]][["marker113_13H_30"]][["pvalue_micron"]])
# temp_vennlist[["A Within vents"]][["marker113_13H_30"]][["pvalue_micron"]]
# 
# full_functions_list %>%
#   dplyr::filter(gene_id %in% setdiff(temp_vennlist[["A Within vents"]][["marker113_13H_30"]][["padj_milli"]], temp_vennlist[["A Within vents"]][["marker113_13H_30"]][["pvalue_micron"]])) %>%
#   distinct(.)
# 
# full_functions_list %>%
#   dplyr::filter(gene_id %in%  temp_vennlist[["A Within vents"]][["marker113_13H_30"]][["pvalue_micron"]]) %>%
#   distinct(.)
# 
# 
# 
#             