# 01_DifferentialExpression.R

# 01_metat_coa_deseq.R
# this script begins DeSeq
# this script follows:
# 00_process_coa_mtgs.R
# 00_process_coa_mtts.R

#this script takes the combined counts from the metatranscriptomes
#and processed them through DESeq
#calculating baseline expressions from either:
#12L+13H (made irrelevant)
#13H only

#and for each group of fractions, using the genes from either
#all temperatures
#or temperature-specific (30, 55, 80)


# Resource allocation time ------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_resource_allocation.R"))){
  cat("Preparing resource allocations.")
  suppressWarnings(
    source(paste0(getwd(), "/", "00_resource_allocation.R"), local = FALSE,
           verbose = getOption("verbose"), prompt.echo = getOption("prompt"),
           # verbose = FALSE, prompt.echo = FALSE,
           echo = FALSE),
    classes = "warnings")
  
} else {
  cat("You have to manually load these packages and set paths.")
  # #packages needed for this script:
  
  library(tidyverse)
  library(data.table)
  library(BiocManager)
  library(BiocParallel)
  library(DESeq2)
  library(cli)
  library(furrr)
  library(progressr)
  
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


no_progress <- function() {} #your placeholder function for progress reporting

coalesce2 <- function(x, y, sep = ".") ifelse(x == y, coalesce(x, y, sep = sep), paste0(x, "_vs_", y))

F_counts_combiner <- function(x){
  p <- progressr::progressor(along = x)
  for(group in seq_along(x)){
    df <- get(x[group], inherits = TRUE) %>%
      tibble::rownames_to_column(., var = "gene_id") %>%
      dplyr::mutate(gene_id = factor(gene_id))
    
    Sys.sleep(1/100)
    if(exists("temp_combo")){
      temp_combo <- full_join(temp_combo,
                              df,
                              # keep = TRUE,
                              by = "gene_id",
                              na_matches = "never")
    }
    if(!exists("temp_combo")){
      temp_combo <- df  
    }
  }
  temp_combo
}

F_deseq_frac <- function(counts, mdata, mdata_filter, deseq_design, prefix,
                         p = no_progress){
  #here is a generalized function to make the comparisons where the baseline is
  #13H + 12L
  #13H only
  #12L only
  
  #counts : your counts table
  #mdata: your metadata
  #mdata_filter: optional selector of your data
  #provide the filter a list of named vectors
  # metadata_filter <- list(`Fraction` = c("12L", "13H"))
  # metadata_filter <- list(`Type` = c("SIP"),
  #                         `Temperature` = c("30; 55"))
  #due to DESeq's design, all samples included via this mdata_filter will be used
  #to calculate a baseline expression for the genes
  #if you don't want that, you will need to modify your filter, or account for it post-hoc
  #deseq_design: the formula
  #`prefix`: name prefix for the output list
  if(!missing(prefix)){
    prefix <- rlang::parse_expr(prefix)
  } else {
    prefix <- "default"
  }
  
  if(exists(paste0(prefix, "_deseq_list"), envir = .GlobalEnv)){
    cli::cli_alert_warning("A DESeq list object already exists with this name: {prefix}_deseq_list", wrap = TRUE)
    cli::cli_abort("Ending DESeq attempt.")
  }
  
  #d
  deseq_design <- rlang::parse_expr(deseq_design)
  
  #m
  if(exists("mdata_filter", inherits = TRUE)){
    metadata <- map2(mdata_filter, names(mdata_filter),
                     ~mdata %>%
                       dplyr::filter(
                         grepl(paste0(.x, collapse = "|"), .data[[.y]], ignore.case = TRUE)
                       )) %>%
      purrr::reduce(inner_join) %>%
      bind_rows %>%
      distinct(sample_name, .keep_all = TRUE) %>%
      droplevels
  } else {
    metadata <- mdata
  }
  
  if(any(grepl("13H", mdata_filter)) & ("13H" == mdata_filter$Fraction[1])){
    # if(any(grepl("13H", mdata_filter))){
    metadata <- metadata %>%
      dplyr::mutate(Fraction = relevel(Fraction, ref = "13H")) %>%
      dplyr::arrange(Fraction, .by_group = TRUE) %>%
      dplyr::mutate(across(contains("frac", ignore.case = TRUE), ~factor(.x, levels = unique(.x)))) %>%
      droplevels
  } else {
    metadata <- metadata %>%
      dplyr::mutate(Fraction = relevel(Fraction, ref = "12L")) %>%
      dplyr::arrange(Fraction, .by_group = TRUE) %>%
      dplyr::mutate(base_12L_frac = ifelse(grepl("A", frac.3), 
                                           as.character(temp_fraction),
                                           as.character(vent_temp_fraction))) %>%
      dplyr::mutate(across(contains("frac", ignore.case = TRUE), ~factor(.x, levels = unique(.x)))) %>%
      droplevels
  } 
  
  
  #check the design formula against the available metadata to make sure we have enough samples
  mod_mat <- model.matrix(eval(deseq_design),
                          data = metadata)
  
  if(length(grep("0", colSums(mod_mat))) > 0){ #none of these should be 0--if so, you need to fix your design/metadata coding
    idx <- names(colSums(mod_mat))[grep("0", colSums(mod_mat))]
    idx <- cli::cli_vec(idx)
    cli::cli_alert_warning("One or more groups in the comparison does not have enough replicates:", wrap = TRUE)
    cli::cli_bullets_raw(idx)
    cli::cli_abort("Ending DESeq attempt.")
  }
  
  #data
  if(!("sample_name" %in% colnames(counts))){
    cli::cli_alert_info("Formatting the counts data to be a long table.")
    
    counts <- counts %>%
      tidyr::pivot_longer(., 
                          cols = !contains("gene"),
                          names_to = "sample_name",
                          values_to = "counts")
    assign("counts", counts, envir = .GlobalEnv, inherits = FALSE)
  }
  matrix <- counts %>%
    dplyr::mutate(across(!contains("counts"), factor)) %>%
    dplyr::filter(sample_name %in% as.character(metadata[["sample_name"]])) %>%
    droplevels %>%
    tidyr::drop_na() %>%
    pivot_wider(., id_cols = "gene_id",
                names_from = "sample_name",
                values_fill = 0,
                values_from = "counts") %>%
    column_to_rownames(., var = "gene_id") %>%
    dplyr::select(metadata[["sample_name"]]) %>%
    dplyr::slice(which(rowSums(.) > 0))
  
  dds <- DESeqDataSetFromMatrix(countData = matrix,
                                tidy = FALSE,
                                colData = metadata,
                                design = eval(deseq_design)) %>%
    DESeq(., parallel = TRUE,
          BPPARAM = bpparam_multi)
  
  # resultsNames(dds)
  # 
  # mod_mat <- model.matrix(design(dds), colData(dds))
  # 
  temp_list <- list(metadata, matrix, dds, mod_mat) %>%
    setNames(., c("metadata", "matrix", "dds", "mod_mat"))
  assign(paste0(prefix, "_deseq_list"), temp_list, envir = .GlobalEnv)
}

# Lookup vectors ----------------------------------------------------------

vents <- c("anemone", "marker113", "marker33")
set.alpha <- 0.05


# Import sample metadata --------------------------------------------------

sample_metadata <- readr::read_delim(paste0(projectpath, "/data/", "sample_metadata.tsv"),
                                     col_names = TRUE,
                                     show_col_types = FALSE,
                                     quote = "",
                                     comment = "#",
                                     delim = "\t",
                                     num_threads = nthreads) %>%
  dplyr::mutate(Temperature = ifelse(is.na(Temperature), "00", Temperature),
                Fraction = ifelse(is.na(Fraction), "W", Fraction)) %>%
  dplyr::mutate(across(everything(), ~as.factor(.x))) %>%
  tidyr::drop_na(Vent) %>%
  dplyr::mutate(Temperature = fct_relevel(Temperature, "00"),
                Fraction = fct_relevel(Fraction, "W"),
                Vent = fct_relevel(Vent, "anemone"),
                Timepoint = fct_relevel(Timepoint, "None"),
                Type = fct_relevel(Type, "Untreated")) %>%
  dplyr::mutate(Vent = gsub("_", "", Vent)) %>%
  dplyr::mutate(vent.1 = dplyr::case_match(Vent,
                                           "anemone" ~ "N",
                                           "marker113" ~ "O",
                                           "marker33" ~ "T")) %>%
  dplyr::mutate(temp.2 = dplyr::case_match(Temperature, 
                                           "00" ~ "0",
                                           "30" ~ "1",
                                           "55" ~ "2",
                                           "80" ~ "3")) %>%
  dplyr::mutate(frac.3 = dplyr::case_match(Fraction,
                                           "W" ~ "W",
                                           "12L" ~ "A",
                                           "13H" ~ "B",
                                           "12H" ~ "C",
                                           "13L" ~ "D")) %>%
  dplyr::mutate(time.4 = dplyr::case_match(Timepoint,
                                           "None" ~ "0",
                                           "TP1" ~ "1",
                                           "TP2" ~ "2")) %>%
  dplyr::mutate(year.5 = dplyr::case_match(year,
                                           "2013" ~ "3",
                                           "2014" ~ "4")) %>%
  dplyr::mutate(vent.1 = fct_relevel(vent.1, "N"),
                temp.2 = fct_relevel(temp.2, "0"),
                frac.3 = fct_relevel(frac.3, "W"),
                time.4 = fct_relevel(time.4, "0"),
                year.5 = fct_relevel(year.5, "3")) %>%
  dplyr::mutate(temp_fraction = interaction(Temperature, frac.3)) %>%
  droplevels %>%
  dplyr::mutate(temp_fraction = factor(temp_fraction, 
                                       levels = levels(.$temp_fraction),
                                       labels = c(LETTERS[1:length(levels(.$temp_fraction))]))) %>%
  droplevels %>%
  dplyr::mutate(vent_temp_fraction = paste0(vent.1, temp_fraction) %>%
                  factor(.)) %>%
  dplyr::mutate(vent_year_temp_fraction = paste0(vent.1, year.5, temp_fraction) %>%
                  factor(.)) %>%
  dplyr::mutate(temp_fraction_timepoint = interaction(temp_fraction, time.4, sep = "")) %>%
  droplevels %>%
  dplyr::mutate(vent_temp_fraction_timepoint = paste0(vent.1, temp_fraction_timepoint) %>%
                  factor(.)) %>%
  dplyr::mutate(vent_year_temp_fraction_timepoint = paste0(vent.1, year.5, temp_fraction_timepoint) %>%
                  factor(.)) %>%
  droplevels %>%
  dplyr::mutate(across(everything(), ~as.factor(.x)))

# Import counts matrices --------------------------------------------------

to_combine <- NULL
to_import_txi <- NULL
to_find_txi <- NULL

if(!exists("coa_counts_combined", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "coa_counts_combined", ".rds"))){
    cli::cli_alert_info("Reading in the combined counts table...")
    progressr::with_progress({
      coa_counts_combined <- readr::read_rds(file = paste0(projectpath, "/data/", "coa_counts_combined", ".rds"))
    })
  } else {
    cli::cli_alert_info("Processing the counts tables and combining them...")
    
    to_combine <- levels(sample_metadata$sample_source) %>%
      paste0("txi_df_", .) %>%
      paste0("coa_", .)
    
    cli::cli_progress_bar("Converting", total = length(to_combine), type = "tasks")
    for(i in 1:length(to_combine)){
      
      txi_df <- to_combine[i]
      Sys.sleep(1/100)
      while(!exists(txi_df, inherits = TRUE)){
        if(file.exists(paste0(projectpath, "/data/", txi_df, ".rds"))){
          progressr::with_progress({
            temp_df <- readr::read_rds(file = paste0(projectpath, "/data/", txi_df, ".rds"))
            assign(txi_df, temp_df, envir = .GlobalEnv)
          })
        } else {
          cli::cli_alert_danger("The tximport file does not exist for this series: {txi_df}")
          to_find_txi <- c(to_find_txi, txi_df)
          break
        }
        to_import_txi <- c(to_import_txi, txi_df)
      # } else {
      #   cli::cli_alert_info("Dataframe for {txi_df} already loaded.")
      }
      # to_import_txi <- to_import_txi[grep(txi_df, to_import_txi, value = FALSE, invert = TRUE)]

      cli::cli_progress_update()
      rm(txi_df)
    }
    cli::cli_alert_success("Imported {length(to_import_txi)} datasets..")
    cli::cli_progress_done() 
    
    if(exists("temp_combo", inherits = TRUE)){
      rm(temp_combo)  
    } else {
      NULL
    }
    
    coa_counts_combined <- F_counts_combiner(to_combine) %>%
      dplyr::mutate(across(everything(), ~tidyr::replace_na(.x, replace = 0)))
    readr::write_rds(coa_counts_combined, file = paste0(projectpath, "/data/", "coa_counts_combined", ".rds"))
  }
  if(length(to_find_txi) > 1){
    cli::cli_alert_danger("You are missing datasets to import: {.val {to_find_txi}}.", wrap = TRUE)
  }
} else {
  cli::cli_alert_info("Combined transcript counts table already imported: coa_counts_combined")
  idx <- grep("gene_id", colnames(coa_counts_combined),invert = TRUE, value = TRUE)
  if(!(idx %in% sample_metadata[["sample_name"]])){
    missing_idx <- setdiff(idx, sample_metadata[["sample_name"]])
    cli::cli_alert_warning("You are missing some samples from the imported combined counts table. {.val {missing_idx}}.")
  } else {
    NULL
  }
}


# Run DESeqDataSetFromMatrix and DESeq ----------------------------------------
# this list provides the DESeq model objects for specific comparisons:
# "axial_temp_13H": all 13C-heavy fractions of Marker113, Marker33, and Anemone vent fluids incubated at 30˚C, 55˚C, and 80˚C
# "axial_30_13H": all 13C-heavy fractions of Marker113, Marker33, and Anemone vent fluids incubated at 30˚C
# "axial_55_13H": all 13C-heavy fractions of Marker113, Marker33, and Anemone vent fluids incubated at 55˚C
# "axial_80_13H": all 13C-heavy fractions of Marker113, Marker33, and Anemone vent fluids incubated at 80˚C
# "axial_30_13H_vent_year": comparison between 2013 and 2014 of vent- and temp- specific 13C-heavy fractions incubated at 30˚C (Marker33)
# "axial_55_13H_vent_year": comparison between 2013 and 2014 of vent- and temp- specific 13C-heavy fractions incubated at 55˚C (Marker113, Marker33)
# "axial_80_13H_vent_year": comparison between 2013 and 2014 of vent- and temp- specific 13C-heavy fractions incubated at 80˚C (Marker113, Marker33, Anemone)
# (supplementary) "axial_temp_bothfrac": all 13C-heavy and 12C-light fractions of Marker113, Marker33, and Anemone vent fluids incubated at 30˚C, 55˚C, and 80˚C
# (supplementary) "axial_30_bothfrac": all 13C-heavy and 12C-light fractions of Marker113, Marker33, and Anemone vent fluids incubated at 30˚C
# (supplementary) "axial_55_bothfac": all 13C-heavy and 12C-light fractions of Marker113, Marker33, and Anemone vent fluids incubated at 55˚C
# (supplementary) "axial_80_bothfrac": all 13C-heavy and 12C-light fractions of Marker113, Marker33, and Anemone vent fluids incubated at 80˚C
# (deprecated) "axial_temp_12L": all 12C-light fractions of Marker113, Marker33, and Anemone vent fluids incubated at 30˚C, 55˚C, and 80˚C
# "axial_eachtemp_13H": extract the temp-specific comparisons and compile in a new R object

to_process_results <- NULL
to_process_deseq <- NULL

for(i in c("axial_temp_13H", "axial_30_13H", "axial_55_13H", "axial_80_13H",
           # "axial_temp_12L","axial_temp_bothfrac",
           # "axial_eachtemp_13H",
           "axial_temp_bothfrac", "axial_30_bothfrac", "axial_55_bothfrac", "axial_80_bothfrac",
           "axial_30_13H_vent_year", "axial_55_13H_vent_year", "axial_80_13H_vent_year"
           )){
  namevar <- paste0("coa_", i)
  if(!file.exists(paste0(projectpath, "/data/", namevar, "_deseq_list", ".rds"))){
    cli::cli_alert_info("DESeq model object does not exist with this name: {namevar}. Please process through DESeq.")
    if(exists("to_process_deseq", envir = .GlobalEnv)){
      to_process_deseq <- c(to_process_deseq, namevar)
    } else {
      to_process_deseq <- c(namevar)
    }
  } 
    if(!file.exists(paste0(projectpath, "/data/", namevar, "_res_list", ".rds"))){
      cli::cli_alert_info("DESeq results list does not exist with this name: {namevar}. Please process through DESeq.")
      if(exists("to_process_results", envir = .GlobalEnv)){
        to_process_results <- c(to_process_results, namevar)
      } else {
        to_process_results <- c(namevar)
      }
    }
  rm(i)
}



if(!is.null(to_process_deseq)){
  
  for(deseq_name in to_process_deseq){
  # for(deseq_name in to_process_deseq[5]){
    ndataset <- match(deseq_name, to_process_deseq)
      cli::cli_alert_info("{ndataset}: Processing through DESeq: {deseq_name}")
      
    metadata_filter <- NULL
    if(any(grepl("13H", deseq_name))){
      metadata_filter <- list(`Fraction` = c("13H"))
    } else {
      metadata_filter <- list(`Fraction` = c("13H", "12L"))
    }
    if(any(grepl("30|55|80", deseq_name))){
      deseq_temp <- stringr::str_extract(deseq_name, "30|55|80") %>%
        as.numeric
      metadata_filter <- append(metadata_filter,
                                list(`Temperature` = deseq_temp))
    } 
    if(any(grepl("vent_year", deseq_name))){ #only working with 13H samples here
      model_design = "~ 0 + vent_year_temp_fraction" 
    } else if(any(grepl("temp_13H", deseq_name))){ #all temps, just 13H
      model_design = "~ 0 + Vent:Temperature"
    } else if(any(grepl("bothfrac", deseq_name))){ #both fractions
      if(any(grepl("temp", deseq_name))){ model_design = "~ 0 + Fraction:Vent:Temperature"  } #all-temps
      else { model_design = "~ 0 + Fraction:Vent" } #temp-specific
    } else {
      model_design = "~ 0 + Vent"
    }
    cli::cli_alert_info("\n using only these sample types: {metadata_filter} \n and this model formula: {model_design}")
    F_deseq_frac(counts = coa_counts_combined,
                 mdata = sample_metadata,
                 mdata_filter = metadata_filter,
                 deseq_design = model_design,
                 # prefix = deseq_name)
                 prefix = "temp")
    assign(paste0(deseq_name, "_deseq_list"), temp_deseq_list, envir = .GlobalEnv)
    readr::write_rds(temp_deseq_list, file = paste0(projectpath, "/data/", deseq_name, "_deseq_list", ".rds"))
    rm(temp_deseq_list)
  }
  #one by one:
  #supplementary:
  {
    # if(any(grepl("coa_axial_temp_bothfrac", to_process_deseq))){
    #   prefix <- "coa_axial_temp_bothfrac"
    #   cli::cli_alert_info("Processing through DESeq: {prefix}")
    #   metadata_filter <- list(`Fraction` = c("13H", "12L"))
    #   F_deseq_frac(counts = coa_counts_combined,
    #                mdata = sample_metadata,
    #                mdata_filter = metadata_filter,
    #                deseq_design = "~ 0 + Fraction:Vent:Temperature",
    #                prefix = "coa_axial_temp_bothfrac")
    #   readr::write_rds(coa_axial_temp_bothfrac_deseq_list, file = paste0(projectpath, "/", "coa_", "axial_temp_bothfrac_deseq_list", ".rds"))
    # 
    # }
    #deprecated:
    # if(any(grepl("coa_axial_temp_12L", to_process_deseq))){
    #   prefix <- "coa_axial_temp_12L"
    #   cli::cli_alert_info("Processing through DESeq: {prefix}")
    #   metadata_filter <- list(`Fraction` = c("12L", "13H"))
    #   F_deseq_frac(counts = coa_counts_combined,
    #                mdata = sample_metadata,
    #                mdata_filter = metadata_filter,
    #                deseq_design = "~ 1 + base_12L_frac",
    #                prefix = "coa_axial_temp_12L")
    #   readr::write_rds(coa_axial_temp_12L_deseq_list, file = paste0(projectpath, "/data/", "coa_", "axial_temp_12L_deseq_list", ".rds"))
    # 
    # }
    # if(any(grepl("coa_axial_30_bothfrac", to_process_deseq))){
    #   deseq_name <- "coa_axial_30_bothfrac"
    #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
    #   metadata_filter <- list(`Fraction` = c("13H", "12L"),
    #                           `Temperature` = 30)
    #   F_deseq_frac(counts = coa_counts_combined,
    #                mdata = sample_metadata,
    #                mdata_filter = metadata_filter,
    #                deseq_design = "~ 0 + Fraction:Vent",
    #                # deseq_design = "~ 0 + Vent:Fraction",
    #                prefix = deseq_name)
    #                # prefix = "coa_axial_30_bothfrac")
    #   readr::write_rds(coa_axial_30_bothfrac_deseq_list, file = paste0(projectpath, "/", "coa_axial_30_bothfrac_deseq_list", ".rds"))
    # }
    # if(any(grepl("coa_axial_55_bothfrac", to_process_deseq))){
    #   deseq_name <- "coa_axial_55_bothfrac"
    #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
    #   metadata_filter <- list(`Fraction` = c("13H", "12L"),
    #                           `Temperature` = 55)
    #   F_deseq_frac(counts = coa_counts_combined,
    #                mdata = sample_metadata,
    #                mdata_filter = metadata_filter,
    #                deseq_design = "~ 0 + Fraction:Vent",
    #                prefix = deseq_name)
    #                # prefix = "coa_axial_55_bothfrac")
    #   readr::write_rds(coa_axial_55_bothfrac_deseq_list, file = paste0(projectpath, "/", "coa_axial_55_bothfrac_deseq_list", ".rds"))
    # }
    # if(any(grepl("coa_axial_80_bothfrac", to_process_deseq))){
    #   deseq_name <- "coa_axial_80_bothfrac"
    #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
    #   metadata_filter <- list(`Fraction` = c("13H", "12L"),
    #                           `Temperature` = 80)
    #   F_deseq_frac(counts = coa_counts_combined,
    #                mdata = sample_metadata,
    #                mdata_filter = metadata_filter,
    #                deseq_design = "~ 0 + Fraction:Vent",
    #                # deseq_design = "~ 0 + Vent:Fraction",
    #                # prefix = "coa_axial_80_bothfrac")
    #                prefix = deseq_name)
    #   readr::write_rds(coa_axial_80_bothfrac_deseq_list, file = paste0(projectpath, "/", "coa_axial_80_bothfrac_deseq_list", ".rds"))
    # }
  #   
  # if(any(grepl("coa_axial_temp_13H", to_process_deseq))){
  #   deseq_name <- "coa_axial_temp_13H"
  #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
  #   metadata_filter <- list(`Fraction` = c("13H"))
  #   F_deseq_frac(counts = coa_counts_combined,
  #                mdata = sample_metadata,
  #                mdata_filter = metadata_filter,
  #                deseq_design = "~ 0 + Vent:Temperature",
  #                prefix = deseq_name)
  #                # prefix = "coa_axial_temp_13H")
  #   readr::write_rds(coa_axial_temp_13H_deseq_list, file = paste0(projectpath, "/data/", "coa_", "axial_temp_13H_deseq_list", ".rds"))
  # }
  # 
  # if(any(grepl("coa_axial_30_13H", to_process_deseq))){
  #   deseq_name <- "coa_axial_30_13H"
  #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
  #   metadata_filter <- list(`Fraction` = c("13H"),
  #                           `Temperature` = 30)
  #   F_deseq_frac(counts = coa_counts_combined,
  #                mdata = sample_metadata,
  #                mdata_filter = metadata_filter,
  #                deseq_design = "~ 0 + Vent",
  #                prefix = deseq_name)
  #                # prefix = "coa_axial_30_13H")
  #   readr::write_rds(coa_axial_30_13H_deseq_list, file = paste0(projectpath, "/data/", "coa_axial_30_13H_deseq_list", ".rds"))
  # }
  # if(any(grepl("coa_axial_55_13H", to_process_deseq))){
  #   deseq_name <- "coa_axial_55_13H"
  #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
  #   metadata_filter <- list(`Fraction` = c("13H"),
  #                           `Temperature` = 55)
  #   F_deseq_frac(counts = coa_counts_combined,
  #                mdata = sample_metadata,
  #                mdata_filter = metadata_filter,
  #                deseq_design = "~ 0 + Vent",
  #                prefix = deseq_name)
  #                # prefix = "coa_axial_55_13H")
  #   readr::write_rds(coa_axial_55_13H_deseq_list, file = paste0(projectpath, "/data/", "coa_axial_55_13H_deseq_list", ".rds"))
  # }
  # if(any(grepl("coa_axial_80_13H", to_process_deseq))){
  #   deseq_name <- "coa_axial_80_13H"
  #   cli::cli_alert_info("Processing through DESeq: {deseq_name}")
  #   metadata_filter <- list(`Fraction` = c("13H"),
  #                           `Temperature` = 80)
  #   F_deseq_frac(counts = coa_counts_combined,
  #                mdata = sample_metadata,
  #                mdata_filter = metadata_filter,
  #                deseq_design = "~ 0 + Vent",
  #                prefix = deseq_name)
  #                # prefix = "coa_axial_80_13H")
  #   readr::write_rds(coa_axial_80_13H_deseq_list, file = paste0(projectpath, "/data/", "coa_axial_80_13H_deseq_list", ".rds"))
  # }
  }
}





# Extract significant results ---------------------------------------------

## for troubleshooting:
# to_process_results <- grep("bothfrac", to_process_results, value = TRUE)

if(!is.null(to_process_results)){

  for(res_name in to_process_results){
    ndataset <- match(res_name, to_process_results)
    condition_filter <- NULL
    baseline_filter <- NULL
    res <- NULL 
    contrast_id <- NULL
    
    cli::cli_alert_info("{ndataset}: Extracting DESeq results: {res_name}")
    
    #supplementary: 13H and 12L
    if(any(grepl("bothfrac", res_name))){ #supplementary: 13H and 12L
      #this gets temp_bothfrac, 30_bothfrac, 55_bothfrac, and 80_bothfrac 
      baseline_filter <- list(`Fraction` = c("13H", "12L"))
      condition_filter <- list(`Vent` = c("anemone", "marker113", "marker33"),
                               `Fraction` = c("13H"))
      
      #all temps, both fractions
      if(any(grepl("temp", res_name))){ #all temps, both fractions
        #identify the conditions you are testing against the baseline, and determine interactions = contrasts
        condition_filter <- append(condition_filter,
                                   list(`Temperature` = c(30, 55, 80)))
        prefix <- "coa_axial_temp_bothfrac"
        list2env(coa_axial_temp_bothfrac_deseq_list, envir = .GlobalEnv)
        baseline_idx <- map2(baseline_filter, names(baseline_filter),
                             ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
          purrr::reduce(inner_join)
        
        baseline_vec <- colMeans(mod_mat[baseline_idx, ])
        
        condition_idx_names <- map(condition_filter,
                                   ~unlist(., recursive = FALSE)) %>%
          purrr::reduce(interaction, sep = "_") %>%
          levels %>%
          as.character
        
        #brute force: list the model vectors
        for(i in vents){
          namevar <- i
          for(j in c(30, 55, 80)){
            tempvar <- j
            if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar,]))){
              temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ]
            } else {
              temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ])
            }
            assign(paste0(namevar, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
            rm(temp_vec)
            rm(tempvar)
          }
          rm(namevar)
        }
        
        contrast_id <- condition_idx_names %>%
          as.data.frame() %>%
          dplyr::mutate(baseline = "baseline_vec") %>%
          dplyr::mutate(across(everything(), as.character))
        
        label_id <- paste0(condition_idx_names, " - axial_temp_", "both")
        
        y <- nrow(contrast_id)
        for(i in seq(y)){
          Sys.sleep(1/100)
          pair1 <- get0(contrast_id[i,1], inherits = TRUE)
          pair2 <- get0(contrast_id[i,2], inherits = TRUE)
          res1 <- results(dds,
                          contrast = pair1 - pair2,
                          parallel = TRUE,
                          BPPARAM = bpparam_multi,
                          alpha = set.alpha) %>%
            as.data.frame() %>%
            rownames_to_column("gene_symbol") %>%
            left_join(., lfcShrink(dds,
                                   contrast = pair1 - pair2,
                                   parallel = TRUE,
                                   BPPARAM = bpparam_multi,
                                   type = "ashr",
                                   quiet = TRUE) %>%
                        as.data.frame() %>%
                        rownames_to_column("gene_symbol") %>%
                        dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                        dplyr::select(gene_symbol, log2FoldChange_MMSE),
                      by = "gene_symbol") %>%
            mutate(contrast = label_id[i]) %>%
            relocate(contrast, .before = "gene_symbol")
          if(exists("res")){
            res <- bind_rows(res, res1)
            rm(res1)}
          if(!exists("res")){
            res <- res1
            rm(res1)}
        }
        
        assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
        readr::write_rds(coa_axial_temp_bothfrac_res_list, file = paste0(projectpath, "/data/", "coa_axial_temp_bothfrac_res_list", ".rds"))
        to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
        rm(prefix)
        
      } #temp-specific, both fractions
      else { #temp-specific, both fractions
        #do the temp-specific 13H+12L result extraction
        baseline_filter <- list(`Fraction` = c("13H", "12L"))
        tempid <- stringr::str_extract(res_name, "(30|55|80)") %>%
          as.numeric(.)
        
        res <- NULL
        prefix <- paste0("coa_axial_", tempid, "_bothfrac")
        condition_filter <- append(condition_filter,
                                   list(`Temperature` = tempid))
        
        deseq_list <- get0(paste0(prefix, "_deseq_list"), envir = .GlobalEnv, inherits = TRUE)
        list2env(deseq_list, envir = .GlobalEnv)
        rm(deseq_list)
        
        baseline_idx <- map2(baseline_filter, names(baseline_filter),
                             ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
          purrr::reduce(inner_join)
        
        baseline_vec <- colMeans(mod_mat[baseline_idx, ])
        
        condition_idx_names <- map(condition_filter,
                                   ~unlist(., recursive = FALSE)) %>%
          purrr::reduce(interaction, sep = "_") %>%
          levels %>%
          as.character
        
        #brute force: list the model vectors
        for(i in vents){
          namevar <- i
          for(j in c(30, 55, 80)){
            tempvar <- j
            if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar,]))){
              temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ]
            } else {
              temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ])
            }
            assign(paste0(namevar, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
            rm(temp_vec)
            rm(tempvar)
          }
          rm(namevar)
        }
        
        contrast_id <- condition_idx_names %>%
          as.data.frame() %>%
          dplyr::mutate(baseline = "baseline_vec") %>%
          dplyr::mutate(across(everything(), as.character))
        
        label_id <- paste0(condition_idx_names, " - axial_", tempid, "_both")
        
        y <- nrow(contrast_id)
        for(i in seq(y)){
          Sys.sleep(1/100)
          pair1 <- get0(contrast_id[i,1], inherits = TRUE)
          pair2 <- get0(contrast_id[i,2], inherits = TRUE)
          res1 <- results(dds,
                          contrast = pair1 - pair2,
                          parallel = TRUE,
                          BPPARAM = bpparam_multi,
                          alpha = set.alpha) %>%
            as.data.frame(.) %>%
            rownames_to_column("gene_symbol") %>%
            left_join(., lfcShrink(dds,
                                   contrast = pair1 - pair2,
                                   parallel = TRUE,
                                   BPPARAM = bpparam_multi,
                                   type = "ashr",
                                   quiet = TRUE) %>%
                        as.data.frame(.) %>%
                        rownames_to_column("gene_symbol") %>%
                        dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                        dplyr::select(gene_symbol, log2FoldChange_MMSE),
                      by = "gene_symbol") %>%
            mutate(contrast = label_id[i]) %>%
            relocate(contrast, .before = "gene_symbol")
          
          if(exists("res")){
            res <- bind_rows(res, res1)
            rm(res1)}
          if(!exists("res")){
            res <- res1
            rm(res1)}
        }
        assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
        readr::write_rds(res, file = paste0(projectpath, "/data/", prefix, "_res_list", ".rds"))
        to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
        rm(prefix)
      }
    }
    
    #main analyses: 13H only
    if(any(grepl("13H", res_name))){ #this gets temp_13H, 30_13H, 55_13H, and 80_13H, and the 13H-vent-temp specific temporal comparisons 
      res <- NULL
      baseline_filter <- list(`Fraction` = c("13H"))
      condition_filter <- list(`Vent` = c("anemone", "marker113", "marker33"),
                               `Fraction` = c("13H"))
      
      #all-temps, 13H only
      if(any(grepl("temp", res_name))){ #all-temps, 13H only
        condition_filter <- append(condition_filter,
                                   list(`Temperature` = c(30, 55, 80)))
        prefix <- "coa_axial_temp_13H"
        list2env(coa_axial_temp_13H_deseq_list, envir = .GlobalEnv)
        
        baseline_idx <- map2(baseline_filter, names(baseline_filter),
                             ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
          purrr::reduce(inner_join)
        
        baseline_vec <- colMeans(mod_mat[baseline_idx, ])
        #identify the conditions you are testing against the baseline, and determine interactions = contrasts
        
        condition_idx_names <- map(condition_filter,
                                   ~unlist(., recursive = FALSE)) %>%
          purrr::reduce(interaction, sep = "_") %>%
          levels %>%
          as.character
        
        #brute force: list the model vectors
        for(i in vents){
          namevar <- i
          for(j in c(30, 55, 80)){
            tempvar <- j
            if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar,]))){
              temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ]
            } else {
              temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ])
            }
            assign(paste0(namevar, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
            rm(temp_vec)
            rm(tempvar)
          }
          rm(namevar)
        }
        
        contrast_id <- condition_idx_names %>%
          as.data.frame() %>%
          dplyr::mutate(baseline = "baseline_vec") %>%
          dplyr::mutate(across(everything(), as.character))
        
        label_id <- paste0(condition_idx_names, " - axial_temp_", "13H")
        
        y <- nrow(contrast_id)
        for(i in seq(y)){
          Sys.sleep(1/100)
          pair1 <- get0(contrast_id[i,1], inherits = TRUE)
          pair2 <- get0(contrast_id[i,2], inherits = TRUE)
          res1 <- results(dds,
                          contrast = pair1 - pair2,
                          parallel = TRUE,
                          BPPARAM = bpparam_multi,
                          alpha = set.alpha) %>%
            as.data.frame() %>%
            rownames_to_column("gene_symbol") %>%
            left_join(., lfcShrink(dds,
                                   contrast = pair1 - pair2,
                                   parallel = TRUE,
                                   BPPARAM = bpparam_multi,
                                   type = "ashr",
                                   quiet = TRUE) %>%
                        as.data.frame() %>%
                        rownames_to_column("gene_symbol") %>%
                        dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                        dplyr::select(gene_symbol, log2FoldChange_MMSE),
                      by = "gene_symbol") %>%
            mutate(contrast = label_id[i]) %>%
            relocate(contrast, .before = "gene_symbol")
          if(exists("res")){
            res <- bind_rows(res, res1)
            rm(res1)}
          if(!exists("res")){
            res <- res1
            rm(res1)}
        }
        assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
        readr::write_rds(res, file = paste0(projectpath, "/data/", prefix, "_res_list", ".rds"))
        to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
        rm(prefix)
        
      } 
      # #vent-specific, temp-specific temporal comparison
      else if(any(grepl("vent_year", res_name))){ # #vent-specific, temp-specific temporal comparison
        
        tempid <- stringr::str_extract(res_name, "(30|55|80)") %>%
          as.numeric(.)
        contrast_id <- NULL
        prefix <- paste0("coa_axial_", tempid, "_13H_vent_year")
        deseq_list <- get0(paste0(prefix, "_deseq_list"), envir = .GlobalEnv, inherits = TRUE)
        list2env(deseq_list, envir = .GlobalEnv)
        rm(deseq_list)
        
        baseline_idx <- map2(baseline_filter, names(baseline_filter),
                             ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
          purrr::reduce(inner_join)
        
        baseline_vec <- colMeans(mod_mat[baseline_idx, ])
        #identify the conditions you are testing against the baseline, and determine interactions = contrasts
        condition_filter <- list(`Vent` = c("anemone", "marker113", "marker33"),
                                 `Temperature` = tempid,
                                 `Fraction` = c("13H"),
                                 `year` = c(2013, 2014))
        
        condition_idx_names <- map(condition_filter,
                                   ~unlist(., recursive = FALSE)) %>%
          purrr::reduce(interaction, sep = "_") %>%
          levels %>%
          as.character
        
        for(i in vents){
          namevar <- i
          for(j in c(2013, 2014)){
            tempvar <- j
            if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$year == tempvar,]))){
              temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$year == tempvar, ]
            } else {
              temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$year == tempvar, ])
            }
            if(any(is.na(temp_vec))){
              contrast_id <- contrast_id
            } else {
              assign(paste0(namevar, "_", tempid, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
              if(exists("contrast_id", inherits = FALSE)){
                contrast_id <- bind_rows(contrast_id, tibble::as_tibble(paste0(namevar, "_", tempid, "_13H_", tempvar)))  
              } else {
                contrast_id <- tibble::as_tibble(paste0(namevar, "_", tempid, "_13H_", tempvar))
              }
            }
            rm(temp_vec)
            rm(tempvar)
          }
          rm(namevar)
        }
        
        if(exists("contrast_id", inherits = FALSE)){
          label_id <- contrast_id %>%
            unlist %>%
            simplify %>%
            as.character %>%
            paste0(., " - axial_", tempid, "_13H")
          contrast_id <- contrast_id %>%
            as.data.frame() %>%
            dplyr::mutate(baseline = "baseline_vec") %>%
            dplyr::mutate(across(everything(), as.character))
        } else {
          contrast_id <- condition_idx_names %>%
            as.data.frame() %>%
            dplyr::mutate(baseline = "baseline_vec") %>%
            dplyr::mutate(across(everything(), as.character))
          label_id <- paste0(condition_idx_names, " - axial_", tempid, "_13H")
        }
        
        y <- nrow(contrast_id)
        for(i in seq(y)){
          Sys.sleep(1/100)
          pair1 <- get0(contrast_id[i,1], inherits = TRUE)
          pair2 <- get0(contrast_id[i,2], inherits = TRUE)
          res1 <- results(dds,
                          contrast = pair1 - pair2,
                          parallel = TRUE,
                          BPPARAM = bpparam_multi,
                          alpha = set.alpha) %>%
            as.data.frame(.) %>%
            rownames_to_column("gene_symbol") %>%
            left_join(., lfcShrink(dds,
                                   contrast = pair1 - pair2,
                                   parallel = TRUE,
                                   BPPARAM = bpparam_multi,
                                   type = "ashr",
                                   quiet = TRUE) %>%
                        as.data.frame(.) %>%
                        rownames_to_column("gene_symbol") %>%
                        dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                        dplyr::select(gene_symbol, log2FoldChange_MMSE),
                      by = "gene_symbol") %>%
            mutate(contrast = label_id[i]) %>%
            relocate(contrast, .before = "gene_symbol")
          
          if(exists("res")){
            res <- bind_rows(res, res1)
            rm(res1)}
          if(!exists("res")){
            res <- res1
            rm(res1)}
        }
        #anemone: only for 80˚C are there enough 13H samples to compare 2013 and 2014.
        #marker113: between 2013 and 2014, 55˚C and 80˚C can be compared
        #marker33: all temps can be be compared btw 2013 and 2014
        #filter out the unwanted comparisons
        
        if(tempid == 30){
          res <- res %>%
            # coa_axial_30_13H_vent_year_res_list <- coa_axial_30_13H_vent_year_res_list %>%
            dplyr::filter(grepl("marker33", contrast)) %>%
            droplevels
          
        } else if(tempid == 55){
          res <- res %>%
            # coa_axial_55_13H_vent_year_res_list <- coa_axial_55_13H_vent_year_res_list %>%
            dplyr::filter(grepl("marker33|marker113", contrast)) %>%
            droplevels
        }
        assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
        readr::write_rds(res, file = paste0(projectpath, "/data/", prefix, "_res_list", ".rds"))
        to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
        rm(prefix)
        
      }
      #temp specific comparisons between 13H vents
      else { #temp specific comparisons between 13H vents
        tempid <- stringr::str_extract(res_name, "(30|55|80)") %>%
          as.numeric(.)
        
        condition_filter <- append(condition_filter,
                                   list(`Temperature` = tempid))
        
        prefix <- paste0("coa_axial_", tempid, "_13H")
        #identify what the baseline should consist of
        deseq_list <- get0(paste0("coa_axial_", tempid, "_13H_deseq_list"), envir = .GlobalEnv, inherits = TRUE)
        list2env(deseq_list, envir = .GlobalEnv)
        rm(deseq_list)
        
        baseline_idx <- map2(baseline_filter, names(baseline_filter),
                             ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
          purrr::reduce(inner_join)
        
        baseline_vec <- colMeans(mod_mat[baseline_idx, ])
        
        condition_idx_names <- map(condition_filter,
                                   ~unlist(., recursive = FALSE)) %>%
          purrr::reduce(interaction, sep = "_") %>%
          levels %>%
          as.character
        
        #brute force: list the model vectors
        for(i in vents){
          namevar <- i
          for(j in c(30, 55, 80)){
            tempvar <- j
            if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar,]))){
              temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ]
            } else {
              temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ])
            }
            assign(paste0(namevar, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
            rm(temp_vec)
            rm(tempvar)
          }
          rm(namevar)
        }
        
        contrast_id <- condition_idx_names %>%
          as.data.frame() %>%
          dplyr::mutate(baseline = "baseline_vec") %>%
          dplyr::mutate(across(everything(), as.character))
        
        # label_id <- paste0(condition_idx_names, " - axial_temp_", "13H")
        label_id <- paste0(condition_idx_names, " - axial_", tempid, "_13H")
        
        y <- nrow(contrast_id)
        for(i in seq(y)){
          Sys.sleep(1/100)
          pair1 <- get0(contrast_id[i,1], inherits = TRUE)
          pair2 <- get0(contrast_id[i,2], inherits = TRUE)
          res1 <- results(dds,
                          contrast = pair1 - pair2,
                          parallel = TRUE,
                          BPPARAM = bpparam_multi,
                          alpha = set.alpha) %>%
            as.data.frame(.) %>%
            rownames_to_column("gene_symbol") %>%
            left_join(., lfcShrink(dds,
                                   contrast = pair1 - pair2,
                                   parallel = TRUE,
                                   BPPARAM = bpparam_multi,
                                   type = "ashr",
                                   quiet = TRUE) %>%
                        as.data.frame(.) %>%
                        rownames_to_column("gene_symbol") %>%
                        dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                        dplyr::select(gene_symbol, log2FoldChange_MMSE),
                      by = "gene_symbol") %>%
            mutate(contrast = label_id[i]) %>%
            relocate(contrast, .before = "gene_symbol")
          
          if(exists("res")){
            res <- bind_rows(res, res1)
            rm(res1)}
          if(!exists("res")){
            res <- res1
            rm(res1)}
        }
        assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
        readr::write_rds(res, file = paste0(projectpath, "/data/", prefix, "_res_list", ".rds"))
        to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
        rm(prefix)
      }
    }
  }
  
  #deprecated:
  {
    # if(any(grepl("coa_axial_temp_12L", to_process_results))){
    #   prefix <- "coa_axial_temp_12L"
    #   cli::cli_alert_info("Extracting DESeq results: {prefix}")
    # 
    #   baseline_filter <- list(`Fraction` = c("12L"))
    #   {
    #     res <- NULL
    #     list2env(coa_axial_temp_12L_deseq_list, envir = .GlobalEnv)
    #     baseline_idx <- map2(baseline_filter, names(baseline_filter),
    #                          ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
    #       purrr::reduce(inner_join)
    # 
    #     colSums(mod_mat[baseline_idx, ]) == colSums(mod_mat[dds$Fraction == "12L", ])
    # 
    #     baseline_vec <- colMeans(mod_mat[baseline_idx, ])
    #     #identify the conditions you are testing against the baseline, and determine interactions = contrasts
    #     condition_filter <- list(`Vent` = c("anemone", "marker113", "marker33"),
    #                              `Fraction` = c("13H"),
    #                              `Temperature` = c(30, 55, 80))
    # 
    #     condition_idx_names <- map(condition_filter,
    #                                ~unlist(., recursive = FALSE)) %>%
    #       purrr::reduce(interaction, sep = "_") %>%
    #       levels %>%
    #       as.character
    # 
    #     #brute force: list the model vectors
    #     for(i in vents){
    #       namevar <- i
    #       for(j in c(30, 55, 80)){
    #         tempvar <- j
    #         if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar,]))){
    #           temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ]
    #         } else {
    #           temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ])
    #         }
    #         assign(paste0(namevar, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
    #         rm(temp_vec)
    #         rm(tempvar)
    #       }
    #       rm(namevar)
    #     }
    # 
    #     contrast_id <- condition_idx_names %>%
    #       as.data.frame() %>%
    #       dplyr::mutate(baseline = "baseline_vec") %>%
    #       dplyr::mutate(across(everything(), as.character))
    # 
    #     label_id <- paste0(condition_idx_names, " - axial_temp_", "12L")
    # 
    #     y <- nrow(contrast_id)
    #     for(i in seq(y)){
    #       Sys.sleep(1/100)
    #       pair1 <- get0(contrast_id[i,1], inherits = TRUE)
    #       pair2 <- get0(contrast_id[i,2], inherits = TRUE)
    #       res1 <- results(dds,
    #                       contrast = pair1 - pair2,
    #                       parallel = TRUE,
    #                       BPPARAM = bpparam_multi,
    #                       alpha = set.alpha) %>%
    #         as.data.frame() %>%
    #         rownames_to_column("gene_symbol") %>%
    #         left_join(., lfcShrink(dds,
    #                                contrast = pair1 - pair2,
    #                                parallel = TRUE,
    #                                BPPARAM = bpparam_multi,
    #                                type = "ashr",
    #                                quiet = TRUE) %>%
    #                     as.data.frame() %>%
    #                     rownames_to_column("gene_symbol") %>%
    #                     dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
    #                     dplyr::select(gene_symbol, log2FoldChange_MMSE),
    #                   by = "gene_symbol") %>%
    #         mutate(contrast = label_id[i]) %>%
    #         relocate(contrast, .before = "gene_symbol")
    #       if(exists("res")){
    #         res <- bind_rows(res, res1)
    #         rm(res1)}
    #       if(!exists("res")){
    #         res <- res1
    #         rm(res1)}
    #     }
    # 
    #     res_12L <- counts(dds[,dds$Fraction == "12L"], normalized = TRUE) %>%
    #       as.data.frame() %>%
    #       tibble::rownames_to_column("gene_symbol") %>%
    #       dplyr::rowwise() %>%
    #       dplyr::mutate(baseMean = mean(c_across(!gene_symbol))) %>%
    #       dplyr::select(gene_symbol, baseMean)
    # 
    #     res <- res %>%
    #       dplyr::rename(priorbaseMean = "baseMean") %>%
    #       left_join(., res_12L,
    #                 relationship = "many-to-many") %>%
    #       dplyr::relocate(baseMean, .after = "gene_symbol")
    #     # assign(paste0("coa_", "axial_temp_12L", "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
    #     assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
    #   }
    #   readr::write_rds(coa_axial_temp_12L_res_list, file = paste0(projectpath, "/data/", "coa_axial_temp_12L_res_list", ".rds"))
    #   to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
    #   rm(prefix)
    # }
  }
  
  #one by one:
  {
    #one by one:
    # if(any(grepl("coa_axial_temp_13H", to_process_results))){
    #   prefix <- "coa_axial_temp_13H"
    #   cli::cli_alert_info("Extracting DESeq results: {prefix}")
    #   
    #   baseline_filter <- list(`Fraction` = c("13H"))
    #   {
    #     res <- NULL
    #     #identify what the baseline should consist of
    #     list2env(coa_axial_temp_13H_deseq_list, envir = .GlobalEnv)
    #     baseline_idx <- map2(baseline_filter, names(baseline_filter),
    #                          ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
    #       purrr::reduce(inner_join)
    #     
    #     baseline_vec <- colMeans(mod_mat[baseline_idx, ])
    #     #identify the conditions you are testing against the baseline, and determine interactions = contrasts
    #     condition_filter <- list(`Vent` = c("anemone", "marker113", "marker33"),
    #                              `Fraction` = c("13H"),
    #                              `Temperature` = c(30, 55, 80))
    #     
    #     condition_idx_names <- map(condition_filter,
    #                                ~unlist(., recursive = FALSE)) %>%
    #       purrr::reduce(interaction, sep = "_") %>%
    #       levels %>%
    #       as.character
    #     
    #     #brute force: list the model vectors
    #     for(i in vents){
    #       namevar <- i
    #       for(j in c(30, 55, 80)){
    #         tempvar <- j
    #         if(is.null(dim(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar,]))){
    #           temp_vec <- mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ]
    #         } else {
    #           temp_vec <- colMeans(mod_mat[dds$Vent == namevar & dds$Fraction == "13H" & dds$Temperature == tempvar, ])
    #         }
    #         assign(paste0(namevar, "_13H_", tempvar), temp_vec, envir = .GlobalEnv)
    #         rm(temp_vec)
    #         rm(tempvar)
    #       }
    #       rm(namevar)
    #     }
    #     
    #     contrast_id <- condition_idx_names %>%
    #       as.data.frame() %>%
    #       dplyr::mutate(baseline = "baseline_vec") %>%
    #       dplyr::mutate(across(everything(), as.character))
    #     
    #     label_id <- paste0(condition_idx_names, " - axial_temp_", "13H")
    #     
    #     y <- nrow(contrast_id)
    #     for(i in seq(y)){
    #       Sys.sleep(1/100)
    #       pair1 <- get0(contrast_id[i,1], inherits = TRUE)
    #       pair2 <- get0(contrast_id[i,2], inherits = TRUE)
    #       res1 <- results(dds,
    #                       contrast = pair1 - pair2,
    #                       parallel = TRUE,
    #                       BPPARAM = bpparam_multi,
    #                       alpha = set.alpha) %>%
    #         as.data.frame() %>%
    #         rownames_to_column("gene_symbol") %>%
    #         left_join(., lfcShrink(dds,
    #                                contrast = pair1 - pair2,
    #                                parallel = TRUE,
    #                                BPPARAM = bpparam_multi,
    #                                type = "ashr",
    #                                quiet = TRUE) %>%
    #                     as.data.frame() %>%
    #                     rownames_to_column("gene_symbol") %>%
    #                     dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
    #                     dplyr::select(gene_symbol, log2FoldChange_MMSE),
    #                   by = "gene_symbol") %>%
    #         mutate(contrast = label_id[i]) %>%
    #         relocate(contrast, .before = "gene_symbol")
    #       if(exists("res")){
    #         res <- bind_rows(res, res1)
    #         rm(res1)}
    #       if(!exists("res")){
    #         res <- res1
    #         rm(res1)}
    #     }
    #     assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
    #   }
    #   readr::write_rds(coa_axial_temp_13H_res_list, file = paste0(projectpath, "/data/", "coa_axial_temp_13H_res_list", ".rds"))
    #   to_process_results <- grep(prefix, to_process_results, invert = TRUE, value = TRUE)
    #   rm(prefix)
    # }
  }
  
}


if(!file.exists(paste0(projectpath, "/data/", "coa_axial_eachtemp_bothfrac_res_list", ".rds"))){
  if(length(apropos(".*_\\d\\d_bothfrac_res_list", mode = "any") == 3)){
    coa_axial_eachtemp_bothfrac_res_list <- bind_rows(coa_axial_30_bothfrac_res_list,
                                                      coa_axial_55_bothfrac_res_list,
                                                      coa_axial_80_bothfrac_res_list)
    readr::write_rds(coa_axial_eachtemp_bothfrac_res_list, file = paste0(projectpath, "/data/", "coa_axial_eachtemp_bothfrac_res_list", ".rds"))
  } else {
    cli::cli_alert_warning("Please parse the temperature-specific comparisons for 13H+12L.")
  }
}

if(!file.exists(paste0(projectpath, "/data/", "coa_axial_eachtemp_13H_res_list", ".rds"))){
  if(length(apropos(".*_\\d\\d_13H_res_list", mode = "any")) == 3){
    coa_axial_eachtemp_13H_res_list <- bind_rows(coa_axial_30_13H_res_list,
                                                 coa_axial_55_13H_res_list,
                                                 coa_axial_80_13H_res_list)
    readr::write_rds(coa_axial_eachtemp_13H_res_list, file = paste0(projectpath, "/data/", "coa_axial_eachtemp_13H_res_list", ".rds"))
  } else {
    cli::cli_alert_warning("Please parse the temperature-specific 13H results and make a new data object.")
  }
} 

if(!file.exists(paste0(projectpath, "/data/", "coa_axial_eachtemp_13H_vent_year_res_list", ".rds"))){
  if(length(apropos(".*_\\d\\d_13H_vent_year_res_list", mode = "any")) == 3){
    coa_axial_eachtemp_13H_vent_year_res_list <- bind_rows(coa_axial_30_13H_vent_year_res_list,
                                                           coa_axial_55_13H_vent_year_res_list,
                                                           coa_axial_80_13H_vent_year_res_list)
    readr::write_rds(coa_axial_eachtemp_13H_vent_year_res_list, file = paste0(projectpath, "/data/", "coa_axial_eachtemp_13H_vent_year_res_list", ".rds"))
  } else {
    cli::cli_alert_warning("Please parse the temperature-specific, vent-specific, temporal 13H results and make a new data object.")
  }
} 

