# 00_resource_allocation.R

#prepare resources for R on your computer to run this project

# Load packages -----------------------------------------------------------

pkgs <- c("cli", "data.table", "rlang", "stringr", "tibble", "tidyr", "viridis", "viridisLite", "parallelly", "future", "furrr", "progressr", "dtplyr", "RVenn", "remotes", "UpSetR", "formattable", "pals", "dendextend", "gplots", "pheatmap", "RColorBrewer", "ggtext", "ggvegan", "ggsankey", "ggVennDiagram", "ggplotify", "patchwork", "multidplyr", 
          "tidyverse")
bioc_pkgs <- c("BiocParallel", "GenomicRanges", "GenomicFeatures", "phyloseq", "DESeq2")
# #for troubleshooting:
# pkgs <- c("cli", "data.table")
# bioc_pkgs <- c("BiocParallel", "DESeq2")

install_pkgs <- NULL

for(temp_lib in pkgs){
  # if (!require(temp_lib, quietly = TRUE, character.only = TRUE)) { install.packages(temp_lib) }
  if (!(temp_lib %in% installed.packages(fields = NULL))){
    install_pkgs <- append(install_pkgs, temp_lib)
  }
}
if(!is.null(install_pkgs)){
  # install_pkgs <- str2lang(deparse1(as.vector(install_pkgs), collapse = ","))
  # install.packages(pkgs = eval(install_pkgs), 
  install.packages(pkgs = eval(str2lang(deparse1(as.vector(install_pkgs), collapse = ","))),
                   dependencies = c("Depends", "Imports"),
                   clean = TRUE)
}
# load_pkgs <- pkgs[1:5]
load_pkgs <- pkgs
for(temp_lib in load_pkgs){
  if(!(temp_lib %in% library())){
    require(temp_lib, character.only = TRUE, quietly = TRUE)  
  }
  load_pkgs <- load_pkgs[grep(temp_lib, load_pkgs, value = FALSE, invert = TRUE)]
}


if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library(BiocManager)
for(temp_lib in bioc_pkgs){
  if (!require(temp_lib, quietly = TRUE, character.only = TRUE)) { BiocManager::install(pkgs = temp_lib) }
  require(temp_lib, character.only = TRUE)
  # bioc_pkgs <- bioc_pkgs[grep(temp_lib, bioc_pkgs, value = FALSE, invert = TRUE)]
}


#this is a loop that checks for each package:
# if (!require("furrr", quietly = TRUE)){
#   install.packages("furrr")
#   install.packages("future")
# }
# library(furrr)


# Plan for resource allocation --------------------------------------------
if(grepl("arch64", Sys.getenv("R_PLATFORM"))){
  print("Detected Mac, using parallelly...")
  # nthreads <- parallelly::availableCores(omit = 1) - 1
  nthreads <- availableCores(omit = 1) - 1
  future::plan(multisession, workers = 1)
  # future::plan(multisession, workers = nthreads)
  options(future.globals.maxSize = 10e9)
} else {
  if(grepl("x86_64", Sys.getenv("R_PLATFORM"))){
    print("Detected Windows")
    nthreads <- parallelly::availableCores(omit = 1) - 1
    future::plan(sequential)
    options(future.globals.maxSize = 10e9)
  } else {
    print("Using data.table")
    nthreads <- data.table::getDTthreads()
    future::plan(sequential)
    options(future.globals.maxSize = 10e9)
  }
}


if(nthreads > 4){
  bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                workers = nthreads - 2,
                                                stop.on.error = TRUE,
                                                RNGseed = 48105,
                                                progressbar = TRUE)
  register(bpparam_multi)
} else {
  bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                workers = nthreads,
                                                stop.on.error = TRUE,
                                                RNGseed = 48105,
                                                progressbar = TRUE)
  register(bpparam_multi)
}

# User input to proceed with code chunk -----------------------------------


f_run_chunk <- function() {
  #run 'try(f_run_chunk())' to set variable 'execchunk' to TRUE or FALSE
  #after running this function, use 'if(execchunk){}' to run a code chunk within the curly braces
  ANSWER <- readline("Do you want to run this chunk? ")
  if (!grepl("T|y", substr("TRUE", 1, 1), ignore.case = TRUE)){
  # if (substr(ANSWER, 1, 1) != "y"){
    assign("execchunk", FALSE, envir = .GlobalEnv, inherits = TRUE)
    stop("Skipping the code chunk.")
  } else {
    #do the code chunk
    cli::cli_alert_info("Running the code chunk... \n") 
    assign("execchunk", TRUE, envir = .GlobalEnv, inherits = TRUE)
  }
}


f_projectpath <- function() {
  ANSWER <- readline("Type in your full project path. ")
  if (substr(ANSWER, 1, 1) == "/"){
    if(file.exists(ANSWER)){
      cli::cli_alert_info("Setting path to {ANSWER} \n") 
      assign("projectpath", ANSWER, inherits = TRUE, envir = .GlobalEnv)  
    } else
      stop("Your provided path does not seem to exist.")
  } else {
    if(ANSWER == "getwd()"){
      assign("projectpath", getwd(), inherits = TRUE, envir = .GlobalEnv)  
    } else {
      stop("Please provide a path.\n")   
    }
  }
}

f_remotepath <- function() {
  ANSWER <- readline("Type in your full remote path. ")
  if (grepl(":", ANSWER, ignore.case = TRUE)){
  # if (substr(ANSWER, 1, 1) == "/"){
      cli::cli_alert_info("Setting remote path to {ANSWER} \n") 
      assign("remotepath", ANSWER, inherits = TRUE, envir = .GlobalEnv)  
  } else {
    stop("Please provide a remote path.\n") 
  }
}

f_serverpath <- function() {
  ANSWER <- readline("Type in your full server path. ")
  if (substr(ANSWER, 1, 1) == "/"){
    cli::cli_alert_info("Setting server path to {ANSWER} \n") 
      assign("serverpath", ANSWER, inherits = TRUE, envir = .GlobalEnv)  
  } else {
    stop("Please provide a path.\n") 
  }
}




