# 09_active_mag_expression.R
# this script follows:
# 00_process_coa_mtgs.R
# 00_process_coa_mtts.R
# 01_metat_coa_deseq.R
# 01_metat_coa_tpm.R
# 02_metat_coa_sde_venn.R
# 03_metat_coa_sde_plots.R
# 04_metat_collate_bins.R
# 05_metat_mag_tx.R
# 06_metat_gene_taxonomy.R
# 07_mag_bin_abundances.R
# 08_mag_bin_expression.R

#this script plots the expression of highly-active, high-quality MAG bins



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
  library(viridis)
  library(BiocManager)
  library(BiocParallel)
  # library(DESeq2)
  library(cli)
  library(furrr)
  library(progressr)
  library(formattable)
  library(pals)
  library(pheatmap)
  library(gplots)
  library(dendextend)
  # library(ggpmisc)
  library(ggplotify)
  
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

# Import sample metadata --------------------------------------------------

vents <- c("marker113", "anemone", "marker33")
set.alpha <- 0.05
vent_lookup <- c("marker113" = "Marker 113",
                 "anemone" = "Anemone",
                 "marker33" = "Marker 33")
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
  dplyr::mutate(Vent = factor(Vent, levels = names(vent_lookup))) %>%
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
  dplyr::mutate(sample_name = dplyr::case_when(grepl("12L10|13H7", sample_name) ~ stringr::str_remove_all(sample_name, "(10|7)"),
                                               .default = sample_name)) %>%
  droplevels %>%
  dplyr::mutate(across(everything(), ~as.factor(.x)))


metadata <- sample_metadata %>%
  dplyr::filter(grepl("(13H|13L|12H|12L)", Fraction)) %>%
  # dplyr::filter(grepl("13H", Fraction)) %>%
  droplevels %>%
  dplyr::arrange(by = Vent, Temperature, sample_source, Timepoint, Fraction, year) %>%
  dplyr::mutate(across(c(Vent,Fraction, sample_source, Timepoint), ~as.factor(.x))) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = unique(.[["sample_name"]]))) %>%
  tibble::rowid_to_column(., var = "order")


# Import binning and expression data --------------------------------------
if(file.exists(paste0(projectpath, "/data/", "binned_contigs_tbl", ".rds"))){
  binned_contigs_tbl <- readr::read_rds(paste0(getwd(), "/", "binned_contigs_tbl", ".rds")) %>%
    dplyr::mutate(contig = factor(contig, levels = unique(.[["contig"]]))) %>%
    dplyr::mutate(bin = factor(bin, levels = unique(.[["bin"]])))
} else {
  cli::cli_alert_info("Please parse the list of binned contigs.")
}
if(!exists("coa_contig_gene_coords_tbl", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "coa_contig_gene_coords_tbl", ".rds"))){
    coa_contig_gene_coords_tbl <- readr::read_rds(paste0(projectpath, "/data/", "coa_contig_gene_coords_tbl", ".rds"))
  } else {
    cli::cli_alert_warning("you need to parse the gene coordinates dataset for contigs in the coassembly.")
  }  
}

# Import high-quality binning metadata ------------------------------------

if(!exists("coa_bins_highq_taxonomy_gtdbtk_df", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "coa_bins_highq_taxonomy_gtdbtk_df", ".tsv"))){
    coa_bins_highq_taxonomy_gtdbtk_df <- readr::read_delim(paste0(projectpath, "/data/", "coa_bins_highq_taxonomy_gtdbtk_df", ".tsv"),
                                                           delim = "\t", 
                                                           show_col_types = FALSE,
                                                           col_names = TRUE, 
                                                           num_threads = nthreads)
  } else {
    cli::cli_alert_warning("You need to finalize the taxonomy of the high-quality MAG bins.")
  }
}

quality_bin_label_lookup <- c(coa_bins_highq_taxonomy_gtdbtk_df[["relabeled"]]) %>%
  setNames(., coa_bins_highq_taxonomy_gtdbtk_df[["bin"]])


if(!exists("coa_highq_newick_pruned_relabeled", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "coa_bins_highq_phylo", ".rds"))){
    coa_bins_highq_phylo <- readr::read_rds((paste0(projectpath, "/data/", "coa_bins_highq_phylo", ".rds")))
    coa_highq_newick_pruned_relabeled <- phyloseq::phy_tree(coa_bins_highq_phylo)
  } else {
    cli::cli_alert_info("You need to finalize the cladogram of MAG bin phylogeny.")
  }
}


if(!exists("mag_bins_highq_annotation_taxa_colors", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "mag_bins_highq_annotation_taxa_colors", ".rds"))){
    mag_bins_highq_annotation_taxa_colors <- readr::read_rds(paste0(projectpath, "/data/", "mag_bins_highq_annotation_taxa_colors", ".rds"))
  }  
}
if(file.exists(paste0(projectpath, "/data/", "mag_bins_highq_dend_rows_taxonomy", ".rds"))){
  dend_rows_taxonomy <- readr::read_rds(paste0(projectpath, "/data/", "mag_bins_highq_dend_rows_taxonomy", ".rds"))
} else {
  cli::cli_alert_warning("You need to prepare the dendrogram of the high-quality MAG bins.")
}

annotation_taxa_colors_bin <- coa_bins_highq_taxonomy_gtdbtk_df %>%
  dplyr::select(relabeled, taxonomy.gtdbtk, bin) %>%
  dplyr::left_join(., (purrr::flatten(mag_bins_highq_annotation_taxa_colors) %>%
                         setNames(., stringr::str_split_i(names(.), ";", -1)) %>%
                         tibble::enframe(., name = "taxonomy", value = "color")), 
                   by = join_by("taxonomy.gtdbtk" == "taxonomy"), multiple = "all", relationship = "many-to-many") %>%
  dplyr::select(relabeled, color) %>%
  tibble::deframe(.) %>%
  simplify

annotation_taxa_colors_taxonomy <- coa_bins_highq_taxonomy_gtdbtk_df %>%
  dplyr::select(relabeled, taxonomy.gtdbtk, bin) %>%
  dplyr::left_join(., (purrr::flatten(mag_bins_highq_annotation_taxa_colors) %>%
                         setNames(., stringr::str_split_i(names(.), ";", -1)) %>%
                         tibble::enframe(., name = "taxonomy", value = "color")), 
                   by = join_by("taxonomy.gtdbtk" == "taxonomy"), multiple = "all", relationship = "many-to-many") %>%
  dplyr::distinct(taxonomy.gtdbtk, color, .keep_all = FALSE) %>%
  tibble::deframe(.) %>%
  simplify

# color_names <- c(annotation_taxa_colors,
color_names <- c(
  list(Vent = c("navy", "pink", "yellow", "lavender", "cornflowerblue") %>%
         setNames(., c("marker113", "anemone", "marker33", "plume", "CTD")),
       bin = annotation_taxa_colors_bin %>% setNames(., names(annotation_taxa_colors_bin)),
       taxonomy.gtdbtk = annotation_taxa_colors_taxonomy %>% setNames(., names(annotation_taxa_colors_taxonomy)),
       Temperature = c("skyblue", "forestgreen", "purple") %>%
         setNames(c("30", "55", "80")),
       year = c("brown", "khaki", "azure") %>%
         setNames(c("2013", "2014", "2015"))))

dend_colorbar_row <- dend_rows_taxonomy %>%
  dplyr::arrange(match(relabeled, labels(coa_highq_newick_pruned_relabeled))) %>%
  dplyr::mutate(Genus = dplyr::case_when(!is.na(Genus) ~ paste(Family, Genus, sep = ";"),
                                         .default = Genus)) %>%
  dplyr::mutate(Family = dplyr::case_when(!is.na(Family) ~ paste(Order, Family, sep = ";"),
                                          .default = Family)) %>%
  dplyr::mutate(Order = dplyr::case_when(!is.na(Order) ~ paste(Class, Order, sep = ";"),
                                         .default = Order)) %>%
  dplyr::mutate(Class = dplyr::case_when(!is.na(Class) ~ paste(Phylum, Class, sep = ";"),
                                         .default = Class)) %>%
  dplyr::mutate(Phylum = dplyr::case_when(!is.na(Phylum) ~ paste(Domain, Phylum, sep = ";"),
                                          .default = Phylum)) %>%
  dplyr::mutate(across(everything(), ~stringr::str_remove_all(.x, ";NA"))) %>%
  droplevels %>%
  dplyr::mutate(taxonomy = grouping) %>%
  dplyr::select(-c(grouping, starts_with("taxonomy"), taxonomy_string))  %>%
  dplyr::relocate(Genus, Family, Order, Class, Phylum, Domain, .after = "bin") %>%
  dplyr::select(-bin) %>%
  tibble::column_to_rownames(var = "relabeled")

mags_mtg_colorbar <- data.frame(relabeled = rownames(dend_colorbar_row)) %>%
  dplyr::left_join(coa_bins_highq_taxonomy_gtdbtk_df, by = join_by("relabeled")) %>%
  dplyr::select(relabeled, taxonomy.gtdbtk) %>%
  # dplyr::mutate(bin = relabeled) %>%
  # dplyr::relocate(bin) %>%
  dplyr::mutate(across(everything(), factor)) %>%
  tibble::column_to_rownames(var = "relabeled") %>%
  droplevels
mtg_contig_colorbar <- binned_contigs_tbl %>%
  dplyr::left_join(coa_bins_highq_taxonomy_gtdbtk_df, by = join_by("bin")) %>%
  droplevels %>%
  dplyr::select(relabeled, contig) %>%
  dplyr::arrange(rev(relabeled)) %>%
  dplyr::mutate(relabeled = factor(relabeled, levels = unique(.[["relabeled"]]))) %>%
  dplyr::rename(bin = "relabeled") %>%
  tibble::column_to_rownames(var = "contig")

dend_colorbar_col <- metadata %>%
  dplyr::select(sample_name, Vent, Temperature, year) %>%
  droplevels %>%
  dplyr::arrange(Temperature, Vent, year) %>%
  tibble::column_to_rownames(var = "sample_name")

x_histogram_labels <- metadata %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(sample_name, Vent, Temperature, .keep_all = FALSE) %>%
  dplyr::mutate(labeled = paste0(Vent, "_", Temperature)) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::distinct(sample_source, labeled, .keep_all = FALSE) %>%
  # dplyr::distinct(sample_name, labeled, .keep_all = FALSE) %>%
  droplevels %>%
  tibble::deframe(.)

# Import TPM tables for contigs across samples ----------------------------


#for just 13H samples:
{
  # if(!exists("coa_genes_tpm_df", envir = .GlobalEnv)){
  #   if(file.exists(paste0(projectpath, "/data/", "coa_genes_tpm_table", ".tsv.gz"))){
  #     cli::cli_alert_info("reading in genes TPM table for 13H samples")
  #     coa_genes_tpm_df <- readr::read_delim(paste0(projectpath, "/data/", "coa_genes_tpm_table", ".tsv.gz"), 
  #                                           col_names = TRUE,
  #                                           show_col_types = FALSE,
  #                                           quote = "",
  #                                           comment = "#",
  #                                           delim = "\t",
  #                                           num_threads = nthreads) %>%
  #       # coa_genes_tpm_df <- coa_genes_tpm_table %>%
  #       tidyr::pivot_longer(., cols = !c(gene_id),
  #                           names_to = "sample_name",
  #                           values_to = "TPM") %>%
  #       dplyr::left_join(., coa_contig_gene_coords_tbl %>%
  #                          dplyr::select(target_id, contig) %>%
  #                          dplyr::rename(gene_id = "target_id"), by = join_by(gene_id)) %>%
  #       dplyr::left_join(., binned_contigs_tbl, by = join_by(contig))
  #   } else {
  #     cli::cli_alert_warning("You need to process the genes TPM table for 13H samples.")
  #   }
  # }
  # 
  # if(!exists("coa_binned_genes_tpm_df", envir = .GlobalEnv)){
  #   if(file.exists(paste0(projectpath, "/data/", "coa_binned_genes_tpm_table", ".tsv"))){
  #     cli::cli_alert_info("reading in binned genes TPM table for 13H samples")
  #     coa_binned_genes_tpm_df <- readr::read_delim(paste0(projectpath, "/data/", "coa_binned_genes_tpm_table", ".tsv"), 
  #                                                  col_names = TRUE,
  #                                                  show_col_types = FALSE,
  #                                                  quote = "",
  #                                                  comment = "#",
  #                                                  delim = "\t",
  #                                                  num_threads = nthreads) %>%
  #       # coa_binned_genes_tpm_df <- coa_binned_genes_tpm_table %>%
  #       tidyr::pivot_longer(., cols = !c(gene_id),
  #                           names_to = "sample_name",
  #                           values_to = "TPM") %>%
  #       dplyr::left_join(., coa_contig_gene_coords_tbl %>%
  #                          dplyr::select(target_id, contig) %>%
  #                          dplyr::rename(gene_id = "target_id"), by = join_by(gene_id)) %>%
  #       dplyr::left_join(., binned_contigs_tbl, by = join_by(contig))
  #     
  #   } else {
  #     cli::cli_alert_warning("You need to process the binned genes TPM table for 13H samples.")
  #   }
  # }  
}


#for all fractions:
if(!exists("coa_genes_tpm_df", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "coa_frac_genes_tpm_table", ".tsv.gz"))){
    cli::cli_alert_info("reading in genes TPM table for all fractions of samples")
    coa_genes_tpm_df <- readr::read_delim(paste0(projectpath, "/data/", "coa_frac_genes_tpm_table", ".tsv.gz"), 
                                          col_names = TRUE,
                                          show_col_types = FALSE,
                                          quote = "",
                                          comment = "#",
                                          delim = "\t",
                                          num_threads = nthreads) %>%
      tidyr::pivot_longer(., cols = !c(gene_id),
                          names_to = "sample_name",
                          values_to = "TPM") %>%
      dplyr::left_join(., coa_contig_gene_coords_tbl %>%
                         dplyr::select(target_id, contig) %>%
                         dplyr::rename(gene_id = "target_id"), by = join_by(gene_id)) %>%
      dplyr::left_join(., binned_contigs_tbl, by = join_by(contig))
  } else {
    cli::cli_alert_warning("You need to process the genes TPM table for all fractions' of samples.")
  }
}
if(!exists("coa_binned_genes_tpm_df", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/data/", "coa_frac_binned_genes_tpm_table", ".tsv.gz"))){
    cli::cli_alert_info("reading in binned genes TPM table for all fractions of samples")
    coa_binned_genes_tpm_df <- readr::read_delim(paste0(projectpath, "/data/", "coa_frac_binned_genes_tpm_table", ".tsv.gz"), 
                                                 col_names = TRUE,
                                                 show_col_types = FALSE,
                                                 quote = "",
                                                 comment = "#",
                                                 delim = "\t",
                                                 num_threads = nthreads) %>%
      tidyr::pivot_longer(., cols = !c(gene_id),
                          names_to = "sample_name",
                          values_to = "TPM") %>%
      dplyr::left_join(., coa_contig_gene_coords_tbl %>%
                         dplyr::select(target_id, contig) %>%
                         dplyr::rename(gene_id = "target_id"), by = join_by(gene_id)) %>%
      dplyr::left_join(., binned_contigs_tbl, by = join_by(contig))
    
  } else {
    cli::cli_alert_warning("You need to process the binned genes TPM table for all fractions' of samples.")
  }
}


#bin summary:
if(!exists("bin_summary_df", envir = .GlobalEnv)){
  if(!exists("bin_summary_tbl", envir = .GlobalEnv)){
    if(file.exists(paste0(projectpath, "/data/", "bin_summary_tbl", ".rds"))){
      bin_summary_tbl <- readr::read_rds(paste0(projectpath, "/data/", "bin_summary_tbl", ".rds"))
    } else {
      bin_summary_tbl <- readr::read_delim(paste0(projectpath, "Axial_coa_bins_quality_report.tsv"),
                                           col_names = TRUE,
                                           show_col_types = FALSE,
                                           quote = "",
                                           comment = "#",
                                           delim = "\t",
                                           num_threads = nthreads)
      
      if(!file.exists(paste0(getwd(), "/", "bin_summary_tbl", ".rds"))){
        readr::write_rds(bin_summary_tbl, paste0(getwd(), "/", "bin_summary_tbl", ".rds"), compress = "gz")
      }
    }
    bin_summary_df <- bin_summary_tbl %>%
      dplyr::rename(bin = "Name") %>%
      dplyr::mutate(bin = stringr::str_replace_all(bin, "MetaBat2_Bins_NEW.", "coa_bin_")) %>%
      dplyr::select(c(bin, Completeness, Contamination, Contig_N50, Genome_Size)) %>%
      dplyr::filter(bin %in% unique(coa_binned_genes_tpm_df[["bin"]])) %>%
      droplevels
  }
}

if(file.exists(paste0(projectpath, "/data/", "binned_contig_counts_vents_frac_df", ".tsv"))){
  binned_contig_counts_vents_df <- readr::read_delim(paste0(projectpath, "/data/", "binned_contig_counts_vents_frac_df", ".tsv"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads, show_col_types = FALSE) %>%
    dplyr::mutate(sample_name = dplyr::case_when(grepl("12L10|13H7", sample_name) ~ stringr::str_remove_all(sample_name, "(10|7)"),
                                                 .default = sample_name)) %>%
    droplevels
} else {
  cli::cli_alert_info("You need to summarize the transcript counts for bins in samples.")
}



# Find the most active bins -----------------------------------------------


#make a histogram showing the number of exprssed genes per bin per sample.

#now do the sus bins:
{
  # namevar_idx <- coa_bins_highq_taxonomy_gtdbtk_df %>%
  #   dplyr::filter(grepl(paste0(c("Desulfurobacteriaceae_bin_324",
  #                                "Methanothermococcus",
  #                                "Aquificales",
  #                                "Endoriftia",
  #                                "Marinobacter",
  #                                "Alloalcanivorax",
  #                                "Thermovibrio"), collapse = "|"), relabeled)) %>%
  #   dplyr::mutate(relabeled = factor(relabeled, levels = unique(.[["relabeled"]]))) %>%
  #   dplyr::arrange(rev(relabeled)) %>%
  #   droplevels %>%
  #   dplyr::distinct(bin, .keep_all = FALSE) %>%
  #   simplify %>%
  #   as.character
  
  # temp_df <- binned_contig_counts_vents_df %>%
  #   dplyr::select(bin, sample_name, num_exp_genes, num_tx) %>%
  #   dplyr::left_join(., metadata %>%
  #                      dplyr::select(sample_name, Temperature, Vent) %>%
  #                      droplevels, by = join_by(sample_name)) %>%
  #   dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
  #                                            .default = "low")) %>%
  #   dplyr::left_join(., (coa_genes_tpm_df %>%
  #                          tidyr::drop_na(bin) %>%
  #                          dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
  #                          dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  #   dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  #   dplyr::filter(bin %in% namevar_idx) %>%
  #   dplyr::mutate(num_exp_genes = num_exp_genes*100) %>% #transform so it's on the same scale as transcripts...
  #   tidyr::pivot_longer(., cols = c(num_exp_genes, num_tx, sumTPM),
  #                       names_to = "metric",
  #                       values_to = "counts") %>%
  #   dplyr::mutate(metric = factor(metric, levels = c("num_exp_genes", "num_tx", "sumTPM"))) %>%
  #   droplevels
}

# #first threshold: bin has to be 0.1% of TPM in at least one sample, to be considered
# #37 bins
# mags_active_1000_idx <- coa_genes_tpm_df %>%
#   tidyr::drop_na(bin) %>%
#   dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::filter(sumTPM > 1000) %>%
#   droplevels %>%
#   dplyr::distinct(bin, .keep_all = FALSE) %>%
#   dplyr::left_join(., coa_bins_highq_taxonomy_gtdbtk_df %>% dplyr::select(bin, relabeled), by = join_by(bin)) %>%
#   # tibble::deframe(.)
#   # simplify
#   droplevels

# #second threshold: mags have to be high quality and at least 0.1% TPM of at least one sample
# #24 bins
# mags_highq_active_1000_idx <- coa_genes_tpm_df %>%
#   tidyr::drop_na(bin) %>%
#   dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
#   dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::filter(sumTPM > 1000) %>%
#   droplevels %>%
#   dplyr::distinct(bin, .keep_all = FALSE) %>%
#   dplyr::left_join(., coa_bins_highq_taxonomy_gtdbtk_df %>% dplyr::select(bin, relabeled), by = join_by(bin)) %>%
#   # tibble::deframe(.)
#   # simplify
#   droplevels

# #third threshold: mags have to be high quality, at least 0.1% TPM of a sample across 1/3 of the samples
# #this yielded 5 bins: coa_bin_340, coa_bin_102, coa_bin_374, coa_bin_320, and coa_bin_104
# mags_highq_majority_1000_idx <- coa_genes_tpm_df %>%
#   tidyr::drop_na(bin) %>%
#   dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
#   dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::mutate(sumTPM = dplyr::case_when((sumTPM > 1000) ~ sumTPM,
#                                           .default = NA)) %>%
#   tidyr::drop_na(sumTPM) %>%
#   dplyr::summarise(ActiveSamples = length(sumTPM), .by = bin) %>%
#   dplyr::filter(ActiveSamples > 6) %>%
#   droplevels %>%
#   # dplyr::distinct(bin, .keep_all = FALSE) %>%
#   dplyr::left_join(., coa_bins_highq_taxonomy_gtdbtk_df %>% dplyr::select(bin, relabeled), by = join_by(bin)) %>%
#   # tibble::deframe(.)
#   # simplify
#   droplevels

# #fourth threshold: retain high-quality bins that had more than 5 genes expressed
# #42 bins
# #fifth threshold: high-quality bins with 5 or more genes expressed, and 0.1% TPM or more in a sample
# #19 bins
# 
# mags_highq_majority_1000_idx <- binned_contig_counts_vents_df %>%
#   dplyr::filter(num_exp_genes > 4) %>%
#   dplyr::select(bin, sample_name) %>%
#   dplyr::inner_join(., coa_genes_tpm_df, by = join_by(bin, sample_name)) %>%
#   droplevels %>%
#   dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
#   dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::filter(sumTPM > 1000) %>%
#   # dplyr::mutate(sumTPM = dplyr::case_when((sumTPM > 1000) ~ sumTPM,
#   #                                         .default = NA)) %>% tidyr::drop_na(sumTPM) %>%
#   # dplyr::summarise(ActiveSamples = length(sumTPM), .by = bin) %>% dplyr::filter(ActiveSamples > 6) %>%
#   droplevels %>%
#   # dplyr::distinct(bin, .keep_all = FALSE) %>%
#   dplyr::left_join(., coa_bins_highq_taxonomy_gtdbtk_df %>% dplyr::select(bin, relabeled), by = join_by(bin)) %>%
#   # tibble::deframe(.)
#   # simplify
#   droplevels
# 
# mags_weak_idx <- binned_contig_counts_vents_df %>%
#   dplyr::select(bin, sample_name) %>%
#   dplyr::inner_join(., coa_genes_tpm_df, by = join_by(bin, sample_name)) %>%
#   droplevels %>%
#   dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
#   dplyr::filter(!(bin %in% unique(mags_active_1000_idx[["bin"]]))) %>%
#   dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#   dplyr::ungroup(.) %>%
#   droplevels %>%
#   dplyr::left_join(., coa_bins_highq_taxonomy_gtdbtk_df %>% dplyr::select(bin, relabeled), by = join_by(bin)) %>%
#   droplevels

# #plotting:
# #first plot the number of expressed genes in each group of bins
# #then plot the TPM and TPM-only-binned for each group of bins
# 
# #let's plot these high-quality bins:
# 
# temp_df <- binned_contig_counts_vents_df %>%
#   dplyr::select(bin, sample_name, num_exp_genes, num_tx) %>%
#   dplyr::left_join(., metadata %>%
#                      dplyr::select(sample_name, Temperature, Vent) %>%
#                      droplevels, by = join_by(sample_name)) %>%
#   dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
#                                            .default = "low")) %>%
#   dplyr::left_join(., (coa_genes_tpm_df %>%
#                          tidyr::drop_na(bin) %>%
#                          dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#                          dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
#   dplyr::left_join(., (coa_binned_genes_tpm_df %>%
#                          tidyr::drop_na(bin) %>%
#                          dplyr::summarise(sumTPM_onlybinned = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
#                          dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
#   dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
#   # dplyr::mutate(sample_name = factor(sample_name, levels = unique(.[["sample_name"]])[match(rownames(dend_colorbar_col), unique(.[["sample_name"]]))])) %>%
#   # dplyr::filter(bin %in% unique(mags_highq_active_1000_idx[["bin"]])) %>%
#   dplyr::arrange(desc(sumTPM)) %>%
#   dplyr::mutate(bin = factor(bin)) %>%
#   # dplyr::mutate(num_exp_genes = num_exp_genes*100) %>% #transform so it's on the same scale as transcripts...
#   tidyr::pivot_longer(., cols = c(num_exp_genes, num_tx, sumTPM, sumTPM_onlybinned),
#                       names_to = "metric",
#                       values_to = "counts") %>%
#   dplyr::mutate(metric = factor(metric, levels = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"))) %>%
#   droplevels
# 
# x_histogram_labels <- temp_df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(sample_name, Vent, Temperature, .keep_all = FALSE) %>%
#   dplyr::mutate(labeled = paste0(Vent, "_", Temperature)) %>%
#   dplyr::distinct(sample_name, labeled, .keep_all = FALSE) %>%
#   droplevels %>%
#   tibble::deframe(.)
# 
# bin_groupings <- list((temp_df %>%
#                          dplyr::filter(bin %in% unique(mags_highq_active_1000_idx[["bin"]])) %>%
#                          dplyr::distinct(bin, .keep_all = FALSE) %>%
#                          droplevels %>%
#                          simplify %>% as.character(.)),
#                       (temp_df %>%
#                          dplyr::filter(bin %in% unique(mags_active_1000_idx[["bin"]])) %>%
#                          dplyr::filter(quality == "low") %>%
#                          dplyr::distinct(bin, .keep_all = FALSE) %>%
#                          droplevels %>%
#                          simplify %>% as.character(.)),
#                       (temp_df %>% #bin was 0.1% of TPM in at least one sample, had 5 or more expressed genes
#                          dplyr::filter(bin %in% unique(mags_highq_majority_1000_idx[["bin"]])) %>%
#                          # dplyr::filter(!(bin %in% unique(mags_active_1000_idx[["bin"]]))) %>%
#                          dplyr::distinct(bin, .keep_all = FALSE) %>%
#                          droplevels %>%
#                          simplify %>% as.character(.)),
#                       (temp_df %>% 
#                          dplyr::filter(bin %in% unique(mags_weak_idx[["bin"]])) %>%
#                          dplyr::distinct(bin, .keep_all = FALSE) %>%
#                          droplevels %>%
#                          simplify %>% as.character(.))) %>%
#   setNames(., c("Active and high quality", "Active and low quality", "Active, high quality, frequent", "Low activity, high quality bins"))

#loop it:
# # for(i in 1:length(bin_groupings)){
# for(i in c(1:3)){
#   namevar <- names(bin_groupings)[i]
#   g1 <- print(
#     ggplot(data = temp_df %>%
#              dplyr::filter(bin %in% bin_groupings[[i]]) %>%
#              dplyr::filter(metric == c("sumTPM")) %>%
#              # dplyr::mutate(counts = dplyr::case_when((metric == "sumTPM_onlybinned") ~ counts/10,
#              #                                         .default = counts)) %>%
#              droplevels,
#            aes(x = sample_name,
#                y = counts,
#                group = interaction(sample_name, metric),
#                fill = metric)) 
#     + geom_col(width = 0.90, show.legend = TRUE, position = "dodge")
#     + geom_col(color = "grey20", width = 0.90, show.legend = FALSE, position = "dodge")
#     + theme_dark()
#     + scale_y_continuous(name = "TPM", 
#                          labels = scales::number,
#                          # sec.axis = sec_axis( ~ log10(.), name = "log10(expression)"),
#                          # sec.axis = sec_axis( ~ . * 10, name = "TPM (among binned contigs)", labels = scales::number),
#                          expand = expansion(mult = c(0,0.1)))
#     + scale_fill_manual(values = viridisLite::mako(n = 4),
#                         breaks = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"),
#                         labels = c("Expressed genes", "Number of recruited transcripts", "TPM", "TPM (among binned contigs)"))
#     + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#     + facet_wrap(.~ bin, shrink = TRUE, drop = TRUE, as.table = TRUE,
#                  labeller = labeller(bin = quality_bin_label_lookup),
#                  scales = "free_y", ncol = 3, strip.position = "top")
#     + theme(panel.border = element_rect(fill = NA, colour = "black"),
#             axis.text.x = element_text(angle = 90, hjust = 0.5),
#             # axis.text.x = element_blank(),
#             panel.grid.major.x = element_blank())
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Metric",  direction = "vertical", override.aes = list(color = "black", shape = NA)),
#              color = "none")
#     # + ggtitle(paste0(namevar))
#   )
#   g2 <- print(
#     ggplot(data = temp_df %>%
#              dplyr::filter(bin %in% bin_groupings[[i]]) %>%
#              dplyr::filter(metric %in% c("sumTPM_onlybinned")) %>%
#              # dplyr::mutate(counts = dplyr::case_when((metric == "sumTPM_onlybinned") ~ counts/10,
#              #                                         .default = counts)) %>%
#              droplevels,
#            aes(x = sample_name,
#                y = counts,
#                group = interaction(sample_name, metric),
#                fill = metric)) 
#     + geom_col(width = 0.90, show.legend = TRUE, position = "dodge")
#     + geom_col(color = "grey20", width = 0.90, show.legend = FALSE, position = "dodge")
#     + theme_dark()
#     + scale_y_continuous(name = "TPM (among binned contigs)", 
#                          labels = scales::number,
#                          # sec.axis = sec_axis( ~ log10(.), name = "log10(expression)"),
#                          # sec.axis = sec_axis( ~ . * 10, name = "TPM (among binned contigs)", labels = scales::number),
#                          expand = expansion(mult = c(0,0.1)))
#     + scale_fill_manual(values = viridisLite::mako(n = 4),
#                         breaks = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"),
#                         labels = c("Expressed genes", "Number of recruited transcripts", "TPM", "TPM (among binned contigs)"))
#     + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#     + facet_wrap(.~ bin, shrink = TRUE, drop = TRUE, as.table = TRUE,
#                  labeller = labeller(bin = quality_bin_label_lookup),
#                  scales = "free_y", ncol = 3, strip.position = "top")
#     + theme(panel.border = element_rect(fill = NA, colour = "black"),
#             axis.text.x = element_text(angle = 90, hjust = 0.5),
#             # axis.text.x = element_blank(),
#             panel.grid.major.x = element_blank())
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Metric",  direction = "vertical", override.aes = list(color = "black", shape = NA)),
#              color = "none")
#     # + ggtitle(paste0(namevar))
#   )
#   
#   g3 <- print(
#     ggplot(data = temp_df %>%
#              dplyr::filter(bin %in% bin_groupings[[i]]) %>%
#              dplyr::filter(metric == "num_exp_genes") %>%
#              droplevels,
#            aes(x = sample_name,
#                y = counts,
#                group = interaction(sample_name, metric),
#                fill = metric)) 
#     + geom_col(width = 0.90, show.legend = TRUE, position = "dodge")
#     + geom_col(color = "grey20", width = 0.90, show.legend = FALSE, position = "dodge")
#     + theme_dark()
#     + scale_y_continuous(name = "Number of genes", 
#                          labels = scales::number,
#                          # sec.axis = sec_axis( ~ . * 0.01, name = "Number of genes"),
#                          expand = expansion(mult = c(0,0.1)))
#     + scale_fill_manual(values = viridisLite::mako(n = 4),
#                         breaks = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"),
#                         labels = c("Expressed genes", "Number of recruited transcripts", "TPM", "TPM (among binned contigs)"))
#     + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#     + facet_wrap(.~ bin, shrink = TRUE, drop = TRUE, as.table = TRUE,
#                  labeller = labeller(bin = quality_bin_label_lookup),
#                  scales = "free_y", ncol = 3, strip.position = "top")
#     + theme(panel.border = element_rect(fill = NA, colour = "black"),
#             axis.text.x = element_text(angle = 90, hjust = 0.5),
#             # axis.text.x = element_blank(),
#             panel.grid.major.x = element_blank())
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Metric",  direction = "vertical", override.aes = list(color = "black", shape = NA)),
#              color = "none")
#     # + ggtitle(paste0(namevar))
#   )
#   # gpatch_layout <- "
#   #   AAAAABB
#   #   AAAAABB
#   # "
#   # 
#   # gpatch <- (g2 | g3) + patchwork::plot_layout(guides = "collect", design = gpatch_layout) + patchwork::plot_annotation(title = paste0(namevar), tag_levels = "A")
#   gpatch_layout <- "
#     AAABBBCC
#     AAABBBCC
#   "
#   
#   gpatch <- (g1 | g2 | g3) + patchwork::plot_layout(guides = "collect", design = gpatch_layout) + patchwork::plot_annotation(title = paste0(namevar), tag_levels = "A")
#   assign(paste0("g_mags_hm_", i), gpatch, envir = .GlobalEnv)
# }


# for(i in 1:length(bin_groupings)){
# for(i in c(1:3)){
#   # namevar <- names(bin_groupings)[i]
#   gpatch <- get0(paste0("g_mags_hm_", i), envir = .GlobalEnv)
#   ggsave(paste0(projectpath, "/figures/", "coa_mags_by_quality_expression_", i, ".png"),
#          gpatch,
#          width = 25, height = 15, units = "in")
# }



# Selected bins to plot ---------------------------------------------------

# keep_bins_idx <- c("335", "370", "340", "41",
#   "319", "102", "104", "374",
#   "320", "324", "285", "303",
#   "221") %>%
#   paste0("coa_bin_", .)

{
 
# g1 <- print(
#   ggplot(data = temp_df %>%
#            dplyr::filter(bin %in% keep_bins_idx) %>%
#            dplyr::filter(metric == c("sumTPM")) %>%
#            droplevels,
#          aes(x = sample_name,
#              y = counts,
#              group = interaction(sample_name, metric),
#              fill = metric)) 
#   + geom_col(width = 0.90, show.legend = TRUE, position = "dodge")
#   + geom_col(color = "grey20", width = 0.90, show.legend = FALSE, position = "dodge")
#   + theme_dark()
#   + scale_y_continuous(name = "TPM", 
#                        labels = scales::number,
#                        expand = expansion(mult = c(0,0.1)))
#   + scale_fill_manual(values = viridisLite::mako(n = 4),
#                       breaks = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"),
#                       labels = c("Expressed genes", "Number of recruited transcripts", "TPM", "TPM (among binned contigs)"))
#   + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#   + facet_wrap(.~ bin, shrink = TRUE, drop = TRUE, as.table = TRUE,
#                labeller = labeller(bin = quality_bin_label_lookup),
#                scales = "free_y", ncol = 3, strip.position = "top")
#   + theme(panel.border = element_rect(fill = NA, colour = "black"),
#           axis.text.x = element_text(angle = 90, hjust = 0.5),
#           # axis.text.x = element_blank(),
#           panel.grid.major.x = element_blank())
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Metric",  direction = "vertical", override.aes = list(color = "black", shape = NA)),
#            color = "none")
# )
# 
# g2 <- print(
#   ggplot(data = temp_df %>%
#            dplyr::filter(bin %in% keep_bins_idx) %>%
#            dplyr::filter(metric %in% c("sumTPM_onlybinned")) %>%
#            droplevels,
#          aes(x = sample_name,
#              y = counts,
#              group = interaction(sample_name, metric),
#              fill = metric)) 
#   + geom_col(width = 0.90, show.legend = TRUE, position = "dodge")
#   + geom_col(color = "grey20", width = 0.90, show.legend = FALSE, position = "dodge")
#   + theme_dark()
#   + scale_y_continuous(name = "TPM (among binned contigs)", 
#                        labels = scales::number,
#                        expand = expansion(mult = c(0,0.1)))
#   + scale_fill_manual(values = viridisLite::mako(n = 4),
#                       breaks = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"),
#                       labels = c("Expressed genes", "Number of recruited transcripts", "TPM", "TPM (among binned contigs)"))
#   + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#   + facet_wrap(.~ bin, shrink = TRUE, drop = TRUE, as.table = TRUE,
#                labeller = labeller(bin = quality_bin_label_lookup),
#                scales = "free_y", ncol = 3, strip.position = "top")
#   + theme(panel.border = element_rect(fill = NA, colour = "black"),
#           axis.text.x = element_text(angle = 90, hjust = 0.5),
#           panel.grid.major.x = element_blank())
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Metric",  direction = "vertical", override.aes = list(color = "black", shape = NA)),
#            color = "none")
# )
# 
# g3 <- print(
#   ggplot(data = temp_df %>%
#            dplyr::filter(bin %in% keep_bins_idx) %>%
#            dplyr::filter(metric == "num_exp_genes") %>%
#            droplevels,
#          aes(x = sample_name,
#              y = counts,
#              group = interaction(sample_name, metric),
#              fill = metric)) 
#   + geom_col(width = 0.90, show.legend = TRUE, position = "dodge")
#   + geom_col(color = "grey20", width = 0.90, show.legend = FALSE, position = "dodge")
#   + theme_dark()
#   + scale_y_continuous(name = "Number of genes", 
#                        labels = scales::number,
#                        # sec.axis = sec_axis( ~ . * 0.01, name = "Number of genes"),
#                        expand = expansion(mult = c(0,0.1)))
#   + scale_fill_manual(values = viridisLite::mako(n = 4),
#                       breaks = c("num_exp_genes", "num_tx", "sumTPM", "sumTPM_onlybinned"),
#                       labels = c("Expressed genes", "Number of recruited transcripts", "TPM", "TPM (among binned contigs)"))
#   + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#   + facet_wrap(.~ bin, shrink = TRUE, drop = TRUE, as.table = TRUE,
#                labeller = labeller(bin = quality_bin_label_lookup),
#                scales = "free_y", ncol = 3, strip.position = "top")
#   + theme(panel.border = element_rect(fill = NA, colour = "black"),
#           axis.text.x = element_text(angle = 90, hjust = 0.5),
#           panel.grid.major.x = element_blank())
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Metric",  direction = "vertical", override.aes = list(color = "black", shape = NA)),
#            color = "none")
# )
# 
# gpatch_layout <- "
#     AAABBBCCC
#     AAABBBCCC
#   "
# 
# gpatch <- (g1 | g2 | g3) + patchwork::plot_layout(guides = "collect", design = gpatch_layout) + patchwork::plot_annotation(title = "Interesting bins?", tag_levels = "A")
# assign(paste0("g_mags_hm_selected"), gpatch, envir = .GlobalEnv)


# ggsave(paste0(projectpath, "/figures/", "coa_mags_by_quality_expression_selected", ".png"),
#        g_mags_hm_selected,
#        width = 25, height = 15, units = "in")


  }

# Other criteria to select MAGs -------------------------------------------

#drop bins that were in the CTD background and/or heterotrophs:
drop_bins_idx <- c(377, 128, 186, 247, 64, 226, 249, 5, 170, 144, 145, 68, 155, 302, 258, 31, 369, 2, 304, 345, 87) %>%
  paste0("coa_bin_", .)

# highq_mtg_grouping_idx <- list(Archaea = c(1:14),
#                                Bacteroidia = c(36:49),
#                                # Campylobacteria = c(61:78),
#                                Aquificota = c(81:104),
#                                Gammaproteobacteria = c(107:119),
#                                Midgroup = c(15:30, 32:35),
#                                Calditrichota = c(31))
#                                # Midgroup = c(15:35))


#make groupings of MAG bins based on taxonomy, phylogeny, and relative importance
highq_mtg_grouping_idx <- list(Archaea = c(1:14),
                               Bacteroidia = c(36:49),
                               Campylobacteria = c(61:78),
                               # Aquificota = c(81:104),
                               # Gammaproteobacteria = c(107:119),
                               Midgroup = c(15:30, 32:35),
                               Calditrichota = c(31))

if(length(coa_bins_highq_taxonomy_gtdbtk_df[["bin"]]) == 122){
  highq_mtg_grouping_idx <- c(highq_mtg_grouping_idx,
                              list(Aquificota = c(81:106), #with the 2 additional MAGs that were observed in the other fractions but not in 13H
                                   Gammaproteobacteria = c(109:121)))
  highq_mtg_grouping_idx <- c(highq_mtg_grouping_idx,
                              list(All_else = setdiff(c(1:122), flatten(highq_mtg_grouping_idx))))
} else {
  highq_mtg_grouping_idx <- c(highq_mtg_grouping_idx,
                              list(Aquificota = c(81:104), #with the 2 additional MAGs that were observed in the other fractions but not in 13H
                                   Gammaproteobacteria = c(107:119)))
  highq_mtg_grouping_idx <- c(highq_mtg_grouping_idx,
                              list(All_else = setdiff(c(1:120), flatten(highq_mtg_grouping_idx))))
}
highq_mtg_grouping_idx <- highq_mtg_grouping_idx[!grepl("Campylobacteria", names(highq_mtg_grouping_idx))]
campylo_mtg_grouping_idx <- dend_rows_taxonomy %>%
                                tibble::rowid_to_column(var = "position") %>%
                                dplyr::filter(if_any(everything(), ~grepl("Campylobacter", .x))) %>%
                                dplyr::mutate(grouping = dplyr::case_when(grepl("Sulfurovaceae", taxonomy_string_gtdbtk) ~ "Sulfurovaceae",
                                                                       grepl("Sulfurimonadaceae", taxonomy_string_gtdbtk) ~ "Sulfurimonadaceae",
                                                                       grepl("Hydrogenimonadaceae", taxonomy_string_gtdbtk) ~ "Hydrogenimonadaceae",
                                                                       .default = "Campylobacteria")) %>%
                                dplyr::select(grouping, position) %>%
                                droplevels %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>% 
        dplyr::select(position) %>%
        tibble::deframe(.) %>%
        unlist(.) %>%
        simplify)
  # tibble::deframe(.)
  # tibble::deframe(.) %>% tibble::enframe(., value = "position", name = "grouping")


#20240903: make plots for just these bins
# namevar <- "selected_20240903"

# Archaea:
#   Methanothermococcus bin 335
# Nanoarchaea bin 370
# All of the Epsilons/Campylobacter (All for now)
# Aquificota:
#   Thermovibrio bin 374
# Thermovibrio bin 320
# Desulfurobacteri bin 324
# Aquificatles bin 285
# Gamma: None
# Others:
#   Caldithrix bin 303

namevar <- "selected_20240917"
keep_bins_idx <- c("335", "370", 
                   "340", "41",
                   "319", "102", "104", 
                   "374", "320", "324", "285", 
                   "303") %>%
  paste0("coa_bin_", .) %>%
  append(., c(dend_rows_taxonomy %>%
                dplyr::filter(if_any(everything(), ~grepl("Campylobacter", .x))) %>%
                dplyr::select(bin) %>%
                droplevels %>%
                simplify %>%
                as.character)) %>%
  unique(.)

# highq_mtg_grouping_list <- c(highq_mtg_grouping_idx, list(All_else = setdiff(c(1:120), flatten(highq_mtg_grouping_idx))))  %>%
highq_mtg_grouping_list <- c(campylo_mtg_grouping_idx,
                             highq_mtg_grouping_idx) %>%
                            # list(All_else = setdiff(c(1:120), flatten(highq_mtg_grouping_idx))))  %>%
  map(., ~names(quality_bin_label_lookup)[.x]) %>%
  # map(., ~.x %>%
  #       setdiff(., drop_bins_idx)) %>%
  map(., ~.x %>% 
        tibble::enframe(., name = NULL, value = "bin")) %>%
  bind_rows(., .id = "grouping")


#plot heatmap of number of expressed genes per bin per sample
temp_df2 <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent, Fraction) %>%
                     droplevels, by = join_by(sample_name)) %>%
  dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
                                           .default = "low")) %>%
  dplyr::left_join(., (coa_genes_tpm_df %>%
                         tidyr::drop_na(bin) %>%
                         dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
                         dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  dplyr::left_join(., (coa_binned_genes_tpm_df %>%
                         tidyr::drop_na(bin) %>%
                         dplyr::summarise(sumTPM_onlybinned = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
                         dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(desc(sumTPM)) %>%
  dplyr::mutate(bin = factor(bin)) %>%
  dplyr::mutate(bin = factor(bin, levels = unique(dend_rows_taxonomy[["bin"]]))) %>%
  tidyr::pivot_longer(., cols = c(num_exp_genes, length_exp_genes, sumTPM, sumTPM_onlybinned),
                      names_to = "metric",
                      values_to = "counts") %>%
  dplyr::mutate(metric = factor(metric, levels = c("num_exp_genes", "length_exp_genes", "sumTPM", "sumTPM_onlybinned"))) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  dplyr::mutate(Fraction = factor(Fraction, levels = c("12H", "13L", "12L", "13H"))) %>%
  droplevels 

temp_df2_breaks <- temp_df2 %>%
  dplyr::filter(metric == "num_exp_genes") %>%
  dplyr::ungroup(.) %>%
  dplyr::reframe(quantiles = quantile(log10(counts), probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)) %>%
  # dplyr::reframe(quantiles = quantile(10^log_norm_TPM, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)) %>% log10(.) %>%
  # round(., digits = 0) %>%
  signif(., digits = 0) %>%
  ceiling(.) %>%
  unique(.) %>%
  deframe(.)

temp_df2_labels <- 10^(temp_df2_breaks) %>% 
  replace_na(., "") %>%
  as.character(.)
tpm_hm_color <- rev(colorRampPalette(pals::cubehelix(n = 20)[5:19])(10))

g4 <- print(
  ggplot(data = temp_df2 %>%
           # dplyr::filter(bin %in% keep_bins_idx) %>%
           dplyr::filter(metric == "num_exp_genes") %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(y = bin, group = bin,
                  # x = sample_name,
                  x = sample_source,
                  fill = log10(counts), color = log10(counts)),
              show.legend = TRUE)
  + scale_fill_gradientn(colors = tpm_hm_color,
                         breaks = temp_df2_breaks,
                         labels = temp_df2_labels,
                         aesthetics = "fill", expand = expansion(1.1,1.1),
                         na.value = "black")
  + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
  + scale_y_discrete(name = "Bin", labels = quality_bin_label_lookup)
  + facet_grid(rows = vars(grouping), 
               # cols = NULL,
               cols = vars(Fraction),
               scales = "free_y", space = "free",
               shrink = TRUE, drop = TRUE)
  + theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90), 
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  + labs(x = "Sample", y = "binned contig")
  + guides(fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = paste0("Number of \n", "expressed genes"),
                                 direction = "vertical"),
           color = "none")
)

#now normalize num of expressed genes by gene lengths
#and genome size, and number of genes in genomes

temp_df3 <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, contains("exp_genes")) %>%
  dplyr::left_join(., (bin_summary_tbl %>%
                         # dplyr::mutate(bin = gsub("MetaBat2_Bins_NEW.", "coa_bin_", Name)) %>%
                         dplyr::select(bin, Genome_Size, Total_Coding_Sequences) %>%
                         droplevels),
                       by = join_by(bin)) %>%
  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent, Fraction) %>%
                     droplevels, by = join_by(sample_name)) %>%
  dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
                                           .default = "low")) %>%
  dplyr::left_join(., (coa_genes_tpm_df %>%
                         tidyr::drop_na(bin) %>%
                         dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
                         dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  dplyr::left_join(., (coa_binned_genes_tpm_df %>%
                         tidyr::drop_na(bin) %>%
                         dplyr::summarise(sumTPM_onlybinned = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
                         dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(desc(sumTPM)) %>%
  dplyr::mutate(bin = factor(bin)) %>%
  dplyr::mutate(bin = factor(bin, levels = unique(dend_rows_taxonomy[["bin"]]))) %>%
  dplyr::mutate(norm_exp_genes_obs = (num_exp_genes / Total_Coding_Sequences),
                norm_exp_length = (length_exp_genes / Genome_Size),
                norm_exp_genes = (num_exp_genes / Total_Coding_Sequences) * (length_exp_genes / Genome_Size)) %>%
  dplyr::select(bin, sample_name, grouping, Temperature, Vent, Fraction, quality, contains("norm_"), contains("TPM")) %>%
  tidyr::pivot_longer(., cols = !c(bin, sample_name, grouping, Temperature, Vent, quality, Fraction),
                      names_to = "metric",
                      values_to = "counts") %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  dplyr::mutate(Fraction = factor(Fraction, levels = c("12H", "13L", "12L", "13H"))) %>%
  # dplyr::mutate(metric = factor(metric, levels = c("norm_exp_genes", "sumTPM", "sumTPM_onlybinned"))) %>%
  droplevels 

#sanity check:
wide_temp_df3 <- temp_df3 %>%
  dplyr::filter(grepl("13H", sample_name)) %>%
  dplyr::select(bin, grouping, sample_name, metric, counts) %>%
  dplyr::filter(grepl("norm_exp_length", metric)) %>%
  dplyr::mutate(counts = 100*counts) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = c(bin, grouping),
                     names_from = "sample_name",
                     values_fill = 0,
                     values_from = "counts")

temp_breaks <- temp_df3 %>%
  dplyr::filter(metric == "norm_exp_length") %>%
  dplyr::select(counts) %>% droplevels

g5_breaks <- temp_df3 %>%
  dplyr::filter(metric == "norm_exp_length") %>%
  dplyr::ungroup(.) %>%
  dplyr::reframe(quantiles = quantile(log10(counts), probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)) %>%
  signif(., digits = 0) %>%
  round(., digits = 0) %>%
  # unique(.) %>% deframe(.) #do you want to make it tidy?
  deframe(.) %>% c(., ceiling(log10(max(temp_breaks))), floor(log10(min(temp_breaks)))) %>% unique(.) %>% sort(.) #do you want to include the min and max on the scale?

g5_labels_basic <- 10^(g5_breaks) %>% 
  signif(., digits = 0)
g5_labels_basic <- 100*g5_labels_basic

# g5_labels <- 10^(g5_breaks) %>% 
#   # scales::scientific(., digits = 2, scale = 1)
#   # scales::number(., scale = 100)
#   scales::percent(., accuracy = NULL, trim = TRUE)

g5_hm_color <- (colorRampPalette(pals::brewer.ylorrd(n = 20)[1:15])(length(g5_breaks)))
g5 <- print(
  ggplot(data = temp_df3 %>%
           dplyr::filter(metric == "norm_exp_length") %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(y = bin, group = bin,
                  # x = sample_name, 
                  x = sample_source, 
                  fill = log10(counts), color = log10(counts)), 
              show.legend = TRUE)
  + scale_fill_gradientn(colors = g5_hm_color,
                         breaks = g5_breaks,
                         # labels = g5_labels,
                         labels = paste0(g5_labels_basic, "%"),
                         aesthetics = "fill", limits = c(min(g5_breaks), max(g5_breaks)),
                         # expand = expansion(1.1,1.1),
                         na.value = "black")
  + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
  + scale_y_discrete(name = "Bin", labels = quality_bin_label_lookup)
  + facet_grid(rows = vars(grouping), 
               # cols = NULL,
               cols = vars(Fraction),
               scales = "free_y", space = "free",
               shrink = TRUE, drop = TRUE)
  + theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90), 
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  + labs(x = "Sample", y = "binned contig")
  + guides(fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = paste0("Coding \n", "proportion"),
                                 direction = "vertical"),
           color = "none")
)


temp_breaks <- temp_df3 %>%
  dplyr::filter(metric == "norm_exp_genes") %>%
  dplyr::select(counts) %>% droplevels
g6_breaks <- 10^c(seq(from = min(temp_breaks, na.rm = TRUE), to = max(temp_breaks, na.rm = TRUE), length = 5)) %>%
  ceiling(.) %>%
  # signif(., digits = 0) %>%
  # unique(.) %>% log10(.)
  c(10^(min(temp_breaks)), .) %>%
  log10(.) %>% signif(., digits = 1) %>% unique(.) %>% sort(.) #do you want to include the min and max on the scale?

g6_labels <- (g6_breaks) %>%
  replace_na(., "") %>%
  # scales::number(., big.mark = ",", decimal.mark = ".", accuracy = 0.1)
  # scales::scientific(digits = 1, scale = 1)
  scales::number(trim = TRUE, accuracy = NULL)
  
# g6_hm_color <- rev(colorRampPalette(pals::cubehelix(n = 20)[1:10])(length(g6_breaks)))
g6_hm_color <- (colorRampPalette(pals::brewer.greys(n = 20)[5:20])(length(g6_breaks)))
g6 <- print(
  ggplot(data = temp_df3 %>%
           dplyr::filter(metric == "norm_exp_genes") %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(y = bin, group = bin,
                  # x = sample_name, 
                  x = sample_source,
                  fill = counts, color = counts), 
                  # fill = log10(counts), color = log10(counts)), 
              show.legend = TRUE)
  + scale_fill_gradientn(colors = g6_hm_color,
                         breaks = g6_breaks,
                         labels = g6_labels,
                         aesthetics = "fill", limits = c(0, 1),
                         # expand = expansion(1.1,1.1),
                         na.value = "black")
  + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
  + scale_y_discrete(name = "Bin", labels = quality_bin_label_lookup)
  + facet_grid(rows = vars(grouping), 
               # cols = NULL,
               cols = vars(Fraction),
               scales = "free_y", space = "free",
               shrink = TRUE, drop = TRUE)
  + theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90), 
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  + labs(x = "Sample", y = "binned contig")
  + guides(fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = paste0("Expression \n", "density"),
                                 direction = "vertical"),
           color = "none")
)


gpatch <- ((g4 + ggtitle("Number of expressed genes in each genome") + theme(strip.text.y = element_blank())) | (g5 + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + theme(strip.text.y = element_blank()) + ggtitle("Proportion of coding genome that is expressed")) | (g6 + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + ggtitle("Expression density"))) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Expression in bins across samples",
                                                                                                                                                                                                                                                                                                                                                                                                                        subtitle = "Coding Proportion = (length of expressed genes)/(genome length); \n Expression Density = (num expressed genes/number of genes in genome)*(Expressed proportion)",
                                                 tag_levels = "A")
gpatch
ggsave(paste0(projectpath, "/figures/", namevar, "_coa_bins_coding_density", ".png"),
       gpatch,
       width = 30, height = 15, units = "in")

# Statistics on the MAG bin activity metrics ------------------------------

#these were the metrics calculated for each MAG bin in each fraction and sample:
# norm_exp_genes = (num_exp_genes / Total_Coding_Sequences),
# norm_length_exp = (length_exp_genes / Genome_Size),
# density = (num_exp_genes / Total_Coding_Sequences) * (length_exp_genes / Genome_Size)


exp_lookup <- data.frame(metric = c("exp_genes", "num_exp_genes", "norm_exp_genes",
                "exp_length", "length_exp_genes",  "norm_length_exp",
                "norm_exp_density"),
                label = c("Total genes expressed in subset", "Expressed genes in sample", "Normalized expression",
                "Total length expressed in subset", "Length of expressed genes in sample", "Normalized expressed length",
                "Expression density")) %>%
  tibble::deframe(.)

temp_df5 <- (coa_binned_genes_tpm_df %>%
    dplyr::filter((TPM > 0)) %>%
    droplevels %>%
    dplyr::distinct(gene_id, bin, sample_name, .keep_all = FALSE) %>%
    droplevels %>%
      dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
      dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
      # dplyr::distinct(gene_id, bin, sample_source, .keep_all = TRUE) %>%
      dplyr::mutate(fraction = dplyr::case_when(grepl("(13L|12H)", sample_name) ~ "12H_13L",
                                                grepl("(13H|12L)", sample_name) ~ "13H_12L",
                                                .default = NA)) %>%
      dplyr::distinct(gene_id, bin, sample_source, fraction, .keep_all = TRUE) %>%
    dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                  gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
    dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
      dplyr::distinct(gene_id, bin, gene_length, fraction, .keep_all = FALSE) %>%
      droplevels %>%
              dplyr::group_by(bin, fraction) %>%
              dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE),
                               exp_genes = length(gene_id), .groups = "keep") %>%
    #   split(., f = .$sample_source) %>%
    #   map(., ~.x %>%
    #         dplyr::ungroup(.) %>%
    #         dplyr::select(bin, gene_length, gene_id) %>%
    #         dplyr::distinct(., .keep_all = TRUE) %>%
    #         droplevels %>%
    #         dplyr::group_by(bin) %>%
    #         dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE),
    #                          exp_genes = length(gene_id), .groups = "keep") %>%
    #           dplyr::ungroup(.) %>%
    #           dplyr::arrange(desc(exp_length), desc(exp_genes)) %>%
    #           # dplyr::select(-sample_source) %>%
    #           droplevels) %>%
    #   dplyr::bind_rows(., .id = "sample_source") %>%
    # # dplyr::mutate(bin = factor(bin)) %>%
    # # dplyr::group_by(bin, sample_source) %>%
    # # dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE),
    # #                  exp_genes = length(gene_id), .groups = "keep") %>%
    # #   dplyr::ungroup(.) %>%
    # #   dplyr::arrange(desc(exp_length), desc(exp_genes)) %>%
    # #   dplyr::select(-sample_source) %>%
    # #   droplevels %>%
    #   # dplyr::distinct(bin, .keep_all = TRUE) %>%
    # dplyr::mutate(fraction = "total") %>%
    droplevels)

binned_mags_exp_summary <- bind_rows(
  # (coa_binned_genes_tpm_df %>%
  #    dplyr::filter((TPM > 0)) %>%
  #    droplevels %>%
  #    dplyr::distinct(gene_id, bin, .keep_all = FALSE) %>%
  #    droplevels %>%
  #    dplyr::group_by(bin) %>%
  #    dplyr::summarise(exp_genes = length(gene_id), .groups = "keep") %>%
  #    dplyr::mutate(fraction = "total") %>%
  #    dplyr::select(bin, fraction, exp_genes) %>%
  #    dplyr::mutate(fraction = factor(fraction),
  #                  bin = factor(bin)) %>%
  #    droplevels),
  # (coa_binned_genes_tpm_df %>% #how many total genes were expressed by bin in 13H+12L and 13L+12H
  #    dplyr::filter((TPM > 0)) %>%
  #    dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
  #    dplyr::mutate(fraction = dplyr::case_when(grepl("(13H|12L)", fraction) ~ "13H_12L",
  #                                              grepl("(13L|12H)", fraction) ~ "12H_13L",
  #                                              .default =  NA)) %>%
  #    dplyr::select(gene_id, bin, fraction) %>%
  #    dplyr::mutate(fraction = factor(fraction),
  #                  bin = factor(bin)) %>%
  #    droplevels %>%
  #    dplyr::distinct(gene_id, bin, fraction, .keep_all = FALSE) %>%
  #    dplyr::group_by(bin, fraction) %>%
  #    dplyr::summarise(exp_genes = length(gene_id), .groups = "keep")),
                      (coa_binned_genes_tpm_df %>% #how many total genes were expressed by fraction and bin
                         dplyr::filter((TPM > 0)) %>%
                         dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
                         dplyr::select(gene_id, bin, fraction) %>%
                         droplevels %>%
                         dplyr::distinct(gene_id, bin, fraction, .keep_all = FALSE) %>%
                         dplyr::mutate(fraction = factor(fraction),
                                       bin = factor(bin)) %>%
                         dplyr::group_by(bin, fraction) %>%
                         dplyr::summarise(exp_genes = length(gene_id), .groups = "keep"))) %>%
  dplyr::arrange(bin) %>%
  dplyr::left_join(., 
                   (bind_rows((coa_binned_genes_tpm_df %>%
                                 dplyr::filter((TPM > 0)) %>%
                                 droplevels %>%
                                 dplyr::distinct(gene_id, bin, .keep_all = FALSE) %>%
                                 droplevels %>%
                                 dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                                               gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
                                 dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
                                 dplyr::rowwise(.) %>%
                                 dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
                                 dplyr::select(bin, gene_length) %>%
                                 dplyr::distinct(., .keep_all = TRUE) %>%
                                 dplyr::mutate(bin = factor(bin)) %>%
                                 dplyr::group_by(bin) %>%
                                 dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE)) %>%
                                 dplyr::mutate(fraction = "total") %>%
                                 droplevels),
                              (coa_binned_genes_tpm_df %>%
                                 dplyr::filter((TPM > 0)) %>%
                                 dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
                                 droplevels %>%
                                 dplyr::distinct(gene_id, bin, fraction, .keep_all = FALSE) %>%
                                 droplevels %>%
                                 dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                                               gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
                                 dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
                                 dplyr::rowwise(.) %>%
                                 dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
                                 dplyr::distinct(bin, fraction, gene_length, gene_id,  .keep_all = FALSE) %>%
                                 dplyr::select(bin, fraction, gene_length) %>%
                                 droplevels %>%
                                 dplyr::mutate(fraction = factor(fraction),
                                               bin = factor(bin)) %>%
                                 dplyr::group_by(bin, fraction) %>%
                                 dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE), .groups = "keep") %>%
                                 droplevels),
                              (coa_binned_genes_tpm_df %>%
                                 dplyr::filter((TPM > 0)) %>%
                                 dplyr::filter(grepl("(13H|12L)", sample_name)) %>%
                                 dplyr::mutate(fraction = "13H_12L") %>%
                                 droplevels %>%
                                 dplyr::distinct(gene_id, bin, fraction, .keep_all = FALSE) %>%
                                 droplevels %>%
                                 dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                                               gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
                                 dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
                                 dplyr::rowwise(.) %>%
                                 dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
                                 dplyr::distinct(bin, fraction, gene_length, gene_id,  .keep_all = FALSE) %>%
                                 dplyr::select(bin, fraction, gene_length) %>%
                                 droplevels %>%
                                 dplyr::mutate(fraction = factor(fraction),
                                               bin = factor(bin)) %>%
                                 dplyr::group_by(bin, fraction) %>%
                                 dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE), .groups = "keep") %>%
                                 droplevels),
                              (coa_binned_genes_tpm_df %>%
                                 dplyr::filter((TPM > 0)) %>%
                                 dplyr::filter(grepl("(13L|12H)", sample_name)) %>%
                                 dplyr::mutate(fraction = "12H_13L") %>%
                                 droplevels %>%
                                 dplyr::distinct(gene_id, bin, fraction, .keep_all = FALSE) %>%
                                 droplevels %>%
                                 dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                                               gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
                                 dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
                                 dplyr::rowwise(.) %>%
                                 dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
                                 dplyr::distinct(bin, fraction, gene_length, gene_id,  .keep_all = FALSE) %>%
                                 dplyr::select(bin, fraction, gene_length) %>%
                                 droplevels %>%
                                 dplyr::mutate(fraction = factor(fraction),
                                               bin = factor(bin)) %>%
                                 dplyr::group_by(bin, fraction) %>%
                                 dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE), .groups = "keep") %>%
                                 droplevels)) %>%
                      dplyr::arrange(bin)),
                   by = join_by(bin, fraction)) %>%
  dplyr::bind_rows(., (coa_binned_genes_tpm_df %>%
                         dplyr::filter((TPM > 0)) %>%
                         droplevels %>%
                         dplyr::distinct(gene_id, bin, sample_name, .keep_all = FALSE) %>%
                         droplevels %>%
                         dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
                         dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
                         dplyr::distinct(gene_id, bin, sample_source, .keep_all = TRUE) %>%
                         dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                                       gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
                         dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
                         dplyr::rowwise(.) %>%
                         dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
                         dplyr::distinct(gene_id, bin, gene_length, .keep_all = FALSE) %>%
                         droplevels %>%
                         dplyr::group_by(bin) %>%
                         dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE),
                                          exp_genes = length(gene_id), .groups = "keep") %>%
                         dplyr::mutate(fraction = "total") %>%
                         droplevels)) %>%
  dplyr::bind_rows(., (coa_binned_genes_tpm_df %>%
                         dplyr::filter((TPM > 0)) %>%
                         droplevels %>%
                         dplyr::distinct(gene_id, bin, sample_name, .keep_all = FALSE) %>%
                         droplevels %>%
                         dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
                         dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
                         dplyr::mutate(fraction = dplyr::case_when(grepl("(13L|12H)", sample_name) ~ "12H_13L",
                                                                   grepl("(13H|12L)", sample_name) ~ "13H_12L",
                                                                   .default = NA)) %>%
                         dplyr::distinct(gene_id, bin, sample_source, fraction, .keep_all = TRUE) %>%
                         dplyr::mutate(gene_start = stringr::str_split_i(gene_id, "_", 3),
                                       gene_stop = stringr::str_split_i(gene_id, "_", 4)) %>%
                         dplyr::mutate(across(c(gene_start, gene_stop), ~as.numeric(.x))) %>%
                         dplyr::rowwise(.) %>%
                         dplyr::mutate(gene_length = (gene_stop - gene_start + 1)) %>%
                         dplyr::distinct(gene_id, bin, gene_length, fraction, .keep_all = FALSE) %>%
                         droplevels %>%
                         dplyr::group_by(bin, fraction) %>%
                         dplyr::summarise(exp_length = sum(gene_length, na.rm = TRUE),
                                          exp_genes = length(gene_id), .groups = "keep") %>%
                         droplevels)) %>%
  dplyr::mutate(fraction = factor(fraction),
                bin = factor(bin)) %>%
  split(., f = .$fraction) %>%
  map(., ~.x %>%
        dplyr::ungroup(.) %>%
        # dplyr::select(-fraction) %>%
        droplevels)
  
binned_genes_tpm_mags_summary <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
  dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::mutate(fraction = factor(fraction),
                bin = factor(bin)) %>%
  bind_rows(., binned_contig_counts_vents_df %>%
              dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
              dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
              dplyr::mutate(fraction = dplyr::case_when(grepl("(13H|12L)", fraction) ~ "13H_12L",
                                                        grepl("(13L|12H)", fraction) ~ "12H_13L",
                                                        .default =  NA)) %>%
              dplyr::mutate(fraction = factor(fraction),
                            bin = factor(bin)) %>%
              droplevels) %>%
  bind_rows(., (binned_contig_counts_vents_df %>%
                  dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
                  dplyr::mutate(fraction = "total") %>%
                  dplyr::mutate(fraction = factor(fraction),
                                bin = factor(bin)))) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::mutate(fraction = factor(fraction, levels = c("13H", "12L", "13L", "12H", "13H_12L", "12H_13L", "total")),
                bin = factor(bin)) %>%
  dplyr::distinct(bin, sample_name, fraction, .keep_all = TRUE) %>%
  droplevels %>%
  split(., f = .$fraction) %>%
  map(., ~.x %>%
        droplevels %>%
        # dplyr::select(-fraction) %>%
        dplyr::distinct(bin, sample_name, fraction, .keep_all = TRUE) %>%
        droplevels)

# temp_df5 <- imap(binned_genes_tpm_mags_summary,
binned_genes_tpm_mags_summary <- imap(binned_genes_tpm_mags_summary,
                 ~.x %>%
                   as.data.frame(.) %>%
                   dplyr::left_join(., binned_mags_exp_summary[[.y]],
# temp_df5 <- names(binned_genes_tpm_mags_summary) %>%
# # binned_genes_tpm_mags_summary <- names(binned_genes_tpm_mags_summary) %>%
#   map(., ~binned_genes_tpm_mags_summary[[.x]] %>%
#         tibble::as_tibble(.) %>%
#         droplevels %>%
#         dplyr::left_join(., binned_mags_exp_summary[[.x]] %>%
#                            tibble::as_tibble(.) %>%
#                            droplevels,
                         by = join_by(bin, fraction),
                         relationship = "many-to-many", multiple = "all") %>%
        dplyr::distinct(bin, sample_name, num_exp_genes, length_exp_genes, .keep_all = TRUE) %>%
        dplyr::rowwise(.) %>%
        dplyr::mutate(norm_exp_genes = (num_exp_genes / exp_genes),
                         norm_length_exp = (length_exp_genes / exp_length)) %>%
        dplyr::select(-c(exp_genes, exp_length)) %>%
#         droplevels) %>%
# setNames(., names(binned_mags_exp_summary))
  droplevels)

{
  temp_list <- list(binned_mags_exp_summary, binned_genes_tpm_mags_summary) %>%
    setNames(., c("binned_mags_exp_summary", "binned_genes_tpm_mags_summary"))
  readr::write_rds(temp_list, paste0(projectpath, "/output/", "binned_mags_exp_summary", ".rds"), compress = "gz")
  rm(temp_list)

  # binned_genes_tpm_mags_summary <- binned_contig_counts_vents_df %>%
#   dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
#   dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
  # dplyr::left_join(., (coa_binned_genes_tpm_df %>%
  #                        dplyr::filter((TPM > 0)) %>%
  #                        dplyr::distinct(gene_id, bin, .keep_all = FALSE) %>%
  #                        dplyr::group_by(bin) %>%
  #                        dplyr::summarise(total_exp_genes = length(gene_id), .groups = "keep")),
  #                  by = join_by(bin), relationship = "many-to-many", multiple = "all") %>%
  # dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
  # dplyr::left_join(., (coa_binned_genes_tpm_df %>% #how many total genes were expressed by fraction and bin
  #                        dplyr::filter((TPM > 0)) %>%
  #                        dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
  #                        dplyr::select(gene_id, bin, fraction) %>%
  #                        droplevels %>%
  #                        dplyr::distinct(gene_id, bin, fraction, .keep_all = TRUE) %>%
  #                        dplyr::group_by(bin, fraction) %>%
  #                        dplyr::summarise(frac_exp_genes = length(gene_id), .groups = "keep")),
  #                  by = join_by(bin, fraction)) %>%
  # dplyr::mutate(heavy_light = dplyr::case_when(grepl("(13H|12L)", fraction) ~ "13H_12L",
  #                                              .default = "12H_13L")) %>%
  # dplyr::left_join(., (coa_binned_genes_tpm_df %>% #how many total genes were expressed by bin and in either 13H+12L or 12H+13L
  #                        dplyr::filter((TPM > 0)) %>%
  #                        dplyr::mutate(heavy_light = dplyr::case_when(grepl("(13H|12L)", sample_name) ~ "13H_12L",
  #                                                                     .default = "12H_13L")) %>%
  #                        dplyr::select(gene_id, bin, heavy_light) %>%
  #                        dplyr::distinct(gene_id, bin,heavy_light, .keep_all = TRUE) %>%
  #                        droplevels %>%
  #                        dplyr::group_by(bin, heavy_light) %>%
  #                        dplyr::summarise(heavy_light_exp_genes = length(gene_id), .groups = "keep")),
  #                  by = join_by(bin, heavy_light), relationship = "many-to-many", multiple = "all") %>%
  # dplyr::left_join(., (binned_contig_counts_vents_df %>%
  #                        dplyr::select(bin, sample_name, length_exp_genes) %>%
  #                        dplyr::mutate(fraction = stringr::str_extract(sample_name, "(13H|13L|12H|12L)")) %>%
  #                        dplyr::mutate(fraction = factor(fraction, levels = c("13H", "12L", "12H", "13L"))) %>%
  #                        dplyr::group_by(bin, sample_name, fraction) %>%
  #                        dplyr::summarise(frac_length_exp_genes = sum(length_exp_genes), .groups = "keep") %>%
  #                        # dplyr::select(-length_exp_genes) %>%
  #                        dplyr::distinct(bin, sample_name, fraction, frac_length_exp_genes, .keep_all = FALSE) %>%
  #                        droplevels), 
  #                  by = join_by(bin, fraction), relationship = "many-to-many", multiple = "all") %>%
  # dplyr::left_join(., (binned_contig_counts_vents_df %>%
  #                        dplyr::select(bin, sample_name, length_exp_genes) %>%
  #                        dplyr::mutate(heavy_light = dplyr::case_when(grepl("(13H|12L)", sample_name) ~ "13H_12L",
  #                                                                     .default = "12H_13L")) %>%
  #                        dplyr::mutate(heavy_light = factor(heavy_light, levels = c("13H_12L", "12H_13L"))) %>%
  #                        dplyr::group_by(bin, sample_name, heavy_light) %>%
  #                        dplyr::summarise(heavy_light_length_exp_genes = sum(length_exp_genes), .groups = "keep") %>%
  #                        # dplyr::select(-length_exp_genes) %>%
  #                        dplyr::distinct(bin, sample_name, heavy_light, heavy_light_length_exp_genes, .keep_all = FALSE) %>%
  #                        droplevels),
  #                  by = join_by(bin, heavy_light), relationship = "many-to-many", multiple = "all") %>%
  # dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes, starts_with("frac"), starts_with("heavy_light"), total_exp_genes) %>%
  # droplevels

  }

binned_genes_tpm_mags_df <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, contains("obs_genes")) %>%
  dplyr::left_join(., (bin_summary_tbl %>%
                         dplyr::mutate(bin = gsub("MetaBat2_Bins_NEW.", "coa_bin_", Name)) %>%
                         dplyr::select(bin, Genome_Size, Total_Coding_Sequences) %>%
                         droplevels),
                   by = join_by(bin)) %>%
  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  # dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent, Fraction) %>%
                     droplevels, by = join_by(sample_name)) %>%
  dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
                                           .default = "low")) %>%
  dplyr::left_join(., (coa_genes_tpm_df %>%
                         tidyr::drop_na(bin) %>%
                         dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
                         dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  dplyr::left_join(., (coa_binned_genes_tpm_df %>%
                         tidyr::drop_na(bin) %>%
                         dplyr::summarise(sumTPM_onlybinned = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
                         dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(desc(sumTPM)) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  dplyr::mutate(Fraction = factor(Fraction, levels = c("12H", "13L", "12L", "13H"))) %>%
  dplyr::mutate(bin = factor(bin)) %>%
  dplyr::mutate(bin = factor(bin, levels = unique(dend_rows_taxonomy[["bin"]])))


# Plot the expression density, proportion across fractions ----------------

#plot the num_exp_genes in each bin across sample_sources by fraction

# temp_df4 <- binned_genes_tpm_mags_summary %>%
#   dplyr::bind_rows(., .id = "subset") %>%
#   dplyr::mutate(subset = factor(subset, levels = unique(.[["subset"]]))) %>%
#   dplyr::mutate(subset = factor(subset, levels = c("12H", "13L", "12H_13L", "12L", "13H", "13H_12L", "total"))) %>%
#   dplyr::select(-c(num_exp_genes, length_exp_genes)) %>%
#   dplyr::distinct(., .keep_all = TRUE) %>%
#   dplyr::filter(!(bin %in% drop_bins_idx)) %>%
#   dplyr::filter((bin %in% keep_bins_idx)) %>%
#   droplevels %>%
#   dplyr::left_join(., highq_mtg_grouping_list,
#                    by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
#   dplyr::left_join(., metadata %>%
#                      dplyr::select(sample_name, Temperature, Vent) %>%
#                      droplevels, by = join_by(sample_name)) %>%
#   tidyr::pivot_longer(., cols = contains("norm"),
#                       names_to = "metric",
#                       values_to = "counts") %>%
#   # dplyr::mutate(Fraction = factor(Fraction, levels = c("13H", "12L", "12H", "13L"))) %>%
#   dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
#   dplyr::arrange(sample_name) %>%
#   dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
#   dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
#   dplyr::relocate(bin, grouping, sample_name, Temperature, Vent, subset, sample_source) %>%
#   dplyr::mutate(metric = factor(metric, levels = c("norm_exp_genes", "norm_length_exp"))) %>%
#   # dplyr::slice_max(., counts, by = c(bin, sample_source, subset, metric)) %>%
#   droplevels
# 
# 
# print(
#   ggplot(data = temp_df4 %>%
#            dplyr::filter(grepl("coa_bin_320", bin)) %>%
#            droplevels,
#          aes(fill = metric, group = interaction(subset, metric, sample_source)))
#   + theme_bw() 
#   + geom_col(aes(x = sample_source, y = counts), alpha = 0.5, color = "grey20", width = 0.90, show.legend = FALSE, position = "identity")
#   + scale_fill_viridis(discrete = TRUE, option = "turbo", labels = exp_lookup, breaks = names(exp_lookup))
#   + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
#   + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Normalized")
#   + facet_grid(metric ~ subset,
#                scales = "free_y", space = "free",
#                labeller = labeller(metric = exp_lookup),
#                shrink = TRUE, drop = TRUE)
#   + theme(strip.text.y = element_text(angle = 90),
#           axis.text.x = element_text(angle = 90), 
#           # axis.ticks.y = element_blank(),
#           # axis.text.y = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank())
# )

temp_df4 <- c("13H", "13H_12L", "total") %>%
  map(., ~binned_genes_tpm_mags_summary %>%
        purrr::pluck(.x)) %>%
  bind_rows(.) %>%
  dplyr::rename(subset = "fraction") %>%
  dplyr::select(-c(num_exp_genes, length_exp_genes)) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  dplyr::filter(grepl("13H", sample_name)) %>%
  droplevels %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent) %>%
                     droplevels, by = join_by(sample_name)) %>%
  # dplyr::group_by(bin, sample_name, subset) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(norm_exp_density = (norm_exp_genes * norm_length_exp)) %>%
  tidyr::pivot_longer(., cols = contains("norm"),
                      names_to = "metric",
                      values_to = "counts") %>%
  # dplyr::mutate(Fraction = factor(Fraction, levels = c("13H", "12L", "12H", "13L"))) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  dplyr::relocate(bin, grouping, sample_name, Temperature, Vent, subset, sample_source) %>%
  dplyr::mutate(metric = factor(metric, levels = names(exp_lookup))) %>%
  # dplyr::mutate(metric = factor(metric, levels = c("norm_exp_genes", "norm_length_exp"))) %>%
  # dplyr::slice_max(., counts, by = c(bin, sample_source, subset, metric)) %>%
  droplevels %>%
  bind_rows(., (binned_contig_counts_vents_df %>%
                  dplyr::select(bin, sample_name, contains("exp_genes")) %>%
                  dplyr::filter(grepl("13H", sample_name)) %>%
                  dplyr::left_join(., (bin_summary_tbl %>%
                                         dplyr::mutate(bin = gsub("MetaBat2_Bins_NEW.", "coa_bin_", Name)) %>%
                                         dplyr::select(bin, Genome_Size, Total_Coding_Sequences) %>%
                                         droplevels),
                                   by = join_by(bin)) %>%
                  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
                  dplyr::left_join(., highq_mtg_grouping_list,
                                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
                  dplyr::filter((bin %in% keep_bins_idx)) %>%
                  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
                  dplyr::left_join(., metadata %>%
                                     dplyr::select(sample_name, Temperature, Vent) %>%
                                     droplevels, by = join_by(sample_name)) %>%
                  dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
                                                           .default = "low")) %>%
                  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
                  dplyr::mutate(bin = factor(bin)) %>%
                  dplyr::mutate(bin = factor(bin, levels = unique(dend_rows_taxonomy[["bin"]]))) %>%
                  dplyr::mutate(norm_exp_genes = (num_exp_genes / Total_Coding_Sequences),
                                norm_length_exp = (length_exp_genes / Genome_Size),
                                norm_exp_density = (num_exp_genes / Total_Coding_Sequences) * (length_exp_genes / Genome_Size)) %>%
                  dplyr::mutate(subset = "whole_genome_potential") %>%
                  dplyr::select(bin, sample_name, grouping, subset, Temperature, Vent, quality, contains("norm_")) %>%
                  tidyr::pivot_longer(., cols = !c(bin, sample_name, grouping, Temperature, Vent, quality, subset),
                                      names_to = "metric",
                                      values_to = "counts") %>%
                  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
                  dplyr::arrange(sample_name) %>%
                  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
                  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
                  # dplyr::mutate(Fraction = factor(Fraction, levels = c("12H", "13L", "12L", "13H"))) %>%
                  # dplyr::mutate(metric = factor(metric, levels = c("norm_exp_genes", "sumTPM", "sumTPM_onlybinned"))) %>%
                  droplevels )) %>%
  dplyr::mutate(bin = factor(bin)) %>%
  dplyr::mutate(bin = factor(bin, levels = unique(dend_rows_taxonomy[["bin"]]))) %>%
  dplyr::arrange(grouping, bin) %>%
  droplevels


#what kind of statistics do we have for norm_exp_density and norm_exp_length
temp_df4 %>%
  dplyr::filter(!grepl("norm_exp_genes", metric)) %>%
  # dplyr::group_by(bin, metric, subset) %>%
  dplyr::group_by(metric, subset) %>%
  dplyr::summarise(across("counts", list(min = min, median = median, mean = mean, max = max)), .groups = "keep") %>%
  # dplyr::filter(grepl("13H", subset)) %>%
  droplevels

g8_parsed <- print(
    ggplot(data = temp_df4 %>%
             # dplyr::filter(grepl("norm_exp", metric)) %>%
             dplyr::filter(grepl("norm_exp_density", metric)) %>%
             droplevels,
           aes(x = subset, y = counts,
               # group = interaction(bin, Temperature)
               # group = bin
               group = interaction(bin, subset)
               ))
    + geom_boxplot(position = position_dodge(0.8), show.legend = FALSE, color = "black", alpha = 1)
    # + geom_violin(position = position_dodge(0.8), show.legend = FALSE, color = "black", alpha = 1)
    + geom_point(aes(fill = Temperature), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8),
                 alpha = 0.7, shape = 21, color = "black", size = 2, show.legend = TRUE)
    + theme_bw()
    # + scale_x_discrete(name = "Bin", labels = quality_bin_label_lookup)
    + scale_y_continuous(name = "normalized expression",
                        labels = scales::number)
    + scale_fill_manual(values = viridisLite::turbo(3),
                        # labels = c("30", "55", "80"),
                        # breaks = c("30", "55", "80"), 
                        drop = FALSE)
    + facet_wrap(bin~., ncol = 3, shrink = TRUE, scales = "free_y", labeller = labeller(bin = quality_bin_label_lookup), drop = TRUE)
    # + facet_grid(bin~., space = "free", labeller = labeller(bin = quality_bin_label_lookup),
    # # + facet_grid(subset~metric, scales = "free_y", space = "free",
    #              shrink = TRUE, drop = TRUE)
    + theme(strip.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 90), 
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
)

ggsave(paste0(projectpath, "/figures/", "fractionated_bins_exp_density_bw", ".png"),
       g8_parsed,
       width = 10, height = 16, units = "in")

temp_breaks <- temp_df4 %>%
  dplyr::filter(grepl("norm_length_exp", metric)) %>%
  dplyr::select(counts) %>% droplevels

g7_breaks <- temp_df4 %>%
  dplyr::filter(grepl("norm_length_exp", metric)) %>%
  dplyr::ungroup(.) %>%
  dplyr::reframe(quantiles = quantile(log10(counts), probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)) %>%
  signif(., digits = 0) %>%
  round(., digits = 0) %>%
  # unique(.) %>% deframe(.) #do you want to make it tidy?
  deframe(.) %>% c(., ceiling(log10(max(temp_breaks))), floor(log10(min(temp_breaks)))) %>% unique(.) %>% sort(.) #do you want to include the min and max on the scale?

g7_labels_basic <- 10^(g7_breaks) %>% 
  signif(., digits = 0)
g7_labels_basic <- 100*g7_labels_basic

g7_hm_color <- (colorRampPalette(pals::brewer.ylorrd(n = 20)[1:15])(length(g7_breaks)))
g7 <- print(
  ggplot(data = temp_df4 %>%
           dplyr::filter(metric == "norm_length_exp") %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(y = bin, group = bin,
                  # x = sample_name, 
                  x = sample_source, 
                  fill = log10(counts), color = log10(counts)), 
              show.legend = TRUE)
  + scale_fill_gradientn(colors = g7_hm_color,
                         breaks = g7_breaks,
                         labels = paste0(g7_labels_basic, "%"),
                         aesthetics = "fill", limits = c(min(g7_breaks), max(g7_breaks)),
                         # expand = expansion(1.1,1.1),
                         na.value = "black")
  + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
  + scale_y_discrete(name = "Bin", labels = quality_bin_label_lookup)
  + facet_grid(rows = vars(grouping), 
               # cols = NULL,
               cols = vars(subset),
               scales = "free_y", space = "free",
               shrink = TRUE, drop = TRUE)
  + theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90), 
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  # + labs(x = "Sample", y = "binned contig")
  + guides(fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = paste0("Coding \n", "proportion"),
                                 direction = "vertical"),
           color = "none")
)



temp_breaks <- temp_df4 %>%
  dplyr::filter(metric == "norm_exp_density") %>%
  dplyr::select(counts) %>% droplevels
g8_breaks <- 10^c(seq(from = min(temp_breaks, na.rm = TRUE), to = max(temp_breaks, na.rm = TRUE), length = 5)) %>%
  ceiling(.) %>%
  # signif(., digits = 0) %>%
  # unique(.) %>% log10(.)
  c(10^(min(temp_breaks)), .) %>%
  log10(.) %>% signif(., digits = 1) %>% unique(.) %>% sort(.) #do you want to include the min and max on the scale?

g8_labels <- (g8_breaks) %>%
  replace_na(., "") %>%
  # scales::number(., big.mark = ",", decimal.mark = ".", accuracy = 0.1)
  # scales::scientific(digits = 1, scale = 1)
  scales::number(trim = TRUE, accuracy = NULL)

# g8_hm_color <- rev(colorRampPalette(pals::cubehelix(n = 20)[1:10])(length(g8_breaks)))
g8_hm_color <- (colorRampPalette(pals::brewer.greys(n = 20)[5:20])(length(g8_breaks)))
g8 <- print(
  ggplot(data = temp_df4 %>%
           dplyr::filter(metric == "norm_exp_density") %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(y = bin, group = bin,
                  # x = sample_name, 
                  x = sample_source,
                  fill = counts, color = counts), 
              # fill = log10(counts), color = log10(counts)), 
              show.legend = TRUE)
  + scale_fill_gradientn(colors = g8_hm_color,
                         breaks = g8_breaks,
                         labels = g8_labels,
                         aesthetics = "fill", limits = c(0, 1),
                         # expand = expansion(1.1,1.1),
                         na.value = "black")
  + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
  + scale_y_discrete(name = "Bin", labels = quality_bin_label_lookup)
  + facet_grid(rows = vars(grouping), 
               # cols = NULL,
               cols = vars(subset),
               scales = "free_y", space = "free",
               shrink = TRUE, drop = TRUE)
  + theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90), 
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  # + labs(x = "Sample", y = "binned contig")
  + guides(fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = paste0("Expression \n", "density"),
                                 direction = "vertical"),
           color = "none")
)

temp_df2_13H <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
  dplyr::filter(grepl("13H", sample_name)) %>%
  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent, Fraction) %>%
                     droplevels, by = join_by(sample_name)) %>%
  dplyr::mutate(quality = dplyr::case_when((bin %in% names(quality_bin_label_lookup)) ~ "high",
                                           .default = "low")) %>%
  # dplyr::left_join(., (coa_genes_tpm_df %>%
  #                        tidyr::drop_na(bin) %>%
  #                        dplyr::summarise(sumTPM = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
  #                        dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  # dplyr::left_join(., (coa_binned_genes_tpm_df %>%
  #                        tidyr::drop_na(bin) %>%
  #                        dplyr::summarise(sumTPM_onlybinned = sum(TPM, na.rm = TRUE), .by = c(bin, sample_name)) %>%
  #                        dplyr::ungroup(.)), by = join_by(bin, sample_name)) %>%
  # dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  # dplyr::arrange(desc(sumTPM)) %>%
  # dplyr::mutate(bin = factor(bin)) %>%
  # dplyr::mutate(bin = factor(bin, levels = unique(dend_rows_taxonomy[["bin"]]))) %>%
  # tidyr::pivot_longer(., cols = c(num_exp_genes, length_exp_genes, sumTPM, sumTPM_onlybinned),
  #                     names_to = "metric",
  #                     values_to = "counts") %>%
  # dplyr::mutate(metric = factor(metric, levels = c("num_exp_genes", "length_exp_genes", "sumTPM", "sumTPM_onlybinned"))) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  # dplyr::mutate(Fraction = factor(Fraction, levels = c("12H", "13L", "12L", "13H"))) %>%
  droplevels 

coa_highq_bins_num_exp_genes_tbl <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
  dplyr::filter(grepl("13H", sample_name)) %>%
  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent, Fraction) %>%
                     droplevels, by = join_by(sample_name)) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  dplyr::left_join(., tibble::enframe(quality_bin_label_lookup, name = "bin", value = "relabel"),
                   by = join_by(bin)) %>%
  dplyr::select(bin, relabel, num_exp_genes, sample_name, grouping) %>%
  tidyr::pivot_wider(., id_cols = c(bin, relabel, grouping),
                     names_from = "sample_name",
                     values_fill = 0,
                     values_from = "num_exp_genes")
if(!file.exists(paste0(projectpath, "/output/", "coa_highq_bins_num_exp_genes_tbl", ".tsv"))){
  readr::write_delim(coa_highq_bins_num_exp_genes_tbl, paste0(projectpath, "/output/", "coa_highq_bins_num_exp_genes_tbl", ".tsv"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
}

coa_highq_bins_length_exp_genes_tbl <- binned_contig_counts_vents_df %>%
  dplyr::select(bin, sample_name, num_exp_genes, length_exp_genes) %>%
  dplyr::filter(grepl("13H", sample_name)) %>%
  dplyr::filter(bin %in% names(quality_bin_label_lookup)) %>%
  dplyr::left_join(., highq_mtg_grouping_list,
                   by = join_by(bin), multiple = "all", relationship = "many-to-many") %>%
  dplyr::filter((bin %in% keep_bins_idx)) %>%
  dplyr::filter(!(bin %in% drop_bins_idx)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_name, Temperature, Vent, Fraction) %>%
                     droplevels, by = join_by(sample_name)) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = rownames(dend_colorbar_col))) %>%
  dplyr::arrange(sample_name) %>%
  dplyr::mutate(sample_source = stringr::str_remove_all(sample_name, "_(13H|13L|12H|12L)")) %>%
  dplyr::mutate(sample_source = factor(sample_source, levels = unique(.[["sample_source"]]))) %>%
  dplyr::left_join(., tibble::enframe(quality_bin_label_lookup, name = "bin", value = "relabel"),
                   by = join_by(bin)) %>%
  dplyr::select(bin, relabel, length_exp_genes, sample_name, grouping) %>%
  tidyr::pivot_wider(., id_cols = c(bin, relabel, grouping),
                     names_from = "sample_name",
                     values_fill = 0,
                     values_from = "length_exp_genes")
if(!file.exists(paste0(projectpath, "/output/", "coa_highq_bins_length_exp_genes_tbl", ".tsv"))){
  readr::write_delim(coa_highq_bins_length_exp_genes_tbl, paste0(projectpath, "/output/", "coa_highq_bins_length_exp_genes_tbl", ".tsv"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
}

temp_df2_breaks <- temp_df2_13H %>%
  dplyr::filter(metric == "num_exp_genes") %>%
  dplyr::ungroup(.) %>%
  dplyr::reframe(quantiles = quantile(log10(counts), probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)) %>%
  # dplyr::reframe(quantiles = quantile(10^log_norm_TPM, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)) %>% log10(.) %>%
  # round(., digits = 0) %>%
  signif(., digits = 0) %>%
  ceiling(.) %>%
  unique(.) %>%
  deframe(.)

temp_df2_labels <- 10^(temp_df2_breaks) %>% 
  replace_na(., "") %>%
  as.character(.)
tpm_hm_color <- rev(colorRampPalette(pals::cubehelix(n = 20)[5:19])(10))

g4_13H <- print(
  ggplot(data = temp_df2_13H %>%
           # dplyr::filter(bin %in% keep_bins_idx) %>%
           dplyr::filter(metric == "num_exp_genes") %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(y = bin, group = bin,
                  # x = sample_name,
                  x = sample_source,
                  fill = log10(counts), color = log10(counts)),
              show.legend = TRUE)
  + scale_fill_gradientn(colors = tpm_hm_color,
                         breaks = temp_df2_breaks,
                         labels = temp_df2_labels,
                         aesthetics = "fill", expand = expansion(1.1,1.1),
                         na.value = "black")
  + scale_x_discrete(name = "Sample name", labels = x_histogram_labels)
  + scale_y_discrete(name = "Bin", labels = quality_bin_label_lookup)
  + facet_grid(rows = vars(grouping), 
               # cols = NULL,
               cols = vars(Fraction),
               scales = "free_y", space = "free",
               shrink = TRUE, drop = TRUE)
  + theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90), 
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  + labs(x = "Sample", y = "binned contig")
  + guides(fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = paste0("Number of \n", "expressed genes"),
                                 direction = "vertical"),
           color = "none")
)
gpatch_layout <- "
    AABBBBCCCCC
    AABBBBCCCCC
    AABBBBCCCCC
  "
gpatch2 <- ((g4_13H + theme(axis.title.y = element_blank(), strip.text.y = element_blank()) + ggtitle("Number of expressed genes in each bin in each sample")) | (g7 + theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + theme(strip.text.y = element_blank()) + ggtitle("Proportion of coding genome that is expressed")) | (g8 + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + ggtitle("Expression density"))) + patchwork::plot_layout(guides = "collect", design = gpatch_layout) + patchwork::plot_annotation(title = "Expression in bins across samples",
                                                                                                                                                                                                                                                                                                                                                     subtitle = "Coding Proportion = (length of expressed genes)/(genome length); \n Expression Density = (num expressed genes/number of genes in genome)*(Expressed proportion)",
                                                                                                                                                                                                                                                                                                                                                     tag_levels = "A")
gpatch2
ggsave(paste0(projectpath, "/figures/", "fractionated_coa_bins_coding_density", ".png"),
       gpatch2,
       width = 30, height = 15, units = "in")


