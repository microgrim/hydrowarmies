## 00_global_labellers.R

#project-wide renamers and lookups

vents <- c("anemone", "marker113", "marker33")

set.alpha <- 0.05

vent_lookup <- data.frame(vent = c("marker113", "anemone",  "marker33"),
                 label = c("Marker 113", "Anemone", "Marker 33")) %>%
  tibble::deframe(.)

contrast_lookup <- data.frame(shortname = c("axial_temp_13H", "axial_30_13H", "axial_55_13H", "axial_80_13H", "marker113", "anemone", "marker33", "axial_30_both", "axial_55_both", "axial_80_both", "axial_temp_both"),
                              label = c("13H samples as baseline", "13H samples as baseline", "13H samples as baseline", "13H samples as baseline", "Marker 113", "Anemone", "Marker 33", "13H and 12L samples as baseline", "13H and 12L samples as baseline", "13H and 12L samples as baseline", "13H and 12L samples as baseline")) %>%
  tibble::deframe(.)

# contrast_lookup <- c(axial_temp_13H = "13H samples as baseline",
#                      axial_30_13H = "13H samples as baseline",
#                      axial_55_13H = "13H samples as baseline",
#                      axial_80_13H = "13H samples as baseline",
#                      "marker113" = "Marker 113",
#                      "anemone" = "Anemone",
#                      "marker33" = "Marker 33",
#                      # axial_temp_12L = "12L samples as baseline",
#                      axial_30_both = "13H and 12L samples as baseline",
#                      axial_55_both = "13H and 12L samples as baseline",
#                      axial_80_both = "13H and 12L samples as baseline",
#                      axial_temp_both = "13H and 12L samples as baseline")
# log2_lookup <- c("log2FC_TPM" = "log2FC(mean TPM)",
#                  "log2FoldChange" = "DESeq log2FC(normalized counts)",
#                  "log2FoldChange_MMSE" = "DESeq log2FC corrected")
# temp_lookup <- c("30" = "13H of vents at 30˚C",
#                  "55" = "13H of vents at 55˚C",
#                  "80" = "13H of vents at 80˚C",
#                  "all" = "13H of vents at all temperatures")
# baseline_lookup <- c("all" = "13H of vents at all temperatures",
#                      "30" = "30˚C",
#                      "55" = "55˚C",
#                      "80" = "80˚C",
#                      # "30" = "13H of vents at 30˚C",
#                      # "55" = "13H of vents at 55˚C",
#                      # "80" = "13H of vents at 80˚C",
#                      "vent_13H_30_2013" = "30˚C in 2013",
#                      "vent_13H_30_2014" = "30˚C in 2014",
#                      "vent_13H_55_2013" = "55˚C in 2013",
#                      "vent_13H_55_2014" = "55˚C in 2014",
#                      "vent_13H_80_2013" = "80˚C in 2013",
#                      "vent_13H_80_2014" = "80˚C in 2014",
#                      "individual" = "13H of vents individual temperatures",
#                      "vent_13H_30" = "30˚C",
#                      "vent_13H_55" = "55˚C",
#                      "vent_13H_80" = "80˚C")

baseline_lookup <- data.frame(short = c("all", "30", "55", "80", "vent_13H_30_2013", "vent_13H_30_2014", "vent_13H_55_2013", "vent_13H_55_2014", "vent_13H_80_2013", "vent_13H_80_2014", "individual", "vent_13H_30", "vent_13H_55", "vent_13H_80"),
                              label = c("13H of vents at all temperatures", "30˚C", "55˚C", "80˚C", "30˚C in 2013", "30˚C in 2014", "55˚C in 2013", "55˚C in 2014", "80˚C in 2013", "80˚C in 2014", "13H of vents individual temperatures", "30˚C", "55˚C", "80˚C" )) %>%
  tibble::deframe(.)


tpm_lookup <- c(axial_temp_13H = "13H samples as baseline",
                # axial_temp_12L = "12L samples as baseline",
                axial_temp_both = "13H and 12L samples as baseline",
                vent_temp_13H = "vent-temp-specific 13H samples")
log2_lookup <- data.frame(short = c("log2FC_TPM", "log2FoldChange", "log2FoldChange_MMSE"),
                          label = c("log2FC(mean TPM)", "DESeq log2FC(normalized counts)", "DESeq log2FC corrected" )) %>%
  tibble::deframe(.)

temp_lookup <- data.frame(temp = c(30, 55, 80, "all"),
                          label = c("13H of vents at 30˚C", "13H of vents at 55˚C", "13H of vents at 80˚C", "13H of vents at all temperatures")) %>%
  tibble::deframe(.)

year_lookup <- data.frame(year = c(2013, 2014),
                          label = c("2013", "2014")) %>%
  tibble::deframe(.)

vent_colors <- c(viridisLite::plasma(n = length(names(vent_lookup)), direction = 1)) %>%
  setNames(., names(vent_lookup))

colors_temp <- c(pals::ocean.phase(n = length(names(temp_lookup)))) %>%
  setNames(., names(temp_lookup))
colors_year <- c(pals::ocean.oxy(n = length(names(year_lookup)))) %>%
  setNames(., names(year_lookup))
colors_log2 <- c(pals::tol(n = length(names(log2_lookup)))) %>%
  setNames(., names(log2_lookup))

