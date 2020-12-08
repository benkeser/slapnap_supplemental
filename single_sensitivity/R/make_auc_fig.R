# create a figure comparing Rawi et al. (2019) and Hake and Pfeifer (2017) to SLAPNAP

# --------------------------------------------------------------------------
# Set up
# --------------------------------------------------------------------------
# load required libraries
library("dplyr")
library("tibble")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("here")

# create a tibble with the data
# note Hake and Pfeifer results are from Table S5, Rawi et al. results are from Table S1
rawi_results <- bind_rows(tibble(epitope = "V1V2", bnab = "PG9", AUC = 0.85, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "PG16", AUC = 0.79, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "2G12", AUC = 0.93, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "PGT145", AUC = 0.86, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "PGDM1400", AUC = 0.83, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "VRC26.08", AUC = 0.89, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "VRC26.25", AUC = 0.89, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "PGT128", AUC = 0.89, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "PGT121", AUC = 0.92, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "10-996", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "10-1074", AUC = 0.95, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "VRC38.01", AUC = 0.87, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "PGT135", AUC = 0.77, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "DH270.1", AUC = 0.92, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "DH270.5", AUC = 0.93, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "DH270.6", AUC = 0.93, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "VRC29.03", AUC = 0.82, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC01", AUC = 0.89, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC-PG04", AUC = 0.78, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "3BNC117", AUC = 0.88, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "NIH45-46", AUC = 0.8, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC03", AUC = 0.83, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC-CH31", AUC = 0.78, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "CH01", AUC = 0.77, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "HJ16", AUC = 0.67, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC07", AUC = 0.78, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "b12", AUC = 0.82, cil = NA, ciu = NA),
                           tibble(epitope = "Fusion peptide", bnab = "PGT151", AUC = 0.78, cil = NA, ciu = NA),
                           tibble(epitope = "Fusion peptide", bnab = "VRC34.01", AUC = 0.78, cil = NA, ciu = NA),
                           tibble(epitope = "Subunit interface", bnab = "8ANC195", AUC = 0.90, cil = NA, ciu = NA),
                           tibble(epitope = "Subunit interface", bnab = "35022", AUC = 0.63, cil = NA, ciu = NA),
                           tibble(epitope = "MPER", bnab = "4E10", AUC = 0.82, cil = NA, ciu = NA),
                           tibble(epitope = "MPER", bnab = "2F5", AUC = 0.97, cil = NA, ciu = NA)) %>%
    mutate(method = "Rawi et al. (2019)")
hake_results <- bind_rows(tibble(epitope = "V1V2", bnab = "PG9", AUC = 0.67, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "PG16", AUC = 0.71, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "2G12", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "PGT145", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "PGDM1400", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "VRC26.08", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V1V2", bnab = "VRC26.25", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "PGT128", AUC = 0.68, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "PGT121", AUC = 0.79, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "10-996", AUC = 0.84, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "10-1074", AUC = 0.81, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "VRC38.01", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "PGT135", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "DH270.1", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "DH270.5", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "DH270.6", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "V3", bnab = "VRC29.03", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC01", AUC = 0.71, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC-PG04", AUC = 0.69, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "3BNC117", AUC = 0.70, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "NIH45-46", AUC = 0.76, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC03", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC-CH31", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "CH01", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "HJ16", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "VRC07", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "CD4bs", bnab = "b12", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "Fusion peptide", bnab = "PGT151", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "Fusion peptide", bnab = "VRC34.01", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "Subunit interface", bnab = "8ANC195", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "Subunit interface", bnab = "35022", AUC = 0.65, cil = NA, ciu = NA),
                           tibble(epitope = "MPER", bnab = "4E10", AUC = NA, cil = NA, ciu = NA),
                           tibble(epitope = "MPER", bnab = "2F5", AUC = NA, cil = NA, ciu = NA)) %>%
    mutate(method = "Hake and Pfeifer (2017)")
slapnap_results_init <- load(here("R_output", "rslt_df.RData"))
# tack on epitope
slapnap_results <- slapnap_results_init %>%
    as_tibble() %>%
    mutate(
        epitope = case_when(
            bnab %in% c("2G12", "PG16", "PG9", "PGDM1400",
            "PGT145", "VRC26.08", "VRC26.25") ~ "V1V2",
            bnab %in% c("10-1074", "10-996", "DH270.1",
            "DH270.5", "DH270.6", "PGT121", "PGT128",
            "PGT135", "VRC29.03", "VRC38.01") ~ "V3",
            bnab %in% c("3BNC117", "b12", "CH01",
            "HJ16", "NIH45-46", "VRC-CH31", "VRC-PG04",
            "VRC01", "VRC03", "VRC07") ~ "CD4bs",
            bnab %in% c("PGT151", "VRC34.01") ~ "Fusion peptide",
            bnab %in% c("35O22", "8ANC195") ~ "Subunit interface",
            bnab %in% c("2F5", "4E10") ~ "MPER"
        ),
    method = "SLAPNAP")

compare_tib <- bind_rows(rawi_results, hake_results, slapnap_results)

# --------------------------------------------------------------------------
# Create the plot
# --------------------------------------------------------------------------
epitope_labs <- unique(compare_tib$epitope)[order(unique(compare_tib$epitope), decreasing = TRUE)]
cd4bs_plot <- compare_tib %>%
    filter(epitope == "CD4bs") %>%
    ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
               y = AUC, ymin = cil, ymax = ciu, shape = method,
               group = paste0(epitope, "_", method))) +
    geom_pointrange(position = position_dodge(width = 0.75, preserve = "total"),
                    size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ylim(c(0.3, 1)) +
    ggtitle("CD4bs") +
    xlab("") +
    labs(y = NULL) +
    guides(x = guide_axis(n.dodge = 2), shape = FALSE, y = "none") +
    theme(plot.title = element_text(hjust = 0.5))
fusion_plot <- compare_tib %>%
    filter(epitope == "Fusion peptide") %>%
    ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
               y = AUC, ymin = cil, ymax = ciu, shape = method,
               group = paste0(epitope, "_", method))) +
    geom_pointrange(position = position_dodge(width = 0.75, preserve = "total"),
                    size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ylim(c(0.3, 1)) +
    labs(y = NULL) +
    ggtitle("Fusion peptide") +
    xlab("") +
    guides(x = guide_axis(n.dodge = 2), shape = FALSE, y = "none") +
    theme(plot.title = element_text(hjust = 0.5))
mper_plot <- compare_tib %>%
    filter(epitope == "MPER") %>%
    ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
               y = AUC, ymin = cil, ymax = ciu, shape = method,
               group = paste0(epitope, "_", method))) +
    geom_pointrange(position = position_dodge(width = 0.75, preserve = "total"),
                    size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ylim(c(0.3, 1)) +
    labs(y = NULL) +
    ggtitle("MPER") +
    xlab("") +
    guides(x = guide_axis(n.dodge = 2), y = "none") +
    theme(plot.title = element_text(hjust = 0.5))
subunit_plot <- compare_tib %>%
    filter(epitope == "Subunit interface") %>%
    ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
               y = AUC, ymin = cil, ymax = ciu, shape = method,
               group = paste0(epitope, "_", method))) +
    geom_pointrange(position = position_dodge(width = 0.75, preserve = "total"),
                    size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ylim(c(0.3, 1)) +
    labs(y = NULL) +
    ggtitle("Subunit interface") +
    xlab("") +
    guides(x = guide_axis(n.dodge = 2), shape = FALSE, y = "none") +
    theme(plot.title = element_text(hjust = 0.5))
v1v2_plot <- compare_tib %>%
    filter(epitope == "V1V2") %>%
    ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
               y = AUC, ymin = cil, ymax = ciu, shape = method,
               group = paste0(epitope, "_", method))) +
    geom_pointrange(position = position_dodge(width = 0.75, preserve = "total"),
                    size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ylim(c(0.3, 1)) +
    ylab("CV-AUC") +
    ggtitle("V1V2") +
    xlab("") +
    guides(x = guide_axis(n.dodge = 2), shape = FALSE) +
    theme(plot.title = element_text(hjust = 0.5))
v3_plot <- compare_tib %>%
    filter(epitope == "V3") %>%
    ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
               y = AUC, ymin = cil, ymax = ciu, shape = method,
               group = paste0(epitope, "_", method))) +
    geom_pointrange(position = position_dodge(width = 0.75, preserve = "total"),
                    size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ylim(c(0.3, 1)) +
    labs(y = NULL) +
    ggtitle("V3") +
    xlab("") +
    guides(x = guide_axis(n.dodge = 2), y = "none", shape = FALSE) +
    theme(plot.title = element_text(hjust = 0.5))


compare_plot <- plot_grid(plot_grid(v1v2_plot, v3_plot, cd4bs_plot, fusion_plot, subunit_plot, mper_plot,
                                    nrow = 1, rel_widths = c(1, 1, 1, 0.475, 0.475, 1)),
                          ggplot() + ggtitle("bnAb") + guides(x = "none", y = "none") + theme(plot.title = element_text(hjust = 0.5, face = "plain")),
                          nrow = 2, rel_heights = c(1, 0.1))


ggsave(filename = here("fig", "auc_fig.tiff"),
       plot = compare_plot,
       width = 45, height = 20, units = "cm")
