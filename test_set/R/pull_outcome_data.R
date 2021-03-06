# pull IC-50 values from Mendoza et al. (2018) and Bar-On et al. (2018)

library("here")
library("tibble")
library("dplyr")
library("tidyr")
library("readr")
library("stringr")

# ------------------------------------------------------------------------------
# Read in the data
# ------------------------------------------------------------------------------

# read in the outcome data
mendoza_outcomes <- readr::read_csv(here("test_set/data", 
                                         "mendoza-etal_outcome_data.csv"))%>% 
    rename(id = participant, seq_id = virus_id, ic50 = IC50, ic80 = IC80) %>% 
    select(id, seq_id, time_point, mab, ic50, ic80) %>% 
    mutate(id = as.character(id))
bar_on_outcomes <- readr::read_csv(here("test_set/data", 
                                        "bar-on-etal_outcome_data.csv")) %>% 
    mutate(time_point = paste0("Week ", time_point))
all_outcomes <- dplyr::bind_rows(mendoza_outcomes %>% mutate(study = "Mendoza"),
                                 bar_on_outcomes %>% mutate(study = "Bar-On"))

# randomly select one per participant; first, widen
wide_outcomes <- all_outcomes %>% 
  mutate(
    ic50.imputed = ifelse(grepl(">", ic50, fixed = TRUE),
                          100, as.numeric(ic50)),
    ic80.imputed = ifelse(grepl(">", ic80, fixed = TRUE),
                          100, as.numeric(ic80))) %>% 
  select(-ic50, -ic80) %>% 
  pivot_wider(names_from = mab, values_from = c(ic50.imputed, ic80.imputed)) %>% 
  rename(nab_3BNC117.ic50.imputed = ic50.imputed_3BNC117, 
         nab_3BNC117.ic80.imputed = ic80.imputed_3BNC117,
         `nab_10-1074.ic50.imputed` = `ic50.imputed_10-1074`,
         `nab_10-1074.ic80.imputed` = `ic80.imputed_10-1074`) %>% 
  select(-ends_with("NA"))

# read in the viral sequence data
viral_sequence_data <- readr::read_csv(
    here("test_set/data", "nussenzweig_gb_env_aa_r2_validation.csv")
    )

# read in the hash table for matching the outcome data with genbank data
hash_table <- readr::read_csv(
    here("test_set/data", "genbank_lookup.csv")
)

# ------------------------------------------------------------------------------
# Compute combination IC50, IC80 and binary outcomes
# ------------------------------------------------------------------------------
# create the final outcomes
multsens_nab <- 1
wagh.additive.method <- function(x) 1 / sum(1 / x)
indiv_sens <- function(x, sens_thresh) sum(as.numeric(x < sens_thresh))
get_iip <- function(x) {
    iip.c <- 10
    iip.m <- log10(4) /(x[4] - x[3])
    iip.f.c <-(iip.c ^ iip.m) /((x[1] ^ iip.m) + (iip.c ^ iip.m))
    iip.f.c[iip.f.c >= 1] <- 1 - .Machine$double.neg.eps
    (-1) * log10(1 - iip.f.c)
}
pc.ic50 <- apply(wide_outcomes[, grepl("ic50", names(wide_outcomes))],
                 1, wagh.additive.method)
pc.ic80 <- apply(wide_outcomes[, grepl("ic80", names(wide_outcomes))],
                 1, wagh.additive.method)
log10.pc.ic50 <- log10(pc.ic50)
log10.pc.ic80 <- log10(pc.ic80)
iip.c <- 10
iip.m <- log10(4) /(log10.pc.ic80 - log10.pc.ic50)
iip.f.c <-(iip.c ^ iip.m) /((pc.ic50 ^ iip.m) + (iip.c ^ iip.m))
iip.f.c[iip.f.c >= 1] <- 1 - .Machine$double.neg.eps
iip <- (-1) * log10(1 - iip.f.c)
estsens_1 <- as.numeric(pc.ic50 < 1)
estsens_2 <- as.numeric(pc.ic50 < 2)
multsens_2 <- as.numeric(apply(wide_outcomes[, grepl("ic50", names(wide_outcomes))], 
                  1, indiv_sens, sens_thresh = 2) > multsens_nab)
multsens_1 <- as.numeric(apply(wide_outcomes[, grepl("ic50", names(wide_outcomes))], 
                               1, indiv_sens, sens_thresh = 1) > multsens_nab)
final_outcomes <- wide_outcomes %>% 
    mutate(pc.ic50 = pc.ic50, pc.ic80 = pc.ic80,
           log10.pc.ic50 = log10.pc.ic50, log10.pc.ic80 = log10.pc.ic80,
           iip = iip, estsens_2 = estsens_2, multsens_2 = multsens_2,
           estsens_1 = estsens_1, multsens_1 = multsens_1)

# ------------------------------------------------------------------------------
# Merge final outcome data with viral sequence data, using hash table
# ------------------------------------------------------------------------------
seq_ids <- unlist(lapply(str_split(final_outcomes$seq_id, pattern = "_", n = 2), 
                         function(x) x[2]))
final_outcomes$seq_id <- seq_ids
final_outcomes$virus_id <- apply(final_outcomes[, c("id", "time_point", "seq_id")],
      1, function(x) {
          paste0(x[1], "_",
                 switch((grepl("-", x[2], fixed = TRUE)) + 1,
                        switch((grepl("R", x[3])) + 1, 
                               paste0("W", gsub("Week ", "", x[2]), "_"),
                               "rebound_"),
                        "pre_"),
                 x[3])
      })
final_outcomes$virus_id <- gsub("W0", "D0", final_outcomes$virus_id)
outcomes_for_slapnap <- final_outcomes %>% 
    select(study, id, virus_id, log10.pc.ic50, log10.pc.ic80, starts_with("estsens"), 
           starts_with("multsens")) %>% 
    rename(ic50 = log10.pc.ic50, ic80 = log10.pc.ic80)


# to match virus_id from outcomes table with accession number, I need to do the following:
# for each outcome, look through and get the accession number; tack it on
# note that the Bar-On et al. ones may not correspond exactly
# note that zeros are okay, since I typically have more than one match per person
# also, some are close, but that's not good enough
outcome_accession_numbers <- vector("character", length = nrow(outcomes_for_slapnap))
id_virus_fuzzymatch <- str_split(outcomes_for_slapnap$virus_id, "_", n = 3)
num_matches <- vector("numeric", length = nrow(outcomes_for_slapnap))
for (i in 1:nrow(outcomes_for_slapnap)) {
  # first, check for an exact match
  this_indx <- which(grepl(
    paste0(outcomes_for_slapnap$virus_id[i], " "), 
    hash_table$definition
  ))
  if (length(this_indx) == 0) {
    this_id <- id_virus_fuzzymatch[[i]][1]
    this_virus <- paste0(id_virus_fuzzymatch[[i]][3], " ")
    this_baron_virus <- paste0(gsub("day_", "", id_virus_fuzzymatch[[i]][3]), "_")
    this_timepoint <- id_virus_fuzzymatch[[i]][2]
    id_match <- grepl(this_id, hash_table$definition)
    virus_match <- grepl(this_virus, hash_table$definition)
    timepoint_match <- grepl(this_timepoint, hash_table$definition) |
      grepl(tolower(this_timepoint), hash_table$definition) |
      grepl("rebound", hash_table$definition) |
      grepl("pre", hash_table$definition)
    this_indx <- which(id_match & virus_match & timepoint_match)
    if (length(this_indx) == 0) {
      virus_match <- grepl(this_baron_virus, hash_table$definition) 
      this_indx <- which(id_match & virus_match & timepoint_match)
      if (length(this_indx) > 1) {
        this_indx <- this_indx[1]
      }
    }
  }
  outcome_accession_numbers[i] <- switch((length(this_indx) == 0) + 1,
                                         hash_table$accnum[this_indx],
                                         NA)
  num_matches[i] <- length(this_indx)
}
check_tib <- outcomes_for_slapnap %>% 
  mutate(seq.id.lanl = outcome_accession_numbers) %>% 
  filter(!is.na(seq.id.lanl))

# check accession numbers 
used_accession_numbers <- check_tib %>% 
  select(seq.id.lanl, virus_id) %>% 
  rename(accnum = seq.id.lanl) %>% 
  left_join(hash_table, by = "accnum")
mean(mapply(grepl, pattern=used_accession_numbers$virus_id, 
       x=used_accession_numbers$definition))
length(which(!mapply(grepl, pattern=used_accession_numbers$virus_id, 
              x=used_accession_numbers$definition)))
hand_check_ids <- unique(unlist(
  lapply(str_split(used_accession_numbers$virus_id, "_", n = 2), 
         function(x) x[1]
  )
))
# hand-verify the rest; 
# the Mendoza et al. ones all check out
# the Bar-On et al. ones all check out
for (i in 1:length(hand_check_ids)) {
  used_accession_numbers[which(!mapply(grepl, pattern=used_accession_numbers$virus_id, 
                                       x=used_accession_numbers$definition)), ] %>%
    filter(grepl(hand_check_ids[i], virus_id)) %>% 
    print(width = Inf)  
}

# add them on
outcome_tib <- outcomes_for_slapnap %>% 
  mutate(seq.id.lanl = outcome_accession_numbers) %>% 
  select(-virus_id) %>% 
  filter(!is.na(seq.id.lanl))

# select one virus per donor
set.seed(4747)
unique_outcomes_tib <- outcome_tib %>% 
  group_by(study, id) %>% 
  sample_n(1) %>% 
  ungroup()

# merge together
analysis_dataset <- dplyr::left_join(
    viral_sequence_data %>% select(-ic50, -ic80, -estsens, -multsens, -iip),
    unique_outcomes_tib %>% select(-study, -id),
    by = "seq.id.lanl"
) %>% 
    select(seq.id.lanl, seq.id.catnap, ic50, ic80, 
           starts_with("estsens"), starts_with("multsens"), everything())
saveRDS(analysis_dataset, here("test_set/data", "analysis_dataset.rds"))
