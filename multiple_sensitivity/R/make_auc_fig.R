library(ggplot2)
library(here)

load(here("R_output", "rslt_df.RData"))
load(here("R_output", "n_sens_df.RData"))

rslt <- merge(rslt_df, n_sens_df)
rslt <- rslt[complete.cases(rslt), ]

# fix case
rslt$Combination <- factor(toupper(gsub("_", " + ", rslt$nab)))
rslt$sens <- as.character(rslt$sens)
rslt$sens[rslt$sens == "est"] <- "Estimated"
rslt$sens[rslt$sens == "mult"] <- "Multiple"
# re-level by number resistant sequences
non_duplicated_res <- rslt$n_res[!duplicated(rslt$Combination)]
non_duplicated_nab <- rslt$Combination[!duplicated(rslt$Combination)]
rslt$Combination <- factor(rslt$Combination, 
                           levels = as.character(non_duplicated_nab[order(non_duplicated_res)]))

tiff(here("fig", "auc_fig.tiff"), 
     width = 6, height = 5, units = "in", res = 300)
ggplot(rslt, aes(x = Combination, y = auc, group = sens, color = sens))+
  ylab("AUC") + labs(color = "Sensitivity") + 
  geom_point(position = position_dodge(width = 0.5)) + 
  scale_y_continuous(limits = c(0.5, 1)) + 
  theme_bw() + 
  scale_color_manual(values = c("gray75", "gray50")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(label = rslt$n_res, 
    y = 0.5, size = 3, 
    color = "black")
dev.off()