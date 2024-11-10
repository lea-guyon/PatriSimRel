library(cowplot)
library(ggplot2)
library(scales)
library(boot)
library(reshape2)
library(tidyr)
library(stringr)
library(dplyr)
library("ggpattern")

fun_median <- function(data, indices, metric) {
  d <- data[indices,]
  median = median(d[[metric]])
  return(median)
}

median_ratio <- function(data, indices, metric1, metric2) {
  d <- data[indices,]
  median_ratio = median(d[[metric1]])/median(d[[metric2]])
  return(median_ratio)
}

reshape_file <- function(path, model, mode, descent, sample=F) {
  setwd(path)
  
  if (mode == "village") {
    file = read.table("relatedness.txt", header = T)
    if (sample == T) {
      var = c("sampleRelFreqM", "sampleRelFreqF", 
              "sampleMeanRelM", "sampleMeanRelF",
              "sampleRelFirstM", "sampleRelFirstF")
    }
    else {
      var = c("relFreqM", "relFreqF",
              "meanRelM", "meanRelF",
              "relFirstM", "relFirstF")
    }
  }
  
  if (mode == "village_father") {
    file = read.table("relatedness_father.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "village_mother") {
    file = read.table("relatedness_mother.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "group") {
    file = read.table("relatedness_group.txt", header = T)
    if (sample == T) {
      var = c("sampleRelGroupFreqM", "sampleRelGroupFreqF", 
              "sampleMeanGroupRelM", "sampleMeanGroupRelF", 
              "sampleRelGroupFirstM", "sampleRelGroupFirstF")
    }
    else {
      var = c("relGroupFreqM", "relGroupFreqF", 
              "meanGroupRelM", "meanGroupRelF", 
              "relGroupFirstM", "relGroupFirstF")
    }
  }
  
  if (mode == "patrilocal_group") {
    file = read.table("relatedness_patriline.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "matrilocal_group") {
    file = read.table("relatedness_matriline.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "local_group") {
    file = read.table("relatedness_local.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "local_group_father") {
    file = read.table("relatedness_patriline_father.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "local_group_mother") {
    file = read.table("relatedness_matriline_mother.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  reshaped_file <- pivot_longer(file, 
                                cols = var,
                                names_to = "variable",
                                values_to = "value")
  reshaped_file$category = str_sub(reshaped_file$variable, end=-2)
  reshaped_file$kinship = rep(model, length(reshaped_file$variable))
  reshaped_file$kinvar = paste(reshaped_file$kinship, reshaped_file$variable)
  reshaped_file$descent = rep(descent, length(reshaped_file$variable))
  return(reshaped_file)
}

bootstrap <- function(data, metric=NaN, metric1=NaN, metric2=NaN, ratio) {
  if (ratio == F) {
    b <- boot(data=data, statistic = fun_median, R = 10000, metric = metric)
  }
  else {
    b <- boot(data=data, statistic=median_ratio, R=10000, metric1=metric1, metric2=metric2)
  }
  ci <- boot.ci(b, conf = 0.95, type="bca")
  median <- b$t0
  ICinf <- ci$bca[4]
  ICsup <- ci$bca[5]
  return(list(median, ICinf, ICsup))
}

CI_metrics <- function(data, metrics, col, kinship) {
  m = metrics[1]
  b = bootstrap(data, metric=m, ratio=F)
  median = b[[1]]
  ICinf = b[[2]]
  ICsup = b[[3]]
  df_metric <- data.frame("metric" = m, "median" = median, "ICinf" = ICinf, "ICsup" = ICsup, "col" = col[1], "kinship" = kinship)
  for (i in seq(2,length(metrics))) {
    m = metrics[i]
    b = bootstrap(data, metric=m, ratio=F)
    median = b[[1]]
    ICinf = b[[2]]
    ICsup = b[[3]]
    df_metric <- merge(df_metric, 
                       data.frame("metric" = m, "median" = median, "ICinf" = ICinf, "ICsup" = ICsup, "col" = col[i], "kinship" = kinship), 
                       all = T)
  }
  return(df_metric)
}

mean_hg <- function(path, kinship, mode, descent) {
  setwd(path)
  
  df = read.table(paste0("hg_freq_", mode, ".txt"), header=T)
  
  df2 <- df %>%
    group_by(Replicat, Generation, Chromosome, Village) %>%
    summarize(Hg_diversity = nIndividuals[1] / (nIndividuals[1] - 1) * (1 - sum(Frequency^2, na.rm=T))) %>%
    ungroup()
  
  mean_df <- df2 %>%
    group_by(Replicat, Generation, Chromosome) %>%
    summarize(mean_hg_div = mean(Hg_diversity, na.rm=T)) %>%
    ungroup()
  
  if (kinship == "SM2B_SM" | kinship == "SP2B_SP") {
    mean_df = mean_df[mean_df$Generation == 20200,]
  }
  else {
    mean_df = mean_df[mean_df$Generation == 20100,]
  }
  mean_df$kinship = rep(kinship, length(mean_df$Replicat))
  mean_df$descent = rep(descent, length(mean_df$Replicat))
  mean_df$descent = factor(mean_df$descent, levels = c("matrilineal", "bilateral", "patrilineal"))
  return(mean_df)
}

rearrange_reshape <- function(reshape2) {
  relabel.variable = as_labeller(c(meanRelM = "male", 
                                   meanRelF = "female", 
                                   relFirstM = "male", relFirstF = "female",
                                   relFreqM = "male", relFreqF = "female"))
  category.label =  c("mean rel", "freq 1st degree rel", "freq 1st, 2nd and 3rd degree rel")
  names(category.label) = c("meanRel", "relFirst", "relFreq")
  reshape_bis = data.frame("variable"=reshape2$variable,"value"=reshape2$value, 
                           "category"=reshape2$category, "kinship"=reshape2$kinship,
                           "descent" = reshape2$descent)
  
  reshape2 <- reshape2 %>% 
    as_tibble() %>%
    rowwise() %>% 
    mutate(sex = factor(relabel.variable(variable))) %>%
    arrange(value) %>%
    mutate(kinship = factor(kinship, levels = c("SM", "LM", "SM_M", "LM_M", 
                                                "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                                "SP", "LP", "SP_F", "LP_F")))
  
  reshape2 = reshape2[with(reshape2, order(sex, kinship)),]
  
  data = reshape2[reshape2$category == "meanRel",]
  data = data[data$Generation == 20100 | data$Generation == 20200,]
  data = data[!((data$kinship == "SM2B_SM" & data$Generation == 20100) | (data$kinship == "SP2B_SP" & data$Generation == 20100)),]
  data <- data %>%
    group_by(descent, kinship, sex, Replicat) %>%
    summarise(value = mean(value))
  data$descent = factor(data$descent, levels = c("matrilineal", "bilateral", "patrilineal"))
  data = data[with(data, order(sex, kinship)),]
  return(data)
}

rearrange_pitab <- function(tab) {
  tab = data.frame("Pi" = tab$Pi, "kinship" = as.factor(tab$kinship), "sex" = tab$sex, 
                   "descent" = factor(tab$descent, levels = c("matrilineal", 
                                                              "bilateral",
                                                           "patrilineal")))
  tab = na.omit(tab)
  tab <- tab %>%
    arrange(Pi) %>%
    mutate(kinship = factor(kinship, levels=c("SM", "LM", "SM_M", "LM_M",
                                              "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                              "SP", "LP", "SP_F", "LP_F")))
  
  tab = tab[with(tab, order(sex, kinship)),]
  return(tab)
}

significance <- function(data_stat, value, param) {
  if (param == "sex") {
    data_stat$var = factor(paste(data_stat$kinship, data_stat[[param]]),
                           levels = c("SM female", "LM female", "SM_M female", "LM_M female",
                                      "B_SM female", "B_LM female", "B_A female",
                                      "B_SP female", "B_LP female", "SP female",
                                      "LP female", "SP_F female", "LP_F female", "SM male",
                                      "LM male", "SM_M male", "LM_M male",
                                      "B_SM male", "B_LM male", "B_A male", "B_SP male", "B_LP male",
                                      "SP male", "LP male", "SP_F male", "LP_F male"))
  }
  else {
    data_stat$var = factor(paste(data_stat$kinship, data_stat[[param]]),
                           levels = c("SM M", "LM M", "SM_M M", "LM_M M",
                                      "B_SM M", "B_LM M", "B_A M",
                                      "B_SP M", "B_LP M", "SP M",
                                      "LP M", "SP_F M", "LP_F M", "SM Y",
                                      "LM Y", "SM_M Y", "LM_M Y",
                                      "B_SM Y", "B_LM Y", "B_A Y", "B_SP Y", "B_LP Y",
                                      "SP Y", "LP Y", "SP_F Y", "LP_F Y"))
  }
  
  pval_matrix = pairwise.wilcox.test(data_stat[[value]], data_stat$var, p.adjust.method="bonferroni")$p.value
  
  # Define the significance levels
  significance_matrix <- matrix("ns", nrow=nrow(pval_matrix), ncol=ncol(pval_matrix)) # Start with "not significant"
  
  # Apply thresholds
  significance_matrix[is.na(pval_matrix)] <- NA
  significance_matrix[pval_matrix < 0.05] <- "*"
  significance_matrix[pval_matrix < 0.01] <- "**"
  significance_matrix[pval_matrix < 0.001] <- "***"
  significance_matrix[pval_matrix < 0.0001] <- "****"
  significance_matrix = rbind("SM_M female"=rep(NA, 25), significance_matrix)
  significance_matrix = cbind(significance_matrix, "L_PF male" = rep(NA, 26))
  
  if (param == "sex") {
    colnames(significance_matrix) = c("SM female", "LM female", "SM_M female", "LM_M female", 
                                      "B_SM female", "B_LM female", "B_A female", 
                                      "B_SP female", "B_LP female", "SP female", 
                                      "LP female", "SP_F female", "LP_F female", "SM male", 
                                      "LM male", "SM_M male", "LM_M male", 
                                      "B_SM male", "B_LM male", "B_A male", "B_SP male", "B_LP male", 
                                      "SP male", "LP male", "SP_F male", "LP_F male")
  }
  else {
    colnames(significance_matrix) = c("SM M", "LM M", "SM_M M", "LM_M M",
                                      "B_SM M", "B_LM M", "B_A M",
                                      "B_SP M", "B_LP M", "SP M",
                                      "LP M", "SP_F M", "LP_F M", "SM Y",
                                      "LM Y", "SM_M Y", "LM_M Y",
                                      "B_SM Y", "B_LM Y", "B_A Y", "B_SP Y", "B_LP Y",
                                      "SP Y", "LP Y", "SP_F Y", "LP_F Y")
  }
  rownames(significance_matrix) = colnames(significance_matrix)
  
  df <- melt(significance_matrix)
  colnames(df) <- c("x", "y", "value")
  return(df)
}

### set working directory
working_dir = "~"
setwd(working_dir)

########## Relatedness ###########
dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
path = paste(dir, "strict_patrilineal_villages", sep = "")
reshape_file_p1 = reshape_file(path, "SP", "village", "patrilineal")
reshape_file_father_1 = reshape_file(path, "SP_F", "village_father", "patrilineal")
reshape_group_p1 = reshape_file(path, "SP", "group", "patrilineal")
reshape_sample_1 = reshape_file(path, "SP", "village", "patrilineal", T)
reshape_group_sample_1 = reshape_file(path, "SP", "group", "patrilineal", T)
reshape_patriline_1 = reshape_file(path, "SP", "patrilocal_group", "patrilineal")
reshape_patriline_father_1 = reshape_file(path, "SP_F", "local_group_father", "patrilineal")
path = paste(dir, "strict_matrilineal_villages", sep = "")
reshape_file_m1 = reshape_file(path, "SM", "village", "matrilineal")
reshape_file_mother_1 = reshape_file(path, "SM_M", "village_mother", "matrilineal")
reshape_group_m1 = reshape_file(path, "SM", "group", "matrilineal")
reshape_matriline_1 = reshape_file(path, "SM", "matrilocal_group", "matrilineal")
reshape_matriline_mother_1 = reshape_file(path, "SM_M", "local_group_mother", "matrilineal")

dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
path = paste(dir, "loose_patrilineal_villages", sep = "")
reshape_file_p2 = reshape_file(path, "LP", "village", "patrilineal")
reshape_file_father_2 = reshape_file(path, "LP_F", "village_father", "patrilineal")
reshape_group_p2 = reshape_file(path, "LP", "group", "patrilineal")
reshape_sample_2 = reshape_file(path, "LP", "village", "patrilineal", T)
reshape_group_sample_2 = reshape_file(path, "LP", "group", "patrilineal", T)
reshape_patriline_2 = reshape_file(path, "LP", "patrilocal_group", "patrilineal")
reshape_patriline_father_2 = reshape_file(path, "LP_F", "local_group_father", "patrilineal")
path = paste(dir, "loose_matrilineal_villages", sep = "")
reshape_file_m2 = reshape_file(path, "LM", "village", "matrilineal")
reshape_file_mother_2 = reshape_file(path, "LM_M", "village_mother", "matrilineal")
reshape_group_m2 = reshape_file(path, "LM", "group", "matrilineal")
reshape_matriline_2 = reshape_file(path, "LM", "matrilocal_group", "matrilineal")
reshape_matriline_mother_2 = reshape_file(path, "LM_M", "local_group_mother", "matrilineal")

dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
path = paste(dir, "strict_patrilocal_villages", sep = "")
reshape_file_p3 = reshape_file(path, "B_SP", "village", "bilateral")
reshape_sample_3 = reshape_file(path, "B_SP", "village", "bilateral", T)
reshape_patriline_3 = reshape_file(path, "B_SP", "patrilocal_group", "bilateral")

path = paste(dir, "strict_matrilocal_villages", sep = "")
reshape_file_m3 = reshape_file(path, "B_SM", "village", "bilateral")
reshape_matriline_3 = reshape_file(path, "B_SM", "matrilocal_group", "bilateral")

path = paste(dir, "ambilocal_villages", sep = "")
reshape_file_a3 = reshape_file(path, "B_A", "village", "bilateral")
reshape_local_3 = reshape_file(path, "B_A", "local_group", "bilateral")

path = paste(dir, "loose_patrilocal_villages", sep = "")
reshape_file_p4 = reshape_file(path, "B_LP", "village", "bilateral")
reshape_sample_4 = reshape_file(path, "B_LP", "village", "bilateral", T)
reshape_patriline_4 = reshape_file(path, "B_LP", "patrilocal_group", "bilateral")

path = paste(dir, "loose_matrilocal_villages", sep = "")
reshape_file_m4 = reshape_file(path, "B_LM", "village", "bilateral")
reshape_matriline_4 = reshape_file(path, "B_LM", "matrilocal_group", "bilateral")

names = c("SM: Strict matrilineal descent, burial in wife's cemetery", 
          "LM: Loose matrilineal descent, burial in wife's cemetery",
          "SM_M: Strict matrilineal descent, burial in mother's cemetery", 
          "LM_M: Loose matrilineal descent, burial in mother's cemetery",
          "B_SM: Bilateral descent and strict matrilocal residence",
          "B_LM: Bilateral descent and loose matrilocal residence",
          "B_A: Bilateral descent and ambilocal residence",
          "B_SP: Bilateral descent and strict patrilocal residence",
          "B_LP: Bilateral descent and loose patrilocal residence",
          "SP: Strict patrilineal descent, burial in husband's cemetery", 
          "LP: Loose patrilineal descent, burial in husband's cemetery",
          "SP_F: Strict patrilineal descent, burial in father's cemetery", 
          "LP_F: Loose patrilineal descent, burial in father's cemetery")

compress_trans <- trans_new(
  name = "compress", 
  transform = function(x) { ifelse(x == 0, 0, sign(x) * log1p(abs(x)^2)) },   # Exacerbate by squaring values
  inverse = function(x) { sign(x) * (exp(abs(x))^(1/2) - 1) }               # Inverse transformation
)
##################################################################
########################### MAIN #################################
##################################################################

########################## VILLAGE ###############################

reshape = merge(reshape_file_p1, reshape_file_p2, all = T)
reshape2 = merge(reshape, reshape_file_p3, all = T)
reshape2 = merge(reshape2, reshape_file_p4, all = T)
reshape2 = merge(reshape2, reshape_file_m1, all = T)
reshape2 = merge(reshape2, reshape_file_m2, all = T)
reshape2 = merge(reshape2, reshape_file_m3, all = T)
reshape2 = merge(reshape2, reshape_file_m4, all = T)
reshape2 = merge(reshape2, reshape_file_a3, all = T)
reshape2 <- reshape2[, -seq(3,8)]
reshape2 = merge(reshape2, reshape_file_father_1, all = T)
reshape2 = merge(reshape2, reshape_file_father_2, all = T)
reshape2 = merge(reshape2, reshape_file_mother_1, all = T)
reshape2 = merge(reshape2, reshape_file_mother_2, all = T)

# Rearrange reshape2
data = rearrange_reshape(reshape2)
data$sex=factor(data$sex)
data_mirror <- data %>%
  mutate(value = ifelse(sex == "male", -value, value))

p1 <- ggplot(data_mirror) +
  facet_wrap(~descent, ncol = 3, scales='free_x') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(kinship, value, col=kinship, group=sex), size=2, alpha=0.05, 
  position = position_identity()) +
  geom_boxplot_pattern(aes(kinship, value, fill = kinship, col = kinship, 
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship, alpha = sex),
                           lwd=0.75, outlier.shape=NA, 
                       position = position_identity()) +
  scale_alpha_manual(values = c(0.2, 0.2)) +
  scale_fill_manual(labels = names,
                    values = c("brown3", "#e18585",
                               "#f0c2c2", "#faebeb",
                               "orchid3", "plum", 
                               "darkcyan",
                               "seagreen", "darkseagreen",
                               "darkgoldenrod1", "#FFD56F",
                               "#FFC55E", "bisque1")) +
  scale_color_manual(labels = names,
                     values = c("brown", "#dc7070",
                                "brown", "#f5d6d6",
                                "plum4", "orchid4",
                                "#027d7d",
                                "seagreen", "darkseagreen",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = names,
                       values= c("crosshatch", "none", "circle", "none", 
                                 "crosshatch", "none", "none", "crosshatch", "none", 
                                 "crosshatch", "none", "circle", "none")) +
  scale_pattern_density_manual(labels = names,
                       values= c(0.4, 0, 0.5, 0, 0.3, 0, 0, 0.3, 0, 0.4, 0, 0.5, 0)) +
  scale_pattern_fill_manual(labels = names,
                            values=c("#faebeb", "#faebeb",
                                     "brown3", "#faebeb",
                                     "#d5caf1", "orchid4",
                                     "#98bbbb",
                                     "#cdeede", "darkseagreen", 
                                     "bisque1", 
                                     "bisque1", "darkgoldenrod1", "bisque1")) +
  scale_y_continuous(trans=compress_trans,
                     breaks = c(-0.005, -0.004, -0.003, 0, 0.003, 0.004, 0.005),
                     labels = c("0.005", "0.004", "0.003", "0", "0.003", "0.004", "0.005")) +
  xlab("Kinship system") +
  ylab("Mean relatedness (village's cemetery)") +
  labs(fill = "Kinship system", col = "Kinship system",
       pattern_fill = "Kinship system", pattern = "Kinship system",
       pattern_density= "Kinship system") +
  guides(fill=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE), 
                           override.aes = list(alpha = 0.2)),
         col=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE)),
         pattern=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE)),
         pattern_fill=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE), 
                                   override.aes = list(pattern_density=0.2)),
         alpha = 'none') +
  theme_bw()+
  theme(text = element_text(size = 12), legend.key.size = unit(1.2, 'cm'),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

df = significance(data, "value", "sex")

h1 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  geom_vline(xintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("#cc1a1a", "navyblue"), na.value="white") +
  ggtitle("mean relatedness (village's cemetery)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

####################### LOCAL ###########################
reshape = merge(reshape_patriline_1, reshape_patriline_2, all = T)
reshape2 = merge(reshape, reshape_patriline_3, all = T)
reshape2 = merge(reshape2, reshape_patriline_4, all = T)
reshape2 = merge(reshape2, reshape_local_3, all = T)
reshape2 = merge(reshape2, reshape_matriline_1, all = T)
reshape2 = merge(reshape2, reshape_matriline_2, all = T)
reshape2 = merge(reshape2, reshape_matriline_3, all = T)
reshape2 = merge(reshape2, reshape_matriline_4, all = T)
reshape2 = merge(reshape2, reshape_patriline_father_1, all = T)
reshape2 = merge(reshape2, reshape_patriline_father_2, all = T)
reshape2 = merge(reshape2, reshape_matriline_mother_1, all = T)
reshape2 = merge(reshape2, reshape_matriline_mother_2, all = T)

# Rearrange reshape2
data = rearrange_reshape(reshape2)
data$sex=factor(data$sex)
data_mirror <- data %>%
  mutate(value = ifelse(sex == "male", -value, value))

p2 <- ggplot(data_mirror) +
  facet_wrap(~descent, ncol = 3, scales='free_x') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(kinship, value, col=kinship, group=sex), size=2, alpha=0.05, 
  position = position_identity()) +
  geom_boxplot_pattern(aes(kinship, value, fill = kinship, col = kinship, 
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship, alpha = sex),
                           lwd=0.75, outlier.shape=NA, 
                       position = position_identity()) +
  scale_alpha_manual(values = c(0.2, 0.2)) +
  scale_fill_manual(labels = names,
                    values = c("brown3", "#e18585",
                               "#f0c2c2", "#faebeb",
                               "orchid3", "plum", 
                               "darkcyan",
                               "seagreen", "darkseagreen",
                               "darkgoldenrod1", "#FFD56F",
                               "#FFC55E", "bisque1")) +
  scale_color_manual(labels = names,
                     values = c("brown", "#dc7070",
                                "brown", "#f5d6d6",
                                "plum4", "orchid4",
                                "#027d7d",
                                "seagreen", "darkseagreen",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = names,
                       values= c("crosshatch", "none", "circle", "none", 
                                 "crosshatch", "none", "none", "crosshatch", "none", 
                                 "crosshatch", "none", "circle", "none")) +
  scale_pattern_density_manual(labels = names,
                               values= c(0.4, 0, 0.5, 0, 0.3, 0, 0, 0.3, 0, 0.4, 0, 0.5, 0)) +
  scale_pattern_fill_manual(labels = names,
                            values=c("#faebeb", "#faebeb",
                                     "brown3", "#faebeb",
                                     "#d5caf1", "orchid4",
                                     "#98bbbb",
                                     "#cdeede", "darkseagreen", 
                                     "bisque1", 
                                     "bisque1", "darkgoldenrod1", "bisque1")) +
  scale_y_continuous(labels = abs) +
  xlab("Kinship systems") +
  ylab("Mean relatedness (local group's cemetery)") +
  theme_bw()+
  theme(text = element_text(size = 12), axis.text.x = element_blank(),
  axis.ticks.x = element_blank())

df = significance(data, "value", "sex")

h2 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  geom_vline(xintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("pink", "pink2", "#d87272", "#cc1a1a", "navyblue"), na.value="white") +
  ggtitle("mean relatedness (local group's cemetery)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

######################## Nucleotide diversity #########################
##### VILLAGE

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY1 = read.table(path, header=TRUE)
tabY1$kinship = rep("SP", length(tabY1$Gen))
tabY1 = tabY1[tabY1$Gen == 100,]
tabY1$sex = rep("male", length(tabY1$Gen))
tabY1$descent = rep("patrilineal", length(tabY1$Gen))
tabY1$Pi = tabY1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM1 = read.table(path, header=TRUE)
tabM1$kinship = rep("SP", length(tabM1$Gen))
tabM1 = tabM1[tabM1$Gen == 100,]
tabM1$sex = rep("female", length(tabM1$Gen))
tabM1$descent = rep("patrilineal", length(tabM1$Gen))
tabM1$Pi = tabM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/village_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF1 = read.table(path, header=TRUE)
tabYF1$kinship = rep("SP_F", length(tabYF1$Gen))
tabYF1 = tabYF1[tabYF1$Gen == 100,]
tabYF1$sex = rep("male", length(tabYF1$Gen))
tabYF1$descent = rep("patrilineal", length(tabYF1$Gen))
tabYF1$Pi = tabYF1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF1 = read.table(path, header=TRUE)
tabMF1$kinship = rep("SP_F", length(tabMF1$Gen))
tabMF1 = tabMF1[tabMF1$Gen == 100,]
tabMF1$sex = rep("female", length(tabMF1$Gen))
tabMF1$descent = rep("patrilineal", length(tabMF1$Gen))
tabMF1$Pi = tabMF1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_matrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM1 = read.table(path, header=TRUE)
tabYM1$kinship = rep("SM", length(tabYM1$Gen))
tabYM1 = tabYM1[tabYM1$Gen == 100,]
tabYM1$sex = rep("male", length(tabYM1$Gen))
tabYM1$descent = rep("matrilineal", length(tabYM1$Gen))
tabYM1$Pi = tabYM1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM1 = read.table(path, header=TRUE)
tabMM1$kinship = rep("SM", length(tabMM1$Gen))
tabMM1 = tabMM1[tabMM1$Gen == 100,]
tabMM1$sex = rep("female", length(tabMM1$Gen))
tabMM1$descent = rep("matrilineal", length(tabMM1$Gen))
tabMM1$Pi = tabMM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_matrilineal_villages/village_mother/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYMM1 = read.table(path, header=TRUE)
tabYMM1$kinship = rep("SM_M", length(tabYMM1$Gen))
tabYMM1 = tabYMM1[tabYMM1$Gen == 100,]
tabYMM1$sex = rep("male", length(tabYMM1$Gen))
tabYMM1$descent = rep("matrilineal", length(tabYMM1$Gen))
tabYMM1$Pi = tabYMM1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMMM1 = read.table(path, header=TRUE)
tabMMM1$kinship = rep("SM_M", length(tabMMM1$Gen))
tabMMM1 = tabMMM1[tabMMM1$Gen == 100,]
tabMMM1$sex = rep("female", length(tabMMM1$Gen))
tabMMM1$descent = rep("matrilineal", length(tabMMM1$Gen))
tabMMM1$Pi = tabMMM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY2 = read.table(path, header=TRUE)
tabY2$kinship = rep("LP", length(tabY2$Gen))
tabY2 = tabY2[tabY2$Gen == 100,]
tabY2$sex = rep("male", length(tabY2$Gen))
tabY2$descent = rep("patrilineal", length(tabY2$Gen))
tabY2$Pi = tabY2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM2 = read.table(path, header=TRUE)
tabM2$kinship = rep("LP", length(tabM2$Gen))
tabM2 = tabM2[tabM2$Gen == 100,]
tabM2$sex = rep("female", length(tabM2$Gen))
tabM2$descent = rep("patrilineal", length(tabM2$Gen))
tabM2$Pi = tabM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/village_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF2 = read.table(path, header=TRUE)
tabYF2$kinship = rep("LP_F", length(tabYF2$Gen))
tabYF2 = tabYF2[tabYF2$Gen == 100,]
tabYF2$sex = rep("male", length(tabYF2$Gen))
tabYF2$descent = rep("patrilineal", length(tabYF2$Gen))
tabYF2$Pi = tabYF2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF2 = read.table(path, header=TRUE)
tabMF2$kinship = rep("LP_F", length(tabMF2$Gen))
tabMF2 = tabMF2[tabMF2$Gen == 100,]
tabMF2$sex = rep("female", length(tabMF2$Gen))
tabMF2$descent = rep("patrilineal", length(tabMF2$Gen))
tabMF2$Pi = tabMF2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_matrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM2 = read.table(path, header=TRUE)
tabYM2$kinship = rep("LM", length(tabYM2$Gen))
tabYM2 = tabYM2[tabYM2$Gen == 100,]
tabYM2$sex = rep("male", length(tabYM2$Gen))
tabYM2$descent = rep("matrilineal", length(tabYM2$Gen))
tabYM2$Pi = tabYM2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM2 = read.table(path, header=TRUE)
tabMM2$kinship = rep("LM", length(tabMM2$Gen))
tabMM2 = tabMM2[tabMM2$Gen == 100,]
tabMM2$sex = rep("female", length(tabMM2$Gen))
tabMM2$descent = rep("matrilineal", length(tabMM2$Gen))
tabMM2$Pi = tabMM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_matrilineal_villages/village_mother/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYMM2 = read.table(path, header=TRUE)
tabYMM2$kinship = rep("LM_M", length(tabYMM2$Gen))
tabYMM2 = tabYMM2[tabYMM2$Gen == 100,]
tabYMM2$sex = rep("male", length(tabYMM2$Gen))
tabYMM2$descent = rep("matrilineal", length(tabYMM2$Gen))
tabYMM2$Pi = tabYMM2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMMM2 = read.table(path, header=TRUE)
tabMMM2$kinship = rep("LM_M", length(tabMMM2$Gen))
tabMMM2 = tabMMM2[tabMMM2$Gen == 100,]
tabMMM2$sex = rep("female", length(tabMMM2$Gen))
tabMMM2$descent = rep("matrilineal", length(tabMMM2$Gen))
tabMMM2$Pi = tabMMM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/strict_patrilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY3 = read.table(path, header=TRUE)
tabY3$kinship = rep("B_SP", length(tabY3$Gen))
tabY3 = tabY3[tabY3$Gen == 100,]
tabY3$sex = rep("male", length(tabY3$Gen))
tabY3$descent = rep("bilateral", length(tabY3$Gen))
tabY3$Pi = tabY3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM3 = read.table(path, header=TRUE)
tabM3$kinship = rep("B_SP", length(tabM3$Gen))
tabM3 = tabM3[tabM3$Gen == 100,]
tabM3$sex = rep("female", length(tabM3$Gen))
tabM3$descent = rep("bilateral", length(tabM3$Gen))
tabM3$Pi = tabM3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/ambilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYA3 = read.table(path, header=TRUE)
tabYA3$kinship = rep("B_A", length(tabYA3$Gen))
tabYA3 = tabYA3[tabYA3$Gen == 100,]
tabYA3$sex = rep("male", length(tabYA3$Gen))
tabYA3$descent = rep("bilateral", length(tabYA3$Gen))
tabYA3$Pi = tabYA3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMA3 = read.table(path, header=TRUE)
tabMA3$kinship = rep("B_A", length(tabMA3$Gen))
tabMA3 = tabMA3[tabMA3$Gen == 100,]
tabMA3$sex = rep("female", length(tabMA3$Gen))
tabMA3$descent = rep("bilateral", length(tabMA3$Gen))
tabMA3$Pi = tabMA3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/strict_matrilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM3 = read.table(path, header=TRUE)
tabYM3$kinship = rep("B_SM", length(tabYM3$Gen))
tabYM3 = tabYM3[tabYM3$Gen == 100,]
tabYM3$sex = rep("male", length(tabYM3$Gen))
tabYM3$descent = rep("bilateral", length(tabYM3$Gen))
tabYM3$Pi = tabYM3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM3 = read.table(path, header=TRUE)
tabMM3$kinship = rep("B_SM", length(tabMM3$Gen))
tabMM3 = tabMM3[tabMM3$Gen == 100,]
tabMM3$sex = rep("female", length(tabMM3$Gen))
tabMM3$descent = rep("bilateral", length(tabMM3$Gen))
tabMM3$Pi = tabMM3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/loose_patrilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY4 = read.table(path, header=TRUE)
tabY4$kinship = rep("B_LP", length(tabY4$Gen))
tabY4 = tabY4[tabY4$Gen == 100,]
tabY4$sex = rep("male", length(tabY4$Gen))
tabY4$descent = rep("bilateral", length(tabY4$Gen))
tabY4$Pi = tabY4$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM4 = read.table(path, header=TRUE)
tabM4$kinship = rep("B_LP", length(tabM4$Gen))
tabM4 = tabM4[tabM4$Gen == 100,]
tabM4$sex = rep("female", length(tabM4$Gen))
tabM4$descent = rep("bilateral", length(tabM4$Gen))
tabM4$Pi = tabM4$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/loose_matrilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM4 = read.table(path, header=TRUE)
tabYM4$kinship = rep("B_LM", length(tabYM4$Gen))
tabYM4 = tabYM4[tabYM4$Gen == 100,]
tabYM4$sex = rep("male", length(tabYM4$Gen))
tabYM4$descent = rep("bilateral", length(tabYM4$Gen))
tabYM4$Pi = tabYM4$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM4 = read.table(path, header=TRUE)
tabMM4$kinship = rep("B_LM", length(tabMM4$Gen))
tabMM4 = tabMM4[tabMM4$Gen == 100,]
tabMM4$sex = rep("female", length(tabMM4$Gen))
tabMM4$descent = rep("bilateral", length(tabMM4$Gen))
tabMM4$Pi = tabMM4$Pi/(2*5.5e-7)

tabY = merge(tabY1, tabY2, all = T)
tabY = merge(tabY, tabY3, all = T)
tabY = merge(tabY, tabY4, all = T)
tabY = merge(tabY, tabYF1, all = T)
tabY = merge(tabY, tabYMM1, all = T)
tabY = merge(tabY, tabYF2, all = T)
tabY = merge(tabY, tabYMM2, all = T)
tabY = merge(tabY, tabYM1, all = T)
tabY = merge(tabY, tabYM2, all = T)
tabY = merge(tabY, tabYM3, all = T)
tabY = merge(tabY, tabYM4, all = T)
tabY = merge(tabY, tabYA3, all = T)

tabM = merge(tabM1, tabM2, all=T)
tabM = merge(tabM, tabM3, all=T)
tabM = merge(tabM, tabM4, all=T)
tabM = merge(tabM, tabMF1, all=T)
tabM = merge(tabM, tabMF2, all=T)
tabM = merge(tabM, tabMMM1, all=T)
tabM = merge(tabM, tabMMM2, all=T)
tabM = merge(tabM, tabMM1, all = T)
tabM = merge(tabM, tabMM2, all = T)
tabM = merge(tabM, tabMM3, all = T)
tabM = merge(tabM, tabMM4, all = T)
tabM = merge(tabM, tabMA3, all = T)

alltab = merge(tabY, tabM, all=T)
alltab$sex = factor(alltab$sex)

alltab = rearrange_pitab(alltab)
data_mirror <- alltab %>%
  mutate(Pi = ifelse(sex == "male", -Pi, Pi))

p3 <- ggplot(data_mirror) +
  facet_wrap(~descent, ncol=3, scales='free_x') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(kinship, Pi, col=kinship), size=2, alpha=0.05, 
            position = position_identity()) +
  geom_boxplot_pattern(aes(kinship, Pi, fill = kinship, col = kinship,
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship, alpha = sex), 
                      lwd=0.75, outlier.shape=NA, position = position_identity()) +
  scale_alpha_manual(values = c(0.2, 0.2)) +
  scale_fill_manual(labels = names,
                    values = c("brown3", "#e18585",
                               "#f0c2c2", "#faebeb",
                               "orchid3", "plum", 
                               "darkcyan",
                               "seagreen", "darkseagreen",
                               "darkgoldenrod1", "#FFD56F",
                               "#FFC55E", "bisque1")) +
  scale_color_manual(labels = names,
                     values = c("brown", "#dc7070",
                                "brown", "#f5d6d6",
                                "plum4", "orchid4",
                                "#027d7d",
                                "seagreen", "darkseagreen",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = names,
                       values= c("crosshatch", "none", "circle", "none", 
                                 "crosshatch", "none", "none", "crosshatch", "none", 
                                 "crosshatch", "none", "circle", "none")) +
  scale_pattern_density_manual(labels = names,
                               values= c(0.4, 0, 0.5, 0, 0.3, 0, 0, 0.3, 0, 0.4, 0, 0.5, 0)) +
  scale_pattern_fill_manual(labels = names,
                            values=c("#faebeb", "#faebeb",
                                     "brown3", "#faebeb",
                                     "#d5caf1", "orchid4",
                                     "#98bbbb",
                                     "#cdeede", "darkseagreen", 
                                     "bisque1", 
                                     "bisque1", "darkgoldenrod1", "bisque1")) +
  scale_y_continuous(labels = abs) +
  xlab("Kinship systems") +
  ylab("Effective population size (village's cemetery)") +
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x = element_blank(),
  axis.ticks.x = element_blank())

df = significance(alltab, "Pi", "sex")

h3 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  geom_vline(xintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("pink3", "#d87272", "#cc1a1a", "navyblue"), na.value="white") +
  ggtitle("mean effective size (village's cemetery)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

##### LOCAL
dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY1 = read.table(path, header=TRUE)
tabY1$kinship = rep("SP", length(tabY1$Gen))
tabY1 = tabY1[tabY1$Gen == 100,]
tabY1$sex = rep("male", length(tabY1$Gen))
tabY1$descent = rep("patrilineal", length(tabY1$Gen))
tabY1$Pi = tabY1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM1 = read.table(path, header=TRUE)
tabM1$kinship = rep("SP", length(tabM1$Gen))
tabM1 = tabM1[tabM1$Gen == 100,]
tabM1$sex = rep("female", length(tabM1$Gen))
tabM1$descent = rep("patrilineal", length(tabM1$Gen))
tabM1$Pi = tabM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/local_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF1 = read.table(path, header=TRUE)
tabYF1$kinship = rep("SP_F", length(tabYF1$Gen))
tabYF1 = tabYF1[tabYF1$Gen == 100,]
tabYF1$sex = rep("male", length(tabYF1$Gen))
tabYF1$descent = rep("patrilineal", length(tabYF1$Gen))
tabYF1$Pi = tabYF1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF1 = read.table(path, header=TRUE)
tabMF1$kinship = rep("SP_F", length(tabMF1$Gen))
tabMF1 = tabMF1[tabMF1$Gen == 100,]
tabMF1$sex = rep("female", length(tabMF1$Gen))
tabMF1$descent = rep("patrilineal", length(tabMF1$Gen))
tabMF1$Pi = tabMF1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_matrilineal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM1 = read.table(path, header=TRUE)
tabYM1$kinship = rep("SM", length(tabYM1$Gen))
tabYM1 = tabYM1[tabYM1$Gen == 100,]
tabYM1$sex = rep("male", length(tabYM1$Gen))
tabYM1$descent = rep("matrilineal", length(tabYM1$Gen))
tabYM1$Pi = tabYM1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM1 = read.table(path, header=TRUE)
tabMM1$kinship = rep("SM", length(tabMM1$Gen))
tabMM1 = tabMM1[tabMM1$Gen == 100,]
tabMM1$sex = rep("female", length(tabMM1$Gen))
tabMM1$descent = rep("matrilineal", length(tabMM1$Gen))
tabMM1$Pi = tabMM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_matrilineal_villages/local_mother/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYMM1 = read.table(path, header=TRUE)
tabYMM1$kinship = rep("SM_M", length(tabYMM1$Gen))
tabYMM1 = tabYMM1[tabYMM1$Gen == 100,]
tabYMM1$sex = rep("male", length(tabYMM1$Gen))
tabYMM1$descent = rep("matrilineal", length(tabYMM1$Gen))
tabYMM1$Pi = tabYMM1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMMM1 = read.table(path, header=TRUE)
tabMMM1$kinship = rep("SM_M", length(tabMMM1$Gen))
tabMMM1 = tabMMM1[tabMMM1$Gen == 100,]
tabMMM1$sex = rep("female", length(tabMMM1$Gen))
tabMMM1$descent = rep("matrilineal", length(tabMMM1$Gen))
tabMMM1$Pi = tabMMM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY2 = read.table(path, header=TRUE)
tabY2$kinship = rep("LP", length(tabY2$Gen))
tabY2 = tabY2[tabY2$Gen == 100,]
tabY2$sex = rep("male", length(tabY2$Gen))
tabY2$descent = rep("patrilineal", length(tabY2$Gen))
tabY2$Pi = tabY2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM2 = read.table(path, header=TRUE)
tabM2$kinship = rep("LP", length(tabM2$Gen))
tabM2 = tabM2[tabM2$Gen == 100,]
tabM2$sex = rep("female", length(tabM2$Gen))
tabM2$descent = rep("patrilineal", length(tabM2$Gen))
tabM2$Pi = tabM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/local_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF2 = read.table(path, header=TRUE)
tabYF2$kinship = rep("LP_F", length(tabYF2$Gen))
tabYF2 = tabYF2[tabYF2$Gen == 100,]
tabYF2$sex = rep("male", length(tabYF2$Gen))
tabYF2$descent = rep("patrilineal", length(tabYF2$Gen))
tabYF2$Pi = tabYF2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF2 = read.table(path, header=TRUE)
tabMF2$kinship = rep("LP_F", length(tabMF2$Gen))
tabMF2 = tabMF2[tabMF2$Gen == 100,]
tabMF2$sex = rep("female", length(tabMF2$Gen))
tabMF2$descent = rep("patrilineal", length(tabMF2$Gen))
tabMF2$Pi = tabMF2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_matrilineal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM2 = read.table(path, header=TRUE)
tabYM2$kinship = rep("LM", length(tabYM2$Gen))
tabYM2 = tabYM2[tabYM2$Gen == 100,]
tabYM2$sex = rep("male", length(tabYM2$Gen))
tabYM2$descent = rep("matrilineal", length(tabYM2$Gen))
tabYM2$Pi = tabYM2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM2 = read.table(path, header=TRUE)
tabMM2$kinship = rep("LM", length(tabMM2$Gen))
tabMM2 = tabMM2[tabMM2$Gen == 100,]
tabMM2$sex = rep("female", length(tabMM2$Gen))
tabMM2$descent = rep("matrilineal", length(tabMM2$Gen))
tabMM2$Pi = tabMM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_matrilineal_villages/local_mother/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYMM2 = read.table(path, header=TRUE)
tabYMM2$kinship = rep("LM_M", length(tabYMM2$Gen))
tabYMM2 = tabYMM2[tabYMM2$Gen == 100,]
tabYMM2$sex = rep("male", length(tabYMM2$Gen))
tabYMM2$descent = rep("matrilineal", length(tabYMM2$Gen))
tabYMM2$Pi = tabYMM2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMMM2 = read.table(path, header=TRUE)
tabMMM2$kinship = rep("LM_M", length(tabMMM2$Gen))
tabMMM2 = tabMMM2[tabMMM2$Gen == 100,]
tabMMM2$sex = rep("female", length(tabMMM2$Gen))
tabMMM2$descent = rep("matrilineal", length(tabMMM2$Gen))
tabMMM2$Pi = tabMMM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/strict_patrilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY3 = read.table(path, header=TRUE)
tabY3$kinship = rep("B_SP", length(tabY3$Gen))
tabY3 = tabY3[tabY3$Gen == 100,]
tabY3$sex = rep("male", length(tabY3$Gen))
tabY3$descent = rep("bilateral", length(tabY3$Gen))
tabY3$Pi = tabY3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM3 = read.table(path, header=TRUE)
tabM3$kinship = rep("B_SP", length(tabM3$Gen))
tabM3 = tabM3[tabM3$Gen == 100,]
tabM3$sex = rep("female", length(tabM3$Gen))
tabM3$descent = rep("bilateral", length(tabM3$Gen))
tabM3$Pi = tabM3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/ambilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYA3 = read.table(path, header=TRUE)
tabYA3$kinship = rep("B_A", length(tabYA3$Gen))
tabYA3 = tabYA3[tabYA3$Gen == 100,]
tabYA3$sex = rep("male", length(tabYA3$Gen))
tabYA3$descent = rep("bilateral", length(tabYA3$Gen))
tabYA3$Pi = tabYA3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMA3 = read.table(path, header=TRUE)
tabMA3$kinship = rep("B_A", length(tabMA3$Gen))
tabMA3 = tabMA3[tabMA3$Gen == 100,]
tabMA3$sex = rep("female", length(tabMA3$Gen))
tabMA3$descent = rep("bilateral", length(tabMA3$Gen))
tabMA3$Pi = tabMA3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/strict_matrilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM3 = read.table(path, header=TRUE)
tabYM3$kinship = rep("B_SM", length(tabYM3$Gen))
tabYM3 = tabYM3[tabYM3$Gen == 100,]
tabYM3$sex = rep("male", length(tabYM3$Gen))
tabYM3$descent = rep("bilateral", length(tabYM3$Gen))
tabYM3$Pi = tabYM3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM3 = read.table(path, header=TRUE)
tabMM3$kinship = rep("B_SM", length(tabMM3$Gen))
tabMM3 = tabMM3[tabMM3$Gen == 100,]
tabMM3$sex = rep("female", length(tabMM3$Gen))
tabMM3$descent = rep("bilateral", length(tabMM3$Gen))
tabMM3$Pi = tabMM3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/loose_patrilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY4 = read.table(path, header=TRUE)
tabY4$kinship = rep("B_LP", length(tabY4$Gen))
tabY4 = tabY4[tabY4$Gen == 100,]
tabY4$sex = rep("male", length(tabY4$Gen))
tabY4$descent = rep("bilateral", length(tabY4$Gen))
tabY4$Pi = tabY4$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM4 = read.table(path, header=TRUE)
tabM4$kinship = rep("B_LP", length(tabM4$Gen))
tabM4 = tabM4[tabM4$Gen == 100,]
tabM4$sex = rep("female", length(tabM4$Gen))
tabM4$descent = rep("bilateral", length(tabM4$Gen))
tabM4$Pi = tabM4$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/loose_matrilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYM4 = read.table(path, header=TRUE)
tabYM4$kinship = rep("B_LM", length(tabYM4$Gen))
tabYM4 = tabYM4[tabYM4$Gen == 100,]
tabYM4$sex = rep("male", length(tabYM4$Gen))
tabYM4$descent = rep("bilateral", length(tabYM4$Gen))
tabYM4$Pi = tabYM4$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMM4 = read.table(path, header=TRUE)
tabMM4$kinship = rep("B_LM", length(tabMM4$Gen))
tabMM4 = tabMM4[tabMM4$Gen == 100,]
tabMM4$sex = rep("female", length(tabMM4$Gen))
tabMM4$descent = rep("bilateral", length(tabMM4$Gen))
tabMM4$Pi = tabMM4$Pi/(2*5.5e-7)

tabY = merge(tabY1, tabY2, all = T)
tabY = merge(tabY, tabY3, all = T)
tabY = merge(tabY, tabY4, all = T)
tabY = merge(tabY, tabYF1, all = T)
tabY = merge(tabY, tabYMM1, all = T)
tabY = merge(tabY, tabYF2, all = T)
tabY = merge(tabY, tabYMM2, all = T)
tabY = merge(tabY, tabYM1, all = T)
tabY = merge(tabY, tabYM2, all = T)
tabY = merge(tabY, tabYM3, all = T)
tabY = merge(tabY, tabYM4, all = T)
tabY = merge(tabY, tabYA3, all = T)

tabM = merge(tabM1, tabM2, all=T)
tabM = merge(tabM, tabM3, all=T)
tabM = merge(tabM, tabM4, all=T)
tabM = merge(tabM, tabMF1, all=T)
tabM = merge(tabM, tabMF2, all=T)
tabM = merge(tabM, tabMMM1, all=T)
tabM = merge(tabM, tabMMM2, all=T)
tabM = merge(tabM, tabMM1, all = T)
tabM = merge(tabM, tabMM2, all = T)
tabM = merge(tabM, tabMM3, all = T)
tabM = merge(tabM, tabMM4, all = T)
tabM = merge(tabM, tabMA3, all = T)

alltab = merge(tabY, tabM, all=T)
alltab$sex = as.factor(alltab$sex)

alltab = rearrange_pitab(alltab)
data_mirror <- alltab %>%
  mutate(Pi = ifelse(sex == "male", -Pi, Pi))

p4 <- ggplot(data_mirror) +
  facet_wrap(~descent, ncol=3, scales='free_x') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(kinship, Pi, col=kinship), size=2, alpha=0.05, 
            position = position_identity()) +
  geom_boxplot_pattern(aes(kinship, Pi, fill = kinship, col = kinship,
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship, alpha = sex), 
                      lwd=0.75, outlier.shape=NA, position = position_identity()) +
  scale_alpha_manual(values = c(0.2, 0.2)) +
  scale_fill_manual(labels = names,
                    values = c("brown3", "#e18585",
                               "#f0c2c2", "#faebeb",
                               "orchid3", "plum", 
                               "darkcyan",
                               "seagreen", "darkseagreen",
                               "darkgoldenrod1", "#FFD56F",
                               "#FFC55E", "bisque1")) +
  scale_color_manual(labels = names,
                     values = c("brown", "#dc7070",
                                "brown", "#f5d6d6",
                                "plum4", "orchid4",
                                "#027d7d",
                                "seagreen", "darkseagreen",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = names,
                       values= c("crosshatch", "none", "circle", "none", 
                                 "crosshatch", "none", "none", "crosshatch", "none", 
                                 "crosshatch", "none", "circle", "none")) +
  scale_pattern_density_manual(labels = names,
                               values= c(0.4, 0, 0.5, 0, 0.3, 0, 0, 0.3, 0, 0.4, 0, 0.5, 0)) +
  scale_pattern_fill_manual(labels = names,
                            values=c("#faebeb", "#faebeb",
                                     "brown3", "#faebeb",
                                     "#d5caf1", "orchid4",
                                     "#98bbbb",
                                     "#cdeede", "darkseagreen", 
                                     "bisque1", 
                                     "bisque1", "darkgoldenrod1", "bisque1")) +
  scale_y_continuous(labels = abs) +
  xlab("Kinship systems") +
  ylab("Effective population size (local group's cemetery)") +
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x = element_blank(),
  axis.ticks.x = element_blank())

df = significance(alltab, "Pi", "sex")

h4 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  geom_vline(xintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("pink", "pink3", "#d87272", "#cc1a1a", "navyblue"), na.value="white") +
  ggtitle("mean effective size (local group's cemetery)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

########### Haplogroup diversity #############
dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
path = paste(dir, "strict_patrilineal_villages", sep = "")
sp = mean_hg(path, "SP", "village", "patrilineal")

path = paste(dir, "strict_patrilineal_villages", sep = "")
spf = mean_hg(path, "SP_F", "village_father", "patrilineal")

path = paste(dir, "strict_matrilineal_villages", sep = "")
sm = mean_hg(path, "SM", "village", "matrilineal")

path = paste(dir, "strict_matrilineal_villages", sep = "")
smm = mean_hg(path, "SM_M", "village_mother", "matrilineal")

dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
path = paste(dir, "loose_patrilineal_villages", sep = "")
lp = mean_hg(path, "LP", "village", "patrilineal")

path = paste(dir, "loose_patrilineal_villages", sep = "")
lpf = mean_hg(path, "LP_F", "village_father", "patrilineal")

path = paste(dir, "loose_matrilineal_villages", sep = "")
lmm = mean_hg(path, "LM_M", "village_mother", "matrilineal")

path = paste(dir, "loose_matrilineal_villages", sep = "")
lm = mean_hg(path, "LM", "village", "matrilineal")

dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
path = paste(dir, "strict_matrilocal_villages", sep = "")
bsm = mean_hg(path, "B_SM", "village", "bilateral")

path = paste(dir, "loose_matrilocal_villages", sep = "")
blm = mean_hg(path, "B_LM", "village", "bilateral")

path = paste(dir, "ambilocal_villages", sep = "")
ba = mean_hg(path, "B_A", "village", "bilateral")

path = paste(dir, "strict_patrilocal_villages", sep = "")
bsp = mean_hg(path, "B_SP", "village", "bilateral")

path = paste(dir, "loose_patrilocal_villages", sep = "")
blp = mean_hg(path, "B_LP", "village", "bilateral")

df = merge(bsm, ba, all=T)
df = merge(df, bsp, all = T)
df = merge(df, blm, all = T)
df = merge(df, blp, all = T)
df = merge(df, sm, all = T)
df = merge(df, smm, all = T)
df = merge(df, sp, all = T)
df = merge(df, spf, all = T)
df = merge(df, lm, all = T)
df = merge(df, lmm, all = T)
df = merge(df, lp, all = T)
df = merge(df, lpf, all = T)

df$kinship = factor(df$kinship, levels = c("SM", "LM", "SM_M", "LM_M", 
                                           "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                           "SP", "LP", "SP_F", "LP_F"))
df$Generation = as.factor(df$Generation)
data_mirror <- df %>%
  mutate(mean_hg_div = ifelse(Chromosome == "Y", -mean_hg_div, mean_hg_div))

p5 <- ggplot(data_mirror) +
  facet_wrap(~descent, ncol=3, scales = 'free_x') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(kinship, mean_hg_div, col=kinship), size=2, alpha=0.05, 
             position = position_identity()) +
  geom_boxplot_pattern(aes(kinship, mean_hg_div, fill = kinship, col = kinship,
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship, alpha = Chromosome), 
                      lwd=0.75, outlier.shape=NA,
                      position = position_identity()) +
  scale_alpha_manual(values = c(0.2, 0.2)) +
  scale_fill_manual(labels = names,
                    values = c("brown3", "#e18585",
                               "#f0c2c2", "#faebeb",
                               "orchid3", "plum", 
                               "darkcyan",
                               "seagreen", "darkseagreen",
                               "darkgoldenrod1", "#FFD56F",
                               "#FFC55E", "bisque1")) +
  scale_color_manual(labels = names,
                     values = c("brown", "#dc7070",
                                "brown", "#f5d6d6",
                                "plum4", "orchid4",
                                "#027d7d",
                                "seagreen", "darkseagreen",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = names,
                       values= c("crosshatch", "none", "circle", "none", 
                                 "crosshatch", "none", "none", "crosshatch", "none", 
                                 "crosshatch", "none", "circle", "none")) +
  scale_pattern_density_manual(labels = names,
                               values= c(0.4, 0, 0.5, 0, 0.3, 0, 0, 0.3, 0, 0.4, 0, 0.5, 0)) +
  scale_pattern_fill_manual(labels = names,
                            values=c("#faebeb", "#faebeb",
                                     "brown3", "#faebeb",
                                     "#d5caf1", "orchid4",
                                     "#98bbbb",
                                     "#cdeede", "darkseagreen", 
                                     "bisque1", 
                                     "bisque1", "darkgoldenrod1", "bisque1")) +
  scale_y_continuous(labels = abs) +
  xlab("Kinship systems") +
  ylab("Mean haplogroup diversity (village's cemetery)") +
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x = element_blank(),
  axis.ticks.x = element_blank())

df = significance(df, "mean_hg_div", "Chromosome")

h5 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  geom_vline(xintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("#cc1a1a", "navyblue"), na.value="white") +
  ggtitle("mean haplogroup diversity (village's cemetery)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

## local group
dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
path = paste(dir, "strict_patrilineal_villages", sep = "")
sp = mean_hg(path, "SP", "patriline", "patrilineal")

path = paste(dir, "strict_patrilineal_villages", sep = "")
spf = mean_hg(path, "SP_F", "patriline_father", "patrilineal")

path = paste(dir, "strict_matrilineal_villages", sep = "")
sm = mean_hg(path, "SM", "matriline", "matrilineal")

path = paste(dir, "strict_matrilineal_villages", sep = "")
smm = mean_hg(path, "SM_M", "matriline_mother", "matrilineal")

dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
path = paste(dir, "loose_patrilineal_villages", sep = "")
lp = mean_hg(path, "LP", "patriline", "patrilineal")

path = paste(dir, "loose_patrilineal_villages", sep = "")
lpf = mean_hg(path, "LP_F", "patriline_father", "patrilineal")

path = paste(dir, "loose_matrilineal_villages", sep = "")
lmm = mean_hg(path, "LM_M", "matriline_mother", "matrilineal")

path = paste(dir, "loose_matrilineal_villages", sep = "")
lm = mean_hg(path, "LM", "matriline", "matrilineal")

dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
path = paste(dir, "strict_matrilocal_villages", sep = "")
bsm = mean_hg(path, "B_SM", "matriline", "bilateral")

path = paste(dir, "loose_matrilocal_villages", sep = "")
blm = mean_hg(path, "B_LM", "matriline", "bilateral")

path = paste(dir, "ambilocal_villages", sep = "")
ba = mean_hg(path, "B_A", "local", "bilateral")

path = paste(dir, "strict_patrilocal_villages", sep = "")
bsp = mean_hg(path, "B_SP", "patriline", "bilateral")

path = paste(dir, "loose_patrilocal_villages", sep = "")
blp = mean_hg(path, "B_LP", "patriline", "bilateral")

df = merge(bsm, ba, all=T)
df = merge(df, bsp, all = T)
df = merge(df, blm, all = T)
df = merge(df, blp, all = T)
df = merge(df, sm, all = T)
df = merge(df, smm, all = T)
df = merge(df, sp, all = T)
df = merge(df, spf, all = T)
df = merge(df, lm, all = T)
df = merge(df, lmm, all = T)
df = merge(df, lp, all = T)
df = merge(df, lpf, all = T)

df$kinship = factor(df$kinship, levels = c("SM", "LM", "SM_M", "LM_M", 
                                           "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                           "SP", "LP", "SP_F", "LP_F"))

df$Generation = as.factor(df$Generation)
data_mirror <- df %>%
  mutate(mean_hg_div = ifelse(Chromosome == "Y", -mean_hg_div, mean_hg_div))

p6 <- ggplot(data_mirror) +
  facet_wrap(~descent, ncol=3, scales = 'free_x') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(kinship, mean_hg_div, col=kinship), size=2, alpha=0.05, 
             position = position_identity()) +
  geom_boxplot_pattern(aes(kinship, mean_hg_div, fill = kinship, col = kinship,
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship, alpha = Chromosome), 
                      lwd=0.75, outlier.shape=NA,
                      position = position_identity()) +
  scale_alpha_manual(values = c(0.2, 0.2)) +
  scale_fill_manual(labels = names,
                    values = c("brown3", "#e18585",
                               "#f0c2c2", "#faebeb",
                               "orchid3", "plum", 
                               "darkcyan",
                               "seagreen", "darkseagreen",
                               "darkgoldenrod1", "#FFD56F",
                               "#FFC55E", "bisque1")) +
  scale_color_manual(labels = names,
                     values = c("brown", "#dc7070",
                                "brown", "#f5d6d6",
                                "plum4", "orchid4",
                                "#027d7d",
                                "seagreen", "darkseagreen",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = names,
                       values= c("crosshatch", "none", "circle", "none", 
                                 "crosshatch", "none", "none", "crosshatch", "none", 
                                 "crosshatch", "none", "circle", "none")) +
  scale_pattern_density_manual(labels = names,
                               values= c(0.4, 0, 0.5, 0, 0.3, 0, 0, 0.3, 0, 0.4, 0, 0.5, 0)) +
  scale_pattern_fill_manual(labels = names,
                            values=c("#faebeb", "#faebeb",
                                     "brown3", "#faebeb",
                                     "#d5caf1", "orchid4",
                                     "#98bbbb",
                                     "#cdeede", "darkseagreen", 
                                     "bisque1", 
                                     "bisque1", "darkgoldenrod1", "bisque1")) +
  scale_y_continuous(labels = abs) +
  xlab("Kinship systems") +
  ylab("Mean haplogroup diversity (local group's cemetery)") +
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x = element_blank(),
  axis.ticks.x = element_blank())

df = significance(df, "mean_hg_div", "Chromosome")

h6 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  geom_vline(xintercept=13.5, color = "black", linetype = "dashed", linewidth=1) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("pink", "pink3","#d87272", "#cc1a1a", "navyblue"), na.value="white") +
  ggtitle("mean haplogroup diversity (local group's cemetery)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

###### Figure 2 (relatedness + haplogroup diversity) ######

p_temp <- plot_grid(p1 + theme(legend.position = "none"),
                    p2 + theme(legend.position = "none"),
                    p5 + theme(legend.position = "none"),
                    p6 + theme(legend.position = "none"),
                    ncol = 2, labels = "auto")
p <- plot_grid(p_temp, legend = get_legend(p1), ncol = 1, rel_heights = c(0.73, 0.27))

png(file = "figure_2.png", width = 3100, height = 4200)
p
dev.off()

###### Figure S5 (effective population size) ######

p_temp <- plot_grid(p3 + theme(legend.position = "none"),
                    p4 + theme(legend.position = "none"),
                    ncol = 2, labels = "auto")
p <- plot_grid(p_temp, legend = get_legend(p1), ncol = 1, rel_heights = c(0.6, 0.4))

png(file = "figure_S5.png", width = 3200, height = 3000)
p
dev.off()

###### Figure S4 (p-values) ######

p_temp <- plot_grid(h1 + theme(legend.position = "none"),
                    h2 + theme(legend.position = "none"),
                    h5 + theme(legend.position = "none"),
                    h6 + theme(legend.position = "none"),
                    h3 + theme(legend.position = "none"),
                    h4 + theme(legend.position = "none"),
                    ncol = 2, labels = "auto")
p <- plot_grid(p_temp, legend = get_legend(h6), ncol = 1, rel_heights = c(0.95, 0.05))

png(file = "Figure_S4.png", width = 3500, height = 5000)
p
dev.off()

###### Figure 3 (effective population size under the patrilineal to bilateral scenario) ######

fun_mean <- function(data, indices, gen) {
  d <- data[indices,]
  d <- d[d$Gen == gen,]
  mean = mean(d$Pi)
  return(mean)
}

prepare_df <- function(model) {
  "
  Description : transforms result tables into data frames with mean pi + 
  compute mean pi X/A with bootstrap for confidence intervals
  model : string indicating the model
  "
  # Load required library
  library(boot)
  
  # Helper function to compute bootstrap statistics and confidence intervals
  compute_bootstrap_ci <- function(file_path, generations) {
    # Read and clean the data
    tab <- read.table(file_path, header = TRUE) %>% na.omit()
    
    # Bootstrap results for each generation
    bootstrap_results <- lapply(generations, function(gen) {
      bootstrap <- boot(data = tab, statistic = fun_mean, R = 10000, gen = gen)
      ci <- boot.ci(bootstrap, conf = 0.95, type = "bca")
      list(mean = bootstrap$t0, ICinf = ci$bca[4], ICsup = ci$bca[5])
    })
    
    # Transform results into a data frame
    data.frame(
      Gen = generations,
      mean = sapply(bootstrap_results, `[[`, "mean"),
      ICinf = sapply(bootstrap_results, `[[`, "ICinf"),
      ICsup = sapply(bootstrap_results, `[[`, "ICsup")
    )
  }
  
  # Get unique generations from Pi_Y file for consistency
  tab_Y <- read.table("Pi_Y_mean_by_rep.txt", header = TRUE)
  generations <- unique(na.omit(tab_Y)$Gen)
  
  # Compute pi_tab_Y and pi_tab_Mito data frames
  pi_tab_Y <- compute_bootstrap_ci("Pi_Y_mean_by_rep.txt", generations)
  pi_tab_Mito <- compute_bootstrap_ci("Pi_Mito_mean_by_rep.txt", generations)
  
  # Add model column to each data frame
  pi_tab_Y$model <- model
  pi_tab_Mito$model <- model
  
  # Return both data frames in a list
  return(list(pi_tab_Y, pi_tab_Mito))
}

pi_tab <- function(paths, name, scale='local') {
  posd <- position_dodge(5)
  
  flag = F
  for (i in seq(1, length(paths))) {
    path = paths[i]
    model = name[i]
    
    if (flag == F) {
      var = c('pi_tab_Y', 'pi_tab_Mito', 'pi_tab_MY')
    }
    else {
      var = c('pi_tab_bis_Y', 'pi_tab_bis_Mito', 'pi_tab_bis_MY')
    }
    
    setwd(path)
    if (scale == 'global') {
      val = prepare_df_global(model)
    }
    else {
      val = prepare_df(model)
    }
    for (i in seq(var)) assign(var[i], val[i])
    
    if (flag == T) {
      # merge tables
      pi_tab_Y <- merge(pi_tab_Y, pi_tab_bis_Y, all = TRUE)
      pi_tab_Mito <- merge(pi_tab_Mito, pi_tab_bis_Mito, all = TRUE)
    }
    flag = T
  }
  return(list(pi_tab_Y, pi_tab_Mito))
}

figure_2_transitions <- function(pi_tab, legend, xmin, ymin, ymax,  segment, text, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace) {
  
  pi_tab$Gen = pi_tab$Gen + xmin
  
  p1 <- ggplot(pi_tab, aes(Gen, mean, shape = sex, col= sex, linetype = sex)) +
    theme_cowplot() +
    theme(text=element_text(size=fontsize1),
          axis.text = element_text(size=fontsize2),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(color="bisque", linewidth=0.3)) +
    annotate("rect", xmin = -200, xmax = -100, ymin = ymin, ymax = ymax,
             alpha = .1,fill = "orange") +
    annotate("rect", xmin = -100, xmax = 0, ymin = ymin, ymax = ymax,
             alpha = .05,fill = "orange") +
    labs(x= 'Time (generations)', y = 'Effective population size',
         color = "Sex", shape = "Sex", linetype = "Sex") +
    scale_x_continuous(breaks = seq(xmin, 0, 20)) +
    scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(ymin, ymax)) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    annotate("segment", x = xmin, xend = xmin, y = ymin, yend = segment, colour = "salmon2", linetype='longdash', linewidth = 0.5) +
    annotate("text", x = xmin, y = textY, label = bquote(''*t[0]*''), col = 'salmon2', size = 6) +
    annotate("segment", x = -100, xend = -100, y = ymin, yend = segment, colour = "salmon2", linetype='longdash', linewidth=0.5) +
    annotate("text", x = -100, y = textY, label = bquote(''*t[1]*''), col = 'salmon2', size = 6) +
    geom_line(aes(linetype = sex), alpha = 0.5) +
    geom_pointrange(aes(ymin=ICinf, ymax=ICsup), size = 0.75, stroke = 0.5, linewidth = 0.5) +
    theme(legend.key.size = unit(1, "cm"),
          legend.key.height = unit(1, "cm"),
          legend.title = element_text(face = 'bold', size = 14),
          legend.text = element_text(size=12),
          legend.spacing.x = unit(2, 'mm'))
  
  return(p1)
}

paths = c("Tables/Pi/unilineal/extended/r=0/k=0.1/FT=150/Patrilineal_2_strict_patrilocal_villages/village")
name = c("4b")
legend = c("Patrilineal descent, patrilocal residence \n then bilateral descent, patrilocal residence")
xmin=-200
ymin = 0
ymax = 990
segment = 890
text = 945
col = c("brown3", "darkgoldenrod2")
sh = c(16, 16)
lt = c("solid", "solid")
fontsize1 = 11
fontsize2 = 9
nCol=1
legendSpace=0.2

pi_tabs = pi_tab(paths, name, 'local')
pi_tab_Y = pi_tabs[[1]][[1]]
pi_tab_M = pi_tabs[[2]][[1]]

pi_tab_Y$sex = rep("male", length(pi_tab_Y$mean))
pi_tab_M$sex = rep("female", length(pi_tab_M$mean))

### Correct nucleotide diversity by the mutation rate
pi_tab_Y$mean = pi_tab_Y$mean/(2*2.5e-8)
pi_tab_Y$ICinf = pi_tab_Y$ICinf/(2*2.5e-8)
pi_tab_Y$ICsup = pi_tab_Y$ICsup/(2*2.5e-8)
pi_tab_M$mean = pi_tab_M$mean/(2*5.5e-7)
pi_tab_M$ICinf = pi_tab_M$ICinf/(2*5.5e-7)
pi_tab_M$ICsup = pi_tab_M$ICsup/(2*5.5e-7)

pi_tab = merge(pi_tab_Y, pi_tab_M, all=T)
p = figure_2_transitions(pi_tab, legend, xmin,
                          ymin, ymax, segment, text, 
                          col, sh, lt, fontsize1, fontsize2, nCol,
                          legendSpace)
png(file = "figure_3.png", width = 2000, height = 1500, res = 300)
p
dev.off()

###### Effect of parameter change in the model ######

fun_median <- function(data, indices, metric) {
  d <- data[indices,]
  median = median(d[[metric]])
  return(median)
}

median_ratio <- function(data, indices, metric1, metric2) {
  d <- data[indices,]
  median_ratio = median(d[[metric1]])/median(d[[metric2]])
  return(median_ratio)
}

reshape_file <- function(path, model, mode, descent, sample=F) {
  setwd(path)
  
  if (mode == "village") {
    file = read.table("relatedness.txt", header = T)
    if (sample == T) {
      var = c("sampleRelFreqM", "sampleRelFreqF", 
              "sampleMeanRelM", "sampleMeanRelF",
              "sampleRelFirstM", "sampleRelFirstF")
    }
    else {
      var = c("relFreqM", "relFreqF",
              "meanRelM", "meanRelF",
              "relFirstM", "relFirstF")
    }
  }
  
  if (mode == "village_father") {
    file = read.table("relatedness_father.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "village_mother") {
    file = read.table("relatedness_mother.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "group") {
    file = read.table("relatedness_group.txt", header = T)
    if (sample == T) {
      var = c("sampleRelGroupFreqM", "sampleRelGroupFreqF", 
              "sampleMeanGroupRelM", "sampleMeanGroupRelF", 
              "sampleRelGroupFirstM", "sampleRelGroupFirstF")
    }
    else {
      var = c("relGroupFreqM", "relGroupFreqF", 
              "meanGroupRelM", "meanGroupRelF", 
              "relGroupFirstM", "relGroupFirstF")
    }
  }
  
  if (mode == "patrilocal_group") {
    file = read.table("relatedness_patriline.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "matrilocal_group") {
    file = read.table("relatedness_matriline.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "local_group") {
    file = read.table("relatedness_local.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "local_group_father") {
    file = read.table("relatedness_patriline_father.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  if (mode == "local_group_mother") {
    file = read.table("relatedness_matriline_mother.txt", header = T)
    var = c("relFreqM", "relFreqF",
            "meanRelM", "meanRelF",
            "relFirstM", "relFirstF")
  }
  
  reshaped_file <- pivot_longer(file, 
                                cols = var,
                                names_to = "variable",
                                values_to = "value")
  reshaped_file$category = str_sub(reshaped_file$variable, end=-2)
  reshaped_file$kinship = rep(model, length(reshaped_file$variable))
  reshaped_file$kinvar = paste(reshaped_file$kinship, reshaped_file$variable)
  reshaped_file$descent = rep(descent, length(reshaped_file$variable))
  return(reshaped_file)
}

bootstrap <- function(data, metric=NaN, metric1=NaN, metric2=NaN, ratio) {
  if (ratio == F) {
    b <- boot(data=data, statistic = fun_median, R = 10000, metric = metric)
  }
  else {
    b <- boot(data=data, statistic=median_ratio, R=10000, metric1=metric1, metric2=metric2)
  }
  ci <- boot.ci(b, conf = 0.95, type="bca")
  median <- b$t0
  ICinf <- ci$bca[4]
  ICsup <- ci$bca[5]
  return(list(median, ICinf, ICsup))
}

CI_metrics <- function(data, metrics, col, kinship) {
  m = metrics[1]
  b = bootstrap(data, metric=m, ratio=F)
  median = b[[1]]
  ICinf = b[[2]]
  ICsup = b[[3]]
  df_metric <- data.frame("metric" = m, "median" = median, "ICinf" = ICinf, "ICsup" = ICsup, "col" = col[1], "kinship" = kinship)
  for (i in seq(2,length(metrics))) {
    m = metrics[i]
    b = bootstrap(data, metric=m, ratio=F)
    median = b[[1]]
    ICinf = b[[2]]
    ICsup = b[[3]]
    df_metric <- merge(df_metric, 
                       data.frame("metric" = m, "median" = median, "ICinf" = ICinf, "ICsup" = ICsup, "col" = col[i], "kinship" = kinship), 
                       all = T)
  }
  return(df_metric)
}

mean_hg <- function(path, kinship, mode, descent) {
  setwd(path)
  
  df = read.table(paste0("hg_freq_", mode, ".txt"), header=T)
  
  df2 <- df %>%
    group_by(Replicat, Generation, Chromosome, Village) %>%
    summarize(Hg_diversity = nIndividuals[1] / (nIndividuals[1] - 1) * (1 - sum(Frequency^2, na.rm=T))) %>%
    ungroup()
  
  mean_df <- df2 %>%
    group_by(Replicat, Generation, Chromosome) %>%
    summarize(mean_hg_div = mean(Hg_diversity, na.rm=T)) %>%
    ungroup()
  
  if (kinship == "SM2B_SM" | kinship == "SP2B_SP") {
    mean_df = mean_df[mean_df$Generation == 20200,]
  }
  else {
    mean_df = mean_df[mean_df$Generation == 20100,]
  }
  mean_df$kinship = rep(kinship, length(mean_df$Replicat))
  mean_df$descent = rep(descent, length(mean_df$Replicat))
  mean_df$descent = factor(mean_df$descent, levels = c("matrilineal", "bilateral", "patrilineal"))
  return(mean_df)
}

rearrange_reshape <- function(reshape2) {
  relabel.variable = as_labeller(c(meanRelM = "male", 
                                   meanRelF = "female", 
                                   relFirstM = "male", relFirstF = "female",
                                   relFreqM = "male", relFreqF = "female"))
  category.label =  c("mean rel", "freq 1st degree rel", "freq 1st, 2nd and 3rd degree rel")
  names(category.label) = c("meanRel", "relFirst", "relFreq")
  reshape_bis = data.frame("variable"=reshape2$variable,"value"=reshape2$value, 
                           "category"=reshape2$category, "kinship"=reshape2$kinship,
                           "descent" = reshape2$descent)
  
  reshape2 <- reshape2 %>% 
    as_tibble() %>%
    rowwise() %>% 
    mutate(sex = factor(relabel.variable(variable))) %>%
    arrange(value) %>%
    mutate(kinship = factor(kinship, levels = c("SM", "LM", "SM_M", "LM_M",
                                                "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                                "SP", "LP", "SP_F", "LP_F")))
  
  reshape2 = reshape2[with(reshape2, order(sex, kinship)),]
  
  data = reshape2[reshape2$category == "meanRel",]
  data = data[data$Generation == 20100 | data$Generation == 20200,]
  data = data[!((data$kinship == "SM2B_SM" & data$Generation == 20100) | (data$kinship == "SP2B_SP" & data$Generation == 20100)),]
  data <- data %>%
    group_by(descent, kinship, sex, Replicat) %>%
    summarise(value = mean(value))
  data$descent = factor(data$descent, levels = c("matrilineal", "bilateral", "patrilineal"))
  data = data[with(data, order(sex, kinship)),]
  return(data)
}

rearrange_pitab <- function(tab) {
  tab = data.frame("Pi" = tab$Pi, "kinship" = as.factor(tab$kinship), "sex" = tab$sex, 
                   "descent" = factor(tab$descent, levels = c("matrilineal", 
                                                              "bilateral",
                                                           "patrilineal")))
  tab = na.omit(tab)
  tab <- tab %>%
    arrange(Pi) %>%
    mutate(kinship = factor(kinship, levels=c("SM", "LM", "SM_M", "LM_M",
                                              "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                              "SP", "LP", "SP_F", "LP_F")))
  
  tab = tab[with(tab, order(sex, kinship)),]
  return(tab)
}

significance <- function(data_stat, value) {
  data_stat$var = factor(paste(data_stat$kinship, data_stat$sex),
                         levels = c("SM female", "LM female", "SM_M female", "LM_M female", 
                                    "B_SM female", "B_LM female", "B_A female", 
                                    "B_SP female", "B_LP female", "SP female", 
                                    "LP female", "SP_F female", "LP_F female", "SM male", 
                                    "LM male", "SM_M male", "LM_M male", 
                                    "B_SM male", "B_LM male", "B_A male", "B_SP male", "B_LP male", 
                                    "SP male", "LP male", "SP_F male", "LP_F male"))
  
  pval_matrix = pairwise.wilcox.test(data_stat[[value]], data_stat$var, p.adjust.method="bonferroni")$p.value
  
  # Define the significance levels
  significance_matrix <- matrix("ns", nrow=nrow(pval_matrix), ncol=ncol(pval_matrix)) # Start with "not significant"
  
  # Apply thresholds
  significance_matrix[is.na(pval_matrix)] <- NA
  significance_matrix[pval_matrix < 0.05] <- "*"
  significance_matrix[pval_matrix < 0.01] <- "**"
  significance_matrix[pval_matrix < 0.001] <- "***"
  significance_matrix[pval_matrix < 0.0001] <- "****"
  significance_matrix = rbind("B_SP female"=rep(NA, 29), significance_matrix)
  significance_matrix = cbind(significance_matrix, "L_PF male" = rep(NA, 30))
  
  colnames(significance_matrix) = c("SM female", "LM female", "SM_M female", "LM_M female", 
                                    "B_SM female", "B_LM female", "B_A female", 
                                    "B_SP female", "B_LP female", "SP female", 
                                    "LP female", "SP_F female", "LP_F female", "SM male", 
                                    "LM male", "SM_M male", "LM_M male",
                                    "B_SM male", "B_LM male", "B_A male", "B_SP male", "B_LP male", 
                                    "SP male", "LP male", "SP_F male", "LP_F male")
  rownames(significance_matrix) = colnames(significance_matrix)
  
  df <- melt(significance_matrix)
  colnames(df) <- c("x", "y", "value")
  return(df)
}

########## Relatedness ###########
shape_files <- function(param) {
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  reshape_file_p1 = reshape_file(path, "SP", "village", "patrilineal")
  reshape_file_father_1 = reshape_file(path, "SP_F", "village_father", "patrilineal")
  reshape_group_p1 = reshape_file(path, "SP", "group", "patrilineal")
  reshape_sample_1 = reshape_file(path, "SP", "village", "patrilineal", T)
  reshape_group_sample_1 = reshape_file(path, "SP", "group", "patrilineal", T)
  reshape_patriline_1 = reshape_file(path, "SP", "patrilocal_group", "patrilineal")
  reshape_patriline_father_1 = reshape_file(path, "SP_F", "local_group_father", "patrilineal")
  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  reshape_file_m1 = reshape_file(path, "SM", "village", "matrilineal")
  reshape_file_mother_1 = reshape_file(path, "SM_M", "village_mother", "matrilineal")
  reshape_group_m1 = reshape_file(path, "SM", "group", "matrilineal")
  reshape_matriline_1 = reshape_file(path, "SM", "matrilocal_group", "matrilineal")
  reshape_matriline_mother_1 = reshape_file(path, "SM_M", "local_group_mother", "matrilineal")

  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  reshape_file_p2 = reshape_file(path, "LP", "village", "patrilineal")
  reshape_file_father_2 = reshape_file(path, "LP_F", "village_father", "patrilineal")
  reshape_group_p2 = reshape_file(path, "LP", "group", "patrilineal")
  reshape_sample_2 = reshape_file(path, "LP", "village", "patrilineal", T)
  reshape_group_sample_2 = reshape_file(path, "LP", "group", "patrilineal", T)
  reshape_patriline_2 = reshape_file(path, "LP", "patrilocal_group", "patrilineal")
  reshape_patriline_father_2 = reshape_file(path, "LP_F", "local_group_father", "patrilineal")
  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  reshape_file_m2 = reshape_file(path, "LM", "village", "matrilineal")
  reshape_file_mother_2 = reshape_file(path, "LM_M", "village_mother", "matrilineal")
  reshape_group_m2 = reshape_file(path, "LM", "group", "matrilineal")
  reshape_matriline_2 = reshape_file(path, "LM", "matrilocal_group", "matrilineal")
  reshape_matriline_mother_2 = reshape_file(path, "LM_M", "local_group_mother", "matrilineal")

  dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
  path = paste(dir, "Patrilocal_", param, "relatedness", sep = "")
  reshape_file_p3 = reshape_file(path, "B_SP", "village", "bilateral")
  reshape_sample_3 = reshape_file(path, "B_SP", "village", "bilateral", T)
  reshape_patriline_3 = reshape_file(path, "B_SP", "patrilocal_group", "bilateral")

  path = paste(dir, "Matrilocal_", param, "relatedness", sep = "")
  reshape_file_m3 = reshape_file(path, "B_SM", "village", "bilateral")
  reshape_matriline_3 = reshape_file(path, "B_SM", "matrilocal_group", "bilateral")

  path = paste(dir, "Ambilocal_", param, "relatedness", sep = "")
  reshape_file_a3 = reshape_file(path, "B_A", "village", "bilateral")
  reshape_local_3 = reshape_file(path, "B_A", "local_group", "bilateral")

  path = paste(dir, "Patrilocal_maleM_", param, "relatedness", sep = "")
  reshape_file_p4 = reshape_file(path, "B_LP", "village", "bilateral")
  reshape_sample_4 = reshape_file(path, "B_LP", "village", "bilateral", T)
  reshape_patriline_4 = reshape_file(path, "B_LP", "patrilocal_group", "bilateral")

  path = paste(dir, "Matrilocal_femaleM_", param, "relatedness", sep = "")
  reshape_file_m4 = reshape_file(path, "B_LM", "village", "bilateral")
  reshape_matriline_4 = reshape_file(path, "B_LM", "matrilocal_group", "bilateral")

  names = c("SM: Strict matrilineal descent, burial in wife's cemetery", 
          "LM: Loose matrilineal descent, burial in wife's cemetery",
          "SM_M: Strict matrilineal descent, burial in mother's cemetery", 
          "LM_M: Loose matrilineal descent, burial in mother's cemetery",
          "B_SM: Bilateral descent and strict matrilocal residence",
          "B_LM: Bilateral descent and loose matrilocal residence",
          "B_A: Bilateral descent and ambilocal residence",
          "B_SP: Bilateral descent and strict patrilocal residence",
          "B_LP: Bilateral descent and loose patrilocal residence",
          "SP: Strict patrilineal descent, burial in husband's cemetery", 
          "LP: Loose patrilineal descent, burial in husband's cemetery",
          "SP_F: Strict patrilineal descent, burial in father's cemetery", 
          "LP_F: Loose patrilineal descent, burial in father's cemetery")

reshape = merge(reshape_file_p1, reshape_file_p2, all = T)
reshape2 = merge(reshape, reshape_file_p3, all = T)
reshape2 = merge(reshape2, reshape_file_p4, all = T)
reshape2 = merge(reshape2, reshape_file_m1, all = T)
reshape2 = merge(reshape2, reshape_file_m2, all = T)
reshape2 = merge(reshape2, reshape_file_m3, all = T)
reshape2 = merge(reshape2, reshape_file_m4, all = T)
reshape2 = merge(reshape2, reshape_file_a3, all = T)
reshape2 <- reshape2[, -seq(3,8)]
reshape2 = merge(reshape2, reshape_file_father_1, all = T)
reshape2 = merge(reshape2, reshape_file_father_2, all = T)
reshape2 = merge(reshape2, reshape_file_mother_1, all = T)
reshape2 = merge(reshape2, reshape_file_mother_2, all = T)

# Rearrange reshape2
data_village = rearrange_reshape(reshape2)
data_village$sex=factor(data_village$sex)

reshape = merge(reshape_patriline_1, reshape_patriline_2, all = T)
reshape2 = merge(reshape, reshape_patriline_3, all = T)
reshape2 = merge(reshape2, reshape_patriline_4, all = T)
reshape2 = merge(reshape2, reshape_local_3, all = T)
reshape2 = merge(reshape2, reshape_matriline_1, all = T)
reshape2 = merge(reshape2, reshape_matriline_2, all = T)
reshape2 = merge(reshape2, reshape_matriline_3, all = T)
reshape2 = merge(reshape2, reshape_matriline_4, all = T)
reshape2 = merge(reshape2, reshape_patriline_father_1, all = T)
reshape2 = merge(reshape2, reshape_patriline_father_2, all = T)
reshape2 = merge(reshape2, reshape_matriline_mother_1, all = T)
reshape2 = merge(reshape2, reshape_matriline_mother_2, all = T)

# Rearrange reshape2
data_local = rearrange_reshape(reshape2)
data_local$sex=factor(data_local$sex)
return(list(data_village, data_local))
}

shape_files_FT <- function(param) {
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=75/"
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  reshape_file_p1 = reshape_file(path, "SP", "village", "patrilineal")
  reshape_file_father_1 = reshape_file(path, "SP_F", "village_father", "patrilineal")
  reshape_group_p1 = reshape_file(path, "SP", "group", "patrilineal")
  reshape_sample_1 = reshape_file(path, "SP", "village", "patrilineal", T)
  reshape_group_sample_1 = reshape_file(path, "SP", "group", "patrilineal", T)
  reshape_patriline_1 = reshape_file(path, "SP", "patrilocal_group", "patrilineal")
  reshape_patriline_father_1 = reshape_file(path, "SP_F", "local_group_father", "patrilineal")
  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  reshape_file_m1 = reshape_file(path, "SM", "village", "matrilineal")
  reshape_file_mother_1 = reshape_file(path, "SM_M", "village_mother", "matrilineal")
  reshape_group_m1 = reshape_file(path, "SM", "group", "matrilineal")
  reshape_matriline_1 = reshape_file(path, "SM", "matrilocal_group", "matrilineal")
  reshape_matriline_mother_1 = reshape_file(path, "SM_M", "local_group_mother", "matrilineal")
  
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=75/"
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  reshape_file_p2 = reshape_file(path, "LP", "village", "patrilineal")
  reshape_file_father_2 = reshape_file(path, "LP_F", "village_father", "patrilineal")
  reshape_group_p2 = reshape_file(path, "LP", "group", "patrilineal")
  reshape_sample_2 = reshape_file(path, "LP", "village", "patrilineal", T)
  reshape_group_sample_2 = reshape_file(path, "LP", "group", "patrilineal", T)
  reshape_patriline_2 = reshape_file(path, "LP", "patrilocal_group", "patrilineal")
  reshape_patriline_father_2 = reshape_file(path, "LP_F", "local_group_father", "patrilineal")
  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  reshape_file_m2 = reshape_file(path, "LM", "village", "matrilineal")
  reshape_file_mother_2 = reshape_file(path, "LM_M", "village_mother", "matrilineal")
  reshape_group_m2 = reshape_file(path, "LM", "group", "matrilineal")
  reshape_matriline_2 = reshape_file(path, "LM", "matrilocal_group", "matrilineal")
  reshape_matriline_mother_2 = reshape_file(path, "LM_M", "local_group_mother", "matrilineal")
  
  dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
  path = paste(dir, "Patrilocal_", param, "relatedness", sep = "")
  reshape_file_p3 = reshape_file(path, "B_SP", "village", "bilateral")
  reshape_sample_3 = reshape_file(path, "B_SP", "village", "bilateral", T)
  reshape_patriline_3 = reshape_file(path, "B_SP", "patrilocal_group", "bilateral")
  
  path = paste(dir, "Matrilocal_", param, "relatedness", sep = "")
  reshape_file_m3 = reshape_file(path, "B_SM", "village", "bilateral")
  reshape_matriline_3 = reshape_file(path, "B_SM", "matrilocal_group", "bilateral")
  
  path = paste(dir, "Ambilocal_", param, "relatedness", sep = "")
  reshape_file_a3 = reshape_file(path, "B_A", "village", "bilateral")
  reshape_local_3 = reshape_file(path, "B_A", "local_group", "bilateral")
  
  path = paste(dir, "Patrilocal_maleM_", param, "relatedness", sep = "")
  reshape_file_p4 = reshape_file(path, "B_LP", "village", "bilateral")
  reshape_sample_4 = reshape_file(path, "B_LP", "village", "bilateral", T)
  reshape_patriline_4 = reshape_file(path, "B_LP", "patrilocal_group", "bilateral")
  
  path = paste(dir, "Matrilocal_femaleM_", param, "relatedness", sep = "")
  reshape_file_m4 = reshape_file(path, "B_LM", "village", "bilateral")
  reshape_matriline_4 = reshape_file(path, "B_LM", "matrilocal_group", "bilateral")
  
  names = c("SM: Strict matrilineal descent, burial in wife's cemetery", 
            "LM: Loose matrilineal descent, burial in wife's cemetery",
            "SM_M: Strict matrilineal descent, burial in mother's cemetery", 
            "LM_M: Loose matrilineal descent, burial in mother's cemetery",
            "B_SM: Bilateral descent and strict matrilocal residence",
            "B_LM: Bilateral descent and loose matrilocal residence",
            "B_A: Bilateral descent and ambilocal residence",
            "B_SP: Bilateral descent and strict patrilocal residence",
            "B_LP: Bilateral descent and loose patrilocal residence",
            "SP: Strict patrilineal descent, burial in husband's cemetery", 
            "LP: Loose patrilineal descent, burial in husband's cemetery",
            "SP_F: Strict patrilineal descent, burial in father's cemetery", 
            "LP_F: Loose patrilineal descent, burial in father's cemetery")
  
  reshape = merge(reshape_file_p1, reshape_file_p2, all = T)
  reshape2 = merge(reshape, reshape_file_p3, all = T)
  reshape2 = merge(reshape2, reshape_file_p4, all = T)
  reshape2 = merge(reshape2, reshape_file_m1, all = T)
  reshape2 = merge(reshape2, reshape_file_m2, all = T)
  reshape2 = merge(reshape2, reshape_file_m3, all = T)
  reshape2 = merge(reshape2, reshape_file_m4, all = T)
  reshape2 = merge(reshape2, reshape_file_a3, all = T)
  reshape2 <- reshape2[, -seq(3,8)]
  reshape2 = merge(reshape2, reshape_file_father_1, all = T)
  reshape2 = merge(reshape2, reshape_file_father_2, all = T)
  reshape2 = merge(reshape2, reshape_file_mother_1, all = T)
  reshape2 = merge(reshape2, reshape_file_mother_2, all = T)
  
  # Rearrange reshape2
  data_village = rearrange_reshape(reshape2)
  data_village$sex=factor(data_village$sex)
  
  reshape = merge(reshape_patriline_1, reshape_patriline_2, all = T)
  reshape2 = merge(reshape, reshape_patriline_3, all = T)
  reshape2 = merge(reshape2, reshape_patriline_4, all = T)
  reshape2 = merge(reshape2, reshape_local_3, all = T)
  reshape2 = merge(reshape2, reshape_matriline_1, all = T)
  reshape2 = merge(reshape2, reshape_matriline_2, all = T)
  reshape2 = merge(reshape2, reshape_matriline_3, all = T)
  reshape2 = merge(reshape2, reshape_matriline_4, all = T)
  reshape2 = merge(reshape2, reshape_patriline_father_1, all = T)
  reshape2 = merge(reshape2, reshape_patriline_father_2, all = T)
  reshape2 = merge(reshape2, reshape_matriline_mother_1, all = T)
  reshape2 = merge(reshape2, reshape_matriline_mother_2, all = T)
  
  # Rearrange reshape2
  data_local = rearrange_reshape(reshape2)
  data_local$sex=factor(data_local$sex)
  return(list(data_village, data_local))
}

median_ratios <- function(data, mod) {
  data <- data %>%
  group_by(descent, kinship, sex) %>%
  summarise(median = median(value, na.rm=T))

  data_ratio <- data %>%
  group_by(descent, kinship) %>%
  summarise(female_median = median[sex == "female"],
            male_median = median[sex == "male"]) %>%
  # Calculate the ratio of female over male median
  mutate(ratio = female_median - male_median) %>%
  # Select only relevant columns: descent, kinship, and the ratio
  select(descent, kinship, male_median, female_median, ratio)

  data_ratio$model = rep(mod, length(data_ratio$descent))
  return(data_ratio)
}

##################################################################
########################### MAIN #################################
##################################################################

########################## VILLAGE ###############################
# Argument in shape_file depends on the name of the file  
df = shape_files("")
data_village = df[[1]]
data_local = df[[2]]

data_village_base <- median_ratios(data_village, "baseline")
data_local_base <- median_ratios(data_local, "baseline")

df = shape_files("small_pop_")
data_village = df[[1]]
data_local = df[[2]]

data_village_200ind <- median_ratios(data_village, "N=200")
data_local_200ind <- median_ratios(data_local, "N=200")

df = shape_files("2villages_")
data_village = df[[1]]
data_local = df[[2]]

data_village_2v <- median_ratios(data_village, "Nv=2")
data_local_2v <- median_ratios(data_local, "Nv=2")

df = shape_files("low_migration_")
data_village = df[[1]]
data_local = df[[2]]

data_village_m0.1 <- median_ratios(data_village, "m=0.1")
data_local_m0.1 <- median_ratios(data_local, "m=0.1")

df = shape_files_FT("")
data_village = df[[1]]
data_local = df[[2]]

data_village_FT <- median_ratios(data_village, "FT=75")
data_local_FT <- median_ratios(data_local, "FT=75")

df_V = merge(data_village_base, data_village_200ind, all=T)
df_V = merge(df_V, data_village_2v, all=T)
df_V = merge(df_V, data_village_m0.1, all=T)
df_V = merge(df_V, data_village_FT, all=T)

p1 <- ggplot(df_V, aes(kinship, ratio)) +
  geom_hline(yintercept=0, linetype = "dashed", col = "black") +
  geom_path(col = "#9c0aa3", linewidth = 0.5, group=1) +
  geom_point(col = "#9c0aa3", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Female - male relatedness difference (village)") +
  theme_bw()
p1

a <- ggplot(df_V) +
  geom_path(aes(kinship, female_median), col = "brown", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, female_median), col = "brown", size = 2) +
  geom_path(aes(kinship, male_median), col = "darkgoldenrod2", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, male_median), col = "darkgoldenrod2", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  scale_y_continuous(trans="log10") +
  labs(x = "Kinship systems", y = "Median female and male relatedness (village)") +
  theme_bw()
a

####################### LOCAL ###########################

df_L = merge(data_local_base, data_local_200ind, all=T)
df_L = merge(df_L, data_local_2v, all=T)
df_L = merge(df_L, data_local_m0.1, all=T)
df_L = merge(df_L, data_local_FT, all=T)

p2 <- ggplot(df_L, aes(kinship, ratio)) +
  geom_hline(yintercept=0, linetype = "dashed", col = "black") +
  geom_path(col = "#9c0aa3", linewidth=0.5, group=1) +
  geom_point(col = "#9c0aa3", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Female - male relatedness difference (local group)") +
  theme_bw()
p2

b <- ggplot(df_L) +
  geom_path(aes(kinship, female_median), col = "brown", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, female_median), col = "brown", size = 2) +
  geom_path(aes(kinship, male_median), col = "darkgoldenrod2", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, male_median), col = "darkgoldenrod2", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Median female and male relatedness (local group)") +
  theme_bw()
b

######################## Nucleotide diversity #########################
##### VILLAGE
shape_Pi_files <- function(param, level) {
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/Patrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY1 = read.table(path, header=TRUE)
  tabY1$kinship = rep("SP", length(tabY1$Gen))
  tabY1 = tabY1[tabY1$Gen == 100,]
  tabY1$sex = rep("male", length(tabY1$Gen))
  tabY1$descent = rep("patrilineal", length(tabY1$Gen))
  tabY1$Pi = tabY1$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM1 = read.table(path, header=TRUE)
  tabM1$kinship = rep("SP", length(tabM1$Gen))
  tabM1 = tabM1[tabM1$Gen == 100,]
  tabM1$sex = rep("female", length(tabM1$Gen))
  tabM1$descent = rep("patrilineal", length(tabM1$Gen))
  tabM1$Pi = tabM1$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/Patrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYF1 = read.table(path, header=TRUE)
  tabYF1$kinship = rep("SP_F", length(tabYF1$Gen))
  tabYF1 = tabYF1[tabYF1$Gen == 100,]
  tabYF1$sex = rep("male", length(tabYF1$Gen))
  tabYF1$descent = rep("patrilineal", length(tabYF1$Gen))
  tabYF1$Pi = tabYF1$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMF1 = read.table(path, header=TRUE)
  tabMF1$kinship = rep("SP_F", length(tabMF1$Gen))
  tabMF1 = tabMF1[tabMF1$Gen == 100,]
  tabMF1$sex = rep("female", length(tabMF1$Gen))
  tabMF1$descent = rep("patrilineal", length(tabMF1$Gen))
  tabMF1$Pi = tabMF1$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/Matrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM1 = read.table(path, header=TRUE)
  tabYM1$kinship = rep("SM", length(tabYM1$Gen))
  tabYM1 = tabYM1[tabYM1$Gen == 100,]
  tabYM1$sex = rep("male", length(tabYM1$Gen))
  tabYM1$descent = rep("matrilineal", length(tabYM1$Gen))
  tabYM1$Pi = tabYM1$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM1 = read.table(path, header=TRUE)
  tabMM1$kinship = rep("SM", length(tabMM1$Gen))
  tabMM1 = tabMM1[tabMM1$Gen == 100,]
  tabMM1$sex = rep("female", length(tabMM1$Gen))
  tabMM1$descent = rep("matrilineal", length(tabMM1$Gen))
  tabMM1$Pi = tabMM1$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/Matrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYMM1 = read.table(path, header=TRUE)
  tabYMM1$kinship = rep("SM_M", length(tabYMM1$Gen))
  tabYMM1 = tabYMM1[tabYMM1$Gen == 100,]
  tabYMM1$sex = rep("male", length(tabYMM1$Gen))
  tabYMM1$descent = rep("matrilineal", length(tabYMM1$Gen))
  tabYMM1$Pi = tabYMM1$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMMM1 = read.table(path, header=TRUE)
  tabMMM1$kinship = rep("SM_M", length(tabMMM1$Gen))
  tabMMM1 = tabMMM1[tabMMM1$Gen == 100,]
  tabMMM1$sex = rep("female", length(tabMMM1$Gen))
  tabMMM1$descent = rep("matrilineal", length(tabMMM1$Gen))
  tabMMM1$Pi = tabMMM1$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=150/Patrilineal_villages_rf_maleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY2 = read.table(path, header=TRUE)
  tabY2$kinship = rep("LP", length(tabY2$Gen))
  tabY2 = tabY2[tabY2$Gen == 100,]
  tabY2$sex = rep("male", length(tabY2$Gen))
  tabY2$descent = rep("patrilineal", length(tabY2$Gen))
  tabY2$Pi = tabY2$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM2 = read.table(path, header=TRUE)
  tabM2$kinship = rep("LP", length(tabM2$Gen))
  tabM2 = tabM2[tabM2$Gen == 100,]
  tabM2$sex = rep("female", length(tabM2$Gen))
  tabM2$descent = rep("patrilineal", length(tabM2$Gen))
  tabM2$Pi = tabM2$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=150/Patrilineal_villages_rf_maleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYF2 = read.table(path, header=TRUE)
  tabYF2$kinship = rep("LP_F", length(tabYF2$Gen))
  tabYF2 = tabYF2[tabYF2$Gen == 100,]
  tabYF2$sex = rep("male", length(tabYF2$Gen))
  tabYF2$descent = rep("patrilineal", length(tabYF2$Gen))
  tabYF2$Pi = tabYF2$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMF2 = read.table(path, header=TRUE)
  tabMF2$kinship = rep("LP_F", length(tabMF2$Gen))
  tabMF2 = tabMF2[tabMF2$Gen == 100,]
  tabMF2$sex = rep("female", length(tabMF2$Gen))
  tabMF2$descent = rep("patrilineal", length(tabMF2$Gen))
  tabMF2$Pi = tabMF2$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=150/Matrilineal_villages_rf_femaleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM2 = read.table(path, header=TRUE)
  tabYM2$kinship = rep("LM", length(tabYM2$Gen))
  tabYM2 = tabYM2[tabYM2$Gen == 100,]
  tabYM2$sex = rep("male", length(tabYM2$Gen))
  tabYM2$descent = rep("matrilineal", length(tabYM2$Gen))
  tabYM2$Pi = tabYM2$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM2 = read.table(path, header=TRUE)
  tabMM2$kinship = rep("LM", length(tabMM2$Gen))
  tabMM2 = tabMM2[tabMM2$Gen == 100,]
  tabMM2$sex = rep("female", length(tabMM2$Gen))
  tabMM2$descent = rep("matrilineal", length(tabMM2$Gen))
  tabMM2$Pi = tabMM2$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=150/Matrilineal_villages_rf_femaleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYMM2 = read.table(path, header=TRUE)
  tabYMM2$kinship = rep("LM_M", length(tabYMM2$Gen))
  tabYMM2 = tabYMM2[tabYMM2$Gen == 100,]
  tabYMM2$sex = rep("male", length(tabYMM2$Gen))
  tabYMM2$descent = rep("matrilineal", length(tabYMM2$Gen))
  tabYMM2$Pi = tabYMM2$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMMM2 = read.table(path, header=TRUE)
  tabMMM2$kinship = rep("LM_M", length(tabMMM2$Gen))
  tabMMM2 = tabMMM2[tabMMM2$Gen == 100,]
  tabMMM2$sex = rep("female", length(tabMMM2$Gen))
  tabMMM2$descent = rep("matrilineal", length(tabMMM2$Gen))
  tabMMM2$Pi = tabMMM2$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/bilateral/regular/r=0/Patrilocal_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY3 = read.table(path, header=TRUE)
  tabY3$kinship = rep("B_SP", length(tabY3$Gen))
  tabY3 = tabY3[tabY3$Gen == 100,]
  tabY3$sex = rep("male", length(tabY3$Gen))
  tabY3$descent = rep("bilateral", length(tabY3$Gen))
  tabY3$Pi = tabY3$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM3 = read.table(path, header=TRUE)
  tabM3$kinship = rep("B_SP", length(tabM3$Gen))
  tabM3 = tabM3[tabM3$Gen == 100,]
  tabM3$sex = rep("female", length(tabM3$Gen))
  tabM3$descent = rep("bilateral", length(tabM3$Gen))
  tabM3$Pi = tabM3$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/bilateral/regular/r=0/Ambilocal_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYA3 = read.table(path, header=TRUE)
  tabYA3$kinship = rep("B_A", length(tabYA3$Gen))
  tabYA3 = tabYA3[tabYA3$Gen == 100,]
  tabYA3$sex = rep("male", length(tabYA3$Gen))
  tabYA3$descent = rep("bilateral", length(tabYA3$Gen))
  tabYA3$Pi = tabYA3$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMA3 = read.table(path, header=TRUE)
  tabMA3$kinship = rep("B_A", length(tabMA3$Gen))
  tabMA3 = tabMA3[tabMA3$Gen == 100,]
  tabMA3$sex = rep("female", length(tabMA3$Gen))
  tabMA3$descent = rep("bilateral", length(tabMA3$Gen))
  tabMA3$Pi = tabMA3$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/bilateral/regular/r=0/Matrilocal_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM3 = read.table(path, header=TRUE)
  tabYM3$kinship = rep("B_SM", length(tabYM3$Gen))
  tabYM3 = tabYM3[tabYM3$Gen == 100,]
  tabYM3$sex = rep("male", length(tabYM3$Gen))
  tabYM3$descent = rep("bilateral", length(tabYM3$Gen))
  tabYM3$Pi = tabYM3$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM3 = read.table(path, header=TRUE)
  tabMM3$kinship = rep("B_SM", length(tabMM3$Gen))
  tabMM3 = tabMM3[tabMM3$Gen == 100,]
  tabMM3$sex = rep("female", length(tabMM3$Gen))
  tabMM3$descent = rep("bilateral", length(tabMM3$Gen))
  tabMM3$Pi = tabMM3$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/bilateral/regular/r=0/Patrilocal_maleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY4 = read.table(path, header=TRUE)
  tabY4$kinship = rep("B_LP", length(tabY4$Gen))
  tabY4 = tabY4[tabY4$Gen == 100,]
  tabY4$sex = rep("male", length(tabY4$Gen))
  tabY4$descent = rep("bilateral", length(tabY4$Gen))
  tabY4$Pi = tabY4$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM4 = read.table(path, header=TRUE)
  tabM4$kinship = rep("B_LP", length(tabM4$Gen))
  tabM4 = tabM4[tabM4$Gen == 100,]
  tabM4$sex = rep("female", length(tabM4$Gen))
  tabM4$descent = rep("bilateral", length(tabM4$Gen))
  tabM4$Pi = tabM4$Pi/(2*5.5e-7)

  dir = paste0("Tables/Pi/bilateral/regular/r=0/Matrilocal_femaleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM4 = read.table(path, header=TRUE)
  tabYM4$kinship = rep("B_LM", length(tabYM4$Gen))
  tabYM4 = tabYM4[tabYM4$Gen == 100,]
  tabYM4$sex = rep("male", length(tabYM4$Gen))
  tabYM4$descent = rep("bilateral", length(tabYM4$Gen))
  tabYM4$Pi = tabYM4$Pi/(2*2.5e-8)

  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM4 = read.table(path, header=TRUE)
  tabMM4$kinship = rep("B_LM", length(tabMM4$Gen))
  tabMM4 = tabMM4[tabMM4$Gen == 100,]
  tabMM4$sex = rep("female", length(tabMM4$Gen))
  tabMM4$descent = rep("bilateral", length(tabMM4$Gen))
  tabMM4$Pi = tabMM4$Pi/(2*5.5e-7)

  tabY = merge(tabY1, tabY2, all = T)
  tabY = merge(tabY, tabY3, all = T)
  tabY = merge(tabY, tabY4, all = T)
  tabY = merge(tabY, tabYF1, all = T)
  tabY = merge(tabY, tabYMM1, all = T)
  tabY = merge(tabY, tabYF2, all = T)
  tabY = merge(tabY, tabYMM2, all = T)
  tabY = merge(tabY, tabYM1, all = T)
  tabY = merge(tabY, tabYM2, all = T)
  tabY = merge(tabY, tabYM3, all = T)
  tabY = merge(tabY, tabYM4, all = T)
  tabY = merge(tabY, tabYA3, all = T)

  tabM = merge(tabM1, tabM2, all=T)
  tabM = merge(tabM, tabM3, all=T)
  tabM = merge(tabM, tabM4, all=T)
  tabM = merge(tabM, tabMF1, all=T)
  tabM = merge(tabM, tabMF2, all=T)
  tabM = merge(tabM, tabMMM1, all=T)
  tabM = merge(tabM, tabMMM2, all=T)
  tabM = merge(tabM, tabMM1, all = T)
  tabM = merge(tabM, tabMM2, all = T)
  tabM = merge(tabM, tabMM3, all = T)
  tabM = merge(tabM, tabMM4, all = T)
  tabM = merge(tabM, tabMA3, all = T)

  alltab = merge(tabY, tabM, all=T)
  alltab$sex = factor(alltab$sex)

  alltab = rearrange_pitab(alltab)
  return(alltab)
}

shape_Pi_files_FT <- function(param, level) {
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=75/Patrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY1 = read.table(path, header=TRUE)
  tabY1$kinship = rep("SP", length(tabY1$Gen))
  tabY1 = tabY1[tabY1$Gen == 100,]
  tabY1$sex = rep("male", length(tabY1$Gen))
  tabY1$descent = rep("patrilineal", length(tabY1$Gen))
  tabY1$Pi = tabY1$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM1 = read.table(path, header=TRUE)
  tabM1$kinship = rep("SP", length(tabM1$Gen))
  tabM1 = tabM1[tabM1$Gen == 100,]
  tabM1$sex = rep("female", length(tabM1$Gen))
  tabM1$descent = rep("patrilineal", length(tabM1$Gen))
  tabM1$Pi = tabM1$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=75/Patrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYF1 = read.table(path, header=TRUE)
  tabYF1$kinship = rep("SP_F", length(tabYF1$Gen))
  tabYF1 = tabYF1[tabYF1$Gen == 100,]
  tabYF1$sex = rep("male", length(tabYF1$Gen))
  tabYF1$descent = rep("patrilineal", length(tabYF1$Gen))
  tabYF1$Pi = tabYF1$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMF1 = read.table(path, header=TRUE)
  tabMF1$kinship = rep("SP_F", length(tabMF1$Gen))
  tabMF1 = tabMF1[tabMF1$Gen == 100,]
  tabMF1$sex = rep("female", length(tabMF1$Gen))
  tabMF1$descent = rep("patrilineal", length(tabMF1$Gen))
  tabMF1$Pi = tabMF1$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=75/Matrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM1 = read.table(path, header=TRUE)
  tabYM1$kinship = rep("SM", length(tabYM1$Gen))
  tabYM1 = tabYM1[tabYM1$Gen == 100,]
  tabYM1$sex = rep("male", length(tabYM1$Gen))
  tabYM1$descent = rep("matrilineal", length(tabYM1$Gen))
  tabYM1$Pi = tabYM1$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM1 = read.table(path, header=TRUE)
  tabMM1$kinship = rep("SM", length(tabMM1$Gen))
  tabMM1 = tabMM1[tabMM1$Gen == 100,]
  tabMM1$sex = rep("female", length(tabMM1$Gen))
  tabMM1$descent = rep("matrilineal", length(tabMM1$Gen))
  tabMM1$Pi = tabMM1$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0.1/FT=75/Matrilineal_villages_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYMM1 = read.table(path, header=TRUE)
  tabYMM1$kinship = rep("SM_M", length(tabYMM1$Gen))
  tabYMM1 = tabYMM1[tabYMM1$Gen == 100,]
  tabYMM1$sex = rep("male", length(tabYMM1$Gen))
  tabYMM1$descent = rep("matrilineal", length(tabYMM1$Gen))
  tabYMM1$Pi = tabYMM1$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMMM1 = read.table(path, header=TRUE)
  tabMMM1$kinship = rep("SM_M", length(tabMMM1$Gen))
  tabMMM1 = tabMMM1[tabMMM1$Gen == 100,]
  tabMMM1$sex = rep("female", length(tabMMM1$Gen))
  tabMMM1$descent = rep("matrilineal", length(tabMMM1$Gen))
  tabMMM1$Pi = tabMMM1$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=75/Patrilineal_villages_rf_maleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY2 = read.table(path, header=TRUE)
  tabY2$kinship = rep("LP", length(tabY2$Gen))
  tabY2 = tabY2[tabY2$Gen == 100,]
  tabY2$sex = rep("male", length(tabY2$Gen))
  tabY2$descent = rep("patrilineal", length(tabY2$Gen))
  tabY2$Pi = tabY2$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM2 = read.table(path, header=TRUE)
  tabM2$kinship = rep("LP", length(tabM2$Gen))
  tabM2 = tabM2[tabM2$Gen == 100,]
  tabM2$sex = rep("female", length(tabM2$Gen))
  tabM2$descent = rep("patrilineal", length(tabM2$Gen))
  tabM2$Pi = tabM2$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=75/Patrilineal_villages_rf_maleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYF2 = read.table(path, header=TRUE)
  tabYF2$kinship = rep("LP_F", length(tabYF2$Gen))
  tabYF2 = tabYF2[tabYF2$Gen == 100,]
  tabYF2$sex = rep("male", length(tabYF2$Gen))
  tabYF2$descent = rep("patrilineal", length(tabYF2$Gen))
  tabYF2$Pi = tabYF2$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMF2 = read.table(path, header=TRUE)
  tabMF2$kinship = rep("LP_F", length(tabMF2$Gen))
  tabMF2 = tabMF2[tabMF2$Gen == 100,]
  tabMF2$sex = rep("female", length(tabMF2$Gen))
  tabMF2$descent = rep("patrilineal", length(tabMF2$Gen))
  tabMF2$Pi = tabMF2$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=75/Matrilineal_villages_rf_femaleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM2 = read.table(path, header=TRUE)
  tabYM2$kinship = rep("LM", length(tabYM2$Gen))
  tabYM2 = tabYM2[tabYM2$Gen == 100,]
  tabYM2$sex = rep("male", length(tabYM2$Gen))
  tabYM2$descent = rep("matrilineal", length(tabYM2$Gen))
  tabYM2$Pi = tabYM2$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM2 = read.table(path, header=TRUE)
  tabMM2$kinship = rep("LM", length(tabMM2$Gen))
  tabMM2 = tabMM2[tabMM2$Gen == 100,]
  tabMM2$sex = rep("female", length(tabMM2$Gen))
  tabMM2$descent = rep("matrilineal", length(tabMM2$Gen))
  tabMM2$Pi = tabMM2$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/unilineal/regular/r=0/k=0/FT=75/Matrilineal_villages_rf_femaleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYMM2 = read.table(path, header=TRUE)
  tabYMM2$kinship = rep("LM_M", length(tabYMM2$Gen))
  tabYMM2 = tabYMM2[tabYMM2$Gen == 100,]
  tabYMM2$sex = rep("male", length(tabYMM2$Gen))
  tabYMM2$descent = rep("matrilineal", length(tabYMM2$Gen))
  tabYMM2$Pi = tabYMM2$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMMM2 = read.table(path, header=TRUE)
  tabMMM2$kinship = rep("LM_M", length(tabMMM2$Gen))
  tabMMM2 = tabMMM2[tabMMM2$Gen == 100,]
  tabMMM2$sex = rep("female", length(tabMMM2$Gen))
  tabMMM2$descent = rep("matrilineal", length(tabMMM2$Gen))
  tabMMM2$Pi = tabMMM2$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/bilateral/regular/r=0/Patrilocal_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY3 = read.table(path, header=TRUE)
  tabY3$kinship = rep("B_SP", length(tabY3$Gen))
  tabY3 = tabY3[tabY3$Gen == 100,]
  tabY3$sex = rep("male", length(tabY3$Gen))
  tabY3$descent = rep("bilateral", length(tabY3$Gen))
  tabY3$Pi = tabY3$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM3 = read.table(path, header=TRUE)
  tabM3$kinship = rep("B_SP", length(tabM3$Gen))
  tabM3 = tabM3[tabM3$Gen == 100,]
  tabM3$sex = rep("female", length(tabM3$Gen))
  tabM3$descent = rep("bilateral", length(tabM3$Gen))
  tabM3$Pi = tabM3$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/bilateral/regular/r=0/Ambilocal_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYA3 = read.table(path, header=TRUE)
  tabYA3$kinship = rep("B_A", length(tabYA3$Gen))
  tabYA3 = tabYA3[tabYA3$Gen == 100,]
  tabYA3$sex = rep("male", length(tabYA3$Gen))
  tabYA3$descent = rep("bilateral", length(tabYA3$Gen))
  tabYA3$Pi = tabYA3$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMA3 = read.table(path, header=TRUE)
  tabMA3$kinship = rep("B_A", length(tabMA3$Gen))
  tabMA3 = tabMA3[tabMA3$Gen == 100,]
  tabMA3$sex = rep("female", length(tabMA3$Gen))
  tabMA3$descent = rep("bilateral", length(tabMA3$Gen))
  tabMA3$Pi = tabMA3$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/bilateral/regular/r=0/Matrilocal_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM3 = read.table(path, header=TRUE)
  tabYM3$kinship = rep("B_SM", length(tabYM3$Gen))
  tabYM3 = tabYM3[tabYM3$Gen == 100,]
  tabYM3$sex = rep("male", length(tabYM3$Gen))
  tabYM3$descent = rep("bilateral", length(tabYM3$Gen))
  tabYM3$Pi = tabYM3$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM3 = read.table(path, header=TRUE)
  tabMM3$kinship = rep("B_SM", length(tabMM3$Gen))
  tabMM3 = tabMM3[tabMM3$Gen == 100,]
  tabMM3$sex = rep("female", length(tabMM3$Gen))
  tabMM3$descent = rep("bilateral", length(tabMM3$Gen))
  tabMM3$Pi = tabMM3$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/bilateral/regular/r=0/Patrilocal_maleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabY4 = read.table(path, header=TRUE)
  tabY4$kinship = rep("B_LP", length(tabY4$Gen))
  tabY4 = tabY4[tabY4$Gen == 100,]
  tabY4$sex = rep("male", length(tabY4$Gen))
  tabY4$descent = rep("bilateral", length(tabY4$Gen))
  tabY4$Pi = tabY4$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabM4 = read.table(path, header=TRUE)
  tabM4$kinship = rep("B_LP", length(tabM4$Gen))
  tabM4 = tabM4[tabM4$Gen == 100,]
  tabM4$sex = rep("female", length(tabM4$Gen))
  tabM4$descent = rep("bilateral", length(tabM4$Gen))
  tabM4$Pi = tabM4$Pi/(2*5.5e-7)
  
  dir = paste0("Tables/Pi/bilateral/regular/r=0/Matrilocal_femaleM_", param, "relatedness/", level, "/")
  path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
  tabYM4 = read.table(path, header=TRUE)
  tabYM4$kinship = rep("B_LM", length(tabYM4$Gen))
  tabYM4 = tabYM4[tabYM4$Gen == 100,]
  tabYM4$sex = rep("male", length(tabYM4$Gen))
  tabYM4$descent = rep("bilateral", length(tabYM4$Gen))
  tabYM4$Pi = tabYM4$Pi/(2*2.5e-8)
  
  path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
  tabMM4 = read.table(path, header=TRUE)
  tabMM4$kinship = rep("B_LM", length(tabMM4$Gen))
  tabMM4 = tabMM4[tabMM4$Gen == 100,]
  tabMM4$sex = rep("female", length(tabMM4$Gen))
  tabMM4$descent = rep("bilateral", length(tabMM4$Gen))
  tabMM4$Pi = tabMM4$Pi/(2*5.5e-7)
  
  tabY = merge(tabY1, tabY2, all = T)
  tabY = merge(tabY, tabY3, all = T)
  tabY = merge(tabY, tabY4, all = T)
  tabY = merge(tabY, tabYF1, all = T)
  tabY = merge(tabY, tabYMM1, all = T)
  tabY = merge(tabY, tabYF2, all = T)
  tabY = merge(tabY, tabYMM2, all = T)
  tabY = merge(tabY, tabYM1, all = T)
  tabY = merge(tabY, tabYM2, all = T)
  tabY = merge(tabY, tabYM3, all = T)
  tabY = merge(tabY, tabYM4, all = T)
  tabY = merge(tabY, tabYA3, all = T)
  
  tabM = merge(tabM1, tabM2, all=T)
  tabM = merge(tabM, tabM3, all=T)
  tabM = merge(tabM, tabM4, all=T)
  tabM = merge(tabM, tabMF1, all=T)
  tabM = merge(tabM, tabMF2, all=T)
  tabM = merge(tabM, tabMMM1, all=T)
  tabM = merge(tabM, tabMMM2, all=T)
  tabM = merge(tabM, tabMM1, all = T)
  tabM = merge(tabM, tabMM2, all = T)
  tabM = merge(tabM, tabMM3, all = T)
  tabM = merge(tabM, tabMM4, all = T)
  tabM = merge(tabM, tabMA3, all = T)
  
  alltab = merge(tabY, tabM, all=T)
  alltab$sex = factor(alltab$sex)
  
  alltab = rearrange_pitab(alltab)
  return(alltab)
}

median_pi_ratios <- function(data, mod) {
  data <- data %>%
  group_by(descent, kinship, sex) %>%
  summarise(median = median(Pi, na.rm=T))

  data_ratio <- data %>%
  group_by(descent, kinship) %>%
  summarise(female_median = median[sex == "female"],
            male_median = median[sex == "male"]) %>%
  # Calculate the ratio of female over male median
  mutate(ratio = female_median - male_median) %>%
  # Select only relevant columns: descent, kinship, and the ratio
  select(descent, kinship, male_median, female_median, ratio)

  data_ratio$model = rep(mod, length(data_ratio$descent))
  return(data_ratio)
}

pi_village_base = shape_Pi_files("", "village")
pi_village_base = median_pi_ratios(pi_village_base, "baseline")

pi_village_200ind = shape_Pi_files("small_pop_", "village")
pi_village_200ind = median_pi_ratios(pi_village_200ind, "N=200")

pi_village_2v = shape_Pi_files("2villages_", "village")
pi_village_2v = median_pi_ratios(pi_village_2v, "Nv=2")

pi_village_0.5m = shape_Pi_files("low_migration_", "village")
pi_village_0.5m = median_pi_ratios(pi_village_0.5m, "m=0.1")

pi_village_FT = shape_Pi_files_FT("", "village")
pi_village_FT = median_pi_ratios(pi_village_FT, "FT=75")

df_V = merge(pi_village_base, pi_village_200ind, all=T)
df_V = merge(df_V, pi_village_2v, all=T)
df_V = merge(df_V, pi_village_0.5m, all=T)
df_V = merge(df_V, pi_village_FT, all=T)

p3 <- ggplot(df_V, aes(kinship, ratio)) +
  geom_hline(yintercept=0, linetype = "dashed", col = "black") +
  geom_path(col = "#9c0aa3", linewidth=0.5, group=1) +
  geom_point(col = "#9c0aa3", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Female - male  effective size difference (village)") +
  theme_bw()
p3

c <- ggplot(df_V) +
  geom_path(aes(kinship, female_median), col = "brown", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, female_median), col = "brown", size = 2) +
  geom_path(aes(kinship, male_median), col = "darkgoldenrod2", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, male_median), col = "darkgoldenrod2", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Median female and male effective size (village)") +
  theme_bw()
c

##### LOCAL
pi_local_base = shape_Pi_files("", "local")
pi_local_base = median_pi_ratios(pi_local_base, "baseline")

pi_local_200ind = shape_Pi_files("small_pop_", "local")
pi_local_200ind = median_pi_ratios(pi_local_200ind, "N=200")

pi_local_2v = shape_Pi_files("2villages_", "local")
pi_local_2v = median_pi_ratios(pi_local_2v, "Nv=2")

pi_local_0.5m = shape_Pi_files("low_migration_", "local")
pi_local_0.5m = median_pi_ratios(pi_local_0.5m, "m=0.1")

pi_local_FT = shape_Pi_files_FT("", "local")
pi_local_FT = median_pi_ratios(pi_local_FT, "FT=75")

df_L = merge(pi_local_base, pi_local_200ind, all=T)
df_L = merge(df_L, pi_local_2v, all=T)
df_L = merge(df_L, pi_local_0.5m, all=T)
df_L = merge(df_L, pi_local_FT, all=T)

p4 <- ggplot(df_L, aes(kinship, ratio)) +
  geom_hline(yintercept=0, linetype = "dashed", col = "black") +
  geom_path(col = "#9c0aa3", linewidth=0.5, group=1) +
  geom_point(col = "#9c0aa3", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Female - male  effective size difference (local group)") +
  theme_bw()
p4

d <- ggplot(df_L) +
  geom_path(aes(kinship, female_median), col = "brown", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, female_median), col = "brown", size = 2) +
  geom_path(aes(kinship, male_median), col = "darkgoldenrod2", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, male_median), col = "darkgoldenrod2", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Median female and male effective size (local group)") +
  theme_bw()
d

########### Haplogroup diversity #############
shape_hg_files <- function(param) {
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  sp = mean_hg(path, "SP", "village", "patrilineal")

  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  spf = mean_hg(path, "SP_F", "village_father", "patrilineal")

  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  sm = mean_hg(path, "SM", "village", "matrilineal")

  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  smm = mean_hg(path, "SM_M", "village_mother", "matrilineal")

  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lp = mean_hg(path, "LP", "village", "patrilineal")

  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lpf = mean_hg(path, "LP_F", "village_father", "patrilineal")

  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lmm = mean_hg(path, "LM_M", "village_mother", "matrilineal")

  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lm = mean_hg(path, "LM", "village", "matrilineal")

  dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
  path = paste(dir, "Matrilocal_", param, "relatedness", sep = "")
  bsm = mean_hg(path, "B_SM", "village", "bilateral")

  path = paste(dir, "Matrilocal_femaleM_", param, "relatedness", sep = "")
  blm = mean_hg(path, "B_LM", "village", "bilateral")

  path = paste(dir, "Ambilocal_", param, "relatedness", sep = "")
  ba = mean_hg(path, "B_A", "village", "bilateral")

  path = paste(dir, "Patrilocal_", param, "relatedness", sep = "")
  bsp = mean_hg(path, "B_SP", "village", "bilateral")

  path = paste(dir, "Patrilocal_maleM_", param, "relatedness", sep = "")
  blp = mean_hg(path, "B_LP", "village", "bilateral")

  df = merge(bsm, ba, all=T)
  df = merge(df, bsp, all = T)
  df = merge(df, blm, all = T)
  df = merge(df, blp, all = T)
  df = merge(df, sm, all = T)
  df = merge(df, smm, all = T)
  df = merge(df, sp, all = T)
  df = merge(df, spf, all = T)
  df = merge(df, lm, all = T)
  df = merge(df, lmm, all = T)
  df = merge(df, lp, all = T)
  df_village = merge(df, lpf, all = T)

  df_village$kinship = factor(df_village$kinship, levels = c("SM", "LM", "SM_M", "LM_M",
                                            "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                            "SP", "LP", "SP_F", "LP_F"))
  df_village$Generation = as.factor(df_village$Generation)

  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  sp = mean_hg(path, "SP", "patriline", "patrilineal")

  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  spf = mean_hg(path, "SP_F", "patriline_father", "patrilineal")

  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  sm = mean_hg(path, "SM", "matriline", "matrilineal")

  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  smm = mean_hg(path, "SM_M", "matriline_mother", "matrilineal")

  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lp = mean_hg(path, "LP", "patriline", "patrilineal")

  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lpf = mean_hg(path, "LP_F", "patriline_father", "patrilineal")

  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lmm = mean_hg(path, "LM_M", "matriline_mother", "matrilineal")

  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lm = mean_hg(path, "LM", "matriline", "matrilineal")

  dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
  path = paste(dir, "Matrilocal_", param, "relatedness", sep = "")
  bsm = mean_hg(path, "B_SM", "matriline", "bilateral")

  path = paste(dir, "Matrilocal_femaleM_", param, "relatedness", sep = "")
  blm = mean_hg(path, "B_LM", "matriline", "bilateral")

  path = paste(dir, "Ambilocal_", param, "relatedness", sep = "")
  ba = mean_hg(path, "B_A", "local", "bilateral")

  path = paste(dir, "Patrilocal_", param, "relatedness", sep = "")
  bsp = mean_hg(path, "B_SP", "patriline", "bilateral")

  path = paste(dir, "Patrilocal_maleM_", param, "relatedness", sep = "")
  blp = mean_hg(path, "B_LP", "patriline", "bilateral")

  df = merge(bsm, ba, all=T)
  df = merge(df, bsp, all = T)
  df = merge(df, blm, all = T)
  df = merge(df, blp, all = T)
  df = merge(df, sm, all = T)
  df = merge(df, smm, all = T)
  df = merge(df, sp, all = T)
  df = merge(df, spf, all = T)
  df = merge(df, lm, all = T)
  df = merge(df, lmm, all = T)
  df = merge(df, lp, all = T)
  df_local = merge(df, lpf, all = T)

  df_local$kinship = factor(df_local$kinship, levels = c("SM", "LM", "SM_M", "LM_M",
                                            "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                            "SP", "LP", "SP_F", "LP_F"))

  df_local$Generation = as.factor(df_local$Generation)
  return(list(df_village, df_local))
}

shape_hg_files_FT <- function(param) {
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=75/"
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  sp = mean_hg(path, "SP", "village", "patrilineal")
  
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  spf = mean_hg(path, "SP_F", "village_father", "patrilineal")
  
  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  sm = mean_hg(path, "SM", "village", "matrilineal")
  
  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  smm = mean_hg(path, "SM_M", "village_mother", "matrilineal")
  
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=75/"
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lp = mean_hg(path, "LP", "village", "patrilineal")
  
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lpf = mean_hg(path, "LP_F", "village_father", "patrilineal")
  
  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lmm = mean_hg(path, "LM_M", "village_mother", "matrilineal")
  
  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lm = mean_hg(path, "LM", "village", "matrilineal")
  
  dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
  path = paste(dir, "Matrilocal_", param, "relatedness", sep = "")
  bsm = mean_hg(path, "B_SM", "village", "bilateral")
  
  path = paste(dir, "Matrilocal_femaleM_", param, "relatedness", sep = "")
  blm = mean_hg(path, "B_LM", "village", "bilateral")
  
  path = paste(dir, "Ambilocal_", param, "relatedness", sep = "")
  ba = mean_hg(path, "B_A", "village", "bilateral")
  
  path = paste(dir, "Patrilocal_", param, "relatedness", sep = "")
  bsp = mean_hg(path, "B_SP", "village", "bilateral")
  
  path = paste(dir, "Patrilocal_maleM_", param, "relatedness", sep = "")
  blp = mean_hg(path, "B_LP", "village", "bilateral")
  
  df = merge(bsm, ba, all=T)
  df = merge(df, bsp, all = T)
  df = merge(df, blm, all = T)
  df = merge(df, blp, all = T)
  df = merge(df, sm, all = T)
  df = merge(df, smm, all = T)
  df = merge(df, sp, all = T)
  df = merge(df, spf, all = T)
  df = merge(df, lm, all = T)
  df = merge(df, lmm, all = T)
  df = merge(df, lp, all = T)
  df_village = merge(df, lpf, all = T)
  
  df_village$kinship = factor(df_village$kinship, levels = c("SM", "LM", "SM_M", "LM_M",
                                                             "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                                             "SP", "LP", "SP_F", "LP_F"))
  df_village$Generation = as.factor(df_village$Generation)
  
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=75/"
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  sp = mean_hg(path, "SP", "patriline", "patrilineal")
  
  path = paste(dir, "Patrilineal_villages_", param, "relatedness", sep = "")
  spf = mean_hg(path, "SP_F", "patriline_father", "patrilineal")
  
  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  sm = mean_hg(path, "SM", "matriline", "matrilineal")
  
  path = paste(dir, "Matrilineal_villages_", param, "relatedness", sep = "")
  smm = mean_hg(path, "SM_M", "matriline_mother", "matrilineal")
  
  dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=75/"
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lp = mean_hg(path, "LP", "patriline", "patrilineal")
  
  path = paste(dir, "Patrilineal_villages_rf_maleM_", param, "relatedness", sep = "")
  lpf = mean_hg(path, "LP_F", "patriline_father", "patrilineal")
  
  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lmm = mean_hg(path, "LM_M", "matriline_mother", "matrilineal")
  
  path = paste(dir, "Matrilineal_villages_rf_femaleM_", param, "relatedness", sep = "")
  lm = mean_hg(path, "LM", "matriline", "matrilineal")
  
  dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
  path = paste(dir, "Matrilocal_", param, "relatedness", sep = "")
  bsm = mean_hg(path, "B_SM", "matriline", "bilateral")
  
  path = paste(dir, "Matrilocal_femaleM_", param, "relatedness", sep = "")
  blm = mean_hg(path, "B_LM", "matriline", "bilateral")
  
  path = paste(dir, "Ambilocal_", param, "relatedness", sep = "")
  ba = mean_hg(path, "B_A", "local", "bilateral")
  
  path = paste(dir, "Patrilocal_", param, "relatedness", sep = "")
  bsp = mean_hg(path, "B_SP", "patriline", "bilateral")
  
  path = paste(dir, "Patrilocal_maleM_", param, "relatedness", sep = "")
  blp = mean_hg(path, "B_LP", "patriline", "bilateral")
  
  df = merge(bsm, ba, all=T)
  df = merge(df, bsp, all = T)
  df = merge(df, blm, all = T)
  df = merge(df, blp, all = T)
  df = merge(df, sm, all = T)
  df = merge(df, smm, all = T)
  df = merge(df, sp, all = T)
  df = merge(df, spf, all = T)
  df = merge(df, lm, all = T)
  df = merge(df, lmm, all = T)
  df = merge(df, lp, all = T)
  df_local = merge(df, lpf, all = T)
  
  df_local$kinship = factor(df_local$kinship, levels = c("SM", "LM", "SM_M", "LM_M",
                                                         "B_SM", "B_LM", "B_A", "B_SP", "B_LP", 
                                                         "SP", "LP", "SP_F", "LP_F"))
  
  df_local$Generation = as.factor(df_local$Generation)
  return(list(df_village, df_local))
}

median_hg_ratios <- function(data, mod) {
  data <- data %>%
  group_by(descent, kinship, Chromosome) %>%
  summarise(median = median(mean_hg_div, na.rm=T))

  data_ratio <- data %>%
  group_by(descent, kinship) %>%
  summarise(female_median = median[Chromosome == "M"],
            male_median = median[Chromosome == "Y"]) %>%
  # Calculate the ratio of female over male median
  mutate(ratio = female_median - male_median) %>%
  # Select only relevant columns: descent, kinship, and the ratio
  select(descent, kinship, male_median, female_median, ratio)

  data_ratio$model = rep(mod, length(data_ratio$descent))
  return(data_ratio)
}

df = shape_hg_files("")
data_village = df[[1]]
data_local = df[[2]]

data_village_base <- median_hg_ratios(data_village, "baseline")
data_local_base <- median_hg_ratios(data_local, "baseline")

df = shape_hg_files("small_pop_")
data_village = df[[1]]
data_local = df[[2]]

data_village_200ind <- median_hg_ratios(data_village, "N=200")
data_local_200ind <- median_hg_ratios(data_local, "N=200")

df = shape_hg_files("2villages_")
data_village = df[[1]]
data_local = df[[2]]

data_village_2v <- median_hg_ratios(data_village, "Nv=2")
data_local_2v <- median_hg_ratios(data_local, "Nv=2")

df = shape_hg_files("low_migration_")
data_village = df[[1]]
data_local = df[[2]]

data_village_m0.1 <- median_hg_ratios(data_village, "m=0.1")
data_local_m0.1 <- median_hg_ratios(data_local, "m=0.1")

df = shape_hg_files_FT("")
data_village = df[[1]]
data_local = df[[2]]

data_village_FT <- median_hg_ratios(data_village, "FT=75")
data_local_FT <- median_hg_ratios(data_local, "FT=75")

df_V = merge(data_village_base, data_village_200ind, all=T)
df_V = merge(df_V, data_village_2v, all=T)
df_V = merge(df_V, data_village_m0.1, all=T)
df_V = merge(df_V, data_village_FT, all=T)

p5 <- ggplot(df_V, aes(kinship, ratio)) +
  geom_hline(yintercept=0, linetype = "dashed", col = "black") +
  geom_path(col = "#9c0aa3", linewidth=0.5, group=1) +
  geom_point(col = "#9c0aa3", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Mitochondrial - Y haplogroup diversity difference (village)") +
  theme_bw()
p5

e <- ggplot(df_V) +
  geom_path(aes(kinship, female_median), col = "brown", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, female_median), col = "brown", size = 2) +
  geom_path(aes(kinship, male_median), col = "darkgoldenrod2", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, male_median), col = "darkgoldenrod2", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Median mitochondrial and Y chromosome haplogroup diversity (village)") +
  theme_bw()
e

####################### LOCAL ###########################

df_L = merge(data_local_base, data_local_200ind, all=T)
df_L = merge(df_L, data_local_2v, all=T)
df_L = merge(df_L, data_local_m0.1, all=T)
df_L = merge(df_L, data_local_FT, all=T)

p6 <- ggplot(df_L, aes(kinship, ratio)) +
  geom_hline(yintercept=0, linetype = "dashed", col = "black") +
  geom_path(col = "#9c0aa3", linewidth=0.5, group=1) +
  geom_point(col = "#9c0aa3", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Mitochondrial - Y haplogroup diversity difference (local group)") +
  theme_bw()
p6

f <- ggplot(df_L) +
  geom_path(aes(kinship, female_median), col = "brown", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, female_median), col = "brown", size = 2) +
  geom_path(aes(kinship, male_median), col = "darkgoldenrod2", linewidth = 0.5, group=1) +
  geom_point(aes(kinship, male_median), col = "darkgoldenrod2", size = 2) +
  facet_grid(model~descent, scales='free_x') +
  labs(x = "Kinship systems", y = "Median mitochondrial and Y chromosome haplogroup diversity (local group)") +
  theme_bw()
f

###### Figure 4 ######

p <- plot_grid(a + theme(legend.position = "none"),
               b + theme(legend.position = "none"),
               e + theme(legend.position = "none"),
               f + theme(legend.position = "none"),
               ncol = 2, labels = "auto")

png(file = "figure_4.png", width = 4000, height = 4000, res = 300)
p
dev.off()

###### Figure S7 ######
p <- plot_grid(c + theme(legend.position = "none"),
               d + theme(legend.position = "none"),
               ncol = 2, labels = "auto")

png(file = "figure_S7.png", width = 4000, height = 2000, res = 300)
p
dev.off()
