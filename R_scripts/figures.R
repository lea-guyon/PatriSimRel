library(cowplot)
library(ggplot2)
library(scales)
library(boot)
library(xtable)
library(reshape2)
library(ggpubr)
library(tidyr)
library(stringr)
library(rstatix)
library(ggprism)
library(dplyr)
library("ggpattern")

working_dir = "~"

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

reshape_file <- function(path, model, mode, sample=F) {
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
  
  if (mode == "local_group") {
    file = read.table("relatedness_patriline.txt", header = T)
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
  
  reshaped_file <- pivot_longer(file, 
                                cols = var,
                                names_to = "variable",
                                values_to = "value")
  reshaped_file$category = str_sub(reshaped_file$variable, end=-2)
  reshaped_file$kinship = rep(model, length(reshaped_file$variable))
  reshaped_file$kinvar = paste(reshaped_file$kinship, reshaped_file$variable)
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

setwd(working_dir)

dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0.1/FT=150/"
path = paste(dir, "strict_patrilineal_villages", sep = "")
reshape_file_1 = reshape_file(path, "PVRS", "village")
reshape_file_father_1 = reshape_file(path, "PVRS_F", "village_father")
reshape_group_1 = reshape_file(path, "PVRS", "group")
reshape_sample_1 = reshape_file(path, "PVRS", "village", T)
reshape_group_sample_1 = reshape_file(path, "PVRS", "group", T)
reshape_patriline_1 = reshape_file(path, "PVRS", "local_group")
reshape_patriline_father_1 = reshape_file(path, "PVRS_F", "local_group_father")

dir = "Tables/Simulations_metrics/unilineal/regular/r=0/k=0/FT=150/"
path = paste(dir, "loose_patrilineal_villages", sep = "")
reshape_file_2 = reshape_file(path, "P", "village")
reshape_file_father_2 = reshape_file(path, "P_F", "village_father")
reshape_group_2 = reshape_file(path, "P", "group")
reshape_sample_2 = reshape_file(path, "P", "village", T)
reshape_group_sample_2 = reshape_file(path, "P", "group", T)
reshape_patriline_2 = reshape_file(path, "P", "local_group")
reshape_patriline_father_2 = reshape_file(path, "P_F", "local_group_father")

dir = "Tables/Simulations_metrics/bilateral/regular/r=0/"
path = paste(dir, "strict_patrilocal_villages", sep = "")
reshape_file_3 = reshape_file(path, "B", "village")
reshape_sample_3 = reshape_file(path, "B", "village", T)
reshape_patriline_3 = reshape_file(path, "B", "local_group")

path = paste(dir, "loose_patrilocal_villages", sep = "")
reshape_file_4 = reshape_file(path, "Bm", "village")
reshape_sample_4 = reshape_file(path, "Bm", "village", T)
reshape_patriline_4 = reshape_file(path, "Bm", "local_group")

##################################################################
########################### MAIN #################################
##################################################################

########################## VILLAGE ###############################

reshape = merge(reshape_file_1, reshape_file_2, all = T)
reshape2 = merge(reshape, reshape_file_3, all = T)
reshape2 = merge(reshape2, reshape_file_4, all = T)
reshape2 <- reshape2[, -seq(3,8)]
reshape2 = merge(reshape2, reshape_file_father_1, all = T)
reshape2 = merge(reshape2, reshape_file_father_2, all = T)

# Rearrange reshape2
relabel.variable = as_labeller(c(meanRelM = "male", 
                                 meanRelF = "female", 
                                 relFirstM = "male", relFirstF = "female",
                                 relFreqM = "male", relFreqF = "female"))
category.label =  c("mean rel", "freq 1st degree rel", "freq 1st, 2nd and 3rd degree rel")
names(category.label) = c("meanRel", "relFirst", "relFreq")
reshape_bis = data.frame("variable"=reshape2$variable,"value"=reshape2$value, 
                         "category"=reshape2$category, "kinship"=reshape2$kinship)

reshape2 <- reshape2 %>% 
  as_tibble() %>%
  rowwise() %>% 
  mutate(variable2 = factor(relabel.variable(variable))) %>%
  arrange(value) %>%
  mutate(kinship = factor(kinship, 
                          levels=c("B", "Bm", "PVRS", "P", "PVRS_F", "P_F")))

reshape2 = reshape2[with(reshape2, order(variable2, kinship)),]
                          
data = reshape2[reshape2$category == "meanRel" & reshape2$Generation == 20100,]
data <- data %>%
  group_by(kinship, variable2, Replicat) %>%
  summarise(value = mean(value))
data = data[with(data, order(variable2, kinship)),]

data_stat = data
data_stat$var = factor(paste(data_stat$kinship, data_stat$variable2),
                levels = c("B female", "Bm female", "PVRS female", "P female", "PVRS_F female", "P_F female",
  "B male", "Bm male", "PVRS male", "P male", "PVRS_F male", "P_F male"),
                labels = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male"))

pval_matrix = pairwise.wilcox.test(data_stat$value, data_stat$var, p.adjust.method="bonferroni")$p.value

# Define the significance levels
significance_matrix <- matrix("ns", nrow=nrow(pval_matrix), ncol=ncol(pval_matrix)) # Start with "not significant"

# Apply thresholds
significance_matrix[is.na(pval_matrix)] <- NA
significance_matrix[pval_matrix < 0.05] <- "*"
significance_matrix[pval_matrix < 0.01] <- "**"
significance_matrix[pval_matrix < 0.001] <- "***"
significance_matrix[pval_matrix < 0.0001] <- "****"
significance_matrix = rbind("B_SP female"=rep(NA, 11), significance_matrix)
significance_matrix = cbind(significance_matrix, "L_PF male" = rep(NA, 12))

colnames(significance_matrix) = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male")
rownames(significance_matrix) = colnames(significance_matrix)

df <- melt(significance_matrix)
colnames(df) <- c("x", "y", "value")

# Create the heatmap
h1 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  geom_vline(xintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("#cc1a1a", "navyblue"), na.value="white") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_text(angle = 90))

p1 <- ggplot(data) +
  geom_point(aes(variable2, value, col=kinship), size=2, alpha=0.1, position = position_jitterdodge()) +
  geom_boxplot_pattern(aes(variable2, value, fill = kinship, col = kinship, 
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship), 
                       alpha = 0.2, lwd=0.75, outlier.shape=NA) +
  #stat_compare_means(aes(variable2, value), label = "p.signif", label.x = 1.45, label.y = 0.0029) +
  scale_fill_manual(labels = c("Bilateral descent and strict patrilocal residence",
                               "Bilateral descent and loose patrilocal residence",
                               "Strict patrilineal descent, burial in husband's cemetery", 
                               "Loose patrilineal descent, burial in husband's cemetery",
                               "Strict patrilineal descent, burial in father's cemetery", 
                               "Loose patrilineal descent, burial in father's cemetery"),
                    values = c("darkcyan", "#98bbbb",
                               "darkgoldenrod1", "#FFD56F",
                               "darkgoldenrod1", "bisque1")) +
  scale_color_manual(labels = c("Bilateral descent and strict patrilocal residence",
                                "Bilateral descent and loose patrilocal residence",
                                "Strict patrilineal descent, burial in husband's cemetery", 
                                "Loose patrilineal descent, burial in husband's cemetery",
                                "Strict patrilineal descent, burial in father's cemetery", 
                                "Loose patrilineal descent, burial in father's cemetery"),
                     values = c("#027d7d", "#7fb8b8",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(labels = c("Bilateral descent and strict patrilocal residence",
                                  "Bilateral descent and loose patrilocal residence",
                                  "Strict patrilineal descent, burial in husband's cemetery", 
                                  "Loose patrilineal descent, burial in husband's cemetery",
                                  "Strict patrilineal descent, burial in father's cemetery", 
                                  "Loose patrilineal descent, burial in father's cemetery"),
                       values= c("stripe", "none", "crosshatch", "none", 
                                 "crosshatch", "none")) +
  scale_pattern_density_manual(labels = c("Bilateral descent and strict patrilocal residence",
                                          "Bilateral descent and loose patrilocal residence",
                                          "Strict patrilineal descent, burial in husband's cemetery", 
                                          "Loose patrilineal descent, burial in husband's cemetery",
                                          "Strict patrilineal descent, burial in father's cemetery", 
                                          "Loose patrilineal descent, burial in father's cemetery"),
                       values= c(0.35, 0, 0.1, 0, 0.35, 0)) +
  scale_pattern_fill_manual(labels = c("Bilateral descent and strict patrilocal residence",
                                       "Bilateral descent and loose patrilocal residence",
                                       "Strict patrilineal descent, burial in husband's cemetery", 
                                       "Loose patrilineal descent, burial in husband's cemetery",
                                       "Strict patrilineal descent, burial in father's cemetery", 
                                       "Loose patrilineal descent, burial in father's cemetery"),
                            values=c("#98bbbb", "#98bbbb", "darkgoldenrod1", 
                                     "bisque1", "bisque1", "bisque1")) +
  xlab("sex") +
  ylab("mean relatedness (village)") +
  labs(fill = "Kinship system", col = "Kinship system",
       pattern_fill = "Kinship system", pattern = "Kinship system",
       pattern_density= "Kinship system") +
  guides(fill=guide_legend(nrow=3, theme = theme(legend.byrow = TRUE)), 
         col=guide_legend(nrow=3, theme = theme(legend.byrow = TRUE)),
         pattern=guide_legend(nrow=3, theme = theme(legend.byrow = TRUE)), 
         pattern_fill=guide_legend(nrow=3, theme = theme(legend.byrow = TRUE)),
         pattern_density=guide_legend(nrow=3, theme = theme(legend.byrow = TRUE))) +
  theme_bw()+
  theme(text = element_text(size = 12), legend.key.size = unit(0.9, 'cm'))

mean = c()
sd=c()
sex = c()
kinship = c()
var2=c("female", "male")
kin=c("B", "Bm", "PVRS", "P", "PVRS_F", "P_F")
for (s in var2) {
  data_s = data[data$variable2 == s,]
  for (k in kin) {
    data_k = data_s[data_s$kinship == k,]
    mean = c(mean, mean(data_k$value))
    sd = c(sd, sd(data_k$value))
    sex = c(sex, s)
    kinship = c(kinship, k)
  }
}
mean_sd <- data.frame(sex=sex , kinship = kinship, mean = mean, sd=format(sd, scientific=F))
matrix(mean_sd$mean, ncol=2)
matrix(mean_sd$sd, ncol=2)

####################### LOCAL ###########################
reshape = merge(reshape_patriline_1, reshape_patriline_2, all = T)
reshape2 = merge(reshape, reshape_patriline_3, all = T)
reshape2 = merge(reshape2, reshape_patriline_4, all = T)
reshape2 = merge(reshape2, reshape_patriline_father_1, all = T)
reshape2 = merge(reshape2, reshape_patriline_father_2, all = T)

# Rearrange reshape2
relabel.variable = as_labeller(c(meanRelM = "male", 
                                 meanRelF = "female", 
                                 relFirstM = "male", relFirstF = "female",
                                 relFreqM = "male", relFreqF = "female"))
category.label =  c("mean rel", "freq 1st degree rel", "freq 1st, 2nd and 3rd degree rel")
names(category.label) = c("meanRel", "relFirst", "relFreq")

reshape2 <- reshape2 %>% 
  as_tibble() %>%
  rowwise() %>% 
  mutate(variable2 = factor(relabel.variable(variable))) %>%
  arrange(value) %>%
  mutate(kinship = factor(kinship, 
                          levels=c("B", "Bm", "PVRS", "P", "PVRS_F", "P_F")))
reshape2 = reshape2[with(reshape2, order(variable2, kinship)),]
                          
data = reshape2[reshape2$category == "meanRel" & reshape2$Generation == 20100,]
data <- data %>%
  group_by(kinship, variable2, Replicat) %>%
  summarise(value = mean(value))
data = data[with(data, order(variable2, kinship)),]

data_stat = data
data_stat$var = factor(paste(data_stat$kinship, data_stat$variable2),
                levels = c("B female", "Bm female", "PVRS female", "P female", "PVRS_F female", "P_F female",
  "B male", "Bm male", "PVRS male", "P male", "PVRS_F male", "P_F male"),
                labels = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male"))

pval_matrix = pairwise.wilcox.test(data_stat$value, data_stat$var, p.adjust.method="bonferroni")$p.value

# Define the significance levels
significance_matrix <- matrix("ns", nrow=nrow(pval_matrix), ncol=ncol(pval_matrix)) # Start with "not significant"

# Apply thresholds
significance_matrix[is.na(pval_matrix)] <- NA
significance_matrix[pval_matrix < 0.05] <- "*"
significance_matrix[pval_matrix < 0.01] <- "**"
significance_matrix[pval_matrix < 0.001] <- "***"
significance_matrix[pval_matrix < 0.0001] <- "****"
significance_matrix = rbind("B_SP female"=rep(NA, 11), significance_matrix)
significance_matrix = cbind(significance_matrix, "L_PF male" = rep(NA, 12))

colnames(significance_matrix) = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male")
rownames(significance_matrix) = colnames(significance_matrix)

df <- melt(significance_matrix)
colnames(df) <- c("x", "y", "value")

h2 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  geom_vline(xintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("#cc1a1a", "navyblue"), na.value="white") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))

p3 <- ggplot(data) +
  geom_point(aes(variable2, value, col=kinship), size=2, alpha=0.1, position = position_jitterdodge()) +
  geom_boxplot_pattern(aes(variable2, value, fill = kinship, col = kinship, 
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship), 
                       alpha = 0.2, lwd=0.75, outlier.shape=NA) +
  #stat_compare_means(aes(variable2, value), label = "p.signif", label.x = 1.45, label.y = 0.85) +
  scale_fill_manual(values = c("darkcyan", "#98bbbb",
                               "darkgoldenrod1", "#FFD56F",
                               "darkgoldenrod1", "bisque1")) +
  scale_color_manual(values = c("#027d7d", "#7fb8b8",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(values= c("stripe", "none", "crosshatch", "none", 
                                 "crosshatch", "none")) +
  scale_pattern_density_manual(values= c(0.35, 0, 0.1, 0, 0.35, 0)) +
  scale_pattern_fill_manual(values=c("#98bbbb", "#98bbbb", "darkgoldenrod1", 
                                     "bisque1", "bisque1", "bisque1")) +
  xlab("sex") +
  ylab("mean relatedness (local group)") +
  theme_bw()+
  theme(text = element_text(size = 12))

mean = c()
sd=c()
sex = c()
kinship = c()
for (s in var2) {
  data_s = data[data$variable2 == s,]
  for (k in kin) {
    data_k = data_s[data_s$kinship == k,]
    mean = c(mean, mean(data_k$value))
    sd = c(sd, sd(data_k$value))
    sex = c(sex, s)
    kinship = c(kinship, k)
  }
}
mean_sd <- data.frame(sex=sex , kinship = kinship, mean = mean, sd=format(sd, scientific=F))
matrix(mean_sd$mean, ncol=2)
matrix(mean_sd$sd, ncol=2)

######################## Nucleotide diversity #########################
##### VILLAGE

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY1 = read.table(path, header=TRUE)
tabY1$kinship = rep("PVRS", length(tabY1$Gen))
tabY1 = tabY1[tabY1$Gen == 100,]
tabY1$sex = rep("male", length(tabY1$Gen))
tabY1$Pi = tabY1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM1 = read.table(path, header=TRUE)
tabM1$kinship = rep("PVRS", length(tabM1$Gen))
tabM1 = tabM1[tabM1$Gen == 100,]
tabM1$sex = rep("female", length(tabM1$Gen))
tabM1$Pi = tabM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/village_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF1 = read.table(path, header=TRUE)
tabYF1$kinship = rep("PVRS_F", length(tabYF1$Gen))
tabYF1 = tabYF1[tabYF1$Gen == 100,]
tabYF1$sex = rep("male", length(tabYF1$Gen))
tabYF1$Pi = tabYF1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF1 = read.table(path, header=TRUE)
tabMF1$kinship = rep("PVRS_F", length(tabMF1$Gen))
tabMF1 = tabMF1[tabMF1$Gen == 100,]
tabMF1$sex = rep("female", length(tabMF1$Gen))
tabMF1$Pi = tabMF1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY2 = read.table(path, header=TRUE)
tabY2$kinship = rep("P", length(tabY2$Gen))
tabY2 = tabY2[tabY2$Gen == 100,]
tabY2$sex = rep("male", length(tabY2$Gen))
tabY2$Pi = tabY2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM2 = read.table(path, header=TRUE)
tabM2$kinship = rep("P", length(tabM2$Gen))
tabM2 = tabM2[tabM2$Gen == 100,]
tabM2$sex = rep("female", length(tabM2$Gen))
tabM2$Pi = tabM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/village_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF2 = read.table(path, header=TRUE)
tabYF2$kinship = rep("P_F", length(tabYF2$Gen))
tabYF2 = tabYF2[tabYF2$Gen == 100,]
tabYF2$sex = rep("male", length(tabYF2$Gen))
tabYF2$Pi = tabYF2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF2 = read.table(path, header=TRUE)
tabMF2$kinship = rep("P_F", length(tabMF2$Gen))
tabMF2 = tabMF2[tabMF2$Gen == 100,]
tabMF2$sex = rep("female", length(tabMF2$Gen))
tabMF2$Pi = tabMF2$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/strict_patrilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY3 = read.table(path, header=TRUE)
tabY3$kinship = rep("B", length(tabY3$Gen))
tabY3 = tabY3[tabY3$Gen == 100,]
tabY3$sex = rep("male", length(tabY3$Gen))
tabY3$Pi = tabY3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM3 = read.table(path, header=TRUE)
tabM3$kinship = rep("B", length(tabM3$Gen))
tabM3 = tabM3[tabM3$Gen == 100,]
tabM3$sex = rep("female", length(tabM3$Gen))
tabM3$Pi = tabM3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/loose_patrilocal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY4 = read.table(path, header=TRUE)
tabY4$kinship = rep("Bm", length(tabY4$Gen))
tabY4 = tabY4[tabY4$Gen == 100,]
tabY4$sex = rep("male", length(tabY4$Gen))
tabY4$Pi = tabY4$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM4 = read.table(path, header=TRUE)
tabM4$kinship = rep("Bm", length(tabM4$Gen))
tabM4 = tabM4[tabM4$Gen == 100,]
tabM4$sex = rep("female", length(tabM4$Gen))
tabM4$Pi = tabM4$Pi/(2*5.5e-7)

tabY = merge(tabY1, tabY2, all = T)
tabY = merge(tabY, tabY3, all = T)
tabY = merge(tabY, tabY4, all = T)
tabY = merge(tabY, tabYF1, all = T)
tabY = merge(tabY, tabYF2, all = T)

tabM = merge(tabM1, tabM2, all=T)
tabM = merge(tabM, tabM3, all=T)
tabM = merge(tabM, tabM4, all=T)
tabM = merge(tabM, tabMF1, all=T)
tabM = merge(tabM, tabMF2, all=T)

alltab = merge(tabY, tabM, all=T)
alltab$sex = as.factor(alltab$sex)

alltab = data.frame("Pi" = alltab$Pi, "kinship" = as.factor(alltab$kinship), "sex" = alltab$sex)
alltab = na.omit(alltab)
alltab <- alltab %>%
  arrange(Pi) %>%
  mutate(kinship = factor(kinship, levels=c("B", "Bm", "PVRS", "P", "PVRS_F", "P_F")))

alltab = alltab[with(alltab, order(sex, kinship)),]

data_stat = alltab
data_stat$var = factor(paste(data_stat$kinship, data_stat$sex),
                levels = c("B female", "Bm female", "PVRS female", "P female", "PVRS_F female", "P_F female",
  "B male", "Bm male", "PVRS male", "P male", "PVRS_F male", "P_F male"),
                labels = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male"))

pval_matrix = pairwise.wilcox.test(data_stat$Pi, data_stat$var, p.adjust.method="bonferroni")$p.value

# Define the significance levels
significance_matrix <- matrix("ns", nrow=nrow(pval_matrix), ncol=ncol(pval_matrix)) # Start with "not significant"

# Apply thresholds
significance_matrix[is.na(pval_matrix)] <- NA
significance_matrix[pval_matrix < 0.05] <- "*"
significance_matrix[pval_matrix < 0.01] <- "**"
significance_matrix[pval_matrix < 0.001] <- "***"
significance_matrix[pval_matrix < 0.0001] <- "****"
significance_matrix = rbind("B_SP female"=rep(NA, 11), significance_matrix)
significance_matrix = cbind(significance_matrix, "L_PF male" = rep(NA, 12))

colnames(significance_matrix) = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male")
rownames(significance_matrix) = colnames(significance_matrix)

df <- melt(significance_matrix)
colnames(df) <- c("x", "y", "value")
df$value = factor(df$value)

h3 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  geom_vline(xintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("#cc1a1a", "navyblue"), na.value="white",
                    labels = c("p < 0.0001", "ns"), na.translate = F) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))

p4 <- ggplot(alltab) +
  geom_point(aes(sex, Pi, col=kinship), size=2, alpha=0.1, position = position_jitterdodge()) +
  geom_boxplot_pattern(aes(sex, Pi, fill = kinship, col = kinship,
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship), 
                       alpha = 0.2, lwd=0.75, outlier.shape=NA) +
  scale_fill_manual(values = c("darkcyan", "#98bbbb",
                               "darkgoldenrod1", "#FFD56F",
                               "darkgoldenrod1", "bisque1")) +
  scale_color_manual(values = c("#027d7d", "#7fb8b8",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(values= c("stripe", "none", "crosshatch", "none", 
                                 "crosshatch", "none")) +
  scale_pattern_density_manual(values= c(0.35, 0, 0.1, 0, 0.35, 0)) +
  scale_pattern_fill_manual(values=c("#98bbbb", "#98bbbb", "darkgoldenrod1", 
                                     "bisque1", "bisque1", "bisque1")) +
  ylab("Effective population size (village)") +
  xlab("sex") +
  theme_bw() +
  theme(text = element_text(size = 12))

mean = c()
sd=c()
sex = c()
kinship = c()
for (s in var2) {
  data_s = alltab[alltab$sex == s,]
  for (k in kin) {
    data_k = data_s[data_s$kinship == k,]
    mean = c(mean, mean(data_k$Pi))
    sd = c(sd, sd(data_k$Pi))
    sex = c(sex, s)
    kinship = c(kinship, k)
  }
}
mean_sd <- data.frame(sex=sex , kinship = kinship, mean = mean, sd=format(sd, scientific=F))
matrix(mean_sd$mean, ncol=2)
matrix(mean_sd$sd, ncol=2)

##### LOCAL
dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY1 = read.table(path, header=TRUE)
tabY1$kinship = rep("PVRS", length(tabY1$Gen))
tabY1 = tabY1[tabY1$Gen == 100,]
tabY1$sex = rep("male", length(tabY1$Gen))
tabY1$Pi = tabY1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM1 = read.table(path, header=TRUE)
tabM1$kinship = rep("PVRS", length(tabM1$Gen))
tabM1 = tabM1[tabM1$Gen == 100,]
tabM1$sex = rep("female", length(tabM1$Gen))
tabM1$Pi = tabM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/local_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF1 = read.table(path, header=TRUE)
tabYF1$kinship = rep("PVRS_F", length(tabYF1$Gen))
tabYF1 = tabYF1[tabYF1$Gen == 100,]
tabYF1$sex = rep("male", length(tabYF1$Gen))
tabYF1$Pi = tabYF1$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF1 = read.table(path, header=TRUE)
tabMF1$kinship = rep("PVRS_F", length(tabMF1$Gen))
tabMF1 = tabMF1[tabMF1$Gen == 100,]
tabMF1$sex = rep("female", length(tabMF1$Gen))
tabMF1$Pi = tabMF1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY2 = read.table(path, header=TRUE)
tabY2$kinship = rep("P", length(tabY2$Gen))
tabY2 = tabY2[tabY2$Gen == 100,]
tabY2$sex = rep("male", length(tabY2$Gen))
tabY2$Pi = tabY2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM2 = read.table(path, header=TRUE)
tabM2$kinship = rep("P", length(tabM2$Gen))
tabM2 = tabM2[tabM2$Gen == 100,]
tabM2$sex = rep("female", length(tabM2$Gen))
tabM2$Pi = tabM2$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/regular/r=0/k=0/FT=150/loose_patrilineal_villages/local_father/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabYF2 = read.table(path, header=TRUE)
tabYF2$kinship = rep("P_F", length(tabYF2$Gen))
tabYF2 = tabYF2[tabYF2$Gen == 100,]
tabYF2$sex = rep("male", length(tabYF2$Gen))
tabYF2$Pi = tabYF2$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabMF2 = read.table(path, header=TRUE)
tabMF2$kinship = rep("P_F", length(tabMF2$Gen))
tabMF2 = tabMF2[tabMF2$Gen == 100,]
tabMF2$sex = rep("female", length(tabMF2$Gen))
tabMF2$Pi = tabMF2$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/strict_patrilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY3 = read.table(path, header=TRUE)
tabY3$kinship = rep("B", length(tabY3$Gen))
tabY3 = tabY3[tabY3$Gen == 100,]
tabY3$sex = rep("male", length(tabY3$Gen))
tabY3$Pi = tabY3$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM3 = read.table(path, header=TRUE)
tabM3$kinship = rep("B", length(tabM3$Gen))
tabM3 = tabM3[tabM3$Gen == 100,]
tabM3$sex = rep("female", length(tabM3$Gen))
tabM3$Pi = tabM3$Pi/(2*5.5e-7)

dir = "Tables/Pi/bilateral/regular/r=0/loose_patrilocal_villages/local/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY4 = read.table(path, header=TRUE)
tabY4$kinship = rep("Bm", length(tabY4$Gen))
tabY4 = tabY4[tabY4$Gen == 100,]
tabY4$sex = rep("male", length(tabY4$Gen))
tabY4$Pi = tabY4$Pi/(2*2.5e-8)

path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM4 = read.table(path, header=TRUE)
tabM4$kinship = rep("Bm", length(tabM4$Gen))
tabM4 = tabM4[tabM4$Gen == 100,]
tabM4$sex = rep("female", length(tabM4$Gen))
tabM4$Pi = tabM4$Pi/(2*5.5e-7)

tabY = merge(tabY1, tabY2, all = T)
tabY = merge(tabY, tabY3, all = T)
tabY = merge(tabY, tabY4, all = T)
tabY = merge(tabY, tabYF1, all = T)
tabY = merge(tabY, tabYF2, all = T)

tabM = merge(tabM1, tabM2, all=T)
tabM = merge(tabM, tabM3, all=T)
tabM = merge(tabM, tabM4, all=T)
tabM = merge(tabM, tabMF1, all=T)
tabM = merge(tabM, tabMF2, all=T)

alltab = merge(tabY, tabM, all=T)
alltab$sex = as.factor(alltab$sex)

alltab = data.frame("Pi" = alltab$Pi, "kinship" = as.factor(alltab$kinship), "sex" = alltab$sex)
alltab = na.omit(alltab)
alltab <- alltab %>%
  arrange(Pi) %>%
  mutate(kinship = factor(kinship, levels=c("B", "Bm", "PVRS", "P", "PVRS_F", "P_F")))

alltab = alltab[with(alltab, order(sex, kinship)),]

data_stat = alltab
data_stat$var = factor(paste(data_stat$kinship, data_stat$sex),
                levels = c("B female", "Bm female", "PVRS female", "P female", "PVRS_F female", "P_F female",
  "B male", "Bm male", "PVRS male", "P male", "PVRS_F male", "P_F male"),
                labels = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male"))

pval_matrix = pairwise.wilcox.test(data_stat$Pi, data_stat$var, p.adjust.method="bonferroni")$p.value

# Define the significance levels
significance_matrix <- matrix("ns", nrow=nrow(pval_matrix), ncol=ncol(pval_matrix)) # Start with "not significant"

# Apply thresholds
significance_matrix[is.na(pval_matrix)] <- NA
significance_matrix[pval_matrix < 0.05] <- "*"
significance_matrix[pval_matrix < 0.01] <- "**"
significance_matrix[pval_matrix < 0.001] <- "***"
significance_matrix[pval_matrix < 0.0001] <- "****"
significance_matrix = rbind("B_SP female"=rep(NA, 11), significance_matrix)
significance_matrix = cbind(significance_matrix, "L_PF male" = rep(NA, 12))

colnames(significance_matrix) = c("B_SP female", "B_LP female", "SP female", "LP female", "SP_F female", "LP_F female",
  "B_SP male", "B_LP male", "SP male", "LP male", "SP_F male", "LP_F male")
rownames(significance_matrix) = colnames(significance_matrix)

df <- melt(significance_matrix)
colnames(df) <- c("x", "y", "value")

h4 <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            alpha = 0.7) +
  geom_hline(yintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  geom_vline(xintercept=6.5, color = "black", linetype = "dashed", linewidth=2) +
  labs(x="", y="", fill = "p-value significance") +
  scale_fill_manual(values = c("#cc1a1a", "navyblue"), na.value="white") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))

p5 <- ggplot(alltab) +
  geom_point(aes(sex, Pi, col=kinship), size=2, alpha=0.1, position = position_jitterdodge()) +
  geom_boxplot_pattern(aes(sex, Pi, fill = kinship, col = kinship,
                           pattern = kinship, pattern_fill = kinship,
                           pattern_density = kinship), 
                       alpha = 0.2, lwd=0.75, outlier.shape=NA) +
  scale_fill_manual(values = c("darkcyan", "#98bbbb",
                               "darkgoldenrod1", "#FFD56F",
                               "darkgoldenrod1", "bisque1")) +
  scale_color_manual(values = c("#027d7d", "#7fb8b8",
                                "darkgoldenrod2", "#fdd16b",
                                "darkgoldenrod2", "#eed7ba")) +
  scale_pattern_manual(values= c("stripe", "none", "crosshatch", "none", 
                                 "crosshatch", "none")) +
  scale_pattern_density_manual(values= c(0.35, 0, 0.1, 0, 0.35, 0)) +
  scale_pattern_fill_manual(values=c("#98bbbb", "#98bbbb", "darkgoldenrod1", 
                                     "bisque1", "bisque1", "bisque1")) +
  ylab("Effective population size (local group)") +
  xlab("sex") +
  theme_bw() +
  theme(text = element_text(size = 12))

mean = c()
sd=c()
sex = c()
kinship = c()
for (s in var2) {
  data_s = alltab[alltab$sex == s,]
  for (k in kin) {
    data_k = data_s[data_s$kinship == k,]
    mean = c(mean, mean(data_k$Pi))
    sd = c(sd, sd(data_k$Pi))
    sex = c(sex, s)
    kinship = c(kinship, k)
  }
}
mean_sd <- data.frame(sex=sex , kinship = kinship, mean = mean, sd=format(sd, scientific=F))
matrix(mean_sd$mean, ncol=2)
matrix(mean_sd$sd, ncol=2)

p_temp <- plot_grid(p1 + theme(legend.position = "none"),
                    p3  + theme(legend.position = "none"),
                    p4  + theme(legend.position = "none"),
                    p5  + theme(legend.position = "none"),
                    ncol = 2, labels = "auto")
p <- plot_grid(p_temp, legend = get_legend(p1), ncol = 1, rel_heights = c(0.85, 0.15))

png(file = "~/Documents/SLiM_model/Figures/Figures_villages/relatedness/figure_1.png", width = 2900, height = 3000, res = 300)
p
dev.off()

p_temp <- plot_grid(h1 + theme(legend.position = "none"),
                    h2  + theme(legend.position = "none"),
                    h3 + theme(legend.position = "none"),
                    h4  + theme(legend.position = "none"),
                    ncol = 2, labels = "auto", label_size = 24)
p_legend <- plot_grid(NULL, get_legend(h3), NULL, rel_widths = c(0.27, 0.33, 0.4))
p <- plot_grid(p_temp, legend = p_legend, ncol = 1, rel_heights = c(0.9, 0.1))

png(file = "~/Documents/SLiM_model/Figures/Figures_villages/relatedness/pval_heatmap.png", width = 6000, height = 6000, res = 300)
p
dev.off()

################## Comparison with B2P ###################
dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY1 = read.table(path, header=TRUE)
tabY1$kinship = rep("SP", length(tabY1$Gen))
tabY1 = tabY1[tabY1$Gen == 100,]
tabY1$sex = rep("male", length(tabY1$Gen))
tabY1$Pi = tabY1$Pi/(2*2.5e-8)

dir = "Tables/Pi/unilineal/extended/r=0/k=0.1/FT=150/patrilineal_2_patrilocal/village/"
path = paste(dir, "Pi_Y_mean_by_rep.txt", sep = "")
tabY2 = read.table(path, header=TRUE)
tabY2$kinship = rep("SP2B", length(tabY2$Gen))
tabY2 = tabY2[tabY2$Gen == 200,]
tabY2$sex = rep("male", length(tabY2$Gen))
tabY2$Pi = tabY2$Pi/(2*2.5e-8)

dir = "Tables/Pi/unilineal/regular/r=0/k=0.1/FT=150/strict_patrilineal_villages/village/"
path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM1 = read.table(path, header=TRUE)
tabM1$kinship = rep("SP", length(tabM1$Gen))
tabM1 = tabM1[tabM1$Gen == 100,]
tabM1$sex = rep("female", length(tabY2$Gen))
tabM1$Pi = tabM1$Pi/(2*5.5e-7)

dir = "Tables/Pi/unilineal/extended/r=0/k=0.1/FT=150/patrilineal_2_patrilocal/village/"
path = paste(dir, "Pi_Mito_mean_by_rep.txt", sep = "")
tabM2 = read.table(path, header=TRUE)
tabM2$kinship = rep("SP2B", length(tabM2$Gen))
tabM2 = tabM2[tabM2$Gen == 200,]
tabM2$sex = rep("female", length(tabM2$Gen))
tabM2$Pi = tabM2$Pi/(2*5.5e-7)

tabY = merge(tabY1, tabY2, all = T)
tabM = merge(tabM1, tabM2, all=T)
alltab = merge(tabY, tabM, all=T)

stat.test <- alltab %>%
  group_by(sex) %>%
  wilcox_test(Pi ~ kinship) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_xy_position(x = "sex", dodge = 0.8)
stat.test$y.position = c(1600, 300)

alltab = na.omit(alltab)

p1 <- ggplot(alltab) +
  scale_fill_manual(labels = c("Strict patrilineal descent", "Strict patrilineal descent followed by bilateral descent"),
  values = c("darkgoldenrod1", "darkcyan")) +
  scale_color_manual(labels = c("Strict patrilineal descent", "Strict patrilineal descent followed by bilateral descent"),
  values = c("darkgoldenrod1", "darkcyan")) +
  ylab("Effective population size") +
  xlab("sex") +
  geom_boxplot(aes(sex, Pi, fill = kinship), alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(sex, Pi, col = kinship), size=2, alpha=0.05, position = position_jitterdodge()) +
  stat_compare_means(aes(sex, Pi), label = "p.signif", label.x = 1.45, label.y = 2500) + 
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0) +
  theme_bw()

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
  descent : character chain indicating the descent rule
  "
  tab_Y <- read.table("Pi_Y_mean_by_rep.txt", header=TRUE)
  tab_Y = na.omit(tab_Y)
  generations = unique(tab_Y$Gen)
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab_Y, statistic = fun_mean, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  pi_tab_Y <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_Y$model <- model
  
  tab_M <- read.table("Pi_Mito_mean_by_rep.txt", header=TRUE)
  tab_M = na.omit(tab_M)
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab_M, statistic = fun_mean, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  pi_tab_Mito <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_Mito$model <- model

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

figure_2_transitions <- function(pi_tab_Y, pi_tab_Mito, legend, xmin, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace) {
  posd <- position_dodge(15)
  # Build legend
  ggp_split_1 <- ggplot(pi_tab_Y,                    # Create ggplot2 plot of first subset
                        aes(Gen,
                            mean,
                            color = model,
                            shape = model,
                            linetype = model)) +
    geom_errorbar(aes(ymin = ICinf,  ymax = ICsup), position = posd, linewidth=0.5) +
    geom_point(size = 2, stroke = 0.75) +
    scale_color_manual(labels = legend,
                       values = col) +
    scale_shape_manual(labels = legend,
                       values = sh) +
    scale_linetype_manual(labels = legend,
                          values = lt) +
    labs(color = "Scenarios", shape = "Scenarios", linetype = "Scenarios") +
    theme(legend.key.size = unit(0.5, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(face = 'bold', size = 9),
          legend.text = element_text(size=7),
          legend.spacing.x = unit(2, 'mm')) +
    guides(color = guide_legend(ncol=nCol, bycol=TRUE))
  ggp_legend_split_1 <- get_legend(ggp_split_1)
  
  pi_tab_Y$Gen = pi_tab_Y$Gen + xmin
  pi_tab_Mito$Gen = pi_tab_Mito$Gen + xmin
  
  p1 <- ggplot(pi_tab_Y, aes(Gen, mean/(2*2.5e-8), shape = model, col= model)) +
    theme_cowplot() +
    theme(text=element_text(size=fontsize1),
          axis.text = element_text(size=fontsize2),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(color="bisque", linewidth=0.3)) +
    annotate("rect", xmin = -200, xmax = -100, ymin = yminY, ymax = ymaxY,
             alpha = .1,fill = "orange") +
    annotate("rect", xmin = -100, xmax = 0, ymin = yminY, ymax = ymaxY,
             alpha = .05,fill = "orange") +
    labs(x= 'Time (generations)', y = 'Male effective population size') +
    scale_x_continuous(breaks = seq(xmin, 0, 20)) +
    scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(yminY, ymaxY)) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    geom_line(aes(linetype = model), position = posd, alpha = 0.5) +
    geom_pointrange(aes(ymin=ICinf/(2*2.5e-8), ymax=ICsup/(2*2.5e-8)), position = posd, size = 0.5, stroke = 0.5, linewidth = 0.3) +
    annotate("segment", x = xmin, xend = xmin, y = yminY, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth = 0.5) +
    annotate("text", x = xmin, y = textY, label = bquote(''*t[0]*''), col = 'salmon2', size = 5) +
    annotate("segment", x = -100, xend = -100, y = yminY, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth=0.5) +
    annotate("text", x = -100, y = textY, label = bquote(''*t[1]*''), col = 'salmon2', size = 5)
  
  p2 <- ggplot(pi_tab_Mito, aes(Gen, mean/(2*5.5e-7), shape = model, col= model)) +
    theme_cowplot() +
    theme(text=element_text(size=fontsize1),
          axis.text = element_text(size=fontsize2),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(color="bisque", linewidth=0.3)) +
    annotate("rect", xmin = -200, xmax = -100, ymin = yminY, ymax = ymaxY,
             alpha = .1,fill = "orange") +
    annotate("rect", xmin = -100, xmax = 0, ymin = yminY, ymax = ymaxY,
             alpha = .05,fill = "orange") +
    labs(x= 'Time (generations)', y = 'Female effective population size') +
    scale_x_continuous(breaks = seq(xmin, 0, 20)) +
    scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(yminY, ymaxY)) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    geom_line(aes(linetype = model), position = posd, alpha = 0.5) +
    geom_pointrange(aes(ymin=ICinf/(2*5.5e-7), ymax=ICsup/(2*5.5e-7)), position = posd, size = 0.5, stroke = 0.5, linewidth = 0.3) +
    annotate("segment", x = xmin, xend = xmin, y = yminY, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth = 0.5) +
    annotate("text", x = xmin, y = textY, label = bquote(''*t[0]*''), col = 'salmon2', size = 5) +
    annotate("segment", x = -100, xend = -100, y = yminY, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth=0.5) +
    annotate("text", x = -100, y = textY, label = bquote(''*t[1]*''), col = 'salmon2', size = 5)

  return(list(p1, p2))
}

paths = c("/Tables/Pi/unilineal/extended/r=0/k=0.1/FT=150/Patrilineal_2_patrilocal_relatedness/village")
name = c("4b")
legend = c("Patrilineal descent, patrilocal residence \n then bilateral descent, patrilocal residence")
xmin=-200
yminY = 0
ymaxY = 990
yminM = 640
ymaxM = 910
yminMY = 0
ymaxMY = 30
segmentY = 880
textY = 940
segmentM = 880
textM = 895
segmentMY = 26
textMY = 28
col = c("#742994")
sh = c(16)
lt = c("solid")
fontsize1 = 11
fontsize2 = 9
nCol=1
legendSpace=0.2

pi_tabs = pi_tab(paths, name, 'local')
pi_tab_Y = pi_tabs[[1]][[1]]
pi_tab_M = pi_tabs[[2]][[1]]
plots = figure_2_transitions(pi_tab_Y, pi_tab_M, legend, xmin, yminY, ymaxY,
                              yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, 
                              segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol,
                              legendSpace)
p2 <- plots[[1]]
p3 <- plots[[2]]

p_temp <- plot_grid(p2 + theme(legend.position = "none"), 
              p3 + theme(legend.position = "none"), ncol = 2,
              labels = c("b", "c"))

p_temp2 <- plot_grid(NULL, p1, NULL, ncol = 3,
              rel_widths = c(0.1, 0.8, 0.2))

p <- plot_grid(p_temp2, p_temp, ncol = 1, 
              labels = c("a", ""))

png(file = "Figures/YM_SP_SP2B.png", width = 3000, height = 2500, res = 300)
p
dev.off()

