library(tidyverse)
library(aorsf)
library(recipes)
library(conflicted)
library(tidymodels)
library(rstatix)
library(EnhancedVolcano)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::select)

# Fig 5: Rapid quantification of fecal butyrate and DCA correlates with 30-day survival better than metagenomic data

# load data ---------------------------------------------------------------

df_meta <- readRDS("data/first_meta.rds")
df_quant <- readRDS("data/first_metabolomics_quant.rds")
df_mp <- readRDS("data/first_mp.rds")

# B: The top 30 negation important features for RSF model -------------------

## organize data before modeling --------------------------------------------

dat <- df_meta |> 
  filter(disease_stage != "Healthy") |> 
  mutate(event_30day = as.numeric(death_transplant_30d),
         time2event_30day = if_else(days_to_death_transplant < 30, 
                                    days_to_death_transplant,
                                    30),
         time2event_30day = if_else(is.na(time2event_30day), 30, time2event_30day),
         time2event_30day = if_else(time2event_30day ==0, 0.5, time2event_30day))

## metagenomics ------------------------------------------------------------

mpa_plt <- df_mp %>%
  filter(!is.na(taxid),
         seq_id %in% dat$seq_id,
         Species != "") |> 
  add_count(seq_id, wt = pctseqs, name = "total") |> 
  mutate(pctseqs = pctseqs/total)

mpa_mat <- mpa_plt |> 
  select(seq_id, Species, pctseqs) |> 
  pivot_wider( names_from = Species, 
               values_from = pctseqs,
               values_fill = 0,
               values_fn = sum)

### add more metabolites ----------------------------------------------------

metab.more <- df_quant |> 
  filter(metabolomicsID_new %in% !!dat$metabolomicsID_new) |> 
  select(metabolomics_id = metabolomicsID_new, hmmf_panel, compound, value_um) |> 
  distinct()

sel.metab <- metab.more |> 
  count(hmmf_panel, compound) |> 
  filter(n == 307) 

metab.more.sel <- metab.more |> 
  filter(compound %in% sel.metab$compound) |> 
  select(-hmmf_panel) |> 
  mutate(value_um = if_else(value_um < 0, 0, value_um)) |> 
  spread(compound, value_um)

metab.more.sel <- metab.more.sel |> 
  select(-c(butyrate,`deoxycholic acid`)) 

### filter low variance + normalize --------------

mpa_event_id <- mpa_mat |> 
  left_join(dat |> 
              select(seq_id, event_30day)) 

rec <- recipe(event_30day ~ ., data = mpa_event_id) |> 
  update_role(seq_id, new_role = "ID") 

nzv_filter_norm <- rec %>%
  step_zv(all_predictors()) |> 
  step_nzv(all_predictors(), unique_cut = 15) |> 
  step_normalize(all_predictors())

filtered_obj <- prep(nzv_filter_norm, training = mpa_event_id)

map_filtered_id <- bake(filtered_obj, mpa_event_id)

### normalize meta ----------------------------------------------------------

dat.sel <- dat |> 
  select(seq_id, time2event_30day, event_30day, butyrate = butyrate_rapid, 
         dca = butyrate_rapid, 
         race, age, sex, meld_na,
         metabolomics_id = metabolomicsID_new) |> 
  left_join(metab.more.sel) |> 
  select(-metabolomics_id)

meta_rec <- recipe(event_30day ~ .,
                   data = dat.sel) |> 
  update_role(seq_id, new_role = "ID") 

meta_norm_trans <- meta_rec |> 
  step_normalize()

norm_obj <- prep(meta_norm_trans, dat.sel)

dat.trans <- bake(norm_obj, dat.sel)

## assemble final data -----------------------------------------------------

metab.mat <- map_filtered_id |> 
  left_join(dat.trans)  

ori_names <- colnames(metab.mat)

colnames(metab.mat) <- str_replace_all(colnames(metab.mat), regex("\\W+"), "_")
colnames(metab.mat) <- gsub("^_","",colnames(metab.mat))

colnames(metab.mat) <- if_else(grepl("^[0-9]", colnames(metab.mat)),
                               paste0("x_", colnames(metab.mat)),
                               colnames(metab.mat))

names(ori_names) <- colnames(metab.mat)


## modeling ----------------------------------------------------------------
set.seed(897253)

mpa_sv_fit_full <- orsf(data = metab.mat,
                        n_tree = 500, 
                        formula = Surv(time2event_30day, event_30day) ~ . -seq_id)

### negation ----------------------------------------------------------------

neg_imp <- orsf_vi_negate(mpa_sv_fit_full) 

neg_30 <- tibble(neg_imp = neg_imp,
                 variable = names(neg_imp)) |> 
  slice_max(neg_imp, n = 30) 

neg_30 |> 
  ggplot(aes(y = reorder(gsub("_", " ", variable), neg_imp), x = neg_imp)) +
  geom_col(color = "black", fill = "gray", alpha = 0.65, width = 0.65) +
  theme_bw() +
  labs( x = "negation importance",
        y = "feature")

ggsave(paste0("results/Fig5B.pdf"), width = 6.5, height = 4.6)

# C: volcano plot of top 30 features ----------------------------------------

neg_feat_df <- mpa_event_id |> 
  left_join(dat.sel) |> 
  select(seq_id, event_30day,
         all_of(ori_names[neg_30$variable])) |> 
  select(-race) |> 
  gather("variable", "value", -c(seq_id, event_30day)) |> 
  mutate(event_30day = factor(event_30day, levels = c(0, 1),
                              labels = c("survival", "death/\ntransplant")),
         variable = factor(variable,
                           levels = neg_30$variable)) 

neg_padj <- neg_feat_df |> 
  group_by(variable) |> 
  wilcox_test(value ~ event_30day)|>
  mutate(p.adj = p.adjust(p, method = "BH")) |> 
  add_significance(p.col = "p.adj")

neg_fc <- neg_feat_df |> 
  group_by(variable, event_30day) |> 
  summarise(med = mean(value)) |> 
  mutate(med = if_else(med == 0, 0.0001, med)) |> 
  spread(event_30day, med) |> 
  mutate(log2fc = log2(`death/\ntransplant`/`survival`))

neg_tot <- left_join(neg_fc, neg_padj) %>%
  column_to_rownames(var = "variable")

set.seed(123456)

xylims = ceiling(max(abs(neg_tot$log2fc)))
plims <- ceiling(max(-log10(neg_tot$p.adj)))

neg_volcano <-
  EnhancedVolcano(neg_tot,
                  lab = rownames(neg_tot),
                  title = 'FDR corrected: Death/Transplant vs No event',
                  y = "p.adj",
                  x = "log2fc",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 2.5,
                  xlim = c(-xylims,xylims),
                  ylim = c(0, plims + 1),
                  col=paletteer::paletteer_d("ggthemes::fivethirtyeight")[c(3,2,5,4)],
                  colAlpha = 0.65,
                  legendPosition = "right",
                  legendLabels = c(expression(p > 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p > 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1)),
                  legendLabSize = 14,
                  drawConnectors = T,
                  widthConnectors = 0.25,
                  maxoverlapsConnectors = Inf,
                  arrowheads = F,
                  gridlines.minor = F,
                  gridlines.major = F) +
  labs(subtitle = "",
       y = expression( -Log[10] ~ P.adj)) +
  annotate("text", x = 0.65*xylims, y = plims + 0.5, label = "death / transplant",
           size = 6, color = ggsci::pal_igv("default")(2)[2]) +
  annotate("rect", xmin = 1.05, xmax = Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = ggsci::pal_igv("default")(2)[2]) +
  annotate("text", x = -0.65*xylims, y = plims + 0.5, label = "survival",
           size = 6, color = ggsci::pal_igv("default")(2)[1]) +
  annotate("rect", xmin = -1.05, xmax = -Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = ggsci::pal_igv("default")(2)[1]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

neg_volcano

ggsave(plot = neg_volcano,
       paste0("results/Fig5C.pdf"),
       width = 12, height = 6.5, units = "in")

# A: ROC curve of RSF ----------------------------------------

mpa_metab_split <- initial_split(metab.mat, strata = event_30day, prop = 4/5)

train_data <- training(mpa_metab_split)
test_data <- testing(mpa_metab_split)

folds <- vfold_cv(train_data, v = 5, repeats = 2)

library(pROC)

auc_roc_train <- map(
  folds$splits,
  
  function(sp){
    df <- as.data.frame(sp)
    
    mpa_metab_split <- initial_split(df, strata = event_30day)
    train_int_data <- training(mpa_metab_split)
    test_int_data <- testing(mpa_metab_split)
    
    sp_sv_fit <- orsf(data = train_int_data,
                      n_tree = 500, 
                      formula = Surv(time2event_30day, event_30day) ~ . -seq_id)
    
    test_full_nolab <- test_int_data |> 
      select(-c(time2event_30day, event_30day, seq_id))
    
    test_full_lab <- test_int_data |> 
      select(c(time2event_30day, event_30day, seq_id))
    
    full_preds <- predict(sp_sv_fit, test_full_nolab, pred_horizon = 30, pred_type="risk") |> 
      bind_cols(test_full_lab) 
    
    train_test_roc <- roc(full_preds$event_30day, full_preds$...1)
    
    test_data_lab <- test_data |> 
      select(c(time2event_30day, event_30day, seq_id))
    
    test_data_nolab <- test_data |> 
      select(-c(time2event_30day, event_30day, seq_id))
    
    test_preds <- predict(sp_sv_fit, test_data_nolab, pred_horizon = 30, pred_type="risk") |> 
      bind_cols(test_data_lab) 
    
    test_ext_roc <- roc(test_preds$event_30day, test_preds$...1)
    
    list(train = train_test_roc, test = test_ext_roc)
    
  }
  
)


parseROC <- function(auc_roc) {
  
  tibble(sensitivity = auc_roc$train$sensitivities,
         specificity = auc_roc$train$specificities,
         auc = as.numeric(auc_roc$train$auc)) |> 
    mutate(source = "training") |> 
    bind_rows(tibble(sensitivity = auc_roc$test$sensitivities,
                     specificity = auc_roc$test$specificities,
                     auc = as.numeric(auc_roc$test$auc),
                     source = "test"))
  
}

dd.roc <- list()

for (i in seq_along(auc_roc_train)){
  temp_roc <- auc_roc_train[[i]]
  
  dd.roc[[i]] <- parseROC(temp_roc) |> 
    mutate(iter = i)
}

dd.roc <- dd.roc |>  
  bind_rows() 

train_auc_roc <- dd.roc |> 
  filter(source == "training") 

train_auc_ave <- train_auc_roc |> 
  distinct(source, iter, auc) |> 
  summarise(mean_auc = mean(auc)) |> 
  pull(mean_auc)

test_auc_roc <- dd.roc |> 
  filter(source == "test") 

test_auc_ave <- test_auc_roc |> 
  distinct(source, iter, auc) |> 
  summarise(mean_auc = mean(auc)) |> 
  pull(mean_auc)

train_auc_roc |> 
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = iter), color = "grey", alpha = 0.45) +
  geom_smooth(data = train_auc_roc, 
              colour = "darkblue", size = 1.25) +
  geom_smooth(data = test_auc_roc,
              colour = "darkgreen", size = 1.25) +
  ylim(0,1) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = 0.99, y = c(0.07, 0.03),              
           label = c(paste("Mean AUC =", signif(train_auc_ave, 2)),                        
                     paste("Test AUC =", signif(test_auc_ave, 2))),              
           color = c("darkblue", "darkgreen"),              
           hjust = 1, # Right justify              
           vjust = 0, # Bottom align              
           size = 3.5) +          # Customize theme     
  theme_classic() +
  labs(x = "False Positive Rate",          
       y = "True Positive Rate",
       title = "Random Survival Forest")

ggsave(paste0("results/Fig5A.pdf"),
       height = 5.6, width = 6.5)
