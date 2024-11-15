function_compute_alpha_div = function(df_mp_species) {
  tmatm = df_mp_species %>%
    select(seq_id, species_mp, pctseqs) %>%
    pivot_wider(
      names_from = species_mp,
      values_from = pctseqs,
      values_fn = sum,
      values_fill = 0
    ) %>%
    column_to_rownames(var = "seq_id")
  
  tmatm.tran <- round(tmatm / rowSums(tmatm) * 1000000)        
  invm <- vegan::diversity(tmatm, index = "invsimpson")   
  
  invdfm <- tibble(seq_id = names(invm), invSimp =  invm)          
  
  if (length(invm) < 2) {
      invdfm <- invdfm %>%         
        mutate(seq_id = rownames(tmatm))
                        }        
  
  alpha.m <- estimateR(tmatm.tran) %>%     
    t() %>%     
    as.data.frame() %>%     
    rownames_to_column(var = "seq_id") %>%     
    left_join(invdfm)   
  
  return(alpha.m)
}

function_plot_km = function(df, 
                            selected_event = "event_30day",
                            selected_time = "time2event_30day",
                            censor_days = 30,
                            selected_metadata = "invSimp",
                            selected_cutpoint = 3.94,
                            xmax = 90,
                            break.x.by = 10) {
  require(survival)
  require(ggsurvfit)
  require(survminer)
  
  df_surv = df %>% 
    select(metabolomicsID, 
           as.name(selected_event), 
           as.name(selected_time), 
           as.name(selected_metadata))
  
  colnames(df_surv) = c("sample_id", "event", "time", "metadata")
  
  df_surv = df_surv %>% 
    mutate(time = ifelse(time > censor_days, censor_days, time),
           time = ifelse(is.na(time), censor_days, time)) %>% 
    mutate(metadata = ifelse(is.na(metadata), mean(df_surv$metadata, na.rm = T), metadata))
  
  temp = df_surv %>% 
    mutate(metadata_binary = ifelse(metadata > selected_cutpoint, "High", "Low"))
  
  name = paste0(str_to_title(gsub("_", " ", selected_event)), ": ", selected_metadata)
  
  s.fit <- survfit(Surv(time, event) ~ metadata_binary, data = temp) 
  
  ggsurvplot(s.fit, data = temp,
             pval = T,
             palette = c( '#0072b5', '#bc3c29'),
             xlab = "Time (days)",
             conf.int = T,
             risk.table = "nrisk_cumevents",
             risk.table.y.text= F,
             risk.table.pos = "in",
             risk.table.fontsize = 5.25,
             xlim = c(0-0.5, xmax+0.5),
             break.x.by = break.x.by) +
    theme_ggsurvfit_default() +
    labs(x = "Time (days)", 
         title = name,
         caption = paste0(selected_metadata, " = ", selected_cutpoint)) %++%
    theme(legend.position = "top") 
  
}

run_coxPH = function(df, 
                     event = "event_90day", 
                     time = "time2event_90day",
                     cut_var = "invSimp",
                     vars = c("age","sex", "meld_na"),
                     cutpoint = 3.94) {
  
  require("survminer")
  require("survival")
  
  dy <- unique(df$day)
  
  df_surv = df %>% 
    select(metabolomicsID, 
           as.name(event), 
           as.name(time), 
           as.name(cut_var),
           all_of(vars))
  
  colnames(df_surv)[1:4] = c("sample_id", "event", "time", "metadata")
  
  temp = df_surv %>% 
    mutate(metadata_binary = ifelse(metadata > cutpoint, "High", "Low"),
           metadata_binary = factor(metadata_binary, 
                                    levels = c("High", "Low")))
  
  my.for = as.formula(paste0("Surv(time, event) ~ ",
                             paste0(vars, collapse = " + "),
                             " +  metadata_binary"))
  
  s.fit <- coxph(my.for, data = temp)
  
  fit.sum <- summary(s.fit)
  
  ggadjustedcurves(s.fit,
                   data = temp %>%
                     as.data.frame(),
                   variable = "metadata_binary",
                   method = "average") +
    annotate(geom = "text", 
             x = 5, y = 0.25, 
             label = paste0("p-value = ", 
                            signif(fit.sum$coefficients["metadata_binaryLow","Pr(>|z|)"], 3)), 
             color = "black")
  
  ggsave(paste0("adjustedCurve.",dy,".",event,".",cut_var, ".pdf"),
         height = 6.5, width = 7)
  
  return(summary(s.fit)$coefficients)
  
}
