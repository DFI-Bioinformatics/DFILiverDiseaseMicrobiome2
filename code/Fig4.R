library(tidyverse)

# Fig 4: Microbiome characteristics at discharge predict rehospitalization 


# Histograms --------------------------------------------------------------

## A: Histogram of last sample collection  ---------------------------------

last.sub <- read_csv("data/readm_Fig4A.csv")

last.sub %>% 
  mutate(daysSinceDischarge = as.numeric(date_collected - dis_date)) %>% 
  ggplot(aes(daysSinceDischarge)) +
  geom_histogram(stat = "count", bins = 30) +
  theme_bw() +
  labs( x = "Days Since Discharge",
        y = "Frequency",
        title = paste0("1st stool sample collection date to 1st hospitalization discharge date\nfor ", nrow(last.sub), " patients")) +
  geom_rect(aes(xmin = -3.5, xmax = 0.5,
                ymin = 0,  ymax = 78),
            color = "red", alpha = 0, linetype = 2) +
  geom_text(aes(y = (..count..),
                label = (..count..),
                # label =  scales::percent((..count..)/sum(..count..))
  ), 
  stat="bin",colour="blue",vjust=-0.5, hjust = 0.5, size = 3, bins = 30) +
  coord_cartesian(xlim = c(-15,0)) 

ggsave(paste0("results/Fig4A_1st_stool_2_discharge.hist.3days.pdf"), 
       height = 4.5, width = 6.05)


## B: Histogram of 30-day re-hospitalizations or death ---------------------

last_readm <- read_csv("data/readm_Fig4B.csv")

second_adm <- last_readm %>% 
  group_by(subject_id) %>% 
  summarise(second_adm_yes = any(second_adm == "Yes")) %>% 
  ungroup()

yes_readm <- second_adm %>% 
  filter(second_adm_yes) %>% 
  left_join(last_readm) %>% 
  filter(adm_date < second_adm_date) %>% 
  group_by(subject_id) %>% 
  slice_min(second_adm_date, n = 1) 

comp_readm <- second_adm %>% 
  filter(!second_adm_yes) %>% 
  left_join(last_readm) %>% 
  bind_rows(yes_readm) %>% 
  add_count(subject_id) %>% 
  filter(n == 1)  # remove patients with prior admissions

yes_readm %>% 
  mutate(days2readm = as.numeric(second_adm_date - dis_date)) %>%
  # filter(days2readm <0) %>% 
  # view
  ggplot(aes(days2readm)) + 
  geom_histogram(bins = 100) +
  theme_bw() +
  labs( x = "Days to re-admission",
        y = "Frequency",
        title = paste0("1st hospitalization discharge date to the next admission date\nfor ", nrow(yes_readm), " patients with re-admissions"),
        caption = "only patients who had a stool sample within 3 days before discharge are included") +
  geom_rect(aes(xmin = -0.5, xmax = 30.5,
                ymin = 0,  ymax = 23),
            color = "red", alpha = 0, linetype = 2) +
  coord_cartesian(xlim = c(0, 365))

ggsave(paste0("results/Fig4B_discharge_2_readm.hist.3day.pdf"), 
       height = 4.65, width = 6.05)

# perform survival analysis -----------------------------------------------

## C: survival curve for inverse simpson -----------------------------------

dat.app <- read_csv("data/readm_Fig4C-D.csv")

source("code/utility.R")

function_plot_km(dat.app,
                 selected_event = "event_30day",
                 selected_time = "time2event_30day",
                 censor_days = 30,
                 selected_metadata = "invSimp",
                 selected_cutpoint = 4,
                 xmax = 30)

ggsave("results/Fig4C_readm_invSimp.survival.pdf", height= 7.5, width = 6.25)

### coxph -------------------------------------------------------------------

inv_coxph <- run_coxPH(dat.app,
          event = "event_30day", 
          time = "time2event_30day",
          cut_var = "invSimp",
          vars = c("age","sex", "meld_na"),
          cutpoint = 4)

inv_cox.res <- inv_coxph |> 
  as.data.frame() |> 
  mutate(var = "invSimp") |> 
  rownames_to_column()

write_csv(inv_cox.res, "results/Fig4C_readm_invSimp.coxPH.csv")

## survival curve for butyrate and DCA -------------------------------------

newdat <- dat.app %>% 
  mutate(bothlow = if_else(butyrate < 2500 & `deoxycholic acid`< 40,
                           1, 0),
         bothhigh = if_else(butyrate >= 2500 & `deoxycholic acid`>= 40,
                            1, 0)) |> 
  filter(bothlow == 1 | bothhigh == 1)

b.fit <- survfit(Surv(time2event_30day, event_30day) ~ bothlow, data = newdat) 

ggsurvplot(b.fit, data = newdat,
           pval = T,
           palette = c( '#0072b5', '#bc3c29'),
           xlab = "Time (days)",
           conf.int = T,
           risk.table = "nrisk_cumevents",
           risk.table.y.text= F,
           risk.table.pos = "in",
           risk.table.fontsize = 3.5,
           xlim = c(0, 30),
           break.x.by = 10) +
  theme_ggsurvfit_default() +
  labs(x = "Time (days)", 
       title = "Event 30 days: Both metabolites low vs. both high") %++%
  theme(legend.position = "top") 

ggsave("results/Fig4D_readm_bothlow.vs.bothhigh.survival.pdf", height= 7.5, width = 6.25)

### coxph -------------------------------------------------------------------

s.fit <- coxph(Surv(time2event_30day , event_30day) ~ age + sex + meld_na + bothlow, 
               data = newdat)

fit.sum <- summary(s.fit)

cox.res <- summary(s.fit)$coefficients |> 
  as.data.frame() |> 
  rownames_to_column()

write_csv(inv_cox.res, "results/Fig4D_readm_bothhigh.vs.bothlow.coxPH.csv")
