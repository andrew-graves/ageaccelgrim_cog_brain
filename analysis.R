# Author: Andrew J. Graves
# Date: 3/24/24

# Load packages
library(tidyverse)
library(lubridate)
library(EGAnet)
library(brms)
library(tidybayes)
library(parallel)
library(wesanderson)
library(patchwork)
library(kableExtra)

# Set plot themes
theme_set(theme_classic())
theme_update(plot.title = element_text(size = 32, 
                                       hjust = .5),
             plot.subtitle = element_text(size = 18, 
                                          hjust = .5),
             plot.caption = element_text(size = 16,
                                         hjust = 0),
             text = element_text(family = "serif", 
                                 size = 20),
             legend.position = "bottom")

# Get the number of available cores
n_cores = detectCores() - 1

# Identify the four cognitive domains
proc_speed <- c("DigSym", "PatCom", "LetCom")
memory <- c("Recall", "LogMem", "PAssoc")
reasoning <- c("MatRea", "Ship", "LetSet")
spatial <- c("SpaRel", "PapFld", "FrmBrd")
vocab <- c("Vocab", "PVoc", "SynVoc", "AntVoc")

# Define helper functions for plotting correlation matrices

get_lower_tri <- function(x){
  x[lower.tri(x, diag = FALSE)] <- NA
  return(x)
}

plot_lower_tri <- function(mat, func){
  cor_matrix <- mat %>%
    get_lower_tri() %>%
    reshape2::melt(na.rm = FALSE)
  cor_matrix %>%
    ggplot(aes(fct_inorder(Var1), 
               fct_rev(fct_inorder(Var2)), 
               fill = value)) +
    geom_tile() +
    geom_text(aes(Var1, Var2,
                  label = broman::myround(value, digits = 2) %>%
                    str_replace("NA", "") %>%
                    str_replace(fixed("0."), ".") %>%
                    str_replace(fixed("1.00"), "1")
    ), color = "black", size = 3) +
    theme(axis.text.x = element_text(size = 14, angle = 45,
                                     hjust = 1),
          axis.text.y = element_text(size = 14),
          legend.justification = c(1, 0),
          legend.position = c(1, .5)
    ) +
    labs(x = "", y = "")
}

# Define helper functions for table formatting

# Convenient standard digit formatting
rnd <- function(x, digits = 2){
  return(format(
    round(x, digits = digits), 
    nsmall = digits))
}

# Set table format
format_table <- function(x, caption = "", 
                         provide_cols = FALSE, col_names = "", 
                         opts = "HOLD_position") {
  
  if(!provide_cols){
    col_names <- colnames(x)
  }
  
  x %>%
    mutate_if(is.numeric, rnd) %>%
    kbl(booktabs = TRUE, caption = caption, 
        col.names = col_names, escape = FALSE) %>% 
    kable_styling(latex_options = opts) %>%
    return()
}

## Generate latent cognitive variables with EGA

# Select data for EGA
ega_dat <- nearest_time_dat %>%
  select(DigSym:AntVoc, -PVoc)

# Run EGA
set.seed(42)
orig_ega <- EGA(data = ega_dat, 
                model = "glasso",
                algorithm = "walktrap"
)

boot_ega <- bootEGA(data = ega_dat, iter = 1000,
                    typicalStructure = TRUE,
                    model = "glasso",
                    algorithm = "louvain",
                    type = "parametric", 
                    ncores = n_cores,
                    seed = 42)

stab <- itemStability(boot_ega)

# Compute network scores for each cognitive dimension
net_scores <- net.scores(nearest_time_dat,
                         A = boot_ega$typicalGraph$graph,
                         wc = boot_ega$typicalGraph$wc,
                         impute = "none")$std.scores

# Add latent cognitive variables to full dataset
final_df <- nearest_time_dat %>%
  bind_cols(net_scores) %>%
  rename(comm_1 = `1`,
         comm_2 = `2`,
         comm_3 = `3`,
         comm_4 = `4`
  ) %>%
  select(-c(Imaging.ID, scan_date, cog_date)) %>%
  mutate_at(vars(time_diff, time_since), as.numeric) %>%
  scale() %>%
  data.frame() %>%
  mutate(Sex = factor(Sex, labels = c("Male", "Female")),
         SubjID = as.character(SubjID))

# Store community names
comm_names <- paste0("comm_", 1:4)
latent_names <- c("Proc. Speed", "Memory",
                  "Spatial/ Reasoning", "Vocabulary")

## Fit Bayesian hierarchichal multivariate model

```{r fit, cache = TRUE}
# Clean up memory by removing temporary datasets
rm(list = ls()[grepl("dat", ls())])

# Specify model formula
brm_formula <- bf(mvbind(comm_1, comm_2, comm_3, comm_4) ~ 
                    AgeAccelGrim + Age + time_since + Session + Sex +
                    NK + Mono + Gran + PlasmaBlast +
                    CD8pCD28nCD45RAn + CD8.naive + CD4.naive + 
                    (Session|p|SubjID))

# Specify Gaussian priors w/ finite variance on betas
# Modify nu prior for Vocab due to skew (allow estimation of 
# larger df than default). Posterior nu will still be small
priors <- c(
  set_prior("normal(0, 10)",
            class = "b",
            resp = "comm1"),
  set_prior("normal(0, 10)", 
            class = "b", 
            resp = "comm2"),
  set_prior("normal(0, 10)",
            class = "b",
            resp = "comm3"),
  set_prior("normal(0, 10)",
            class = "b",
            resp = "comm4"),
  set_prior("gamma(12, 0.1)",
            class = "nu",
            resp = "comm4")
)

# Fit the model
mod <- brm(brm_formula + set_rescor(FALSE),
           data = final_df, 
           family = "student",
           prior = priors,
           seed = 42,
           chains = n_cores,
           cores = n_cores,
           iter = 20000,
           algorithm = "sampling",
           sample_prior = "yes",
           save_all_pars = FALSE,
           control = list(adapt_delta = .99,
                          max_treedepth = 12)
)

# Set up the hypothesis testing framework
h <- c("comm1_AgeAccelGrim < 0", "comm2_AgeAccelGrim < 0",
       "comm3_AgeAccelGrim < 0", "comm4_AgeAccelGrim < 0")

# Test AgeAccelGrim parameters
hyp_test <- hypothesis(mod, h)

# Computer posterior probabilities
bayes_p_vals <- (1 - hyp_test$hypothesis$Post.Prob)*2
# Adjust using FDR
p_adj <- p.adjust(bayes_p_vals, method = "fdr")
```

## Extract posterior information from the models

# Extract conditional effects from the model
cond_ef <- conditional_effects(mod, effects = "AgeAccelGrim",
                               plot = FALSE) %>%
  map(data.frame)
# Name the elements in the list for identification
names(cond_ef) <- comm_names

# Flatten the list into a single data frame
resp_ef <- cond_ef %>%
  bind_rows(.id = "resp")

# Get posterior samples for GrimAge population parameters
post_samps <- mod %>% 
  posterior_samples(pars = "AgeAccelGrim") 

names(post_samps) <- comm_names

## Plot EGA results (Figure 1)

# EGA plot
ega_1 <- boot_ega$plot.typical.ega + 
  scale_color_manual(values = 
                       c(wes_palette("Darjeeling1", 
                                     n = length(comm_names))[1:3],
                         wes_palette("BottleRocket2")[3]),
                     labels = latent_names) + 
  labs(x = "", y = "", color = "Latent Factor") +
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 20))

# GLASSO adjacency matrix for cognitive data
ega_2 <- boot_ega$typicalGraph$graph %>% 
  na_if(0) %>%
  plot_lower_tri() + 
  scale_fill_gradientn(
    name = "Edge Strength",
    colors = wes_palette(
      "Royal1", n = 10000, type = "continuous"), 
    na.value = "white") + 
  labs(title = "")

# Eigenvalue-eigenvector decomposition
eig <- eigen(cor(ega_dat, 
                 use = "pairwise.complete.obs"))$values
# Scree plot of cognitive data
ega_3 <- data.frame(comp = 1:length(eig), vals = eig) %>%
  ggplot(aes(x = comp, y = vals)) +
  geom_line(color = wes_palette("Rushmore1")[4], size = 1) +
  geom_point(color = wes_palette("Rushmore1")[3], size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  labs(x = "Sorted Eigenvector Index", y = "Eigenvalue",
       title = "")

# Combine the EGA plots w/ patchwork
(ega_1 | (ega_2 / ega_3)) + 
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("fig1.png", dpi = 300, 
       width = 12, height = 8)

## Aggregate plots of AgeAccelGrim associations with cognition (Figure 2)

# Aggregate the data
agg_dat <- final_df %>%
  group_by(SubjID) %>%
  summarize_at(vars(AgeAccelGrim, Age,
                    DigSym:AntVoc, -PVoc, contains("comm")), mean, 
               na.rm = TRUE) %>%
  rename(`Proc. Speed` = comm_1, Memory = comm_2,
         `Spatial/ Reasoning` = comm_3, Vocabulary = comm_4)

# Raw correlation matrix between cognition and age metrics
corr_plot <- agg_dat %>%
  select(-SubjID) %>%
  cor(method = "pearson", use = "pairwise.complete.obs") %>%
  plot_lower_tri() + 
  scale_fill_gradientn(name = expression(paste(
    "Pearson's ", italic("r"))),
    colors = wes_palette("Zissou1", 
                         n = 10000, 
                         type = "continuous"), 
    na.value = "white")

ggsave("fig2.png", dpi = 300, 
       width = 12, height = 8)

## Plot the results for AgeAccelGrim posteriors

# Number of sessions
n_sess <- 3

# Axis label objects
grim_lab <- expression(paste(
  "AgeAccelGrim (", italic("Z"), ")"))
cog_lab <- expression(paste(
  "Cognitive Scores (", italic("Z"), ")"))

plot_dat <- final_df %>%
  select(AgeAccelGrim, Session, comm_1:comm_4) %>%
  pivot_longer(c(-AgeAccelGrim, -Session), 
               names_to = "resp") %>%
  mutate(Session = factor(Session, labels = 1:n_sess))

# Produce labels for facets
resp_labs <- as_labeller(c(
  "comm_1" = latent_names[1], 
  "comm_2" = latent_names[2],
  "comm_3" = latent_names[3], 
  "comm_4" = latent_names[4]))

# Required for adding these plots together via patchwork
blank_labs <- as_labeller(c(
  "comm_1" = "", 
  "comm_2" = "",
  "comm_3" = "", 
  "comm_4" = ""))

# Build annotation text data frame
ann_text <- data.frame(AgeAccelGrim = rep(0, 2), 
                       estimate__ = c(-4, 3.2), 
                       lab = 
                         c("Worse Cognitive Performance", 
                           "Better Cognitive Performance"),
                       resp = factor(rep("comm_1", 2), levels = comm_names))

# Scatterplots for GrimAge posteriors
p1 <- resp_ef %>%
  ggplot(aes(x = AgeAccelGrim)) + 
  geom_point(data = plot_dat, aes(y = value,
                                  color = Session
  ), alpha = 0.95) + 
  geom_line(aes(y = estimate__),
            color = "black", size = 1) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              alpha = .6) +
  facet_wrap(~resp, nrow = 1, labeller = resp_labs
  ) + 
  geom_text(data = ann_text, aes(x = AgeAccelGrim,y = estimate__, label = lab),
            family = "serif") + 
  scale_color_manual(values = wes_palette("GrandBudapest2", 
                                          n = n_sess, 
                                          type = "continuous")) +
  labs(x = grim_lab,
       y = expression(paste(
         "Latent Scores (", italic("Z"), ")")),
       color = "Session") + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme(legend.position = "right")

# Posterior density plots for GrimAge parameters
p2 <- post_samps %>%
  pivot_longer(everything(), names_to = "resp") %>%
  ggplot(aes(x = value)) +
  tidybayes::stat_halfeye(.width = c(0.66, 0.95)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~resp, nrow = 1, labeller = blank_labs) +
  labs(x = expression(paste(beta, 
                            " coefficient values w/ priors N ~ (",
                            mu, "= 0, ", sigma, "= 10)")),
       y = "Posterior Density") +
  theme(strip.background = element_blank())

# Combine the GrimAge posterior plots w/ patchwork
(p1 / p2)  +
  plot_layout(heights = c(2, .6)) +
  plot_annotation(tag_levels = "A")

ggsave("fig3.png", dpi = 300, 
       width = 12, height = 8)

## Posterior predictive check plots

pp_plots <- list()

set.seed(42)
# Iterating manually for plotting convenience
for(i in seq_along(mods)){
  # Generate samples from posterior
  pp_plots[[i]] <- mods[[i]] %>%
    pp_check(nsamples = 100) + 
    scale_color_manual(values = wes_palette("Darjeeling1", n = 2), 
                       labels = c("Observed", "Generated Samples")) +
    labs(color = "", title = latent_names[i])
  
  # Arrange the labeling of the grid
  if(i == 3){
    pp_plots[[i]] <- pp_plots[[i]] + 
      labs(x = "Response Variable Distribution",
           y = "Density")
    
  } else {
    pp_plots[[i]] <- pp_plots[[i]] + 
      guides(color = "none")
  }
}

# Combine the posterior predictive check plots w/ patchwork
(pp_plots[[1]] + pp_plots[[2]] + pp_plots[[3]] + pp_plots[[4]]) + 
  plot_annotation(
    title = "Posterior Predictive Distribution Checks",
    tag_levels = "I")

ggsave("pp_check.png", dpi = 300, 
       width = 12, height = 8)

## Table for full models

# Get table data
tab_dat <- mods %>% 
  map(summary) %>%
  # Extract the fixed effects
  map(~ .x$fixed   %>%
        data.frame() %>%
        select(-c(Est.Error, Bulk_ESS, Tail_ESS, Rhat))
  )
# Name the elements in the list for identification
names(tab_dat) <- latent_names

# For table rows
term_names <- c(
  "Intercept", "GrimAge Acceleration", 
  "Age", "Occasion", "Session", 
  "Sex", "Education"
)

# Tabulate full model results
tab_dat %>%
  bind_rows() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = rep(term_names, length(latent_names))) %>%
  format_table("Population Parameter Estimates and 
               Uncertainty for Full Models", 
               provide_cols = TRUE, 
               col_names = 
                 c("Parameter", 
                   "$\\beta$ coefficient", 
                   "Lower $95\\%$ CI", 
                   "Upper $95\\%$ CI")) %>%
  pack_rows(latent_names[1], 1, 7, latex_gap_space = "1em") %>%
  pack_rows(latent_names[2], 8, 14, latex_gap_space = "1em") %>%
  pack_rows(latent_names[3], 15, 21, latex_gap_space = "1em") %>%
  pack_rows(latent_names[4], 22, 28, latex_gap_space = "1em")