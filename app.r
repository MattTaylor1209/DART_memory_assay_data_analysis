# app.R

# ---- Packages ----
library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(janitor)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(ggtext)
library(scales)
library(ggprism)
library(shinythemes)


# ---- Theme (adapted from your script) ----
myfacettheme <- theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.title.y = element_text(colour = "black", face = "bold", size = 14),
  axis.text.y  = element_text(colour = "grey30", face = "bold", size = 14),
  axis.title.x = element_text(colour = "black", face = "bold", size = 14),
  axis.text.x  = element_text(colour = "grey30", face = "bold", size = 14),
  axis.ticks   = element_line(linewidth = 1, colour = "black"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
  panel.grid   = element_blank(),
  axis.line = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  legend.text  = element_text(colour = "black", size = 10),
  strip.background = element_rect(colour = "black", linewidth = 2, fill = "grey30"),
  strip.text.x = element_text(size = 12, face = "bold", colour = "white")
)

# ---- Helpers ----
required_cols <- c("group", "trial", "stim_number", "gof")
pi_compute_cols <- c("max_fitted_amplitude")

GGPRISM_PALS <- (function() {
  # ggprism_data is available as a dataset in ggprism; depending on versions,
  # it may be directly visible or only via namespace.
  gd <- tryCatch(ggprism::ggprism_data, error = function(e) NULL)
  if (is.null(gd)) {
    gd <- tryCatch(getFromNamespace("ggprism_data", "ggprism"), error = function(e) NULL)
  }
  
  if (!is.null(gd) && is.list(gd$fill_palettes) && length(gd$fill_palettes) > 0) {
    return(names(gd$fill_palettes))
  }
  
  # safe fallback
  return(c("colors"))
})()


make_fill_scale <- function(palette_source, ggprism_palette = "colors") {
  if (identical(palette_source, "ggprism")) {
    sc <- tryCatch(
      ggprism::scale_fill_prism(palette = ggprism_palette),
      error = function(e) NULL
    )
    if (!is.null(sc)) return(sc)
    return(ggplot2::scale_fill_hue())
  }
  ggplot2::scale_fill_hue()
}

pdf_device_fun <- function() {
  # Prefer cairo if available (better text rendering), else pdf()
  if (capabilities("cairo")) {
    return(grDevices::cairo_pdf)
  }
  grDevices::pdf
}

compute_pi_if_needed <- function(df) {
  if ("pi" %in% names(df)) return(df)
  
  if (!all(pi_compute_cols %in% names(df))) {
    stop(
      "Your file does not contain a 'pi' column, and PI cannot be computed because it's missing: ",
      paste(setdiff(pi_compute_cols, names(df)), collapse = ", ")
    )
  }
  
  baseline <- df %>%
    group_by(group, trial, stim_number) %>%
    summarise(mean_response = mean(max_fitted_amplitude, na.rm = TRUE), .groups = "drop") %>%
    # IMPORTANT: compare as characters to avoid 1 vs "1" issues
    filter(as.character(trial) == "1", as.character(stim_number) == "1") %>%
    transmute(group, scale_factor = mean_response)
  
  df %>%
    left_join(baseline, by = "group") %>%
    mutate(pi = max_fitted_amplitude / scale_factor * 100)
}

prep_data <- function(df, group_levels, zeroed) {
  df <- df %>%
    mutate(across(where(is.numeric), ~ pmax(.x, 0))) %>%
    compute_pi_if_needed()
  
  if (isTRUE(zeroed)) {
    df <- df %>% mutate(pi = if_else(gof < 0.2, 0, pi))
  }
  
  df <- df %>%
    mutate(
      group = factor(group, levels = group_levels),
      stim_number = factor(stim_number),
      trial = factor(trial)
    ) %>%
    arrange(group, trial, stim_number) %>%
    group_by(group, trial, stim_number) %>%
    mutate(fly_id = row_number()) %>%
    ungroup()
  
  df
}

# Raw pairwise tests with no p-value adjustment
run_pairwise_raw <- function(data, y, x, test_kind = c("wilcox", "t")) {
  test_kind <- match.arg(test_kind)
  f <- reformulate(x, response = y)
  
  data %>%
    group_by(group) %>%
    {
      if (test_kind == "t") {
        rstatix::pairwise_t_test(., f, paired = TRUE, id = "fly_id", p.adjust.method = "none")
      } else {
        rstatix::pairwise_wilcox_test(., f, paired = TRUE, id = "fly_id", p.adjust.method = "none")
      }
    } %>%
    ungroup()
}


compute_all_stats <- function(df,
                              test_kind = c("wilcox", "t"),
                              trial_within,
                              stim_within,
                              stim_between,
                              trial_between_a,
                              trial_between_b,
                              p_adjust_method = "holm",
                              correct_across = c("within_endpoint", "across_endpoints")) {
  
  test_kind <- match.arg(test_kind)
  correct_across <- match.arg(correct_across)
  
  trial_within  <- as.character(trial_within)
  stim_within   <- as.character(stim_within)
  stim_between  <- as.character(stim_between)
  trial_between_a <- as.character(trial_between_a)
  trial_between_b <- as.character(trial_between_b)
  
  # Planned hypothesis count for the endpoint family:
  # - within-trial: all pairs among selected stimuli
  # - between-trial: trial 1 vs 2 at selected stimulus (1 comparison)
  n_within  <- choose(length(stim_within), 2)
  n_between <- 1L
  n_planned <- n_within + n_between
  
  # --- PI: within trial (stim compare) ---
  pi_within <- df %>%
    filter(trial == trial_within, stim_number %in% stim_within) %>%
    mutate(stim_number = factor(stim_number, levels = stim_within)) %>%
    run_pairwise_raw(y = "pi", x = "stim_number", test_kind = test_kind) %>%
    mutate(
      endpoint  = "pi",
      contrast  = "within_trial",
      trial_sel = trial_within,
      stim_sel  = NA_character_,
      n_planned = n_planned
    )
  
  # --- PI: between trials (trial compare) ---
  pi_between <- df %>%
    filter(stim_number == stim_between, trial %in% c("1", "2")) %>%
    mutate(trial = factor(trial, levels = c(trial_between_a, trial_between_b))) %>%
    run_pairwise_raw(y = "pi", x = "trial", test_kind = test_kind) %>%
    mutate(
      endpoint  = "pi",
      contrast  = "between_trial",
      trial_sel = paste0(trial_between_a, "_vs_", trial_between_b),
      stim_sel  = stim_between,
      n_planned = n_planned
    )
  
  # --- Speed: within trial (stim compare) ---
  sp_within <- df %>%
    filter(trial == trial_within, stim_number %in% stim_within) %>%
    mutate(stim_number = factor(stim_number, levels = stim_within)) %>%
    run_pairwise_raw(y = "mean_pre_stim_speed", x = "stim_number", test_kind = test_kind) %>%
    mutate(
      endpoint  = "speed",
      contrast  = "within_trial",
      trial_sel = trial_within,
      stim_sel  = NA_character_,
      n_planned = n_planned
    )
  
  # --- Speed: between trials (trial compare) ---
  sp_between <- df %>%
    filter(stim_number == stim_between, trial %in% c("1", "2")) %>%
    mutate(trial = factor(trial, levels = c(trial_between_a, trial_between_b))) %>%
    run_pairwise_raw(y = "mean_pre_stim_speed", x = "trial", test_kind = test_kind) %>%
    mutate(
      endpoint  = "speed",
      contrast  = "between_trial",
      trial_sel = paste0(trial_between_a, "_vs_", trial_between_b),
      stim_sel  = stim_between,
      n_planned = n_planned
    )
  
  all <- bind_rows(pi_within, pi_between, sp_within, sp_between)
  
  # Adjust within the intended family.
  # Option A (typical): within each endpoint (PI separate from speed) per genotype
  # Option B (more conservative): across endpoints too (PI+speed together) per genotype
  if (correct_across == "within_endpoint") {
    all <- all %>%
      group_by(group, endpoint) %>%
      mutate(p.adj = p.adjust(p, method = p_adjust_method, n = first(n_planned))) %>%
      ungroup()
  } else {
    # across_endpoints: planned n doubles (PI has n_planned, speed has n_planned)
    all <- all %>%
      group_by(group) %>%
      mutate(p.adj = p.adjust(p, method = p_adjust_method, n = first(n_planned) * 2L)) %>%
      ungroup()
  }
  
  all %>%
    add_significance("p.adj")
}



make_facet_labels <- function(df, denom_len, group_levels, n_label = TRUE) {
  
  lab <- df %>%
    group_by(group) %>%
    summarise(n = n() / denom_len, .groups = "drop") %>%
    mutate(n = as.integer(round(n)))  # avoids "12.0" etc
  
  lab <- if (isTRUE(n_label)) {
    lab %>% mutate(label = paste0(as.character(group), " (n = ", n, ")"))
  } else {
    lab %>% mutate(label = as.character(group))
  }
  
  out <- setNames(lab$label, as.character(lab$group))
  out[group_levels[group_levels %in% names(out)]]
}


make_plot <- function(df, plot_kind, test_kind, zeroed, alpha,
                      stim_within, trial_within, stim_between, n_label,
                      col_by_group = TRUE, ylimit,
                      trial_between_a, trial_between_b,
                      palette_source = c("ggprism", "ggplot2 (hue)"),
                      ggprism_palette = "colors",
                      stats_all) {
  
  test_kind <- match.arg(test_kind, c("wilcox", "t"))
  
  df <- df %>% mutate(response = gof > 0.2)
  
  group_levels <- levels(df$group)
  
  # Points (response) colouring
  resp_cols <- c(`TRUE` = "black", `FALSE` = if (isTRUE(zeroed)) "red" else "black")
  
  # Fill palette scale for group boxes
  fill_scale <- make_fill_scale(palette_source, ggprism_palette)
  
  # Variable to set fill colour by
  fill_var <- if (isTRUE(col_by_group)) "group" else "stim_number"
  
  if (plot_kind == "PI: within trial (stim compare)") {
    fdat <- df %>%
      filter(trial %in% as.character(trial_within),
             stim_number %in% as.character(stim_within)) %>%
      mutate(stim_number = factor(stim_number, levels = as.character(stim_within)))
    
    stats <- stats_all %>%
      filter(endpoint == "pi",
             contrast == "within_trial",
             trial_sel == as.character(trial_within))
    
    ymax <- max(fdat$pi, na.rm = TRUE)
    #ylim_upper <- ceiling(ymax * 1.1 * 10) / 10
    stats$y.position <- seq(ylimit * 0.95, by = 0, length.out = nrow(stats))
    
    facet_labeller <- make_facet_labels(fdat, denom_len = length(stim_within), 
                                        group_levels = group_levels,
                                        n_label = n_label)
    
    p <- ggplot(fdat, aes(x = stim_number, y = pi)) +
      geom_boxplot(aes(fill = .data[[fill_var]]), 
                   colour = "black", linewidth = 1, alpha = alpha) +
      ggbeeswarm::geom_beeswarm(aes(colour = response), alpha = 0.6, size = 2, priority = "ascending") +
      fill_scale +
      scale_colour_manual(values = resp_cols) +
      facet_wrap(~ group, nrow = 1, labeller = labeller(group = facet_labeller)) +
      ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
                                 tip.length = 0, size = 8, bracket.size = 1, bracket.nudge.y = 0) +
      labs(title = paste("Within trial:", trial_within), x = "Stimulus", y = "Performance index (%)") +
      coord_cartesian(ylim = c(0, ylimit)) +
      theme_classic() +
      myfacettheme
    
    return(list(plot = p, stats = stats, data = fdat))
  }
  
  if (plot_kind == "PI: between trials (trial compare)") {
    fdat <- df %>%
      filter(stim_number %in% as.character(stim_between),
             trial %in% c(as.character(trial_between_a), as.character(trial_between_b))) %>%
      mutate(trial = factor(trial, levels = c(as.character(trial_between_a), as.character(trial_between_b))))
    
    stats <- stats_all %>%
      filter(endpoint == "pi",
             contrast == "between_trial",
             stim_sel == as.character(stim_between),
             trial_sel == paste0(as.character(trial_between_a), "_vs_", as.character(trial_between_b)))
    
    ymax <- max(fdat$pi, na.rm = TRUE)
    #ylim_upper <- ceiling(ymax * 1.1 * 10) / 10
    stats$y.position <- seq(ylimit * 0.95, by = 0, length.out = nrow(stats))
    
    facet_labeller <- make_facet_labels(fdat, denom_len = 2, 
                                        group_levels = group_levels,
                                        n_label = n_label)
    
    p <- ggplot(fdat, aes(x = trial, y = pi)) +
      geom_boxplot(aes(fill = .data[[fill_var]]), 
                   colour = "black", linewidth = 1, alpha = alpha) +
      ggbeeswarm::geom_beeswarm(aes(colour = response), alpha = 0.6, size = 2, priority = "ascending") +
      fill_scale +
      scale_colour_manual(values = resp_cols) +
      facet_wrap(~ group, nrow = 1, labeller = labeller(group = facet_labeller)) +
      ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
                                 tip.length = 0, size = 8, bracket.size = 1, bracket.nudge.y = 0) +
      labs(title = paste0("Between trials (", trial_between_a, " vs ", trial_between_b, "), stimulus ", stim_between), 
           x = "Trial", y = "Performance index (%)") +
      coord_cartesian(ylim = c(0, ylimit)) +
      theme_classic() +
      myfacettheme
    
    return(list(plot = p, stats = stats, data = fdat))
  }
  
  if (plot_kind == "Speed: within trial (stim compare)") {
    if (!("mean_pre_stim_speed" %in% names(df))) {
      stop("Missing column 'mean_pre_stim_speed' for speed plots.")
    }
    
    fdat <- df %>%
      filter(trial %in% as.character(trial_within),
             stim_number %in% as.character(stim_within)) %>%
      mutate(stim_number = factor(stim_number, levels = as.character(stim_within)))
    
    stats <- stats_all %>%
      filter(endpoint == "speed",
             contrast == "within_trial",
             trial_sel == as.character(trial_within))
    
    ymax <- max(fdat$mean_pre_stim_speed, na.rm = TRUE)
    #ylim_upper <- ceiling(ymax * 1.1 * 10) / 10
    stats$y.position <- seq(ylimit * 0.95, by = 0, length.out = nrow(stats))
    
    facet_labeller <- make_facet_labels(fdat, denom_len = length(stim_within), 
                                        group_levels = group_levels,
                                        n_label = n_label)
    
    p <- ggplot(fdat, aes(x = stim_number, y = mean_pre_stim_speed)) +
      geom_boxplot(aes(fill = .data[[fill_var]]), 
                   colour = "black", linewidth = 1, alpha = alpha) +
      ggbeeswarm::geom_beeswarm(aes(colour = response), alpha = 0.6, size = 2, priority = "ascending") +
      fill_scale +
      scale_colour_manual(values = resp_cols) +
      facet_wrap(~ group, nrow = 1, labeller = labeller(group = facet_labeller)) +
      ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
                                 tip.length = 0, size = 8, bracket.size = 1, bracket.nudge.y = 0) +
      labs(title = paste("Within trial:", trial_within), x = "Stimulus", y = "Pre stimulus speed (mm/s)") +
      coord_cartesian(ylim = c(0, ylimit)) +
      theme_classic() +
      myfacettheme
    
    return(list(plot = p, stats = stats, data = fdat))
  }
  
  if (plot_kind == "Speed: between trials (trial compare)") {
    if (!("mean_pre_stim_speed" %in% names(df))) {
      stop("Missing column 'mean_pre_stim_speed' for speed plots.")
    }
    
    fdat <- df %>%
      filter(stim_number %in% as.character(stim_between),
             trial %in% c(as.character(trial_between_a), as.character(trial_between_b))) %>%
      mutate(trial = factor(trial, levels = c(as.character(trial_between_a), as.character(trial_between_b))))
    
    stats <- stats_all %>%
      filter(endpoint == "speed",
             contrast == "between_trial",
             stim_sel == as.character(stim_between),
             trial_sel == paste0(as.character(trial_between_a), "_vs_", as.character(trial_between_b)))
    
    ymax <- max(fdat$mean_pre_stim_speed, na.rm = TRUE)
    #ylim_upper <- ceiling(ymax * 1.1 * 10) / 10
    stats$y.position <- seq(ylimit * 0.95, by = 0, length.out = nrow(stats))
    
    facet_labeller <- make_facet_labels(fdat, denom_len = 2, 
                                        group_levels = group_levels,
                                        n_label = n_label)
    
    p <- ggplot(fdat, aes(x = trial, y = mean_pre_stim_speed)) +
      geom_boxplot(aes(fill = .data[[fill_var]]), 
                   colour = "black", linewidth = 1, alpha = alpha) +
      ggbeeswarm::geom_beeswarm(aes(colour = response), alpha = 0.6, size = 2, priority = "ascending") +
      fill_scale +
      scale_colour_manual(values = resp_cols) +
      facet_wrap(~ group, nrow = 1, labeller = labeller(group = facet_labeller)) +
      ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
                                 tip.length = 0, size = 8, bracket.size = 1, bracket.nudge.y = 0) +
      labs(title = paste0("Between trials (", trial_between_a, " vs ", trial_between_b, "), stimulus ", stim_between),
           x = "Trial", y = "Pre stimulus speed (mm/s)") +
      coord_cartesian(ylim = c(0, ylimit)) +
      theme_classic() +
      myfacettheme
    
    return(list(plot = p, stats = stats, data = fdat))
  }
  
  stop("Unknown plot_kind")
}

# ---- Normality helpers ----
compute_within_trial_diffs <- function(res, trial_sel, stim_a, stim_b) {
  res %>%
    filter(as.character(trial) == as.character(trial_sel),
           as.character(stim_number) %in% c(as.character(stim_a), as.character(stim_b))) %>%
    select(group, trial, fly_id, stim_number, pi) %>%
    tidyr::pivot_wider(
      id_cols = c(group, trial, fly_id),
      names_from = stim_number,
      values_from = pi,
      names_prefix = "stim_"
    ) %>%
    {
      col_a <- paste0("stim_", as.character(stim_a))
      col_b <- paste0("stim_", as.character(stim_b))
      if (!all(c(col_a, col_b) %in% names(.))) {
        stop("Missing stim columns after pivot: ", paste(setdiff(c(col_a, col_b), names(.)), collapse = ", "))
      }
      mutate(., delta = .data[[col_b]] - .data[[col_a]])
    } %>%
    filter(!is.na(delta))
}

compute_between_trial_diffs <- function(res, stim_sel, trial_a, trial_b) {
  res %>%
    filter(as.character(stim_number) == as.character(stim_sel),
           as.character(trial) %in% c(as.character(trial_a), as.character(trial_b))) %>%
    select(group, trial, fly_id, stim_number, pi) %>%
    tidyr::pivot_wider(
      id_cols = c(group, stim_number, fly_id),
      names_from = trial,
      values_from = pi,
      names_prefix = "trial_"
    ) %>%
    {
      col_a <- paste0("trial_", as.character(trial_a))
      col_b <- paste0("trial_", as.character(trial_b))
      if (!all(c(col_a, col_b) %in% names(.))) {
        stop("Missing trial columns after pivot: ", paste(setdiff(c(col_a, col_b), names(.)), collapse = ", "))
      }
      mutate(., delta = .data[[col_b]] - .data[[col_a]])
    } %>%
    filter(!is.na(delta))
}

make_delta_plots <- function(diffs, palette_source, ggprism_palette, title_prefix = "") {
  fill_scale <- make_fill_scale(palette_source, ggprism_palette)
  
  p_hist <- ggplot(diffs, aes(x = delta)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(title = paste0(title_prefix, "Histogram of paired differences"), x = "Delta", y = "Count")
  
  p_qq <- ggplot(diffs, aes(sample = delta)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ group) +
    theme_classic() + myfacettheme +
    labs(title = paste0(title_prefix, "Q-Q plot of paired differences"),
         x = "Theoretical Quantiles", y = "Sample Quantiles")
  
  p_density <- ggplot(diffs, aes(x = delta)) +
    geom_density() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ group) +
    theme_classic() + myfacettheme +
    labs(title = paste0(title_prefix, "Density of paired differences"), x = "Delta", y = "Density")
  
  p_violin <- ggplot(diffs, aes(x = group, y = delta, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_violin(alpha = 0.85) +
    geom_jitter(width = 0.12, alpha = 0.6) +
    fill_scale +
    theme_classic() +
    labs(title = paste0(title_prefix, "Paired differences by group"), x = "Group", y = "Delta")
  
  list(hist = p_hist, qq = p_qq, density = p_density, violin = p_violin)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("DART analysis"),
  theme = shinytheme("simplex"),
  sidebarLayout(
    sidebarPanel(
      h4("1) Load data"),
      fileInput("xlsx", "Upload Excel (.xlsx)", accept = c(".xlsx")),
      uiOutput("sheet_ui"),
      actionButton("load_btn", "Load file", class = "btn-primary"),
      tags$hr(),
      
      h4("2) Options"),
      checkboxInput("use_ttest", "Use paired t-test (otherwise paired Wilcoxon)", value = FALSE),
      checkboxInput("zeroed", "Zero PI when GOF < 0.2", value = FALSE),
      checkboxInput("n_label", "Show n numbers in facet labels?", value = TRUE),
      tags$hr(),
      
      h4("3) Group order"),
      helpText("Select groups to include (unselected groups are excluded). Drag to reorder."),
      selectizeInput(
        "group_order",
        "Groups (drag to reorder)",
        choices = NULL,
        selected = NULL,
        multiple = TRUE,
        options = list(
          plugins = list("remove_button", "drag_drop"),
          persist = FALSE
        )
      ),
      tags$hr(),
      
      h4("Colours"),
      checkboxInput("col_by_group", "Colour by group (or stimulus if unchecked)?", value = TRUE),
      selectInput(
        "palette_source",
        "Group palette source",
        choices = c("ggprism", "ggplot2 (hue)"),
        selected = "ggprism"
      ),
      uiOutput("ggprism_palette_ui"),
      sliderInput("alpha", "Fill alpha level", min = 0, max = 1, value = 0.8, step = 0.025),
      tags$hr(),
      
      h4("4) Plot controls"),
      numericInput("ylimit", "y-axis upper limit", value = 400 ,min = 50, max = 700, 
                   step = 10, updateOn = "blur"),
      selectInput(
        "plot_kind", "Plot type",
        choices = c(
          "PI: within trial (stim compare)",
          "PI: between trials (trial compare)",
          "Speed: within trial (stim compare)",
          "Speed: between trials (trial compare)"
        )
      ),
      uiOutput("trial_within_ui"),
      uiOutput("trial_between_ui"),
      uiOutput("stim_within_ui"),
      uiOutput("stim_between_ui"),
      
      tags$hr(),
      h4("Normality checks"),
      uiOutput("normality_controls_ui"),
      
      tags$hr(),
      h4("5) Export PDF"),
      textInput("pdf_name", "Filename", value = "DART_plot.pdf"),
      numericInput("pdf_w", "PDF width (in)", value = 7, min = 3, max = 30, step = 0.5),
      numericInput("pdf_h", "PDF height (in)", value = 5, min = 3, max = 30, step = 0.5),
      sliderInput("preview_scale", "Preview scale", min = 0.5, max = 2, value = 1, step = 0.1),
      downloadButton("download_pdf", "Download PDF")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Preview and Plotting",
          uiOutput("status_ui"),
          plotOutput("plot", height = "auto")
        ),
        tabPanel(
          "Stats",
          h4("Pairwise results (holm-adjusted)"),
          tableOutput("stats_tbl")
        ),
        tabPanel(
          "Normality",
          tabsetPanel(
            tabPanel(
              "Summary stats",
              h4("Summary stats (PI)"),
              tableOutput("summary_stats_tbl")
            ),
            tabPanel(
              "Within-trial diffs",
              h4("Within-trial paired differences"),
              tableOutput("shapiro_within_tbl"),
              plotOutput("norm_hist_within", height = "300px"),
              plotOutput("norm_qq_within", height = "350px"),
              plotOutput("norm_density_within", height = "350px"),
              plotOutput("norm_violin_within", height = "350px")
            ),
            tabPanel(
              "Between-trial diffs",
              h4("Between-trial paired differences"),
              tableOutput("shapiro_between_tbl"),
              plotOutput("norm_hist_between", height = "300px"),
              plotOutput("norm_qq_between", height = "350px"),
              plotOutput("norm_density_between", height = "350px"),
              plotOutput("norm_violin_between", height = "350px")
            )
          )
        ),
        tabPanel(
          "Data peek",
          h4("First rows (processed)"),
          tableOutput("data_head")
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  uploaded_path <- reactive({
    req(input$xlsx)
    input$xlsx$datapath
  })
  
  output$sheet_ui <- renderUI({
    req(input$xlsx)
    sheets <- tryCatch(readxl::excel_sheets(uploaded_path()),
                       error = function(e) character(0))
    if (length(sheets) == 0) return(helpText("Could not read sheets from this file."))
    selectInput("sheet", "Sheet", choices = sheets, selected = sheets[[1]])
  })
  
  # Loaded data lives here
  raw_data <- reactiveVal(NULL)
  
  observeEvent(input$load_btn, {
    req(input$xlsx)
    req(input$sheet)
    
    df <- tryCatch(
      {
        readxl::read_excel(uploaded_path(), sheet = input$sheet) |>
          tibble::as_tibble() |>
          janitor::clean_names()
      },
      error = function(e) {
        showNotification(paste("Load failed:", e$message), type = "error", duration = 10)
        NULL
      }
    )
    
    req(!is.null(df))
    raw_data(df)
    showNotification(paste("Loaded", nrow(df), "rows from sheet:", input$sheet),
                     type = "message", duration = 4)
  }, ignoreInit = TRUE)
  
  available_groups <- reactive({
    df <- raw_data()
    req(!is.null(df))
    if (!("group" %in% names(df))) return(character(0))
    sort(unique(as.character(df$group)))
  })
  
  observeEvent(raw_data(), {
    grps <- available_groups()
    if (length(grps) == 0) {
      showNotification("No 'group' column found after loading.", type = "error", duration = 8)
      return()
    }
    updateSelectizeInput(
      session, "group_order",
      choices = grps,
      selected = grps,
      server = TRUE
    )
  }, ignoreInit = TRUE)
  
  output$ggprism_palette_ui <- renderUI({
    req(input$palette_source)
    if (input$palette_source != "ggprism") return(NULL)
    
    selectInput(
      "ggprism_palette",
      "ggprism fill palette",
      choices = GGPRISM_PALS,
      selected = if ("colors" %in% GGPRISM_PALS) "colors" else GGPRISM_PALS[[1]]
    )
  })
  
  output$trial_within_ui <- renderUI({
    df <- raw_data()
    req(!is.null(df))
    if (!("trial" %in% names(df))) return(NULL)
    trials <- sort(unique(as.character(df$trial)))
    selectInput("trial_within", "Within-trial: trial", choices = trials, selected = trials[[1]])
  })
  
  output$trial_between_ui <- renderUI({
    df <- raw_data()
    req(!is.null(df))
    if (!("trial" %in% names(df))) return(NULL)
    
    trials <- sort(unique(as.character(df$trial)))
    if (length(trials) < 2) return(helpText("Need at least 2 trials for between-trial comparisons."))
    
    a_def <- if ("1" %in% trials) "1" else trials[[1]]
    b_def <- if ("2" %in% trials) "2" else trials[[2]]
    
    tagList(
      tags$strong("Between-trials: compare"),
      fluidRow(
        column(6, selectInput("trial_between_a", "Trial A", choices = trials, selected = a_def)),
        column(6, selectInput("trial_between_b", "Trial B", choices = trials, selected = b_def))
      )
    )
  })
  
  
  output$stim_within_ui <- renderUI({
    df <- raw_data()
    req(!is.null(df))
    if (!all(c("stim_number","trial") %in% names(df))) return(NULL)
    
    # if trial input not ready yet, default to first trial available
    trials <- sort(unique(as.character(df$trial)))
    trial_sel <- if (!is.null(input$trial_within)) as.character(input$trial_within) else trials[[1]]
    
    stims <- df %>%
      filter(as.character(trial) == trial_sel) %>%
      pull(stim_number) %>%
      as.character() %>%
      unique() %>%
      sort()
    
    # sensible default selection
    default <- intersect(c("1","6"), stims)
    if (length(default) < 2 && length(stims) >= 2) default <- stims[1:2]
    
    selectizeInput(
      "stim_within", "Within-trial: stim numbers (choose 2+)",
      choices = stims,
      selected = default,
      multiple = TRUE
    )
  })
  
  output$stim_between_ui <- renderUI({
    df <- raw_data()
    req(!is.null(df))
    if (!("stim_number" %in% names(df))) return(NULL)
    stims <- sort(unique(as.character(df$stim_number)))
    selectInput(
      "stim_between", "Between-trials: stimulus",
      choices = stims,
      selected = stims[[1]]
    )
  })
  
  output$normality_controls_ui <- renderUI({
    df <- raw_data()
    req(!is.null(df))
    
    trials <- if ("trial" %in% names(df)) sort(unique(as.character(df$trial))) else character(0)
    stims  <- if ("stim_number" %in% names(df)) sort(unique(as.character(df$stim_number))) else character(0)
    
    if (length(trials) == 0 || length(stims) == 0) return(helpText("Load data to enable normality controls."))
    
    # sensible defaults
    stim_a_def <- if ("1" %in% stims) "1" else stims[[1]]
    stim_b_def <- {
      if ("6" %in% stims) "6"
      else if (length(stims) >= 2) stims[[2]] else stims[[1]]
    }
    trial_def <- if ("1" %in% trials) "1" else trials[[1]]
    trial_a_def <- if ("1" %in% trials) "1" else trials[[1]]
    trial_b_def <- if ("2" %in% trials) "2" else if (length(trials) >= 2) trials[[2]] else trials[[1]]
    
    tagList(
      tags$strong("Within-trial delta (stim B − stim A)"),
      selectInput("norm_trial", "Trial", choices = trials, selected = trial_def),
      fluidRow(
        column(6, selectInput("norm_stim_a", "Stim A", choices = stims, selected = stim_a_def)),
        column(6, selectInput("norm_stim_b", "Stim B", choices = stims, selected = stim_b_def))
      ),
      tags$hr(),
      tags$strong("Between-trial delta (trial B − trial A)"),
      selectInput("norm_stim_between", "Stimulus", choices = stims, selected = stim_a_def),
      fluidRow(
        column(6, selectInput("norm_trial_a", "Trial A", choices = trials, selected = trial_a_def)),
        column(6, selectInput("norm_trial_b", "Trial B", choices = trials, selected = trial_b_def))
      )
    )
  })
  
  group_levels <- reactive({
    df <- raw_data()
    req(!is.null(df))
    
    grps <- available_groups()
    if (is.null(input$group_order)) return(grps)
    
    chosen <- input$group_order
    chosen <- chosen[chosen %in% grps]
    
    validate(need(length(chosen) > 0, "Select at least one group to include."))
    chosen
  })
  
  processed <- reactive({
    df <- raw_data()
    req(!is.null(df))
    
    miss <- setdiff(required_cols, names(df))
    validate(need(length(miss) == 0, paste("Missing required columns:", paste(miss, collapse = ", "))))
    
    # filter out groups not selected
    df <- df %>% filter(group %in% group_levels())
    
    out <- tryCatch(
      prep_data(df, group_levels = group_levels(), zeroed = input$zeroed),
      error = function(e) {
        showNotification(paste("Processing failed:", e$message), type = "error", duration = 10)
        NULL
      }
    )
    out
  })
  
  # ---- Summary stats ----
  summary_stats <- reactive({
    res <- processed()
    req(!is.null(res))
    res %>%
      group_by(group, stim_number, trial) %>%
      rstatix::get_summary_stats(pi)
  })
  
  output$summary_stats_tbl <- renderTable({
    summary_stats()
  })
  
  # ---- Normality: within-trial diffs ----
  diffs_within <- reactive({
    res <- processed()
    req(!is.null(res), !is.null(input$norm_trial), !is.null(input$norm_stim_a), !is.null(input$norm_stim_b))
    
    validate(need(input$norm_stim_a != input$norm_stim_b, "Stim A and Stim B must be different."))
    
    tryCatch(
      compute_within_trial_diffs(res, input$norm_trial, input$norm_stim_a, input$norm_stim_b),
      error = function(e) {
        showNotification(paste("Within-trial diffs failed:", e$message), type = "error", duration = 10)
        NULL
      }
    )
  })
  
  shapiro_within <- reactive({
    d <- diffs_within()
    req(!is.null(d))
    tryCatch(
      d %>% group_by(group) %>% shapiro_test(delta) %>% dplyr::mutate(
        dplyr::across(where(is.numeric), ~ formatC(.x, format = "g", digits = 3))
      ),
      error = function(e) tibble(group = NA, variable = "delta", statistic = NA_real_, p = NA_real_)
    )
  })
  
  within_plots <- reactive({
    d <- diffs_within()
    req(!is.null(d))
    make_delta_plots(
      d,
      palette_source = input$palette_source,
      ggprism_palette = if (!is.null(input$ggprism_palette)) input$ggprism_palette else "colors",
      title_prefix = paste0("Trial ", input$norm_trial, ", Δ(", input$norm_stim_b, "−", input$norm_stim_a, "): ")
    )
  })
  
  output$shapiro_within_tbl <- renderTable({
    shapiro_within()
  })
  
  output$norm_hist_within <- renderPlot({ within_plots()$hist })
  output$norm_qq_within <- renderPlot({ within_plots()$qq })
  output$norm_density_within <- renderPlot({ within_plots()$density })
  output$norm_violin_within <- renderPlot({ within_plots()$violin })
  
  # ---- Normality: between-trial diffs ----
  diffs_between <- reactive({
    res <- processed()
    req(!is.null(res), !is.null(input$norm_stim_between), !is.null(input$norm_trial_a), !is.null(input$norm_trial_b))
    
    validate(need(input$norm_trial_a != input$norm_trial_b, "Trial A and Trial B must be different."))
    
    tryCatch(
      compute_between_trial_diffs(res, input$norm_stim_between, input$norm_trial_a, input$norm_trial_b),
      error = function(e) {
        showNotification(paste("Between-trial diffs failed:", e$message), type = "error", duration = 10)
        NULL
      }
    )
  })
  
  shapiro_between <- reactive({
    d <- diffs_between()
    req(!is.null(d))
    tryCatch(
      d %>% group_by(group) %>% shapiro_test(delta) %>% dplyr::mutate(
        dplyr::across(where(is.numeric), ~ formatC(.x, format = "g", digits = 3))
      ),
      error = function(e) tibble(group = NA, variable = "delta", statistic = NA_real_, p = NA_real_)
    )
  })
  
  between_plots <- reactive({
    d <- diffs_between()
    req(!is.null(d))
    make_delta_plots(
      d,
      palette_source = input$palette_source,
      ggprism_palette = if (!is.null(input$ggprism_palette)) input$ggprism_palette else "colors",
      title_prefix = paste0("Stim ", input$norm_stim_between, ", Δ(trial ", input$norm_trial_b, "−", input$norm_trial_a, "): ")
    )
  })
  
  output$shapiro_between_tbl <- renderTable({
    shapiro_between()
  })
  
  output$norm_hist_between <- renderPlot({ between_plots()$hist })
  output$norm_qq_between <- renderPlot({ between_plots()$qq })
  output$norm_density_between <- renderPlot({ between_plots()$density })
  output$norm_violin_between <- renderPlot({ between_plots()$violin })
  
  
  all_stats <- reactive({
    df <- processed()
    req(!is.null(df))
    
    req(input$trial_within, input$stim_within, input$stim_between)
    
    validate(need(length(input$stim_within) >= 2, "Select at least 2 stim numbers for within-trial comparison."))
    
    test_kind <- if (isTRUE(input$use_ttest)) "t" else "wilcox"
    req(input$trial_between_a, input$trial_between_b)
    validate(need(input$trial_between_a != input$trial_between_b, "Trial A and Trial B must be different."))
    
    compute_all_stats(
      df = df,
      test_kind = test_kind,
      trial_within = input$trial_within,
      stim_within  = input$stim_within,
      stim_between = input$stim_between,
      trial_between_a = input$trial_between_a,
      trial_between_b = input$trial_between_b,
      p_adjust_method = "holm", # or bonferroni
      correct_across  = "within_endpoint" # or across_endpoints
    )
  })
  
  
  # ---- Main plot bundle ----
  plot_bundle <- reactive({
    df <- processed()
    req(!is.null(df))
    
    test_kind <- if (isTRUE(input$use_ttest)) "t" else "wilcox"
    
    make_plot(
      df = df,
      alpha = input$alpha,
      n_label = input$n_label,
      ylimit = input$ylimit,
      col_by_group = input$col_by_group,
      plot_kind = input$plot_kind,
      test_kind = test_kind,
      zeroed = input$zeroed,
      stim_within = input$stim_within,
      trial_within = input$trial_within,
      stim_between = input$stim_between,
      trial_between_a = input$trial_between_a,
      trial_between_b = input$trial_between_b,
      palette_source = input$palette_source,
      ggprism_palette = if (!is.null(input$ggprism_palette)) input$ggprism_palette else "colors",
      stats_all = all_stats()
    )
  })
  
  
  output$status_ui <- renderUI({
    df <- raw_data()
    req(!is.null(df))
    tags$div(
      tags$p(tags$b("Loaded:"), input$xlsx$name),
      tags$p(tags$b("Included groups:"), paste(group_levels(), collapse = ", "))
    )
  })
  
  output$plot <- renderPlot({
    pb <- plot_bundle()
    pb$plot
  },
  width  = function() round(input$pdf_w * 96 * input$preview_scale),
  height = function() round(input$pdf_h * 96 * input$preview_scale),
  res = 96
  )
  
  output$stats_tbl <- renderTable({
    pb <- plot_bundle()
    
    pb$stats %>%
      as.data.frame() %>%
      dplyr::select(group, everything()) %>%
      dplyr::mutate(
        dplyr::across(where(is.numeric), ~ formatC(.x, format = "g", digits = 3))
      )
  }, striped = TRUE, na = "")
  
  
  output$data_head <- renderTable({
    df <- processed()
    req(!is.null(df))
    head(df, 12)
  })
  
  output$download_pdf <- downloadHandler(
    filename = function() {
      nm <- input$pdf_name
      if (!grepl("\\.pdf$", nm, ignore.case = TRUE)) nm <- paste0(nm, ".pdf")
      nm
    },
    content = function(file) {
      pb <- plot_bundle()
      ggsave(
        filename = file,
        plot = pb$plot,
        device = pdf_device_fun(),
        width = input$pdf_w,
        height = input$pdf_h,
        units = "in"
      )
    }
  )
}

shinyApp(ui, server)
