## Figure 1bcd: PMI-related Metabolite Cluster Heatmap with
##           Representative Time-Series Line Plots
## Five-tissue postmortem metabolomics analysis


# =========================================================
# Section 0: Output directory setup
# =========================================================
desktop_path <- ifelse(Sys.info()["sysname"] == "Windows",
                       file.path(Sys.getenv("USERPROFILE"), "Desktop"),
                       file.path(Sys.getenv("HOME"), "Desktop"))
output_dir <- file.path(desktop_path,
                        paste0("PMI_Project_Output_",
                               format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Section 1: Load required packages
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(mgcv)
  library(ggplot2)
  library(VennDiagram)
  library(RColorBrewer)
  library(cowplot)
  library(grid)
  library(gtable)
  library(ragg)      # Arial font support via system font engine
})

# Explicitly bind to ggplot2 namespace to avoid conflicts
margin       <- ggplot2::margin
element_text <- ggplot2::element_text

# =========================================================
# Section 2: Load input data
# =========================================================
# Select: (1) metabolite intensity matrix CSV
# Select: (2) sample metadata CSV (columns: SampleID, Tissue, AnimalID, PMI_day)
metab <- read.csv(file.choose(), check.names = FALSE)
meta  <- read.csv(file.choose(), check.names = FALSE)

metab <- metab %>% rename(SampleID = 1)
common_ids <- intersect(metab$SampleID, meta$SampleID)
metab2 <- metab %>% filter(SampleID %in% common_ids)
meta2  <- meta  %>% filter(SampleID %in% common_ids) %>% arrange(SampleID)
metab2 <- metab2[match(meta2$SampleID, metab2$SampleID), , drop = FALSE]

X_all   <- metab2 %>% select(-SampleID)
tissues <- sort(unique(meta2$Tissue))
stopifnot(all(c("SampleID", "Tissue", "AnimalID", "PMI_day") %in% colnames(meta2)))

# Select: (3) metabolite class annotation CSV
#   Column 1: metabolite ID | Column 2: chemical class | Column 3 (optional): display name
metab_class_df <- read.csv(file.choose(), check.names = FALSE)
colnames(metab_class_df)[1:2] <- c("metabolite", "metab_class")

if (ncol(metab_class_df) >= 3) {
  colnames(metab_class_df)[3] <- "metab_name"
} else {
  metab_class_df$metab_name <- NA_character_
}

metab_class_df <- metab_class_df %>%
  select(metabolite, metab_class, metab_name) %>%
  distinct() %>%
  mutate(
    metab_class = ifelse(is.na(metab_class) | metab_class == "",
                         "Unknown", metab_class),
    metab_name  = ifelse(is.na(metab_name) | metab_name == "",
                         metabolite, metab_name)
  )

# =========================================================
# Section 3: Analysis parameters
# =========================================================
gam_k           <- 4
global_pct_cut  <- 2.0
cluster_pct_cut <- 3.0
min_count_cut   <- 3L
PLOT_HEATMAP    <- TRUE
N_REP           <- 1
ERROR_TYPE      <- "se"   # "se" = standard error; "sd" = standard deviation

# =========================================================
# Section 4: Function definitions
# =========================================================

# ----------------------------------------------------------
# 4.1 Variance-based feature filtering
# ----------------------------------------------------------
filter_features <- function(X_mat) {
  keep <- apply(X_mat, 2, function(v) sd(v, na.rm = TRUE) > 0)
  colnames(X_mat)[keep]
}

# ----------------------------------------------------------
# 4.2 Spearman correlation + GAM association test
# ----------------------------------------------------------
assoc_full_spearman_gam <- function(meta_df, X_df, gam_k = 4) {
  PMI  <- meta_df$PMI_day
  mets <- colnames(X_df)
  
  sp <- lapply(mets, function(m) {
    v  <- X_df[[m]]
    ct <- suppressWarnings(
      cor.test(v, PMI, method = "spearman", exact = FALSE))
    data.frame(metabolite = m,
               rho        = unname(ct$estimate),
               p_spearman = ct$p.value)
  }) %>% bind_rows() %>%
    mutate(fdr_spearman = p.adjust(p_spearman, "BH"))
  
  gm <- lapply(mets, function(m) {
    y   <- X_df[[m]]
    fit <- tryCatch(
      mgcv::gam(y ~ s(PMI, k = gam_k), method = "REML"),
      error = function(e) NULL)
    if (is.null(fit))
      return(data.frame(metabolite = m,
                        p_gam = NA_real_, r2_gam = NA_real_))
    sm   <- summary(fit)
    pval <- tryCatch(sm$s.table[1, "p-value"],
                     error = function(e) NA_real_)
    data.frame(metabolite = m, p_gam = pval, r2_gam = sm$r.sq)
  }) %>% bind_rows() %>%
    mutate(fdr_gam = p.adjust(p_gam, "BH"))
  
  assoc <- sp %>% left_join(gm, by = "metabolite") %>%
    mutate(sig_spearman = !is.na(fdr_spearman) & fdr_spearman < 0.05,
           sig_gam      = !is.na(fdr_gam)      & fdr_gam      < 0.05,
           sig_any      = sig_spearman | sig_gam,
           fdr_min      = pmin(fdr_spearman, fdr_gam, na.rm = TRUE))
  
  feats_PMI_related_full <- assoc %>%
    filter(sig_any) %>% pull(metabolite) %>% unique()
  list(assoc = assoc, feats_PMI_related_full = feats_PMI_related_full)
}

# ----------------------------------------------------------
# 4.3 Three-cluster classification
#     Cluster 1: Monotone Increase
#     Cluster 2: Monotone Decrease
#     Cluster 3: Non-Monotone
# ----------------------------------------------------------
classify_three_clusters <- function(assoc_tbl, meta_df, X_df,
                                    gam_k = 4) {
  PMI      <- meta_df$PMI_day
  sig_mets <- assoc_tbl %>% filter(sig_any) %>% pull(metabolite)
  
  classify_one <- function(m) {
    row        <- assoc_tbl %>% filter(metabolite == m)
    is_sig_sp  <- isTRUE(row$sig_spearman[1])
    is_sig_gam <- isTRUE(row$sig_gam[1])
    rho        <- row$rho[1]
    
    if (is_sig_sp) {
      cluster <- if (!is.na(rho) && rho > 0) "Monotone_Increase"
      else "Monotone_Decrease"
      if (is_sig_gam && m %in% colnames(X_df)) {
        y   <- X_df[[m]]
        fit <- tryCatch(
          mgcv::gam(y ~ s(PMI, k = gam_k), method = "REML"),
          error = function(e) NULL)
        if (!is.null(fit)) {
          pmi_seq   <- seq(min(PMI, na.rm = TRUE),
                           max(PMI, na.rm = TRUE), length.out = 100)
          pred      <- as.numeric(
            predict(fit, newdata = data.frame(PMI = pmi_seq)))
          diffs     <- diff(pred)
          pos_ratio <- sum(diffs > 0) / length(diffs)
          neg_ratio <- sum(diffs < 0) / length(diffs)
          if ((rho > 0 && neg_ratio > 0.6) ||
              (rho < 0 && pos_ratio > 0.6))
            cluster <- "Non_Monotone"
        }
      }
    } else if (is_sig_gam) {
      cluster <- "Non_Monotone"
    } else {
      cluster <- "Non_Monotone"
    }
    
    data.frame(metabolite   = m, cluster = cluster, rho = rho,
               sig_spearman = is_sig_sp, sig_gam = is_sig_gam,
               fdr_spearman = row$fdr_spearman[1],
               fdr_gam      = row$fdr_gam[1],
               stringsAsFactors = FALSE)
  }
  
  result <- bind_rows(lapply(sig_mets, classify_one))
  result$cluster <- factor(result$cluster,
                           levels = c("Monotone_Increase",
                                      "Monotone_Decrease",
                                      "Non_Monotone"))
  result
}

# ----------------------------------------------------------
# 4.4 Select representative metabolites per cluster
# ----------------------------------------------------------
select_representative_metabolites <- function(cluster_df, assoc_tbl,
                                              metab_class_df,
                                              n_per_cluster = 1) {
  scored <- cluster_df %>%
    left_join(assoc_tbl %>% select(metabolite, r2_gam),
              by = "metabolite") %>%
    left_join(metab_class_df %>%
                select(metabolite, metab_class, metab_name),
              by = "metabolite") %>%
    mutate(
      metab_class = ifelse(is.na(metab_class) | metab_class == "",
                           "Unknown", metab_class),
      metab_name  = ifelse(is.na(metab_name)  | metab_name  == "",
                           metabolite, metab_name),
      r2_gam      = ifelse(is.na(r2_gam), 0, r2_gam),
      abs_rho     = abs(rho),
      sort_key    = case_when(
        cluster == "Non_Monotone" ~ r2_gam,
        TRUE                     ~ abs_rho
      )
    )
  
  scored %>%
    group_by(cluster) %>%
    arrange(desc(sort_key), .by_group = TRUE) %>%
    slice_head(n = n_per_cluster) %>%
    ungroup() %>%
    select(metabolite, metab_name, cluster,
           abs_rho, r2_gam, metab_class, sort_key)
}

# =========================================================
# 4.5 Figure 1b: Heatmap + bridge + time-series composite
# =========================================================
plot_heatmap_with_timeseries <- function(meta_t, X_t, cluster_df,
                                         assoc_tbl, metab_class_df,
                                         tissue_name,
                                         n_rep      = 1,
                                         error_type = "se",
                                         save_path  = NULL,
                                         fig_width  = 19,
                                         fig_height = 9,
                                         res_dpi    = 300) {
  
  CLUSTER_ORDER  <- c("Monotone_Increase", "Monotone_Decrease", "Non_Monotone")
  CLUSTER_COLORS <- c(Monotone_Increase = "#C0392B",
                      Monotone_Decrease = "#2980B9",
                      Non_Monotone      = "#27AE60")
  CLUSTER_LABELS <- c(Monotone_Increase = "Monotone Increase",
                      Monotone_Decrease = "Monotone Decrease",
                      Non_Monotone      = "Non-Monotone")
  PMI_PALETTE <- c("#B39DFF", "#FF9A8B", "#FF59D1",
                   "#7CFC00", "#00C8FF", "#00D7A7",
                   "#FFA500", "#FF4500")
  
  # Step 1: Representative metabolites
  rep_mets_df <- select_representative_metabolites(
    cluster_df, assoc_tbl, metab_class_df, n_per_cluster = n_rep
  )
  if (!is.null(save_path))
    write.csv(rep_mets_df,
              sub("\\.png$", "_RepresentativeMetabolites.csv", save_path),
              row.names = FALSE)
  
  # Step 2: Heatmap matrix (Z-scored, ordered by PMI)
  cluster_df_plot <- cluster_df %>%
    filter(metabolite %in% colnames(X_t)) %>%
    mutate(cluster = factor(as.character(cluster),
                            levels = CLUSTER_ORDER)) %>%
    arrange(cluster)
  
  feats    <- cluster_df_plot$metabolite
  ord      <- order(meta_t$PMI_day, meta_t$AnimalID)
  meta_ord <- meta_t[ord, , drop = FALSE]
  mat      <- as.matrix(X_t[ord, feats, drop = FALSE])
  rownames(mat) <- meta_ord$SampleID
  
  sds  <- apply(mat, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  mat  <- mat[, keep, drop = FALSE]
  cluster_df_plot <- cluster_df_plot %>% filter(metabolite %in% colnames(mat))
  if (ncol(mat) < 2) { return(invisible(NULL)) }
  
  mat_t <- t(mat)
  mat_z <- t(apply(mat_t, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (s == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / s
  }))
  mat_z[mat_z >  3] <-  3
  mat_z[mat_z < -3] <- -3
  
  met_order    <- cluster_df_plot$metabolite
  met_order    <- met_order[met_order %in% rownames(mat_z)]
  n_mets       <- length(met_order)
  sample_order <- meta_ord$SampleID
  n_samples    <- length(sample_order)
  
  cl_sizes <- cluster_df_plot %>%
    filter(metabolite %in% met_order) %>%
    mutate(cluster = factor(cluster, levels = CLUSTER_ORDER)) %>%
    group_by(cluster, .drop = FALSE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 0)
  
  met_idx <- setNames(seq_along(met_order), met_order)
  
  tile_df <- expand.grid(
    metabolite = met_order,
    SampleID   = sample_order,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      row_idx = met_idx[metabolite],
      col_idx = match(SampleID, sample_order),
      z       = mapply(function(m, s) mat_z[m, s],
                       metabolite, SampleID)
    )
  
  pmi_levels <- sort(unique(meta_ord$PMI_day))
  pmi_labels <- paste0("Day", pmi_levels)
  pmi_colors <- setNames(PMI_PALETTE[seq_along(pmi_levels)], pmi_labels)
  
  ann_bar_df <- data.frame(
    col_idx = seq_len(n_samples),
    PMI_day = factor(paste0("Day", meta_ord$PMI_day), levels = pmi_labels)
  )
  
  cluster_ann_df <- cluster_df_plot %>%
    filter(metabolite %in% met_order) %>%
    mutate(row_idx = met_idx[metabolite])
  
  cum_sizes     <- cumsum(cl_sizes$n)
  gap_positions <- cum_sizes[-length(cum_sizes)] + 0.5
  
  # Step 3: Heatmap ggplot
  p_heatmap <- ggplot() +
    geom_tile(data = tile_df,
              aes(x = col_idx, y = row_idx, fill = z),
              width = 1, height = 1) +
    scale_fill_gradientn(
      colors = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
      limits = c(-3, 3), guide = "none"
    ) +
    geom_tile(data = ann_bar_df,
              aes(x = col_idx, y = -1.5, fill = NULL),
              fill = pmi_colors[as.character(ann_bar_df$PMI_day)],
              width = 1, height = 2, color = NA) +
    geom_tile(data = cluster_ann_df,
              aes(x = -2, y = row_idx, fill = NULL),
              fill = CLUSTER_COLORS[as.character(cluster_ann_df$cluster)],
              width = 3, height = 1, color = NA) +
    {if (length(gap_positions) > 0)
      geom_hline(yintercept = gap_positions,
                 color = "white", linewidth = 1.5)} +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(-4.5, n_samples + 0.5)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(-3.5, n_mets + 0.5),
                       trans  = "reverse") +
    labs(title = paste0(tissue_name,
                        ": PMI-related metabolites (3 clusters)")) +
    theme_void(base_size = 10, base_family = "Arial") +
    theme(
      text            = element_text(family = "Arial"),
      plot.title      = element_text(face = "bold", size = 11,
                                     hjust = 0.5, family = "Arial",
                                     margin = margin(b = 4)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = ggplot2::margin(t = 8, r = 0, b = 4, l = 4)
    )
  
  # Step 4: Time-series line plots (equal height, enlarged fonts)
  cl_nz   <- cl_sizes %>% filter(n > 0)
  n_plots <- nrow(cl_nz)
  avail_mets  <- intersect(rep_mets_df$metabolite, colnames(X_t))
  rep_mets_ok <- rep_mets_df %>% filter(metabolite %in% avail_mets)
  TS_MARGIN   <- ggplot2::margin(t = 6, r = 8, b = 6, l = 6)
  
  ts_plots <- lapply(seq_len(n_plots), function(i) {
    cl_name  <- as.character(cl_nz$cluster[i])
    sub_rep  <- rep_mets_ok %>% filter(cluster == cl_name)
    cl_color <- CLUSTER_COLORS[cl_name]
    
    if (nrow(sub_rep) == 0) {
      return(
        ggplot() +
          theme_void(base_family = "Arial") +
          theme(plot.background = element_rect(fill = "white", color = NA),
                plot.margin     = TS_MARGIN)
      )
    }
    
    m_id      <- sub_rep$metabolite[1]
    met_name  <- sub_rep$metab_name[1]
    met_class <- sub_rep$metab_class[1]
    
    ts_data <- data.frame(
      PMI_day = meta_t$PMI_day,
      value   = as.numeric(X_t[, m_id])
    ) %>%
      group_by(PMI_day) %>%
      summarise(mean_val = mean(value, na.rm = TRUE),
                sd_val   = sd(value, na.rm = TRUE),
                n        = sum(!is.na(value)),
                .groups  = "drop") %>%
      mutate(
        se_val = sd_val / sqrt(n),
        err    = if (error_type == "se") se_val else sd_val,
        ymin   = mean_val - err,
        ymax   = mean_val + err
      )
    
    ggplot(ts_data, aes(x = PMI_day, y = mean_val)) +
      geom_ribbon(aes(ymin = ymin, ymax = ymax),
                  fill = cl_color, alpha = 0.18) +
      geom_line(color = cl_color, linewidth = 0.9) +
      geom_point(shape = 23, size = 2.5,
                 color = cl_color, fill = "white", stroke = 1.1) +
      scale_x_continuous(breaks = sort(unique(ts_data$PMI_day)),
                         expand = c(0.08, 0.08)) +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      labs(title    = met_name,
           subtitle = paste0("(", met_class, ")"),
           x        = "PMI (days)",
           y        = "Intensity") +
      theme_classic(base_size = 14, base_family = "Arial") +
      theme(
        text          = element_text(family = "Arial"),
        panel.border  = element_rect(color = "grey30",
                                     linewidth = 0.6, fill = NA),
        plot.title    = element_text(face = "bold.italic", size = 13,
                                     color = cl_color, hjust = 0.5,
                                     family = "Arial",
                                     margin = ggplot2::margin(b = 1)),
        plot.subtitle = element_text(size = 11.5, color = "grey45",
                                     hjust = 0.5, family = "Arial",
                                     margin = ggplot2::margin(b = 2)),
        axis.title    = element_text(size = 12, color = "grey20",
                                     family = "Arial"),
        axis.text     = element_text(size = 11.5, color = "grey20",
                                     family = "Arial"),
        axis.line     = element_blank(),
        axis.ticks    = element_line(linewidth = 0.35, color = "grey50"),
        plot.margin   = TS_MARGIN
      )
  })
  
  ts_panel <- cowplot::plot_grid(
    plotlist    = ts_plots,
    ncol        = 1,
    align       = "v",
    axis        = "lr",
    rel_heights = rep(1, n_plots)
  )
  
  # Step 5: Bridge connector
  cum_n        <- c(0, cumsum(cl_nz$n))
  left_tops    <- 1 - cum_n[-length(cum_n)] / n_mets
  left_bottoms <- 1 - cum_n[-1]             / n_mets
  right_tops    <- 1 - (seq_len(n_plots) - 1) / n_plots
  right_bottoms <- 1 - seq_len(n_plots)        / n_plots
  
  poly_df <- bind_rows(lapply(seq_len(n_plots), function(i) {
    data.frame(
      x   = c(0, 0, 1, 1),
      y   = c(left_bottoms[i], left_tops[i],
              right_tops[i],   right_bottoms[i]),
      grp = as.character(cl_nz$cluster[i]),
      stringsAsFactors = FALSE
    )
  }))
  
  side_df <- bind_rows(lapply(seq_len(n_plots), function(i) {
    data.frame(
      grp     = as.character(cl_nz$cluster[i]),
      y_top_l = left_tops[i],    y_bot_l = left_bottoms[i],
      y_top_r = right_tops[i],   y_bot_r = right_bottoms[i],
      stringsAsFactors = FALSE
    )
  }))
  
  p_bridge <- ggplot() +
    geom_polygon(data = poly_df,
                 aes(x = x, y = y, group = grp, fill = grp),
                 alpha = 0.18, color = NA) +
    geom_segment(data = side_df,
                 aes(x = 0, xend = 0,
                     y = y_bot_l, yend = y_top_l, color = grp),
                 linewidth = 2.5) +
    geom_segment(data = side_df,
                 aes(x = 1, xend = 1,
                     y = y_bot_r, yend = y_top_r, color = grp),
                 linewidth = 1.2) +
    scale_fill_manual(values  = CLUSTER_COLORS, guide = "none") +
    scale_color_manual(values = CLUSTER_COLORS, guide = "none") +
    scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1),       expand = c(0, 0)) +
    theme_void(base_family = "Arial") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin     = margin(0, 0, 0, 0))
  
  # Step 6: Legend panel
  legend_panel <- build_legend_panel_v2(
    pmi_levels     = pmi_levels,
    pmi_palette    = PMI_PALETTE,
    cluster_labels = CLUSTER_LABELS,
    cluster_colors = CLUSTER_COLORS
  )
  
  # Step 7: Assemble and save with ragg device
  bridge_ts <- cowplot::plot_grid(
    p_bridge, ts_panel,
    ncol       = 2,
    rel_widths = c(0.12, 0.88),
    align      = "h",
    axis       = "tb"
  )
  
  combined <- cowplot::plot_grid(
    p_heatmap, bridge_ts, legend_panel,
    ncol       = 3,
    rel_widths = c(0.62, 0.28, 0.10),
    align      = "h",
    axis       = "tb"
  )
  
  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot     = combined,
      width    = fig_width,
      height   = fig_height,
      units    = "in",
      dpi      = res_dpi,
      device   = ragg::agg_png
    )
  }
  
  invisible(rep_mets_df)
}

# ----------------------------------------------------------
# Legend panel builder
# ----------------------------------------------------------
build_legend_panel_v2 <- function(pmi_levels, pmi_palette,
                                  cluster_labels, cluster_colors) {
  
  grad_df <- data.frame(y = seq(-3, 3, length.out = 100), x = 1)
  p_zscore <- ggplot(grad_df, aes(x = x, y = y, fill = y)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
      limits = c(-3, 3),
      breaks = c(-3, -1.5, 0, 1.5, 3),
      name   = "Z-score"
    ) +
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth  = unit(0.5, "cm"),
      barheight = unit(2.5, "cm"),
      ticks.linewidth = 0.5
    )) +
    theme_void(base_family = "Arial") +
    theme(
      text            = element_text(family = "Arial"),
      legend.position = "right",
      legend.title    = element_text(size = 12, face = "bold",
                                     family = "Arial"),
      legend.text     = element_text(size = 11, family = "Arial")
    )
  
  cl_df <- data.frame(
    label = unname(cluster_labels),
    color = unname(cluster_colors),
    y     = rev(seq_along(cluster_labels)),
    stringsAsFactors = FALSE
  )
  p_cluster <- ggplot(cl_df) +
    geom_tile(aes(x = 1, y = y), fill = cl_df$color,
              width = 0.6, height = 0.7) +
    geom_text(aes(x = 1.55, y = y, label = label),
              hjust = 0, size = 3.8, color = "grey15",
              family = "Arial") +
    scale_x_continuous(limits = c(0.5, 5.5)) +
    scale_y_continuous(limits = c(0.2, length(cluster_labels) + 0.8)) +
    labs(title = "Cluster") +
    theme_void(base_family = "Arial") +
    theme(
      text       = element_text(family = "Arial"),
      plot.title = element_text(size = 10, face = "bold",
                                hjust = 0, family = "Arial",
                                margin = ggplot2::margin(b = 3))
    )
  
  pmi_labels_vec <- paste0("Day", pmi_levels)
  pmi_df <- data.frame(
    label = pmi_labels_vec,
    color = pmi_palette[seq_along(pmi_levels)],
    y     = rev(seq_along(pmi_levels)),
    stringsAsFactors = FALSE
  )
  p_pmi <- ggplot(pmi_df) +
    geom_tile(aes(x = 1, y = y), fill = pmi_df$color,
              width = 0.6, height = 0.7) +
    geom_text(aes(x = 1.55, y = y, label = label),
              hjust = 0, size = 2.5, color = "grey15",
              family = "Arial") +
    scale_x_continuous(limits = c(0.5, 4.5)) +
    scale_y_continuous(limits = c(0.2, length(pmi_levels) + 0.8)) +
    labs(title = "PMI (days)") +
    theme_void(base_family = "Arial") +
    theme(
      text       = element_text(family = "Arial"),
      plot.title = element_text(size = 8, face = "bold",
                                hjust = 0, family = "Arial",
                                margin = margin(b = 3))
    )
  
  legend_zscore <- cowplot::get_plot_component(
    p_zscore, "guide-box-right", return_all = FALSE
  )
  cowplot::plot_grid(
    legend_zscore, p_cluster, p_pmi,
    ncol        = 1,
    rel_heights = c(0.28, 0.22, 0.50),
    align       = "v",
    axis        = "lr"
  )
}

# ----------------------------------------------------------
# 4.6 Per-tissue metabolite class composition analysis
# ----------------------------------------------------------
analyze_class_composition_csv <- function(cluster_df, metab_class_df,
                                          tissue_name, output_dir) {
  CLUSTER_ORDER <- c("Monotone_Increase", "Monotone_Decrease", "Non_Monotone")
  
  df <- cluster_df %>%
    left_join(metab_class_df %>% select(metabolite, metab_class),
              by = "metabolite") %>%
    mutate(
      metab_class = ifelse(is.na(metab_class) | metab_class == "",
                           "Unknown", metab_class),
      cluster     = factor(cluster, levels = CLUSTER_ORDER)
    )
  
  write.csv(df %>% arrange(cluster, metab_class, metabolite),
            file.path(output_dir,
                      paste0("ClusterAssignment_", tissue_name, ".csv")),
            row.names = FALSE)
  
  comp_tbl <- df %>%
    group_by(cluster, metab_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(pct = n / sum(n) * 100, total = sum(n)) %>%
    ungroup()
  
  write.csv(comp_tbl,
            file.path(output_dir,
                      paste0("ClusterClassCounts_", tissue_name, ".csv")),
            row.names = FALSE)
  
  fisher_res <- expand.grid(cluster     = CLUSTER_ORDER,
                            metab_class = unique(df$metab_class),
                            stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      a = sum(df$cluster == cluster & df$metab_class == metab_class),
      b = sum(df$cluster == cluster & df$metab_class != metab_class),
      c = sum(df$cluster != cluster & df$metab_class == metab_class),
      d = sum(df$cluster != cluster & df$metab_class != metab_class),
      p_fisher = tryCatch(
        fisher.test(matrix(c(a, b, c, d), 2, 2),
                    alternative = "greater")$p.value,
        error = function(e) NA_real_),
      odds_ratio = tryCatch(
        fisher.test(matrix(c(a, b, c, d), 2, 2))$estimate,
        error = function(e) NA_real_)
    ) %>%
    ungroup() %>%
    mutate(fdr_fisher = p.adjust(p_fisher, "BH")) %>%
    arrange(cluster, fdr_fisher)
  
  write.csv(fisher_res,
            file.path(output_dir,
                      paste0("ClusterClass_FisherTest_",
                             tissue_name, ".csv")),
            row.names = FALSE)
  
  invisible(comp_tbl)
}

# ----------------------------------------------------------
# 4.7 Figure 1c: Five-tissue integrated class composition bar plot
# ----------------------------------------------------------
plot_combined_class_composition <- function(all_comp_list, output_dir) {
  CLUSTER_ORDER  <- c("Monotone_Increase", "Monotone_Decrease", "Non_Monotone")
  CLUSTER_LABELS <- c(Monotone_Increase = "Monotone Increase",
                      Monotone_Decrease = "Monotone Decrease",
                      Non_Monotone      = "Non-Monotone")
  TISSUE_ORDER   <- c("Heart", "Liver", "Lung", "Muscle", "Spleen")
  
  combined <- bind_rows(lapply(names(all_comp_list), function(tis) {
    all_comp_list[[tis]] %>% mutate(tissue = tis)
  }))
  
  global_comp <- combined %>%
    group_by(metab_class) %>%
    summarise(global_n = sum(n), .groups = "drop") %>%
    mutate(global_pct = global_n / sum(global_n) * 100)
  
  keep_global  <- global_comp %>%
    filter(global_pct >= global_pct_cut,
           global_n   >= min_count_cut) %>%
    pull(metab_class)
  keep_cluster <- combined %>%
    filter(pct >= cluster_pct_cut) %>%
    pull(metab_class) %>% unique()
  keep_classes <- union(keep_global, keep_cluster)
  
  plot_df <- combined %>%
    mutate(class_plot = ifelse(metab_class %in% keep_classes,
                               metab_class, "Other")) %>%
    group_by(tissue, cluster, class_plot) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(tissue, cluster) %>%
    mutate(pct = n / sum(n) * 100, total = sum(n)) %>%
    ungroup() %>%
    mutate(
      cluster       = factor(cluster, levels = CLUSTER_ORDER),
      cluster_label = CLUSTER_LABELS[as.character(cluster)],
      tissue        = factor(tissue,
                             levels = intersect(TISSUE_ORDER,
                                                unique(tissue)))
    )
  
  pal_use <- GLOBAL_CLASS_COLORS[unique(plot_df$class_plot)]
  pal_use[is.na(pal_use)] <- "#CCCCCC"
  names(pal_use) <- unique(plot_df$class_plot)
  
  p <- ggplot(plot_df,
              aes(x = cluster_label, y = pct, fill = class_plot)) +
    geom_bar(stat = "identity", position = "stack",
             width = 0.7, color = "white", linewidth = 0.15) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_fill_manual(values = pal_use, name = "Metabolite Class") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    labs(title = "Metabolite class composition across tissues and clusters",
         x = NULL, y = "Percentage (%)") +
    theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      text               = element_text(family = "Arial"),
      legend.position    = "right",
      legend.key.size    = unit(0.4, "cm"),
      legend.text        = element_text(size = 9, family = "Arial"),
      axis.text.x        = element_text(size = 9, angle = 20, hjust = 1,
                                        family = "Arial"),
      strip.text         = element_text(face = "bold", size = 11,
                                        family = "Arial"),
      strip.background   = element_rect(fill = "#F0F0F0", color = NA),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(face = "bold", size = 13,
                                        family = "Arial")
    )
  
  ggsave(
    filename = file.path(output_dir,
                         "FigureS_ClassComposition_AllTissues.png"),
    plot   = p,
    width  = 16,
    height = 6,
    units  = "in",
    dpi    = 300,
    device = ragg::agg_png
  )
  
  invisible(p)
}

# =========================================================
# Section 3b: Global class color palette
# =========================================================
all_classes_global <- metab_class_df %>%
  pull(metab_class) %>% unique() %>% sort()
n_cls <- length(all_classes_global)

soft_palette <- c(
  "#6FA8DC", "#E78B6C", "#6FBF73", "#C27BCB", "#D9B44A",
  "#7FBF9A", "#D98C4A", "#6E9FD6", "#C79A5A", "#8C8CD9",
  "#D97A9A", "#5FC9B2", "#BFAF4A", "#B06BB0", "#68B08A",
  "#D46A6A", "#7CAEE6", "#D9A066", "#7BC67B", "#C88AB8",
  "#A9C95A", "#8FA1C8", "#F09A62", "#69C7C7", "#C9A36B"
)
if (n_cls > length(soft_palette))
  soft_palette <- colorRampPalette(soft_palette)(n_cls)

GLOBAL_CLASS_COLORS <- setNames(soft_palette[seq_len(n_cls)],
                                all_classes_global)
GLOBAL_CLASS_COLORS["Other"]   <- "#D8D8D8"
GLOBAL_CLASS_COLORS["Unknown"] <- "#E8E8E8"

# =========================================================
# Section 5: Main analysis loop (per tissue)
# =========================================================
all_out       <- setNames(vector("list", length(tissues)), tissues)
all_comp_list <- list()

for (tis in tissues) {
  cat("\n=============================\n")
  cat("Tissue:", tis, "\n")
  cat("=============================\n")
  
  idx_t  <- meta2$Tissue == tis
  meta_t <- meta2[idx_t, , drop = FALSE]
  X_t0   <- X_all[idx_t, , drop = FALSE]
  
  keep_full <- filter_features(X_t0)
  X_full    <- X_t0[, keep_full, drop = FALSE]
  
  out_full   <- assoc_full_spearman_gam(meta_t,
                                        as.data.frame(X_full),
                                        gam_k = gam_k)
  assoc_tbl  <- out_full$assoc %>% arrange(fdr_min)
  feats_full <- out_full$feats_PMI_related_full
  
  write.csv(assoc_tbl,
            file.path(output_dir,
                      paste0("Assoc_SpearmanGAM_", tis, ".csv")),
            row.names = FALSE)
  write.csv(data.frame(metabolite = feats_full),
            file.path(output_dir,
                      paste0("PMI_related_", tis, ".csv")),
            row.names = FALSE)
  
  cluster_df <- classify_three_clusters(assoc_tbl, meta_t,
                                        as.data.frame(X_full),
                                        gam_k = gam_k)
  write.csv(cluster_df,
            file.path(output_dir,
                      paste0("ClusterAssignment_raw_", tis, ".csv")),
            row.names = FALSE)
  
  if (PLOT_HEATMAP) {
    plot_heatmap_with_timeseries(
      meta_t         = meta_t,
      X_t            = X_full,
      cluster_df     = cluster_df,
      assoc_tbl      = assoc_tbl,
      metab_class_df = metab_class_df,
      tissue_name    = tis,
      n_rep          = N_REP,
      error_type     = ERROR_TYPE,
      save_path      = file.path(output_dir,
                                 paste0("Figure2_Heatmap_Timeseries_",
                                        tis, ".TIFF")),
      fig_width      = 19,
      fig_height     = 9,
      res_dpi        = 300
    )
  }
  
  comp_tbl <- analyze_class_composition_csv(cluster_df, metab_class_df,
                                            tis, output_dir)
  all_comp_list[[tis]] <- comp_tbl
  all_out[[tis]] <- list(PMI_related_full = feats_full,
                         cluster_df       = cluster_df)
}

# =========================================================
# Section 5b: Five-tissue class composition plot
# =========================================================
plot_combined_class_composition(all_comp_list, output_dir)

# =========================================================
# Section 6: Figure 1d - Five-set Venn diagram (ragg device)
# =========================================================
PMI_related_full_sets <- lapply(tissues,
                                function(tis) all_out[[tis]]$PMI_related_full)
names(PMI_related_full_sets) <- tissues

H  <- PMI_related_full_sets[["Heart"]]
L  <- PMI_related_full_sets[["Liver"]]
S  <- PMI_related_full_sets[["Spleen"]]
Lu <- PMI_related_full_sets[["Lung"]]
M  <- PMI_related_full_sets[["Muscle"]]

venn_file <- file.path(output_dir, "FigureS_Venn_PMIrelated_5tissues.TIFF")

ragg::agg_png(venn_file,
              width  = 2200,
              height = 1800,
              units  = "px",
              res    = 220)
grid::grid.newpage()
VennDiagram::draw.quintuple.venn(
  area1 = length(H),  area2 = length(L),  area3 = length(S),
  area4 = length(Lu), area5 = length(M),
  n12 = length(intersect(H, L)),   n13 = length(intersect(H, S)),
  n14 = length(intersect(H, Lu)),  n15 = length(intersect(H, M)),
  n23 = length(intersect(L, S)),   n24 = length(intersect(L, Lu)),
  n25 = length(intersect(L, M)),   n34 = length(intersect(S, Lu)),
  n35 = length(intersect(S, M)),   n45 = length(intersect(Lu, M)),
  n123  = length(Reduce(intersect, list(H, L, S))),
  n124  = length(Reduce(intersect, list(H, L, Lu))),
  n125  = length(Reduce(intersect, list(H, L, M))),
  n134  = length(Reduce(intersect, list(H, S, Lu))),
  n135  = length(Reduce(intersect, list(H, S, M))),
  n145  = length(Reduce(intersect, list(H, Lu, M))),
  n234  = length(Reduce(intersect, list(L, S, Lu))),
  n235  = length(Reduce(intersect, list(L, S, M))),
  n245  = length(Reduce(intersect, list(L, Lu, M))),
  n345  = length(Reduce(intersect, list(S, Lu, M))),
  n1234 = length(Reduce(intersect, list(H, L, S, Lu))),
  n1235 = length(Reduce(intersect, list(H, L, S, M))),
  n1245 = length(Reduce(intersect, list(H, L, Lu, M))),
  n1345 = length(Reduce(intersect, list(H, S, Lu, M))),
  n2345 = length(Reduce(intersect, list(L, S, Lu, M))),
  n12345 = length(Reduce(intersect, list(H, L, S, Lu, M))),
  category = c("Heart", "Liver", "Spleen", "Lung", "Muscle"),
  fill     = c("#B0D8F6", "#F6D2B2", "#BAECC7", "#FBC9C3", "#BBCBEF"),
  alpha    = rep(0.45, 5),
  cat.col  = c("#5B9BD5", "#E07B39", "#5BAD72", "#D45F5F", "#7B7BBF"),
  cat.cex  = 1.2, cex = 1.8, lwd = 0, lty = "blank"
)
dev.off()

# Export shared and tissue-specific metabolite lists
shared_full_5of5 <- Reduce(intersect, PMI_related_full_sets)
write.csv(data.frame(metabolite = shared_full_5of5),
          file.path(output_dir, "Shared_PMIrelated_5of5tissues.csv"),
          row.names = FALSE)

specific_full <- lapply(tissues, function(tis) {
  others_union <- unique(unlist(
    PMI_related_full_sets[setdiff(tissues, tis)]))
  setdiff(PMI_related_full_sets[[tis]], others_union)
})
names(specific_full) <- tissues

for (tis in tissues) {
  write.csv(
    data.frame(metabolite = specific_full[[tis]]),
    file.path(output_dir,
              paste0("TissueSpecific_PMIrelated_", tis, ".csv")),
    row.names = FALSE
  )
}

cat("\nAnalysis complete. All outputs saved to:", output_dir, "\n")

