############################################################
## Figure 4: Multi-Omics Joint Analysis
## Analyses included:
##   - Procrustes Analysis (metabolome vs microbiome)
##   - Mantel & Partial Mantel Tests
##   - Metabolite-Microbiome Bipartite Correlation Network
## Inputs:
##   - Metabolite feature table (samples x metabolites)
##   - Species-level relative abundance table (samples x species)
##   - Previous metabolomics output directory
##     (containing assoc_spearman_gam_full_*.csv files)
## Output: timestamped directory on Desktop
############################################################

# =========================================================
# Section 0: Output directory setup
# =========================================================
desktop_path <- ifelse(
  Sys.info()["sysname"] == "Windows",
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("HOME"), "Desktop")
)
output_dir <- file.path(
  desktop_path,
  paste0("MultiOmics_Joint_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Section 1: Load required packages
# ragg is used for all figure output to support Arial font
# via the system font engine (bypasses PostScript database)
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(RColorBrewer)
  library(patchwork)
  library(ragg)
})

margin       <- ggplot2::margin
element_text <- ggplot2::element_text

# =========================================================
# Section 2: Analysis parameters
# =========================================================
PMI_MIN_JOINT    <- 2      # Minimum PMI for joint analysis
TOP_N_METAB      <- 30     # Top metabolites per tissue for joint analysis
SPEARMAN_RHO     <- 0.40   # |rho| threshold for network edges
SPEARMAN_FDR     <- 0.05   # FDR threshold for network edges
MANTEL_PERM      <- 999    # Permutations for Mantel tests
CLR_PSEUDO       <- 1e-6   # Pseudocount for CLR transformation
SPECIES_RHO_CUT  <- 0.35   # |rho| threshold for species-PMI association
SPECIES_FDR_CUT  <- 0.05   # FDR threshold for species-PMI association
SPECIES_PREV_CUT <- 0.20   # Prevalence filter: species present in >= 20% samples

# =========================================================
# Section 3: Load input data
# =========================================================
# Locate previous metabolomics output directory
prev_output_dir <- dirname(file.choose())

# Select: (1) metabolite feature table CSV
#   Required columns: SampleID, Organ, PMI, AnimalID + metabolite columns
metab_raw <- read.csv(file.choose(), check.names = FALSE)

# Select: (2) species-level relative abundance table CSV
#   Required columns: SampleID, Organ, PMI, AnimalID + species columns
micro_raw <- read.csv(file.choose(), check.names = FALSE)

# Select: (3) metabolite ID-to-name mapping table CSV
#   Column 1: metabolite_id | Column 2: metabolite_name
metab_name_map <- read.csv(file.choose(), check.names = FALSE) %>%
  rename(metabolite_id   = 1,
         metabolite_name = 2) %>%
  distinct(metabolite_id, .keep_all = TRUE)

# Standardize first column as SampleID
metab_raw <- metab_raw %>% rename(SampleID = 1)
micro_raw <- micro_raw %>% rename(SampleID = 1)

# Validate required columns
stopifnot(all(c("SampleID", "Organ", "PMI", "AnimalID") %in% colnames(metab_raw)))
stopifnot(all(c("SampleID", "Organ", "PMI", "AnimalID") %in% colnames(micro_raw)))

# Separate metadata from feature matrices
meta_cols  <- c("SampleID", "Organ", "PMI", "AnimalID")
meta_metab <- metab_raw %>% select(all_of(meta_cols))
X_metab    <- metab_raw %>% select(-all_of(meta_cols))
meta_micro <- micro_raw %>% select(all_of(meta_cols))
X_micro    <- micro_raw %>% select(-all_of(meta_cols))

cat("Metabolomics: n =", nrow(meta_metab), "samples |",
    ncol(X_metab), "metabolites\n")
cat("Microbiomics: n =", nrow(meta_micro), "samples |",
    ncol(X_micro), "species\n")
cat("Name mapping: n =", nrow(metab_name_map), "entries\n")

# =========================================================
# Section 4: Function definitions
# =========================================================

# ----------------------------------------------------------
# 4.1 CLR transformation (row-wise, with pseudocount)
# ----------------------------------------------------------
clr_transform <- function(X, pseudo = CLR_PSEUDO) {
  X_ps  <- as.matrix(X) + pseudo
  log_X <- log(X_ps)
  sweep(log_X, 1, rowMeans(log_X), "-")
}

# ----------------------------------------------------------
# 4.2 Load top metabolites from previous association results
# ----------------------------------------------------------
get_top_metabolites <- function(assoc_file, top_n = TOP_N_METAB) {
  if (!file.exists(assoc_file)) {
    warning("Association file not found: ", assoc_file)
    return(character(0))
  }
  assoc <- read.csv(assoc_file, check.names = FALSE)
  
  needed <- c("metabolite", "rho", "fdr_spearman", "sig_any", "fdr_min")
  if (!all(needed %in% colnames(assoc))) {
    warning("Missing required columns in association file.")
    return(character(0))
  }
  
  sig <- assoc %>% filter(sig_any == TRUE)
  
  if ("r2_gam" %in% colnames(sig)) {
    sig <- sig %>%
      mutate(r2_gam = ifelse(is.na(r2_gam), 0, r2_gam)) %>%
      arrange(fdr_min, desc(abs(rho) + r2_gam))
  } else {
    sig <- sig %>% arrange(fdr_min, desc(abs(rho)))
  }
  
  top_mets <- sig %>% slice_head(n = top_n) %>%
    pull(metabolite) %>% unique()
  cat("  Top metabolites n =", length(top_mets), "\n")
  top_mets
}

# ----------------------------------------------------------
# 4.3 Species-PMI Spearman association (per tissue)
# ----------------------------------------------------------
get_pmi_related_species <- function(meta_t, X_clr_t,
                                    rho_cut = SPECIES_RHO_CUT,
                                    fdr_cut = SPECIES_FDR_CUT) {
  PMI  <- meta_t$PMI
  sps  <- colnames(X_clr_t)
  
  res <- lapply(sps, function(g) {
    ct <- suppressWarnings(
      cor.test(X_clr_t[, g], PMI, method = "spearman", exact = FALSE))
    data.frame(species = g,
               rho     = unname(ct$estimate),
               p       = ct$p.value)
  }) %>% bind_rows() %>%
    mutate(fdr = p.adjust(p, "BH"),
           sig = fdr < fdr_cut & abs(rho) >= rho_cut)
  res
}

# ----------------------------------------------------------
# 4.4 Procrustes + Mantel + Partial Mantel (per tissue)
#
# Partial Mantel interpretation:
#   pmantel_micro_indep: microbiome independent contribution
#     to PMI variation, controlling for metabolome distance
#   pmantel_metab_indep: metabolome independent contribution
#     to PMI variation, controlling for microbiome distance
# ----------------------------------------------------------
run_procrustes_mantel <- function(meta_t, X_met_t, X_mic_t,
                                  top_mets, top_species,
                                  perm        = MANTEL_PERM,
                                  tissue_name = "") {
  use_mets    <- intersect(top_mets,    colnames(X_met_t))
  use_species <- intersect(top_species, colnames(X_mic_t))
  
  if (length(use_mets) < 3 || length(use_species) < 3 ||
      nrow(meta_t) < 6) {
    cat("  [skip]", tissue_name, "- insufficient samples or features\n")
    return(NULL)
  }
  
  dist_met <- dist(scale(X_met_t[, use_mets,    drop = FALSE]))
  dist_mic <- dist(      X_mic_t[, use_species,  drop = FALSE])
  dist_pmi <- dist(meta_t$PMI)
  
  pco_met <- cmdscale(dist_met, k = 2)
  pco_mic <- cmdscale(dist_mic, k = 2)
  
  proc <- procrustes(pco_met, pco_mic, symmetric = TRUE)
  prot <- protest(pco_met, pco_mic, permutations = perm)
  
  mantel_res <- mantel(dist_met, dist_mic,
                       method = "spearman", permutations = perm)
  
  pmantel_micro_indep <- mantel.partial(
    dist_mic, dist_pmi, dist_met,
    method = "spearman", permutations = perm
  )
  pmantel_metab_indep <- mantel.partial(
    dist_met, dist_pmi, dist_mic,
    method = "spearman", permutations = perm
  )
  
  list(
    proc                = proc,
    prot                = prot,
    mantel              = mantel_res,
    pmantel_micro_indep = pmantel_micro_indep,
    pmantel_metab_indep = pmantel_metab_indep,
    meta_t              = meta_t,
    pco_met             = pco_met,
    pco_mic             = pco_mic,
    tissue              = tissue_name
  )
}

# ----------------------------------------------------------
# 4.5 Plot Procrustes overlay (per tissue panel)
# ----------------------------------------------------------
plot_procrustes <- function(proc_res, tissue_name) {
  if (is.null(proc_res))
    return(ggplot() + theme_void(base_family = "Arial") +
             labs(title = tissue_name))
  
  proc   <- proc_res$proc
  meta_t <- proc_res$meta_t
  
  X <- as.data.frame(proc$X[,    1:2, drop = FALSE])
  Y <- as.data.frame(proc$Yrot[, 1:2, drop = FALSE])
  colnames(X) <- c("Dim1", "Dim2")
  colnames(Y) <- c("Dim1", "Dim2")
  
  df_plot <- bind_rows(
    X %>% mutate(type = "Metabolome", PMI = meta_t$PMI),
    Y %>% mutate(type = "Microbiome", PMI = meta_t$PMI)
  )
  df_seg <- data.frame(
    x1 = X$Dim1, y1 = X$Dim2,
    x2 = Y$Dim1, y2 = Y$Dim2,
    PMI = meta_t$PMI
  )
  
  m2_val <- round(proc_res$prot$ss,           3)
  r_val  <- round(sqrt(1 - proc_res$prot$ss), 3)
  p_val  <- round(proc_res$prot$signif,        3)
  
  ggplot() +
    geom_segment(data = df_seg,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 color = "grey70", linewidth = 0.4) +
    geom_point(data = df_plot,
               aes(x = Dim1, y = Dim2,
                   color = factor(PMI), shape = type),
               size = 2.8, alpha = 0.9) +
    scale_color_brewer(palette = "Set1", name = "PMI (d)") +
    scale_shape_manual(
      values = c(Metabolome = 16, Microbiome = 17),
      name   = "Omics"
    ) +
    labs(
      title    = tissue_name,
      subtitle = bquote(
        italic(M)^2 == .(m2_val) ~ "| r =" ~ .(r_val) ~ "| p =" ~ .(p_val)
      ),
      x = "PC1", y = "PC2"
    ) +
    theme_bw(base_size = 10, base_family = "Arial") +
    theme(
      text             = element_text(family = "Arial"),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
}

# ----------------------------------------------------------
# 4.6 Combined Procrustes plot (all tissues in one figure)
# ----------------------------------------------------------
plot_procrustes_combined <- function(proc_mantel_list, tissues_joint) {
  df_plot_all <- list()
  df_seg_all  <- list()
  
  for (tis in tissues_joint) {
    pm <- proc_mantel_list[[tis]]
    if (is.null(pm)) next
    
    proc   <- pm$proc
    meta_t <- pm$meta_t
    
    X <- as.data.frame(proc$X[,    1:2, drop = FALSE])
    Y <- as.data.frame(proc$Yrot[, 1:2, drop = FALSE])
    colnames(X) <- c("Dim1", "Dim2")
    colnames(Y) <- c("Dim1", "Dim2")
    
    df_plot_all[[tis]] <- bind_rows(
      X %>% mutate(type = "Metabolome", tissue = tis, PMI = meta_t$PMI),
      Y %>% mutate(type = "Microbiome", tissue = tis, PMI = meta_t$PMI)
    )
    df_seg_all[[tis]] <- data.frame(
      x1 = X$Dim1, y1 = X$Dim2,
      x2 = Y$Dim1, y2 = Y$Dim2,
      tissue = tis
    )
  }
  
  df_plot <- bind_rows(df_plot_all)
  df_seg  <- bind_rows(df_seg_all)
  
  stats_df <- bind_rows(lapply(tissues_joint, function(tis) {
    pm <- proc_mantel_list[[tis]]
    if (is.null(pm)) return(NULL)
    data.frame(
      tissue = tis,
      r      = round(sqrt(1 - pm$prot$ss), 3),
      p      = pm$prot$signif
    )
  })) %>% mutate(label = paste0(tissue, " (r=", r, ", p=", p, ")"))
  
  tissue_colors <- setNames(
    brewer.pal(max(3, length(tissues_joint)), "Set2")[seq_along(tissues_joint)],
    tissues_joint
  )
  
  ggplot() +
    geom_segment(
      data = df_seg,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = tissue),
      linewidth = 0.35, alpha = 0.4
    ) +
    geom_point(
      data = df_plot,
      aes(x = Dim1, y = Dim2, color = tissue, shape = type),
      size = 2.5, alpha = 0.85
    ) +
    scale_color_manual(
      values = tissue_colors, name = "Tissue",
      labels = stats_df$label
    ) +
    scale_shape_manual(
      values = c(Metabolome = 16, Microbiome = 17),
      name   = "Omics layer"
    ) +
    labs(
      title    = "Procrustes Analysis: Metabolome vs Microbiome",
      subtitle = paste0("All tissues combined | PMI ", PMI_MIN_JOINT, "\u20137 d"),
      x = "PC1", y = "PC2"
    ) +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(
      text             = element_text(family = "Arial"),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      legend.text      = element_text(size = 8, family = "Arial"),
      panel.grid.minor = element_blank()
    )
}

# ----------------------------------------------------------
# 4.7 Summarise Procrustes + Mantel results into a table
# ----------------------------------------------------------
summarise_procrustes_mantel <- function(pm_list) {
  bind_rows(lapply(pm_list, function(x) {
    if (is.null(x)) return(NULL)
    prot <- x$prot
    man  <- x$mantel
    pm_m <- x$pmantel_micro_indep
    pm_e <- x$pmantel_metab_indep
    data.frame(
      tissue                = x$tissue,
      procrustes_r          = round(sqrt(1 - prot$ss), 4),
      procrustes_M2         = round(prot$ss, 4),
      procrustes_p          = prot$signif,
      mantel_r              = round(man$statistic,  4),
      mantel_p              = man$signif,
      pMantel_micro_indep_r = round(pm_m$statistic, 4),
      pMantel_micro_indep_p = pm_m$signif,
      pMantel_metab_indep_r = round(pm_e$statistic, 4),
      pMantel_metab_indep_p = pm_e$signif
    )
  }))
}

# ----------------------------------------------------------
# 4.8 Build bipartite metabolite-microbiome network (per tissue)
# ----------------------------------------------------------
run_bipartite_network <- function(meta_t, X_met_t, X_mic_t,
                                  top_mets, pmi_species_df,
                                  metab_name_map = NULL,
                                  top_n_met   = 10,
                                  top_n_mic   = 5,
                                  rho_cut     = SPEARMAN_RHO,
                                  fdr_cut     = SPEARMAN_FDR,
                                  tissue_name = "") {
  
  # Select top metabolites (already ranked by FDR + |rho| + GAM R²)
  use_mets <- head(intersect(top_mets, colnames(X_met_t)), top_n_met)
  
  # Select top species by PMI Spearman |rho|
  use_species <- pmi_species_df %>%
    filter(sig == TRUE) %>%
    arrange(desc(abs(rho))) %>%
    slice_head(n = top_n_mic) %>%
    pull(species) %>%
    intersect(colnames(X_mic_t))
  
  if (length(use_mets) < 2 || length(use_species) < 1) {
    cat("  [skip network]", tissue_name, "- insufficient features\n")
    return(NULL)
  }
  
  # Compute all pairwise Spearman correlations
  cor_results <- expand.grid(
    metabolite = use_mets,
    species    = use_species,
    stringsAsFactors = FALSE
  )
  
  cor_results$rho <- mapply(function(m, g)
    suppressWarnings(cor(X_met_t[, m], X_mic_t[, g],
                         method = "spearman", use = "complete.obs")),
    cor_results$metabolite, cor_results$species)
  
  cor_results$p <- mapply(function(m, g)
    suppressWarnings(cor.test(X_met_t[, m], X_mic_t[, g],
                              method = "spearman",
                              exact = FALSE)$p.value),
    cor_results$metabolite, cor_results$species)
  
  cor_results$fdr <- p.adjust(cor_results$p, "BH")
  
  cor_results <- cor_results %>%
    mutate(
      sig       = abs(rho) >= rho_cut & fdr < fdr_cut,
      direction = ifelse(rho > 0, "positive", "negative")
    )
  
  # Metabolite ID to display name mapping
  id_to_name <- function(ids) {
    if (is.null(metab_name_map)) return(ids)
    idx <- match(ids, metab_name_map$metabolite_id)
    ifelse(!is.na(idx), metab_name_map$metabolite_name[idx], ids)
  }
  
  cor_results <- cor_results %>%
    mutate(metabolite_name = id_to_name(metabolite))
  
  met_names <- id_to_name(use_mets)
  mic_names <- use_species
  
  edges_sig <- cor_results %>% filter(sig == TRUE)
  
  cat("  [", tissue_name, "] Significant pairs:",
      nrow(edges_sig), "/", nrow(cor_results), "\n")
  
  list(
    all_cors      = cor_results,
    edges         = edges_sig,
    use_mets      = use_mets,
    use_mets_name = met_names,
    use_species   = use_species,
    tissue        = tissue_name
  )
}

# ----------------------------------------------------------
# 4.9 Plot bipartite network (dual-column layout)
# ----------------------------------------------------------
plot_bipartite_network <- function(net_res, tissue_name) {
  
  if (is.null(net_res)) {
    return(
      ggplot() + theme_void(base_family = "Arial") +
        labs(title = paste0(tissue_name, ": No significant pairs")) +
        theme(plot.background = element_rect(fill = "white",
                                             color = "grey88",
                                             linewidth = 0.5),
              text = element_text(family = "Arial"))
    )
  }
  
  met_ord <- net_res$use_mets_name
  mic_ord <- net_res$use_species
  n_met   <- length(met_ord)
  n_mic   <- length(mic_ord)
  
  # Node positions: metabolites on left (x=0), species on right (x=1)
  node_pos <- bind_rows(
    data.frame(name = met_ord, x = 0,
               y = seq(1, 0, length.out = n_met),
               type = "Metabolite", stringsAsFactors = FALSE),
    data.frame(name = mic_ord, x = 1,
               y = seq(1, 0, length.out = n_mic),
               type = "Microbe", stringsAsFactors = FALSE)
  )
  
  # Scale microbe node size by mean |rho| with all metabolites
  mic_rho <- net_res$all_cors %>%
    group_by(species) %>%
    summarise(mean_abs_rho = mean(abs(rho)), .groups = "drop")
  
  node_pos <- node_pos %>%
    left_join(mic_rho %>% rename(name = species), by = "name") %>%
    mutate(
      mean_abs_rho = ifelse(is.na(mean_abs_rho), 0, mean_abs_rho),
      node_size = case_when(
        type == "Microbe"    ~ scales::rescale(mean_abs_rho, to = c(5, 13)),
        type == "Metabolite" ~ 3.2
      ),
      display_name = ifelse(type == "Microbe",
                            gsub("^s__", "", name), name)
    )
  
  # Edge coordinates (significant edges only)
  edges_sig <- net_res$edges
  edge_df   <- data.frame()
  
  if (!is.null(edges_sig) && nrow(edges_sig) > 0) {
    edge_df <- edges_sig %>%
      left_join(node_pos %>% filter(type == "Metabolite") %>%
                  select(name, y) %>% rename(y0 = y),
                by = c("metabolite_name" = "name")) %>%
      left_join(node_pos %>% filter(type == "Microbe") %>%
                  select(name, y) %>% rename(y1 = y),
                by = c("species" = "name")) %>%
      mutate(x0 = 0, x1 = 1, abs_r = abs(rho)) %>%
      filter(!is.na(y0), !is.na(y1))
  }
  
  n_sig <- nrow(edge_df)
  
  p <- ggplot() +
    # Negative correlation edges (blue)
    {if (nrow(edge_df) > 0 && any(edge_df$rho < 0))
      geom_segment(
        data = edge_df %>% filter(rho < 0),
        aes(x = x0, y = y0, xend = x1, yend = y1,
            alpha = abs_r, linewidth = abs_r),
        color = "#3A7FC1", lineend = "round", show.legend = FALSE
      )} +
    # Positive correlation edges (orange-red)
    {if (nrow(edge_df) > 0 && any(edge_df$rho >= 0))
      geom_segment(
        data = edge_df %>% filter(rho >= 0),
        aes(x = x0, y = y0, xend = x1, yend = y1,
            alpha = abs_r, linewidth = abs_r),
        color = "#E05A3A", lineend = "round", show.legend = FALSE
      )} +
    scale_alpha_continuous(range = c(0.22, 0.68), guide = "none") +
    scale_linewidth_continuous(range = c(0.5, 2.2), guide = "none") +
    # Microbe nodes (size proportional to mean |rho|)
    geom_point(
      data = node_pos %>% filter(type == "Microbe"),
      aes(x = x, y = y, size = node_size),
      shape = 21, fill = "#5BA8D4", color = "#2C6EA0",
      stroke = 0.8, alpha = 0.90, show.legend = FALSE
    ) +
    # Metabolite nodes
    geom_point(
      data = node_pos %>% filter(type == "Metabolite"),
      aes(x = x, y = y, size = node_size),
      shape = 21, fill = "#D94F3D", color = "#A02820",
      stroke = 0.6, alpha = 0.92, show.legend = FALSE
    ) +
    scale_size_identity() +
    # Species labels (right, italic)
    geom_text(
      data = node_pos %>% filter(type == "Microbe"),
      aes(x = x + 0.04, y = y, label = display_name),
      hjust = 0, size = 3.6, color = "grey10",
      fontface = "italic", family = "Arial"
    ) +
    # Metabolite labels (left)
    geom_text(
      data = node_pos %>% filter(type == "Metabolite"),
      aes(x = x - 0.04, y = y, label = display_name),
      hjust = 1, size = 3.6, color = "grey10", family = "Arial"
    ) +
    # Column header labels
    annotate("text", x = 0,  y = -0.13, label = "Metabolites",
             hjust = 0.5, size = 4.2, fontface = "bold",
             color = "#A02820", family = "Arial") +
    annotate("text", x = 1, y = -0.13, label = "Microbiome",
             hjust = 0.5, size = 4.2, fontface = "bold",
             color = "#2C6EA0", family = "Arial") +
    coord_cartesian(xlim = c(-1.6, 2.4), ylim = c(-0.20, 1.08),
                    clip = "off") +
    labs(
      title    = tissue_name,
      subtitle = paste0("Top ", n_met, " metabolites \u00d7 Top ", n_mic,
                        " species  |  Nodes: ", n_met + n_mic,
                        "  Edges: ", n_sig)
    ) +
    theme_void(base_size = 12, base_family = "Arial") +
    theme(
      text            = element_text(family = "Arial"),
      plot.title      = element_text(size = 12, face = "bold",
                                     color = "grey10", family = "Arial",
                                     margin = margin(b = 3)),
      plot.subtitle   = element_text(size = 8, color = "grey45",
                                     family = "Arial",
                                     margin = margin(b = 6)),
      plot.background = element_rect(fill = "white", color = "grey82",
                                     linewidth = 0.6),
      plot.margin     = ggplot2::margin(12, 6, 24, 6)
    )
  
  p
}

# =========================================================
# Section 5: Main analysis pipeline
# =========================================================

# ----------------------------------------------------------
# 5.1 Build joint analysis sample set (PMI >= PMI_MIN_JOINT)
# ----------------------------------------------------------
common_ids_joint <- intersect(
  meta_metab %>% filter(PMI >= PMI_MIN_JOINT) %>% pull(SampleID),
  meta_micro %>% filter(PMI >= PMI_MIN_JOINT) %>% pull(SampleID)
)
cat("\nJoint analysis sample count (PMI >=", PMI_MIN_JOINT, "d):",
    length(common_ids_joint), "\n")

meta_joint <- meta_metab %>%
  filter(SampleID %in% common_ids_joint) %>%
  arrange(SampleID)

X_met_joint <- X_metab[meta_metab$SampleID %in% common_ids_joint, ]
X_met_joint <- X_met_joint[
  match(meta_joint$SampleID,
        meta_metab$SampleID[meta_metab$SampleID %in% common_ids_joint]), ]

meta_mic_joint  <- meta_micro %>% filter(SampleID %in% common_ids_joint)
X_mic_joint_raw <- X_micro[meta_micro$SampleID %in% common_ids_joint, ]
X_mic_joint_raw <- X_mic_joint_raw[
  match(meta_joint$SampleID, meta_mic_joint$SampleID), ]

# CLR transformation
X_mic_clr <- as.data.frame(clr_transform(X_mic_joint_raw))

# Prevalence filter: retain species present in >= 20% of samples
prev_count   <- apply(X_mic_joint_raw > 0, 2, sum)
keep_species <- names(prev_count)[
  prev_count >= SPECIES_PREV_CUT * nrow(X_mic_joint_raw)]
X_mic_clr <- X_mic_clr[, keep_species, drop = FALSE]
cat("Species retained after prevalence filter:", ncol(X_mic_clr), "\n")

write.csv(meta_joint,
          file.path(output_dir, "joint_sample_info.csv"),
          row.names = FALSE)

# Identify tissues with >= 6 samples
tissue_n      <- meta_joint %>% count(Organ) %>% filter(n >= 6)
tissues_joint <- sort(tissue_n$Organ)
cat("Tissues included:", paste(tissues_joint, collapse = ", "), "\n")

# ----------------------------------------------------------
# 5.2 Load top metabolites from previous association files
# ----------------------------------------------------------
cat("\n--- Loading top metabolites per tissue ---\n")
top_mets_per_tissue <- list()

for (tis in tissues_joint) {
  cat("[", tis, "]\n")
  assoc_f <- file.path(prev_output_dir,
                       paste0("Assoc_SpearmanGAM_", tis, ".csv"))
  top_mets_per_tissue[[tis]] <- get_top_metabolites(assoc_f,
                                                    top_n = TOP_N_METAB)
}

top_met_df <- bind_rows(lapply(names(top_mets_per_tissue), function(tis) {
  mets <- top_mets_per_tissue[[tis]]
  if (length(mets) == 0) return(NULL)
  data.frame(tissue = tis, metabolite = mets, rank = seq_along(mets))
}))
write.csv(top_met_df,
          file.path(output_dir, "top_metabolites_per_tissue.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 5.3 Species-PMI Spearman association (per tissue)
# ----------------------------------------------------------
cat("\n--- Species-PMI Spearman association ---\n")
species_pmi_list <- list()

for (tis in tissues_joint) {
  idx     <- meta_joint$Organ == tis
  meta_t  <- meta_joint[idx, , drop = FALSE]
  X_mic_t <- X_mic_clr[idx, , drop = FALSE]
  
  sp_res <- get_pmi_related_species(meta_t, X_mic_t)
  species_pmi_list[[tis]] <- sp_res
  
  write.csv(sp_res,
            file.path(output_dir,
                      paste0("Species_PMI_spearman_", tis, ".csv")),
            row.names = FALSE)
  cat("[", tis, "] PMI-related species n =",
      sum(sp_res$sig, na.rm = TRUE), "\n")
}

# ----------------------------------------------------------
# 5.4 Procrustes + Mantel + Partial Mantel (per tissue)
# ----------------------------------------------------------
cat("\n--- Procrustes / Mantel analysis ---\n")
proc_mantel_list <- list()

for (tis in tissues_joint) {
  cat("[", tis, "]\n")
  idx     <- meta_joint$Organ == tis
  meta_t  <- meta_joint[idx, , drop = FALSE]
  X_met_t <- X_met_joint[idx, , drop = FALSE]
  X_mic_t <- X_mic_clr[idx, , drop = FALSE]
  
  proc_mantel_list[[tis]] <- run_procrustes_mantel(
    meta_t, X_met_t, X_mic_t,
    top_mets_per_tissue[[tis]],
    colnames(X_mic_t),
    perm        = MANTEL_PERM,
    tissue_name = tis
  )
}

pm_summary <- summarise_procrustes_mantel(proc_mantel_list)
write.csv(pm_summary,
          file.path(output_dir, "Procrustes_Mantel_summary.csv"),
          row.names = FALSE)
cat("\n[Procrustes + Mantel summary]\n")
print(pm_summary)

# ----------------------------------------------------------
# 5.5 Figure 4a: Procrustes plots
# ----------------------------------------------------------
cat("\n--- Plotting Procrustes panels ---\n")
proc_plots <- lapply(tissues_joint, function(tis)
  plot_procrustes(proc_mantel_list[[tis]], tis))
names(proc_plots) <- tissues_joint

n_tis      <- length(tissues_joint)
ncols_proc <- min(n_tis, 3)

fig_proc <- wrap_plots(proc_plots, ncol = ncols_proc) +
  plot_annotation(
    title    = "Procrustes Analysis: Metabolome vs Microbiome",
    subtitle = paste0("Per tissue | PMI ", PMI_MIN_JOINT, "\u20137 d")
  ) &
  theme(
    text          = element_text(family = "Arial"),
    plot.title    = element_text(size = 13, face = "bold",
                                 family = "Arial"),
    plot.subtitle = element_text(size = 10, family = "Arial")
  )

ggsave(
  file.path(output_dir, "Fig4A_Procrustes_per_tissue.png"),
  fig_proc,
  width  = 5 * ncols_proc,
  height = 5 * ceiling(n_tis / ncols_proc),
  units  = "in", dpi = 300,
  device = ragg::agg_png
)
cat("[Saved] Fig4A_Procrustes_per_tissue.png\n")

# Combined all-tissue Procrustes
fig_proc_combined <- plot_procrustes_combined(proc_mantel_list,
                                              tissues_joint)
ggsave(
  file.path(output_dir, "Fig4A_Procrustes_combined.png"),
  fig_proc_combined,
  width = 9, height = 6,
  units = "in", dpi = 300,
  device = ragg::agg_png
)
cat("[Saved] Fig4A_Procrustes_combined.png\n")

# ----------------------------------------------------------
# 5.6 Figure 4b: Mantel & Partial Mantel bubble plot
# ----------------------------------------------------------
cat("\n--- Plotting Mantel bubble chart ---\n")

if (!is.null(pm_summary) && nrow(pm_summary) > 0) {
  
  pm_long <- pm_summary %>%
    select(tissue,
           mantel_r,              mantel_p,
           pMantel_micro_indep_r, pMantel_micro_indep_p,
           pMantel_metab_indep_r, pMantel_metab_indep_p) %>%
    pivot_longer(
      cols          = -tissue,
      names_to      = c("test", ".value"),
      names_pattern = "(.+)_(r|p)$"
    ) %>%
    mutate(
      test_label = factor(
        dplyr::recode(
          test,
          "mantel"              = "Mantel\n(Overall)",
          "pMantel_micro_indep" = "pMantel\n(Microbiome)",
          "pMantel_metab_indep" = "pMantel\n(Metabolome)"
        ),
        levels = c("Mantel\n(Overall)",
                   "pMantel\n(Microbiome)",
                   "pMantel\n(Metabolome)")
      ),
      color_group = dplyr::recode(
        test,
        "mantel"              = "Overall",
        "pMantel_micro_indep" = "Microbiome",
        "pMantel_metab_indep" = "Metabolome"
      ),
      sig_label = case_when(
        as.numeric(p) < 0.001 ~ "***",
        as.numeric(p) < 0.01  ~ "**",
        as.numeric(p) < 0.05  ~ "*",
        TRUE                   ~ ""
      ),
      r = as.numeric(r),
      p = as.numeric(p)
    )
  
  color_vals <- c(
    Overall    = "#7B6FA0",
    Microbiome = "#00A896",
    Metabolome = "#E8541A"
  )
  
  fig_mantel <- ggplot(pm_long, aes(x = tissue, y = test_label)) +
    geom_tile(aes(fill = color_group),
              alpha = 0.06, color = NA, width = Inf, height = 1) +
    geom_point(
      aes(size = abs(r), color = color_group,
          alpha = ifelse(p < 0.05, 1, 0.38)),
      shape = 16
    ) +
    geom_text(
      aes(label = sig_label, color = color_group),
      size = 4.5, vjust = -1.1, fontface = "bold",
      family = "Arial", show.legend = FALSE
    ) +
    geom_text(
      data = pm_long %>% filter(p < 0.05),
      aes(label = sprintf("%.2f", r)),
      size = 2.8, vjust = 2.2, color = "grey20",
      family = "Arial", show.legend = FALSE
    ) +
    scale_size_continuous(
      range  = c(2.5, 11),
      breaks = c(0.1, 0.3, 0.5, 0.7),
      name   = "|r|"
    ) +
    scale_color_manual(values = color_vals, guide = "none") +
    scale_fill_manual(values  = color_vals, guide = "none") +
    scale_alpha_identity() +
    geom_hline(yintercept = c(1.5, 2.5),
               color = "grey80", linewidth = 0.4,
               linetype = "dashed") +
    labs(
      title    = "Mantel & Partial Mantel Tests",
      subtitle = paste0(
        "* p<0.05  ** p<0.01  *** p<0.001  (",
        MANTEL_PERM, " permutations)\n",
        "\u25CF Microbiome (teal) vs Metabolome (orange)",
        " independent contributions"
      ),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(
      text             = element_text(family = "Arial"),
      axis.text.x      = element_text(angle = 30, hjust = 1,
                                      size = 10, family = "Arial"),
      axis.text.y      = element_text(size = 10, family = "Arial"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "grey70"),
      plot.title       = element_text(face = "bold", size = 12,
                                      family = "Arial"),
      plot.subtitle    = element_text(size = 9, color = "grey40",
                                      lineheight = 1.3,
                                      family = "Arial"),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      legend.text      = element_text(family = "Arial"),
      plot.margin      = ggplot2::margin(8, 10, 8, 8)
    )
  
  ggsave(
    file.path(output_dir, "Fig4C_Mantel_bubble.png"),
    fig_mantel,
    width = 7, height = 3.8,
    units = "in", dpi = 300,
    device = ragg::agg_png
  )
  cat("[Saved] Fig4C_Mantel_bubble.png\n")
}

# ----------------------------------------------------------
# 5.7 Bipartite metabolite-microbiome network (per tissue)
# ----------------------------------------------------------
cat("\n--- Building bipartite networks ---\n")
network_list <- list()

for (tis in tissues_joint) {
  cat("[", tis, "]\n")
  idx     <- meta_joint$Organ == tis
  meta_t  <- meta_joint[idx, , drop = FALSE]
  X_met_t <- X_met_joint[idx, , drop = FALSE]
  X_mic_t <- X_mic_clr[idx, , drop = FALSE]
  
  network_list[[tis]] <- run_bipartite_network(
    meta_t, X_met_t, X_mic_t,
    top_mets_per_tissue[[tis]],
    species_pmi_list[[tis]],
    metab_name_map = metab_name_map,
    top_n_met      = 10,
    top_n_mic      = 5,
    rho_cut        = SPEARMAN_RHO,
    fdr_cut        = SPEARMAN_FDR,
    tissue_name    = tis
  )
  
  if (!is.null(network_list[[tis]]) &&
      nrow(network_list[[tis]]$edges) > 0) {
    write.csv(
      network_list[[tis]]$edges,
      file.path(output_dir, paste0("network_edges_", tis, ".csv")),
      row.names = FALSE
    )
  }
}

# ----------------------------------------------------------
# 5.8 Figure 4c: Bipartite network multi-panel plot
# ----------------------------------------------------------
cat("\n--- Plotting bipartite network panels ---\n")
net_plots <- lapply(tissues_joint, function(tis)
  plot_bipartite_network(network_list[[tis]], tis))
names(net_plots) <- tissues_joint

if (n_tis == 5) {
  fig_net <-
    (net_plots[[1]] | net_plots[[2]] | net_plots[[3]]) /
    (patchwork::plot_spacer() | net_plots[[4]] |
       net_plots[[5]] | patchwork::plot_spacer()) +
    plot_layout(heights = c(1, 1))
} else {
  fig_net <- wrap_plots(net_plots, ncol = min(n_tis, 3))
}

fig_net <- fig_net +
  plot_annotation(
    title    = "Metabolite\u2013Microbiome Bipartite Correlation Network",
    subtitle = paste0(
      "|rho| \u2265 ", SPEARMAN_RHO, ", FDR < ", SPEARMAN_FDR,
      "  \u2502  \u25AC red = positive  \u25AC blue = negative"
    ),
    theme = theme(
      text          = element_text(family = "Arial"),
      plot.title    = element_text(size = 15, face = "bold",
                                   family = "Arial",
                                   margin = ggplot2::margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "grey40",
                                   family = "Arial")
    )
  )

ggsave(
  file.path(output_dir, "Fig4B_Bipartite_Network_all_tissues.png"),
  fig_net,
  width = 22, height = 15,
  units = "in", dpi = 300,
  device = ragg::agg_png
)
cat("[Saved] Fig4B_Bipartite_Network_all_tissues.png\n")

# ----------------------------------------------------------
# 5.9 Network summary statistics
# ----------------------------------------------------------
net_stats <- bind_rows(lapply(tissues_joint, function(tis) {
  net <- network_list[[tis]]
  if (is.null(net))
    return(data.frame(tissue = tis, n_met_nodes = 0,
                      n_mic_nodes = 0, n_edges = 0,
                      n_pos = 0, n_neg = 0))
  edges <- net$edges
  data.frame(
    tissue      = tis,
    n_met_nodes = length(unique(edges$metabolite)),
    n_mic_nodes = length(unique(edges$species)),
    n_edges     = nrow(edges),
    n_pos       = sum(edges$direction == "positive"),
    n_neg       = sum(edges$direction == "negative")
  )
}))
write.csv(net_stats,
          file.path(output_dir, "network_statistics_all_tissues.csv"),
          row.names = FALSE)
cat("\n[Network statistics]\n")
print(net_stats)

# ----------------------------------------------------------
# 5.10 Cross-tissue shared metabolite-species pairs (>= 2 tissues)
# ----------------------------------------------------------
all_edges <- bind_rows(lapply(tissues_joint, function(tis) {
  if (is.null(network_list[[tis]])) return(NULL)
  network_list[[tis]]$edges %>% mutate(tissue = tis)
}))

if (nrow(all_edges) > 0) {
  shared_pairs <- all_edges %>%
    group_by(metabolite, species) %>%
    summarise(
      n_tissues = n_distinct(tissue),
      tissues   = paste(tissue, collapse = "/"),
      rho_mean  = round(mean(rho), 3),
      direction = ifelse(mean(rho) > 0, "positive", "negative"),
      .groups   = "drop"
    ) %>%
    filter(n_tissues >= 2) %>%
    arrange(desc(n_tissues), desc(abs(rho_mean)))
  
  write.csv(shared_pairs,
            file.path(output_dir,
                      "shared_pairs_ge2tissues.csv"),
            row.names = FALSE)
  cat("\nCross-tissue shared pairs (>= 2 tissues):",
      nrow(shared_pairs), "\n")
  if (nrow(shared_pairs) > 0) print(head(shared_pairs, 10))
}

# ----------------------------------------------------------
# 5.11 PMI-related species cross-tissue overlap
# ----------------------------------------------------------
sig_species_sets <- lapply(tissues_joint, function(tis)
  species_pmi_list[[tis]] %>% filter(sig == TRUE) %>% pull(species))
names(sig_species_sets) <- tissues_joint

species_overlap <- stack(sig_species_sets) %>%
  rename(species = values, tissue = ind) %>%
  distinct() %>%
  count(species, name = "n_tissues") %>%
  arrange(desc(n_tissues))

write.csv(species_overlap,
          file.path(output_dir,
                    "PMI_related_species_tissue_overlap.csv"),
          row.names = FALSE)
cat("\nPMI-related species cross-tissue distribution (Top 10):\n")
print(head(species_overlap, 10))

# =========================================================
# Analysis complete
# =========================================================
cat("\n============================\n")
cat("Multi-omics joint analysis complete.\n")
cat("Output directory:\n", output_dir, "\n")
cat("============================\n")

cat("\nOutput file list:\n")
out_files <- list.files(output_dir, full.names = FALSE)
cat(paste(" -", out_files, collapse = "\n"), "\n")

