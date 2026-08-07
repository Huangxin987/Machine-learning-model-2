setwd("D:/PMI_Project/Figure S1")
output_dir <- file.path(getwd(), "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)


rlsc_span <- 0.75
min_qc_per_batch <- 5

tech_rsd_cutoff <- 30
expected_tech_pairs <- 5


group_missing_cutoff <- 0.50


pca_scale <- TRUE


feat <- read.csv(
  "Feature table.csv",
  row.names = 1,
  check.names = FALSE)

info <- read.csv(
  "Sample_info.csv",
  check.names = FALSE,
  na.strings = c("", "NA")
)

run_order_info <- read.csv(
  "Sample order.csv",
  check.names = FALSE,
  na.strings = c("", "NA")
)


required_feat_meta <- c(
  "mz",
  "RT",
  "Annotation",
  "Is_Endogenous",
  "Mass Error (ppm)"
)

required_info_cols <- c(
  "Sample_ID",
  "Sample_Type",
  "Tissue",
  "PMI_days",
  "Tech_Pair",
  "Tech_Replicate"
)

required_order_cols <- c(
  "Sample_ID",
  "Run_Order"
)


colnames(feat) <- trimws(colnames(feat))

info <- info %>%
  mutate(
    Sample_ID = trimws(as.character(Sample_ID)),
    Sample_Type = case_when(
      tolower(trimws(as.character(Sample_Type))) == "qc" ~ "QC",
      tolower(trimws(as.character(Sample_Type))) %in%
        c("biological", "bio", "sample") ~ "Biological",
      TRUE ~ trimws(as.character(Sample_Type))
    ),
    Tissue = trimws(as.character(Tissue)),
    PMI_days = as.character(PMI_days),
    Tech_Pair = na_if(
      trimws(as.character(Tech_Pair)),
      ""
    ),
    Tech_Replicate = suppressWarnings(
      as.numeric(as.character(Tech_Replicate))
    )
  )

input_sample_order <- colnames(feat)[
  colnames(feat) %in% info$Sample_ID
]

run_order_info <- run_order_info %>%
  transmute(
    Sample_ID = trimws(as.character(Sample_ID)),
    Run_Order = as.numeric(as.character(Run_Order)),
    Batch = ifelse(
      is.na(Batch) | trimws(as.character(Batch)) == "",
      "Batch1",
      trimws(as.character(Batch))
    )
  )

sample_info <- info %>%
  select(
    -any_of(c("Run_Order", "Batch"))
  ) %>%
  left_join(
    run_order_info,
    by = "Sample_ID"
  ) %>%
  filter(
    Sample_ID %in% colnames(feat)
  ) %>%
  arrange(
    Batch,
    Run_Order
  )



sample_cols <- sample_info$Sample_ID


qc_cols <- sample_info %>%
  filter(Sample_Type == "QC") %>%
  pull(Sample_ID)

bio_cols <- sample_info %>%
  filter(Sample_Type == "Biological") %>%
  pull(Sample_ID)
input_bio_order <- input_sample_order[
  input_sample_order %in% bio_cols
]



tech_info <- sample_info %>%
  filter(
    Sample_Type == "Biological",
    !is.na(Tech_Pair)
  ) %>%
  arrange(
    Tech_Pair,
    Tech_Replicate
  )

tech_pair_summary <- tech_info %>%
  group_by(Tech_Pair) %>%
  summarise(
    N = n(),
    N_Tissue = n_distinct(Tissue),
    N_PMI = n_distinct(PMI_days),
    Samples = paste(Sample_ID, collapse = "; "),
    .groups = "drop"
  )


feat[, sample_cols] <- lapply(
  feat[, sample_cols, drop = FALSE],
  function(x) {
    x <- suppressWarnings(
      as.numeric(as.character(x))
    )
    x[
      x == 0 |
        !is.finite(x)
    ] <- NA_real_
    x
  }
)

feat_meta <- feat[
  ,
  required_feat_meta,
  drop = FALSE
]

feat_mat <- as.matrix(
  feat[
    ,
    sample_cols,
    drop = FALSE
  ]
)

storage.mode(feat_mat) <- "numeric"

endo_numeric <- suppressWarnings(
  as.numeric(
    as.character(feat_meta$Is_Endogenous)
  )
)

endo_idx <- !is.na(endo_numeric) &
  endo_numeric == 1

feat_mat_endo_raw <- feat_mat[
  endo_idx,
  ,
  drop = FALSE
]



bio_info <- sample_info %>%
  filter(
    Sample_Type == "Biological"
  ) %>%
  mutate(
    Group = paste(
      Tissue,
      PMI_days,
      sep = "_"
    ),
    PMI_days_num = suppressWarnings(
      as.numeric(PMI_days)
    )
  )

groups <- unique(bio_info$Group)



calc_qc_rsd <- function(x) {

  x <- as.numeric(
    x[
      !is.na(x) &
        is.finite(x) &
        x > 0
    ]
  )

  if (length(x) < 3) {
    return(NA_real_)
  }

  sd(x) / mean(x) * 100
}

calc_tech_rsd <- function(x) {

  x <- as.numeric(
    x[
      !is.na(x) &
        is.finite(x) &
        x > 0
    ]
  )

  if (length(x) != 2) {
    return(NA_real_)
  }

  sd(x) / mean(x) * 100
}

get_fill_value <- function(mat) {

  positive_values <- mat[
    !is.na(mat) &
      is.finite(mat) &
      mat > 0
  ]

  min(positive_values) / 2
}

qc_rlsc_correct_batchwise <- function(
    mat,
    sample_info,
    span = 0.75,
    min_qc = 5
) {

  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  corrected <- mat
  diagnostic_list <- list()
  diagnostic_index <- 1

  for (batch_name in unique(sample_info$Batch)) {

    batch_info <- sample_info %>%
      filter(
        Batch == batch_name
      ) %>%
      arrange(
        Run_Order
      )

    batch_ids <- batch_info$Sample_ID

    batch_qc_ids <- batch_info %>%
      filter(
        Sample_Type == "QC"
      ) %>%
      pull(
        Sample_ID
      )

    x_all <- batch_info$Run_Order

    x_qc <- batch_info %>%
      filter(
        Sample_Type == "QC"
      ) %>%
      pull(
        Run_Order
      )

    for (i in seq_len(nrow(mat))) {

      y_all <- as.numeric(
        mat[
          i,
          batch_ids
        ]
      )

      y_qc <- as.numeric(
        mat[
          i,
          batch_qc_ids
        ]
      )

      valid_qc <- !is.na(y_qc) &
        is.finite(y_qc) &
        y_qc > 0

      n_valid_qc <- sum(valid_qc)

      rsd_before <- calc_qc_rsd(y_qc)
      rsd_after <- rsd_before
      drift_range_fold <- NA_real_
      status <- "Not corrected"

      if (n_valid_qc >= min_qc) {

        fit_df <- data.frame(
          Run_Order = x_qc[valid_qc],
          Log2_Intensity = log2(y_qc[valid_qc])
        )

        fit_degree <- ifelse(
          n_valid_qc >= 6,
          2,
          1
        )

        loess_fit <- loess(
          Log2_Intensity ~ Run_Order,
          data = fit_df,
          span = span,
          degree = fit_degree,
          family = "symmetric",
          normalize = TRUE,
          na.action = na.exclude,
          control = loess.control(
            surface = "direct"
          )
        )

        predicted_all <- as.numeric(
          predict(
            loess_fit,
            newdata = data.frame(
              Run_Order = x_all
            )
          )
        )

        finite_prediction <- is.finite(predicted_all)

        if (sum(finite_prediction) >= 2) {

          predicted_all <- approx(
            x = x_all[finite_prediction],
            y = predicted_all[finite_prediction],
            xout = x_all,
            rule = 2,
            ties = mean
          )$y

          qc_position_in_batch <- match(
            batch_qc_ids[valid_qc],
            batch_ids
          )

          predicted_qc <- predicted_all[
            qc_position_in_batch
          ]

          reference_log2 <- median(
            predicted_qc,
            na.rm = TRUE
          )

          valid_all <- !is.na(y_all) &
            is.finite(y_all) &
            y_all > 0

          corrected[
            i,
            batch_ids[valid_all]
          ] <- 2^(
            log2(y_all[valid_all]) -
              predicted_all[valid_all] +
              reference_log2
          )

          rsd_after <- calc_qc_rsd(
            corrected[
              i,
              batch_qc_ids
            ]
          )

          drift_range_fold <- 2^(
            max(predicted_qc) -
              min(predicted_qc)
          )

          status <- "Corrected"
        }
      }

      diagnostic_list[[diagnostic_index]] <- data.frame(
        Feature = rownames(mat)[i],
        Batch = batch_name,
        N_QC_Valid = n_valid_qc,
        QC_RSD_Before = rsd_before,
        QC_RSD_After = rsd_after,
        Drift_Range_Fold = drift_range_fold,
        Status = status,
        stringsAsFactors = FALSE
      )

      diagnostic_index <- diagnostic_index + 1
    }
  }

  list(
    corrected = corrected,
    diagnostics = bind_rows(diagnostic_list)
  )
}

calc_tech_rsd_matrix <- function(
    mat,
    tech_info
) {

  pair_levels <- unique(tech_info$Tech_Pair)

  rsd_mat <- sapply(
    pair_levels,
    function(pair_id) {

      pair_ids <- tech_info %>%
        filter(
          Tech_Pair == pair_id
        ) %>%
        arrange(
          Tech_Replicate
        ) %>%
        pull(
          Sample_ID
        )

      apply(
        mat[
          ,
          pair_ids,
          drop = FALSE
        ],
        1,
        calc_tech_rsd
      )
    }
  )

  if (is.null(dim(rsd_mat))) {
    rsd_mat <- matrix(
      rsd_mat,
      ncol = 1
    )
  }

  rownames(rsd_mat) <- rownames(mat)
  colnames(rsd_mat) <- paste0(
    "RSD_",
    pair_levels
  )

  rsd_mat
}

run_pca <- function(
    mat,
    stage_label,
    sample_metadata,
    fill_value_override = NULL
) {

  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  mat[
    !is.finite(mat)
  ] <- NA_real_

  fill_value <- if (
    is.null(fill_value_override)
  ) {
    get_fill_value(mat)
  } else {
    fill_value_override
  }

  mat_imputed <- mat
  mat_imputed[
    is.na(mat_imputed)
  ] <- fill_value

  mat_log2 <- log2(mat_imputed)

  feature_sd <- apply(
    mat_log2,
    1,
    sd
  )

  keep_feature <- is.finite(feature_sd) &
    feature_sd > 0

  pca_mat <- mat_log2[
    keep_feature,
    ,
    drop = FALSE
  ]

  pca_model <- prcomp(
    t(pca_mat),
    center = TRUE,
    scale. = pca_scale
  )

  variance_pct <- 100 *
    pca_model$sdev^2 /
    sum(pca_model$sdev^2)

  scores <- as.data.frame(
    pca_model$x[
      ,
      1:2,
      drop = FALSE
    ]
  ) %>%
    rownames_to_column(
      "Sample_ID"
    ) %>%
    left_join(
      sample_metadata,
      by = "Sample_ID"
    ) %>%
    mutate(
      Stage = stage_label
    )

  variance_table <- data.frame(
    PC = paste0(
      "PC",
      seq_along(variance_pct)
    ),
    Variance = variance_pct,
    Stage = stage_label
  )

  list(
    stage = stage_label,
    scores = scores,
    variance = variance_pct,
    variance_table = variance_table,
    fill_value = fill_value,
    n_feature_input = nrow(mat),
    n_feature_used = nrow(pca_mat),
    n_sample = ncol(pca_mat)
  )
}

make_qc_pca_plot <- function(
    pca_result,
    run_order_limits
) {

  ggplot(
    pca_result$scores,
    aes(
      x = PC1,
      y = PC2
    )
  ) +
    geom_point(
      aes(
        fill = Run_Order
      ),
      shape = 21,
      color = "black",
      stroke = 0.55,
      size = 4.3,
      alpha = 1
    ) +
    geom_text(
      aes(
        label = Sample_ID
      ),
      vjust = -1,
      size = 3,
      color = "black",
      show.legend = FALSE
    ) +
    scale_fill_viridis_c(
      name = "Run order",
      option = "C",
      direction = 1,
      limits = run_order_limits
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0.12, 0.12)
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.12, 0.20)
      )
    ) +
    labs(
      title = pca_result$stage,
      subtitle = sprintf(
        "%d QC samples; %d input features; %d PCA features",
        pca_result$n_sample,
        pca_result$n_feature_input,
        pca_result$n_feature_used
      ),
      x = sprintf(
        "PC1 (%.1f%%)",
        pca_result$variance[1]
      ),
      y = sprintf(
        "PC2 (%.1f%%)",
        pca_result$variance[2]
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(
        10,
        15,
        10,
        10
      )
    )
}

make_bio_pca_plot <- function(
    pca_result,
    pmi_limits
) {

  ggplot(
    pca_result$scores,
    aes(
      x = PC1,
      y = PC2
    )
  ) +
    geom_point(
      aes(
        fill = PMI_days_num
      ),
      shape = 21,
      color = "black",
      stroke = 0.45,
      size = 3.4,
      alpha = 0.95
    ) +
    scale_fill_viridis_c(
      name = "PMI (days)",
      option = "C",
      direction = 1,
      limits = pmi_limits
    ) +
    facet_wrap(
      ~ Tissue
    ) +
    labs(
      title = pca_result$stage,
      subtitle = sprintf(
        "%d biological samples; %d PCA features",
        pca_result$n_sample,
        pca_result$n_feature_used
      ),
      x = sprintf(
        "PC1 (%.1f%%)",
        pca_result$variance[1]
      ),
      y = sprintf(
        "PC2 (%.1f%%)",
        pca_result$variance[2]
      )
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      legend.position = "bottom"
    )
}



rlsc_result <- qc_rlsc_correct_batchwise(
  mat = feat_mat_endo_raw,
  sample_info = sample_info,
  span = rlsc_span,
  min_qc = min_qc_per_batch
)

feat_mat_endo_corrected <- rlsc_result$corrected
rlsc_diagnostics <- rlsc_result$diagnostics

write.csv(
  rlsc_diagnostics,
  file.path(
    output_dir,
    "QC_RLSC_feature_batch_diagnostics.csv"
  ),
  row.names = FALSE
)

qc_rsd_diagnostic <- rlsc_diagnostics %>%
  filter(
    is.finite(QC_RSD_Before),
    is.finite(QC_RSD_After)
  ) %>%
  mutate(
    RSD_Decreased = QC_RSD_After <
      QC_RSD_Before,
    Below_30_Before = QC_RSD_Before <
      30,
    Below_30_After = QC_RSD_After <
      30
  )

write.csv(
  qc_rsd_diagnostic,
  file.path(
    output_dir,
    "QC_RSD_before_after_diagnostic.csv"
  ),
  row.names = FALSE
)






rsd_mat_tech_before <- calc_tech_rsd_matrix(
  mat = feat_mat_endo_raw,
  tech_info = tech_info
)

rsd_mat_tech_after <- calc_tech_rsd_matrix(
  mat = feat_mat_endo_corrected,
  tech_info = tech_info
)

n_valid_tech_pairs <- rowSums(
  is.finite(rsd_mat_tech_after)
)

max_tech_rsd_after <- apply(
  rsd_mat_tech_after,
  1,
  function(x) {

    valid_x <- x[
      is.finite(x)
    ]

    if (length(valid_x) == 0) {
      return(NA_real_)
    }

    max(valid_x)
  }
)


all_tech_rsd_missing <- n_valid_tech_pairs == 0


pass_tech <- apply(
  rsd_mat_tech_after,
  1,
  function(x) {
    
    valid_x <- x[
      !is.na(x) &
        is.finite(x)
    ]
    
    !any(
      valid_x > tech_rsd_cutoff
    )
  }
)

rsd_before_out <- as.data.frame(
  rsd_mat_tech_before,
  check.names = FALSE
)

colnames(rsd_before_out) <- paste0(
  colnames(rsd_before_out),
  "_Before"
)

rsd_after_out <- as.data.frame(
  rsd_mat_tech_after,
  check.names = FALSE
)

colnames(rsd_after_out) <- paste0(
  colnames(rsd_after_out),
  "_After"
)

tech_filter_result <- bind_cols(
  data.frame(
    Feature = rownames(
      feat_mat_endo_corrected
    ),
    stringsAsFactors = FALSE
  ),
  rsd_before_out,
  rsd_after_out,
  data.frame(
    N_Valid_Technical_Pairs = n_valid_tech_pairs,
    All_Technical_RSD_Missing = all_tech_rsd_missing,
    Max_Technical_RSD_After = max_tech_rsd_after,
    Pass_Technical_Filter = pass_tech,
    
    Filter_Reason = case_when(
      
      !pass_tech ~ paste0(
        "Removed: at least one corrected technical RSD > ",
        tech_rsd_cutoff,
        "%"
      ),
      
      all_tech_rsd_missing ~
        "Retained, but all corrected technical RSD values are missing",
      
      TRUE ~
        "Retained"
    ),
    
    stringsAsFactors = FALSE
  )
  
)

write.csv(
  tech_filter_result,
  file.path(
    output_dir,
    "technical_replicate_RSD_before_after_and_filter.csv"
  ),
  row.names = FALSE
)

feat_mat_stable_raw <- feat_mat_endo_raw[
  pass_tech,
  ,
  drop = FALSE
]

feat_mat_stable_corrected <- feat_mat_endo_corrected[
  pass_tech,
  ,
  drop = FALSE
]



bio_mat_corrected <- feat_mat_stable_corrected[
  ,
  bio_cols,
  drop = FALSE
]

bio_mat_raw_comparison <- feat_mat_stable_raw[
  ,
  bio_cols,
  drop = FALSE
]

group_missing <- sapply(
  groups,
  function(g) {

    g_samples <- bio_info %>%
      filter(
        Group == g
      ) %>%
      pull(
        Sample_ID
      )

    rowMeans(
      is.na(
        bio_mat_corrected[
          ,
          g_samples,
          drop = FALSE
        ]
      )
    )
  }
)

if (is.null(dim(group_missing))) {
  group_missing <- matrix(
    group_missing,
    ncol = 1,
    dimnames = list(
      rownames(bio_mat_corrected),
      groups[1]
    )
  )
}

rownames(group_missing) <- rownames(
  bio_mat_corrected
)

colnames(group_missing) <- groups

pass_missing <- apply(
  group_missing,
  1,
  function(x) {
    any(
      x <= group_missing_cutoff,
      na.rm = TRUE
    )
  }
)

bio_mat_filtered_corrected <- bio_mat_corrected[
  pass_missing,
  ,
  drop = FALSE
]

bio_mat_filtered_raw <- bio_mat_raw_comparison[
  pass_missing,
  ,
  drop = FALSE
]

cat(
  sprintf(
    "缺失率筛选：%d → %d 个特征\n",
    nrow(bio_mat_corrected),
    nrow(bio_mat_filtered_corrected)
  )
)


ppm_numeric <- suppressWarnings(
  as.numeric(
    as.character(
      feat_meta[
        rownames(bio_mat_filtered_corrected),
        "Mass Error (ppm)"
      ]
    )
  )
)

pass_ppm <- is.na(ppm_numeric) |
  !is.finite(ppm_numeric) |
  abs(ppm_numeric) <= 5

ppm_filter_result <- data.frame(
  Feature = rownames(
    bio_mat_filtered_corrected
  ),
  Mass_Error_ppm = ppm_numeric,
  Pass_PPM_Filter = pass_ppm,
  Filter_Reason = ifelse(
    !pass_ppm,
    "Removed: absolute mass error > 5 ppm",
    "Retained"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  ppm_filter_result,
  file.path(
    output_dir,
    "mass_error_ppm_filter_result.csv"
  ),
  row.names = FALSE
)

n_before_ppm <- nrow(
  bio_mat_filtered_corrected
)

bio_mat_filtered_corrected <-
  bio_mat_filtered_corrected[
    pass_ppm,
    ,
    drop = FALSE
  ]

bio_mat_filtered_raw <-
  bio_mat_filtered_raw[
    pass_ppm,
    ,
    drop = FALSE
  ]





fill_value <- get_fill_value(
  bio_mat_filtered_corrected
)

bio_mat_imputed <- bio_mat_filtered_corrected

bio_mat_imputed[
  is.na(bio_mat_imputed)
] <- fill_value

bio_log2 <- log2(
  bio_mat_imputed
)

final_features <- rownames(
  bio_log2
)

final_meta <- feat_meta[
  final_features,
  ,
  drop = FALSE
]



bio_log2_output <- bio_log2[
  ,
  input_bio_order,
  drop = FALSE
]

write.csv(
  cbind(
    final_meta,
    bio_log2_output
  ),
  file.path(
    output_dir,
    "feature_table_processed_QC_RLSC_log2.csv"
  )
)

bio_corrected_output <- bio_mat_filtered_corrected[
  ,
  input_bio_order,
  drop = FALSE
]

write.csv(
  cbind(
    final_meta,
    bio_corrected_output
  ),
  file.path(
    output_dir,
    "feature_table_QC_RLSC_corrected_before_imputation.csv"
  )
)

write.csv(
  cbind(
    feat_meta[
      rownames(feat_mat_endo_corrected),
      ,
      drop = FALSE
    ],
    feat_mat_endo_corrected
  ),
  file.path(
    output_dir,
    "feature_table_endogenous_QC_RLSC_corrected_all_samples.csv"
  )
)

write.csv(
  sample_info,
  file.path(
    output_dir,
    "sample_information_with_run_order.csv"
  ),
  row.names = FALSE
)




n_raw <- nrow(feat_mat)
n_endo <- nrow(feat_mat_endo_raw)
n_after_tech <- nrow(feat_mat_stable_corrected)
n_after_missing <- sum(pass_missing)
n_final <- nrow(bio_mat_filtered_corrected)

max_valid_rsd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

max_tech_rsd_before <- apply(
  rsd_mat_tech_before,
  1,
  max_valid_rsd
)

max_tech_rsd_after <- apply(
  rsd_mat_tech_after,
  1,
  max_valid_rsd
)

qc_rsd_feature <- rlsc_diagnostics %>%
  group_by(Feature) %>%
  summarise(
    QC_RSD_Before = median(QC_RSD_Before, na.rm = TRUE),
    QC_RSD_After = median(QC_RSD_After, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    QC_RSD_Before = ifelse(
      is.finite(QC_RSD_Before),
      QC_RSD_Before,
      NA_real_
    ),
    QC_RSD_After = ifelse(
      is.finite(QC_RSD_After),
      QC_RSD_After,
      NA_real_
    )
  )

stability_matrix <- apply(
  rsd_mat_tech_after,
  2,
  function(col) {
    case_when(
      !is.finite(col) ~ "absent",
      col < tech_rsd_cutoff ~ "stable",
      TRUE ~ "unstable"
    )
  }
)

if (is.null(dim(stability_matrix))) {
  stability_matrix <- matrix(
    stability_matrix,
    ncol = 1,
    dimnames = list(
      rownames(rsd_mat_tech_after),
      colnames(rsd_mat_tech_after)[1]
    )
  )
}

rownames(stability_matrix) <- rownames(rsd_mat_tech_after)
colnames(stability_matrix) <- colnames(rsd_mat_tech_after)


# FigS


steps <- tibble(
  Step = c(
    "Raw features",
    "Endogenous\nfeatures",
    "QC-RLSC\ncorrected",
    "Technical replicate\nstability filtered",
    "Missing value\nfiltered",
    "Mass-error\nfiltered"
  ),
  Count = c(
    n_raw,
    n_endo,
    n_endo,
    n_after_tech,
    n_after_missing,
    n_final
  )
) %>%
  mutate(
    Lost = c(
      0,
      head(Count, -1) - tail(Count, -1)
    ),
    Step = factor(
      Step,
      levels = Step
    )
  )

p_waterfall <- ggplot(
  steps,
  aes(
    x = Step,
    y = Count
  )
) +
  geom_col(
    fill = c(
      "#4393C3",
      "#74C476",
      "#9ECAE1",
      "#6A51A3",
      "#FD8D3C",
      "#D73027"
    ),
    width = 0.62
  ) +
  geom_text(
    aes(label = Count),
    vjust = -0.45,
    size = 4.3,
    fontface = "bold"
  ) +
  geom_text(
    aes(
      label = ifelse(
        Lost > 0,
        paste0("\u2212", Lost),
        ""
      ),
      y = Count + max(Count) * 0.035
    ),
    vjust = -1.4,
    size = 3.3,
    color = "red"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.20)
    )
  ) +
  labs(
    title = "Feature count at each preprocessing step",
    subtitle = "QC-RLSC corrects signal drift without directly removing features",
    x = NULL,
    y = "Feature count"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

ggsave(
  file.path(
    output_dir,
    "FigS_A_feature_filtering_waterfall.TIFF"
  ),
  p_waterfall,
  width = 9,
  height = 5.5,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


#  FigS：RSD


rsd_df <- bind_rows(
  qc_rsd_feature %>%
    transmute(
      RSD = QC_RSD_Before,
      Type = "Pooled QC RSD",
      Status = "Before QC-RLSC"
    ),
  qc_rsd_feature %>%
    transmute(
      RSD = QC_RSD_After,
      Type = "Pooled QC RSD",
      Status = "After QC-RLSC"
    ),
  tibble(
    RSD = max_tech_rsd_before,
    Type = "Technical-replicate RSD (maximum)",
    Status = "Before QC-RLSC"
  ),
  tibble(
    RSD = max_tech_rsd_after,
    Type = "Technical-replicate RSD (maximum)",
    Status = "After QC-RLSC"
  ),
  tibble(
    RSD = max_tech_rsd_after[pass_tech],
    Type = "Technical-replicate RSD (maximum)",
    Status = "Retained after technical filtering"
  )
) %>%
  filter(
    is.finite(RSD),
    RSD >= 0
  ) %>%
  mutate(
    Status = factor(
      Status,
      levels = c(
        "Before QC-RLSC",
        "After QC-RLSC",
        "Retained after technical filtering"
      )
    )
  )

p_rsd <- ggplot(
  rsd_df,
  aes(
    x = RSD,
    fill = Status
  )
) +
  geom_density(
    alpha = 0.52,
    color = NA,
    adjust = 1
  ) +
  geom_vline(
    xintercept = tech_rsd_cutoff,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = tech_rsd_cutoff + 2,
    y = Inf,
    label = "RSD = 30%",
    hjust = 0,
    vjust = 1.5,
    color = "red",
    size = 3.4
  ) +
  scale_fill_manual(
    values = c(
      "Before QC-RLSC" = "#BDBDBD",
      "After QC-RLSC" = "#FC8D59",
      "Retained after technical filtering" = "#2171B5"
    ),
    drop = FALSE
  ) +
  facet_wrap(
    ~Type,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    title = "RSD distributions before and after QC-RLSC correction",
    x = "RSD (%)",
    y = "Density",
    fill = "Processing stage"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    output_dir,
    "FigS_B_RSD_distribution.TIFF"
  ),
  p_rsd,
  width = 7.5,
  height = 7.5,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

# 14. FigS_
scatter_df <- qc_rsd_feature %>%
  select(
    Feature,
    RSD_QC = QC_RSD_After
  ) %>%
  left_join(
    tibble(
      Feature = rownames(rsd_mat_tech_after),
      RSD_Tech = max_tech_rsd_after,
      Pass = ifelse(
        pass_tech,
        "Technical filter passed",
        "Technical filter failed"
      )
    ),
    by = "Feature"
  ) %>%
  filter(
    is.finite(RSD_QC),
    is.finite(RSD_Tech)
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = RSD_QC,
    y = RSD_Tech,
    color = Pass
  )
) +
  geom_point(
    alpha = 0.52,
    size = 1.6
  ) +
  geom_vline(
    xintercept = 30,
    linetype = "dashed",
    color = "red"
  ) +
  geom_hline(
    yintercept = tech_rsd_cutoff,
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 30,
    ymin = -Inf,
    ymax = tech_rsd_cutoff,
    fill = "#2171B5",
    alpha = 0.045
  ) +
  annotate(
    "text",
    x = 1,
    y = tech_rsd_cutoff - 2,
    hjust = 0,
    vjust = 1,
    label = "Both RSDs <30% (diagnostic region)",
    color = "#2171B5",
    size = 3.2
  ) +
  scale_color_manual(
    values = c(
      "Technical filter passed" = "#2171B5",
      "Technical filter failed" = "#FC8D59"
    )
  ) +
  labs(
    title = "Post-correction QC and technical-replicate RSDs",
    subtitle = "Only technical-replicate RSD was used as a hard stability filter",
    x = "Pooled QC RSD after QC-RLSC (%)",
    y = "Maximum technical-replicate RSD after QC-RLSC (%)",
    color = NULL
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    output_dir,
    "FigS_B2_postcorrection_RSD_scatter.TIFF"
  ),
  p_scatter,
  width = 6.5,
  height = 6.2,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)



set.seed(42)

failed_ids <- rownames(stability_matrix)[
  !pass_tech
]
passed_ids <- rownames(stability_matrix)[
  pass_tech
]

show_failed <- if (
  length(failed_ids) > 0L
) {
  sample(
    failed_ids,
    min(40L, length(failed_ids))
  )
} else {
  character()
}

show_passed <- if (
  length(passed_ids) > 0L
) {
  sample(
    passed_ids,
    min(40L, length(passed_ids))
  )
} else {
  character()
}

show_ids <- unique(
  c(
    show_failed,
    show_passed
  )
)

if (length(show_ids) < min(80L, nrow(stability_matrix))) {
  remaining_ids <- setdiff(
    rownames(stability_matrix),
    show_ids
  )

  add_n <- min(
    80L - length(show_ids),
    length(remaining_ids)
  )

  if (add_n > 0L) {
    show_ids <- c(
      show_ids,
      sample(
        remaining_ids,
        add_n
      )
    )
  }
}

if (length(show_ids) > 0L) {
  status_num <- stability_matrix[
    show_ids,
    ,
    drop = FALSE
  ]

  status_num_plot <- matrix(
    ifelse(
      status_num == "stable",
      1,
      ifelse(
        status_num == "absent",
        0,
        ifelse(
          status_num == "unstable",
          -1,
          NA_real_
        )
      )
    ),
    nrow = nrow(status_num),
    dimnames = dimnames(status_num)
  )

  tiff(
    file.path(
      output_dir,
      "FigS_B3_technical_replicate_stability_heatmap.TIFF"
    ),
    width = 7,
    height = 9,
    units = "in",
    res = 600,
    compression = "lzw"
  )

  pheatmap(
    status_num_plot,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = c(
      "#FC8D59",
      "#FFFFBF",
      "#2171B5"
    ),
    breaks = c(
      -1.5,
      -0.5,
      0.5,
      1.5
    ),
    legend_breaks = c(
      -1,
      0,
      1
    ),
    legend_labels = c(
      "Unstable",
      "Absent",
      "Stable"
    ),
    main = "Technical-replicate stability after QC-RLSC",
    show_rownames = FALSE,
    fontsize_col = 10,
    border_color = "white"
  )

  dev.off()
}


p_qc_rsd_scatter <- ggplot(
  qc_rsd_diagnostic,
  aes(
    x = QC_RSD_Before,
    y = QC_RSD_After
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_point(
    shape = 21,
    fill = "#6A00A8",
    color = "black",
    size = 1.8,
    alpha = 0.45,
    stroke = 0.25
  ) +
  geom_hline(
    yintercept = 30,
    linetype = "dashed",
    color = "#D73027"
  ) +
  geom_vline(
    xintercept = 30,
    linetype = "dashed",
    color = "#D73027"
  ) +
  labs(
    title = "Pooled QC RSD before and after QC-RLSC correction",
    subtitle = "Features below the diagonal showed reduced QC RSD after correction",
    x = "QC RSD before correction (%)",
    y = "QC RSD after correction (%)"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    output_dir,
    "FigS_C_QC_RSD_before_after_scatter.TIFF"
  ),
  p_qc_rsd_scatter,
  width = 7.5,
  height = 5.8,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)



make_boxplot_df <- function(
  mat,
  stage_label
) {
  as.data.frame(
    mat,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "Feature"
    ) %>%
    pivot_longer(
      -Feature,
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    mutate(
      Stage = stage_label
    )
}

box_df <- bind_rows(
  make_boxplot_df(
    bio_mat_imputed,
    "Corrected peak area (imputed)"
  ),
  make_boxplot_df(
    bio_log2,
    "After log2 transformation"
  )
) %>%
  left_join(
    bio_info %>%
      select(
        Sample_ID,
        PMI_days,
        Tissue
      ),
    by = c(
      "Sample" = "Sample_ID"
    )
  )

p_box <- ggplot(
  box_df,
  aes(
    x = Sample,
    y = Intensity,
    fill = factor(PMI_days)
  )
) +
  geom_boxplot(
    outlier.size = 0.3,
    outlier.alpha = 0.3
  ) +
  facet_wrap(
    ~Stage,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_fill_brewer(
    palette = "RdYlBu",
    name = "PMI (days)"
  ) +
  labs(
    title = "Sample intensity distributions before and after log2 transformation",
    x = NULL,
    y = "Intensity"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    legend.position = "right"
  )

ggsave(
  file.path(
    output_dir,
    "FigS_D_intensity_boxplot.TIFF"
  ),
  p_box,
  width = 13,
  height = 8,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)



sample_medians <- apply(
  bio_log2,
  2,
  median,
  na.rm = TRUE
)

median_df <- tibble(
  Sample = names(sample_medians),
  Median = as.numeric(sample_medians)
) %>%
  left_join(
    bio_info %>%
      select(
        Sample_ID,
        Tissue,
        PMI_days
      ),
    by = c(
      "Sample" = "Sample_ID"
    )
  )

p_median <- ggplot(
  median_df,
  aes(
    x = factor(PMI_days),
    y = Median,
    fill = factor(PMI_days)
  )
) +
  geom_boxplot(
    outlier.size = 1
  ) +
  facet_wrap(
    ~Tissue,
    ncol = 5
  ) +
  scale_fill_brewer(
    palette = "RdYlBu",
    name = "PMI (days)"
  ) +
  labs(
    title = "Median log2 intensity across PMI groups",
    subtitle = "No between-sample normalization was applied",
    x = "PMI (days)",
    y = "Median log2 intensity"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    output_dir,
    "FigS_E_sample_median_by_PMI.TIFF"
  ),
  p_median,
  width = 14,
  height = 5.5,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


# 19. FigS_A


qc_info_pca <- sample_info %>%
  filter(
    Sample_Type == "QC"
  )

qc_mat_before <- feat_mat_endo_raw[
  ,
  qc_info_pca$Sample_ID,
  drop = FALSE
]

qc_mat_after <- feat_mat_endo_corrected[
  ,
  qc_info_pca$Sample_ID,
  drop = FALSE
]

qc_common_fill <- get_fill_value(
  c(
    qc_mat_before,
    qc_mat_after
  )
)

qc_pca_before <- run_pca(
  mat = qc_mat_before,
  stage_label = "QC before QC-RLSC correction",
  sample_metadata = qc_info_pca,
  fill_value_override = qc_common_fill
)

qc_pca_after <- run_pca(
  mat = qc_mat_after,
  stage_label = "QC after QC-RLSC correction",
  sample_metadata = qc_info_pca,
  fill_value_override = qc_common_fill
)

qc_run_limits <- range(
  qc_info_pca$Run_Order,
  na.rm = TRUE
)

p_qc_before <- make_qc_pca_plot(
  qc_pca_before,
  qc_run_limits
)

p_qc_after <- make_qc_pca_plot(
  qc_pca_after,
  qc_run_limits
)

p_qc_comparison <- (
  p_qc_before +
    p_qc_after +
    plot_layout(
      ncol = 2,
      guides = "collect"
    )
) +
  plot_annotation(
    title = "QC PCA before and after QC-RLSC correction",
    subtitle = paste0(
      "The same QC injections and endogenous features are shown; ",
      "point fill indicates run order."
    )
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    output_dir,
    "FigS_F_QC_PCA_before_after_QC_RLSC.TIFF"
  ),
  p_qc_comparison,
  width = 14,
  height = 6.5,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)




pattern_df <- tibble(
  Feature = rownames(group_missing),
  N_Detected_Groups = apply(
    group_missing,
    1,
    function(x) {
      sum(
        is.finite(x) &
          x <= group_missing_cutoff
      )
    }
  ),
  N_Evaluable_Groups = apply(
    group_missing,
    1,
    function(x) {
      sum(
        is.finite(x)
      )
    }
  )
) %>%
  mutate(
    Detection_Ratio = ifelse(
      N_Evaluable_Groups > 0,
      N_Detected_Groups /
        N_Evaluable_Groups,
      NA_real_
    )
  ) %>%
  filter(
    Feature %in% final_features
  )

write.csv(
  pattern_df,
  file.path(
    output_dir,
    "feature_detection_pattern.csv"
  ),
  row.names = FALSE
)

qc_scores_before_after <- bind_rows(
  qc_pca_before$scores,
  qc_pca_after$scores
)

qc_variance_before_after <- bind_rows(
  qc_pca_before$variance_table,
  qc_pca_after$variance_table
)

bio_scores_before_after <- bind_rows(
  bio_pca_before$scores,
  bio_pca_after$scores
)

bio_variance_before_after <- bind_rows(
  bio_pca_before$variance_table,
  bio_pca_after$variance_table
)

write.csv(
  qc_scores_before_after,
  file.path(
    output_dir,
    "QC_PCA_scores_before_after.csv"
  ),
  row.names = FALSE
)

write.csv(
  qc_variance_before_after,
  file.path(
    output_dir,
    "QC_PCA_variance_before_after.csv"
  ),
  row.names = FALSE
)

write.csv(
  bio_scores_before_after,
  file.path(
    output_dir,
    "biological_PCA_scores_before_after.csv"
  ),
  row.names = FALSE
)

write.csv(
  bio_variance_before_after,
  file.path(
    output_dir,
    "biological_PCA_variance_before_after.csv"
  ),
  row.names = FALSE
)


