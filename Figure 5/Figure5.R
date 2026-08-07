setwd("D:/PMI_Project/Figure 5")
output_dir <- file.path(getwd(), "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(randomForest)
  library(patchwork)
  library(ggtext)
})

TISSUES <- c("Heart", "Liver", "Spleen", "Lung", "Muscle")
COL_ID <- "SampleID"
COL_ANIMAL <- "AnimalID"
COL_TISSUE <- "Organ"
COL_TARGET <- "PMI_day"

META_REQUIRED <- c(COL_ID, COL_ANIMAL, COL_TISSUE, COL_TARGET)
META_EXCLUDE <- unique(c(META_REQUIRED, "PMI_bin", "PMI", "Day"))

PREV_CUT <- 0.20
ABUND_CUT <- 0.001
PSEUDO <- 1e-6

METAB_ALGORITHM <- "ElasticNet"
MICRO_ALGORITHM <- "RandomForest"
EF_ALGORITHM <- "ElasticNet"

ALPHA <- 0.5
K_METAB_PRIMARY <- 10L
K_MICRO_PRIMARY <- 5L

INNER_FOLDS <- 3L
INNER_REPEATS <- 3L
STACK_FOLDS <- 3L

RF_NTREE <- 500L
RF_NODESIZE <- 5L

SEED_BASE <- 1000L
BOOT_B <- 2000L


SINGLE_MULTI_COMPARISON_B <- 9999L
SINGLE_MULTI_ALPHA <- 0.05
SINGLE_MULTI_P_ADJUST <- "holm"

MULTI_ORGAN_COMPARISON_B <- 1000L


MULTI_ORGAN_OOF_INNER_REPEATS <- 1L

PRIMARY_MODELS <- c(
  "Metab_FixedK",
  "Micro_FixedK",
  "EF_FixedK",
  "Stacking_FixedK"
)

ALL_BENCHMARK_MODELS <- c(
  "Metab_ALL",
  "Micro_ALL",
  "EF_ALL",
  "Stacking_ALL"
)

CORE_MODELS <- c(PRIMARY_MODELS, ALL_BENCHMARK_MODELS)

MODEL_LABELS <- c(
  "Metab_FixedK" = "Metabolomics: fixed K",
  "Micro_FixedK" = "Microbiome: fixed K",
  "EF_FixedK" = "Early fusion: fixed K",
  "Stacking_FixedK" = "Stacking: fixed K",
  "Metab_ALL" = "Metabolomics: all retained features",
  "Micro_ALL" = "Microbiome: all retained features",
  "EF_ALL" = "Early fusion: all retained features",
  "Stacking_ALL" = "Stacking: all retained features"
)

MODEL_COLORS <- c(
  "Metab_FixedK" = "#F4A460",
  "Micro_FixedK" = "#9370DB",
  "EF_FixedK" = "#2C7BB6",
  "Stacking_FixedK" = "#D7191C",
  "Metab_ALL" = "#F4A460",
  "Micro_ALL" = "#9370DB",
  "EF_ALL" = "#2C7BB6",
  "Stacking_ALL" = "#D7191C"
)

calc_r2 <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]

  if (length(observed) < 2L) return(NA_real_)

  den <- sum((observed - mean(observed))^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)

  1 - sum((observed - predicted)^2) / den
}

calc_rmse <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((observed[ok] - predicted[ok])^2))
}

calc_mae <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  mean(abs(observed[ok] - predicted[ok]))
}

safe_cor <- function(x, y) {
  out <- suppressWarnings(
    cor(x, y, use = "pairwise.complete.obs")
  )
  ifelse(is.finite(out), out, 0)
}

stable_seed <- function(...) {
  txt <- enc2utf8(
    paste(..., collapse = "__")
  )
  values <- utf8ToInt(txt)

  
  hash <- as.double(SEED_BASE)

  for (value in values) {
    hash <- (hash * 131 + as.double(value)) %% 2147483646
  }

  as.integer(hash) + 1L
}

id_set_key <- function(ids) {
  ids <- sort(unique(as.character(ids)))
  paste(ids, collapse = "|")
}

model_partition_seed <- function(
  stage,
  tissue,
  train_ids,
  test_ids
) {
  stable_seed(
    stage,
    tissue,
    id_set_key(train_ids),
    id_set_key(test_ids)
  )
}

select_min_rmse <- function(
  summary_df,
  rmse_col = "PooledRMSE",
  tie_cols = character()
) {
  valid <- summary_df %>%
    filter(is.finite(.data[[rmse_col]]))

  order_cols <- c(rmse_col, tie_cols)

  valid %>%
    arrange(!!!rlang::syms(order_cols)) %>%
    slice(1)
}


# LOAO and grouped inner-fold tools
make_outer_splits <- function(meta) {
  animals <- sort(
    unique(
      as.character(meta[[COL_ANIMAL]])
    )
  )

  lapply(seq_along(animals), function(i) {
    test_id <- animals[i]

    list(
      OuterFold = i,
      Seed = stable_seed(
        "OuterLOAO",
        test_id
      ),
      TrainIDs = animals[-i],
      TestIDs = test_id
    )
  })
}

make_group_folds <- function(
  groups,
  y,
  kfold = 3L,
  seed = 42L
) {
  groups <- as.character(groups)
  y <- as.numeric(y)

  split_y <- split(y, groups)
  n_unique_y <- vapply(
    split_y,
    function(z) length(unique(round(z, 8))),
    integer(1)
  )

   group_y <- vapply(
    split_y,
    function(z) z[1],
    numeric(1)
  )

  n_groups <- length(group_y)
  kfold <- min(as.integer(kfold), n_groups)

  

  set.seed(seed)
  ord <- order(group_y, runif(n_groups))
  group_fold <- integer(n_groups)

  blocks <- split(
    ord,
    ceiling(seq_along(ord) / kfold)
  )

  for (b in blocks) {
    group_fold[b] <- sample(
      seq_len(kfold),
      length(b),
      replace = FALSE
    )
  }

  names(group_fold) <- names(group_y)
  foldid <- unname(group_fold[groups])

  check <- tapply(
    foldid,
    groups,
    function(z) length(unique(z))
  )

  if (any(check != 1L)) {
    stop(
      "An AnimalID was assigned to more than one inner fold."
    )
  }

  foldid
}

make_repeated_group_folds <- function(
  groups,
  y,
  kfold = INNER_FOLDS,
  repeats = INNER_REPEATS,
  seed = 42L
) {
  groups <- as.character(groups)
  y <- as.numeric(y)

  plans <- lapply(
    seq_len(as.integer(repeats)),
    function(r) {
      foldid <- make_group_folds(
        groups = groups,
        y = y,
        kfold = kfold,
        seed = seed + r * 1009L
      )

      tibble(
        RowIndex = seq_along(groups),
        AnimalID = groups,
        Repeat = as.integer(r),
        InnerFold = as.integer(foldid)
      )
    }
  )

  bind_rows(plans)
}


#  Leakage-safe preprocessing
impute_numeric_matrix <- function(X, medians) {
  X <- as.data.frame(X)

  for (nm in names(medians)) {
    bad <- !is.finite(X[[nm]])
    if (any(bad)) X[[nm]][bad] <- medians[[nm]]
  }

  X
}

preprocess_metab <- function(X_train_raw) {
  X_train_raw <- as.data.frame(X_train_raw)

  medians <- vapply(
    X_train_raw,
    function(z) {
      z <- z[is.finite(z)]
      if (length(z) == 0L) 0 else median(z)
    },
    numeric(1)
  )

  X_imp <- impute_numeric_matrix(
    X_train_raw,
    medians
  )

  vars <- vapply(X_imp, var, numeric(1))

  keep <- names(vars)[
    is.finite(vars) & vars > 0
  ]

  if (length(keep) == 0L) {
    stop(
      "No non-zero-variance metabolite features remain."
    )
  }

  X_keep <- X_imp[, keep, drop = FALSE]
  medians <- medians[keep]

  mu <- vapply(X_keep, mean, numeric(1))
  sdv <- vapply(X_keep, sd, numeric(1))
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  X_scaled <- sweep(
    sweep(
      as.matrix(X_keep),
      2,
      mu,
      "-"
    ),
    2,
    sdv,
    "/"
  )

  list(
    feats = keep,
    medians = medians,
    mu = mu,
    sdv = sdv,
    X_train = as.data.frame(
      X_scaled,
      check.names = FALSE
    )
  )
}

apply_metab <- function(prep, X_test_raw) {
  X_test_raw <- as.data.frame(X_test_raw)
  X_test <- X_test_raw[
    , prep$feats, drop = FALSE
  ]

  X_test <- impute_numeric_matrix(
    X_test,
    prep$medians
  )

  as.data.frame(
    sweep(
      sweep(
        as.matrix(X_test),
        2,
        prep$mu,
        "-"
      ),
      2,
      prep$sdv,
      "/"
    ),
    check.names = FALSE
  )
}

clr_mat <- function(X) {
  mat <- as.matrix(X)
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  mat <- mat + PSEUDO

  log_mat <- log(mat)
  sweep(
    log_mat,
    1,
    rowMeans(log_mat),
    "-"
  )
}

preprocess_micro <- function(X_train_raw) {
  X_train_raw <- as.data.frame(X_train_raw)

  X_train_raw[] <- lapply(
    X_train_raw,
    function(z) {
      z[!is.finite(z)] <- 0
      z
    }
  )

  prevalence <- vapply(
    X_train_raw,
    function(z) mean(z > 0),
    numeric(1)
  )

  abundance <- vapply(
    X_train_raw,
    mean,
    numeric(1)
  )

  keep <- names(prevalence)[
    prevalence >= PREV_CUT &
      abundance > ABUND_CUT
  ]

  if (length(keep) == 0L) {
    ord <- order(
      prevalence,
      abundance,
      decreasing = TRUE
    )
    if (length(ord) > 0L) {
      keep <- names(prevalence)[ord[1]]
    }
  }

  if (length(keep) == 0L) {
    stop(
      "No microbiome features remain after filtering."
    )
  }

  X_keep <- X_train_raw[, keep, drop = FALSE]
  X_clr <- as.data.frame(
    clr_mat(X_keep),
    check.names = FALSE
  )

  mu <- vapply(X_clr, mean, numeric(1))
  sdv <- vapply(X_clr, sd, numeric(1))
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  X_scaled <- sweep(
    sweep(
      as.matrix(X_clr),
      2,
      mu,
      "-"
    ),
    2,
    sdv,
    "/"
  )

  list(
    feats = keep,
    mu = mu,
    sdv = sdv,
    X_train = as.data.frame(
      X_scaled,
      check.names = FALSE
    )
  )
}

apply_micro <- function(prep, X_test_raw) {
  X_test_raw <- as.data.frame(X_test_raw)
  X_test <- X_test_raw[
    , prep$feats, drop = FALSE
  ]

  X_test[] <- lapply(
    X_test,
    function(z) {
      z[!is.finite(z)] <- 0
      z
    }
  )

  X_clr <- as.data.frame(
    clr_mat(X_test),
    check.names = FALSE
  )

  as.data.frame(
    sweep(
      sweep(
        as.matrix(X_clr),
        2,
        prep$mu,
        "-"
      ),
      2,
      prep$sdv,
      "/"
    ),
    check.names = FALSE
  )
}

prepare_omics <- function(
  X_train_raw,
  X_test_raw,
  omics_type
) {
  prep <- if (omics_type == "metab") {
    preprocess_metab(X_train_raw)
  } else {
    preprocess_micro(X_train_raw)
  }

  X_test <- if (omics_type == "metab") {
    apply_metab(prep, X_test_raw)
  } else {
    apply_micro(prep, X_test_raw)
  }

  list(
    prep = prep,
    X_train = prep$X_train,
    X_test = X_test
  )
}


#  Fixed algorithms
fit_enet_processed <- function(
  X_train,
  y_train,
  X_test,
  groups = NULL,
  seed = 42L,
  alpha = ALPHA
) {
  X_train <- as.data.frame(X_train)
  X_test <- as.data.frame(X_test)
  y_train <- as.numeric(y_train)

  set.seed(seed)

  if (
    !is.null(groups) &&
      length(groups) == nrow(X_train) &&
      length(unique(groups)) >= 3L
  ) {
    foldid <- make_group_folds(
      groups,
      y_train,
      kfold = min(
        INNER_FOLDS,
        length(unique(groups))
      ),
      seed = seed
    )

    cvfit <- cv.glmnet(
      as.matrix(X_train),
      y_train,
      alpha = alpha,
      foldid = foldid,
      standardize = FALSE
    )
  } else {
    nfolds <- min(
      5L,
      nrow(X_train)
    )

    cvfit <- cv.glmnet(
      as.matrix(X_train),
      y_train,
      alpha = alpha,
      nfolds = nfolds,
      standardize = FALSE
    )
  }

  pred_train <- drop(
    predict(
      cvfit,
      as.matrix(X_train),
      s = "lambda.1se"
    )
  )

  pred_test <- drop(
    predict(
      cvfit,
      as.matrix(X_test),
      s = "lambda.1se"
    )
  )

  cm <- coef(
    cvfit,
    s = "lambda.1se"
  )

  coef_values <- setNames(
    as.numeric(cm)[-1],
    rownames(cm)[-1]
  )

  marginal <- vapply(
    X_train,
    function(z) abs(safe_cor(z, y_train)),
    numeric(1)
  )

  importance <- abs(coef_values)

  if (
    all(!is.finite(importance)) ||
      max(importance, na.rm = TRUE) <= 0
  ) {
    importance <- marginal[names(coef_values)]
  } else {
    marginal <- marginal[names(coef_values)]
    marginal[!is.finite(marginal)] <- 0
    importance[!is.finite(importance)] <- 0

    importance <- importance +
      1e-9 * marginal / max(
        1,
        max(marginal)
      )
  }

  list(
    pred_train = pred_train,
    pred_test = pred_test,
    n_features = as.integer(
      sum(abs(coef_values) > 0)
    ),
    importance = importance,
    lambda_1se = cvfit$lambda.1se,
    lambda_min = cvfit$lambda.min,
    alpha = alpha
  )
}

fit_rf_processed <- function(
  X_train,
  y_train,
  X_test,
  mtry,
  seed = 42L
) {
  X_train <- as.data.frame(X_train)
  X_test <- as.data.frame(X_test)
  y_train <- as.numeric(y_train)

  p <- ncol(X_train)
  if (p < 1L) {
    stop("Random Forest received no predictors.")
  }

  mtry <- max(
    1L,
    min(as.integer(mtry), p)
  )

  set.seed(seed)
  fit <- suppressWarnings(
    randomForest(
      x = as.matrix(X_train),
      y = y_train,
      ntree = RF_NTREE,
      mtry = mtry,
      nodesize = RF_NODESIZE,
      importance = TRUE
    )
  )

  imp_obj <- importance(
    fit,
    type = 1
  )

  importance_values <- if (is.matrix(imp_obj)) {
    imp_obj[, 1]
  } else {
    as.numeric(imp_obj)
  }

  importance_names <- if (is.matrix(imp_obj)) {
    rownames(imp_obj)
  } else {
    names(imp_obj)
  }

  if (
    is.null(importance_names) ||
      length(importance_names) !=
        length(importance_values)
  ) {
    importance_names <- colnames(X_train)
  }

  names(importance_values) <- importance_names

  marginal <- vapply(
    X_train,
    function(z) abs(safe_cor(z, y_train)),
    numeric(1)
  )

  importance_values[
    !is.finite(importance_values) |
      importance_values < 0
  ] <- 0

  marginal <- marginal[
    names(importance_values)
  ]
  marginal[!is.finite(marginal)] <- 0

  if (max(importance_values, na.rm = TRUE) <= 0) {
    importance_values <- marginal
  } else {
    importance_values <- importance_values +
      1e-9 * marginal /
      max(1, max(marginal))
  }

  list(
    pred_train = as.numeric(
      predict(
        fit,
        as.matrix(X_train)
      )
    ),
    pred_test = as.numeric(
      predict(
        fit,
        as.matrix(X_test)
      )
    ),
    n_features = p,
    importance = importance_values,
    mtry = mtry
  )
}

rank_features_processed <- function(
  algorithm,
  X_train,
  y_train,
  groups,
  seed
) {
  if (algorithm == "ElasticNet") {
    rank_fit <- fit_enet_processed(
      X_train,
      y_train,
      X_train,
      groups = groups,
      seed = seed,
      alpha = ALPHA
    )
  } else if (algorithm == "RandomForest") {
    rank_fit <- fit_rf_processed(
      X_train,
      y_train,
      X_train,
      mtry = max(
        1L,
        floor(sqrt(ncol(X_train)))
      ),
      seed = seed
    )
  } else {
    stop("Unknown algorithm: ", algorithm)
  }

  rank_fit$importance
}

fallback_rf_ranking <- function(
  X_train,
  y_train
) {
  importance <- vapply(
    as.data.frame(X_train),
    function(z) {
      abs(
        safe_cor(
          as.numeric(z),
          as.numeric(y_train)
        )
      )
    },
    numeric(1)
  )

  importance[!is.finite(importance)] <- 0
  importance
}

select_top_features <- function(importance, k) {
  importance <- importance[
    is.finite(importance)
  ]

  if (length(importance) == 0L) {
    stop(
      "No valid feature-importance values are available."
    )
  }

  k <- min(
    as.integer(k),
    length(importance)
  )

  names(
    sort(
      importance,
      decreasing = TRUE
    )
  )[seq_len(k)]
}

# RF mtry tuning and fixed/all omics pipelines
RF_MTRY_RULES <- c(
  "one",
  "sqrt",
  "one_third",
  "half",
  "all"
)

rf_mtry_from_rule <- function(rule, p) {
  p <- max(1L, as.integer(p))

  value <- switch(
    as.character(rule),
    one = 1,
    sqrt = sqrt(p),
    one_third = p / 3,
    half = p / 2,
    all = p,
    stop("Unknown RF mtry rule: ", rule)
  )

  max(
    1L,
    min(
      p,
      as.integer(round(value))
    )
  )
}

select_rf_mtry_inner <- function(
  meta_train,
  X_train_raw,
  omics_type,
  inner_plan,
  feature_mode,
  fixed_k = NA_integer_,
  seed,
  stage
) {
  groups <- as.character(
    meta_train[[COL_ANIMAL]]
  )
  y <- as.numeric(
    meta_train[[COL_TARGET]]
  )
  n <- nrow(meta_train)

  split_keys <- inner_plan %>%
    distinct(Repeat, InnerFold) %>%
    arrange(Repeat, InnerFold)

  rows <- vector(
    "list",
    nrow(split_keys) *
      length(RF_MTRY_RULES)
  )
  pos <- 1L

  for (s in seq_len(nrow(split_keys))) {
    repeat_id <- split_keys$Repeat[s]
    fold_id <- split_keys$InnerFold[s]

    val_index <- inner_plan %>%
      filter(
        Repeat == repeat_id,
        InnerFold == fold_id
      ) %>%
      pull(RowIndex)

    is_val <- seq_len(n) %in% val_index
    is_train <- !is_val

    prepared <- tryCatch(
      prepare_omics(
        X_train_raw[
          is_train, , drop = FALSE
        ],
        X_train_raw[
          is_val, , drop = FALSE
        ],
        omics_type
      ),
      error = identity
    )

    ranking <- NULL

    if (
      !inherits(prepared, "error") &&
        feature_mode == "FixedK"
    ) {
      ranking <- tryCatch(
        rank_features_processed(
          "RandomForest",
          prepared$X_train,
          y[is_train],
          groups[is_train],
          seed = seed +
            repeat_id * 100000L +
            fold_id * 1000L
        ),
        error = function(e) {
          fallback_rf_ranking(
            prepared$X_train,
            y[is_train]
          )
        }
      )
    }

    for (rule_index in seq_along(RF_MTRY_RULES)) {
      rule <- RF_MTRY_RULES[rule_index]

      rows[[pos]] <- tryCatch({
        if (inherits(prepared, "error")) {
          stop(conditionMessage(prepared))
        }

        selected_features <- if (
          feature_mode == "ALL"
        ) {
          colnames(prepared$X_train)
        } else {
          select_top_features(
            ranking,
            fixed_k
          )
        }

        mtry_effective <- rf_mtry_from_rule(
          rule,
          length(selected_features)
        )

        rf_fit_seed <-
  seed +
  repeat_id * 1000000L +
  fold_id * 10000L +
  mtry_effective

fit <- fit_rf_processed(
  prepared$X_train[
    , selected_features, drop = FALSE
  ],
  y[is_train],
  prepared$X_test[
    , selected_features, drop = FALSE
  ],
  mtry = mtry_effective,
  seed = rf_fit_seed
)


        residuals <- y[is_val] -
          fit$pred_test

        tibble(
          Stage = stage,
          FeatureMode = feature_mode,
          Repeat = repeat_id,
          InnerFold = fold_id,
          K = ifelse(
            feature_mode == "ALL",
            NA_integer_,
            fixed_k
          ),
          MtryRule = rule,
          MtryRuleOrder = rule_index,
          MtryEffective = mtry_effective,
          NValidation = length(residuals),
          SSE = sum(residuals^2),
          RMSE = sqrt(mean(residuals^2)),
          Error = NA_character_
        )
      }, error = function(e) {
        tibble(
          Stage = stage,
          FeatureMode = feature_mode,
          Repeat = repeat_id,
          InnerFold = fold_id,
          K = ifelse(
            feature_mode == "ALL",
            NA_integer_,
            fixed_k
          ),
          MtryRule = rule,
          MtryRuleOrder = rule_index,
          MtryEffective = NA_integer_,
          NValidation = sum(is_val),
          SSE = Inf,
          RMSE = Inf,
          Error = conditionMessage(e)
        )
      })

      pos <- pos + 1L
    }
  }

  folds <- bind_rows(rows)
  expected_folds <- nrow(split_keys)

  summary <- folds %>%
    group_by(
      Stage,
      FeatureMode,
      K,
      MtryRule,
      MtryRuleOrder
    ) %>%
    summarise(
      PooledRMSE = {
        ok <- is.finite(SSE) &
          NValidation > 0

        if (any(ok)) {
          sqrt(
            sum(SSE[ok]) /
              sum(NValidation[ok])
          )
        } else {
          Inf
        }
      },
      MeanFoldRMSE = mean(
        RMSE[is.finite(RMSE)]
      ),
      SDFoldRMSE = sd(
        RMSE[is.finite(RMSE)]
      ),
      ValidFolds = sum(
        is.finite(RMSE)
      ),
      ExpectedFolds = expected_folds,
      .groups = "drop"
    )

  selected <- summary %>%
    filter(
      ValidFolds > 0L,
      is.finite(PooledRMSE)
    ) %>%
    arrange(
      desc(ValidFolds),
      PooledRMSE,
      MtryRuleOrder
    ) %>%
    slice(1)

  if (nrow(selected) == 0L) {
    selected <- tibble(
      Stage = stage,
      FeatureMode = feature_mode,
      K = ifelse(
        feature_mode == "ALL",
        NA_integer_,
        fixed_k
      ),
      MtryRule = "sqrt",
      MtryRuleOrder = match(
        "sqrt",
        RF_MTRY_RULES
      ),
      PooledRMSE = NA_real_,
      MeanFoldRMSE = NA_real_,
      SDFoldRMSE = NA_real_,
      ValidFolds = 0L,
      ExpectedFolds = expected_folds
    )
    
    warning(
      stage,
      ": all inner RF mtry evaluations failed; ",
      "the prespecified sqrt rule was used.",
      call. = FALSE
    )
  } else if (
    selected$ValidFolds[1] <
    expected_folds
  ) {
    warning(
      stage,
      ": mtry was selected from ",
      selected$ValidFolds[1],
      " of ",
      expected_folds,
      " valid inner folds.",
      call. = FALSE
    )
  }
  
  list(
    selected = selected,
    summary = summary,
    folds = folds
  )
}

fit_omics_pipeline <- function(
  meta_train,
  X_train_raw,
  X_test_raw,
  omics_type,
  mode = c("FixedK", "ALL"),
  inner_plan,
  seed,
  stage
) {
  mode <- match.arg(mode)

  algorithm <- if (omics_type == "metab") {
    METAB_ALGORITHM
  } else {
    MICRO_ALGORITHM
  }

  fixed_k <- if (omics_type == "metab") {
    K_METAB_PRIMARY
  } else {
    K_MICRO_PRIMARY
  }

  groups <- as.character(
    meta_train[[COL_ANIMAL]]
  )
  y <- as.numeric(
    meta_train[[COL_TARGET]]
  )

  selected_mtry_rule <- NA_character_
  selected_inner_rmse <- NA_real_
  selection_rule <- NA_character_
  inner_summary <- tibble()
  inner_folds <- tibble()

  if (algorithm == "RandomForest") {
    tuning <- select_rf_mtry_inner(
      meta_train = meta_train,
      X_train_raw = X_train_raw,
      omics_type = omics_type,
      inner_plan = inner_plan,
      feature_mode = mode,
      fixed_k = if (
        mode == "FixedK"
      ) fixed_k else NA_integer_,
      seed = seed + 1000L,
      stage = stage
    )

    selected_mtry_rule <-
      tuning$selected$MtryRule[1]

    selected_inner_rmse <-
      tuning$selected$PooledRMSE[1]

    inner_summary <- tuning$summary
    inner_folds <- tuning$folds
  }

  prepared <- prepare_omics(
    X_train_raw,
    X_test_raw,
    omics_type
  )

  if (mode == "ALL") {
    selected_features <- colnames(
      prepared$X_train
    )
  } else {
    ranking <- tryCatch(
      rank_features_processed(
        algorithm,
        prepared$X_train,
        y,
        groups,
        seed = seed + 2000L
      ),
      error = function(e) {
        if (algorithm != "RandomForest") {
          stop(e)
        }

        warning(
          stage,
          ": RF importance was unavailable; ",
          "absolute training-set correlation was used.",
          call. = FALSE
        )

        fallback_rf_ranking(
          prepared$X_train,
          y
        )
      }
    )

    selected_features <- select_top_features(
      ranking,
      fixed_k
    )
  }

  selected_mtry <- if (
    algorithm == "RandomForest"
  ) {
    rf_mtry_from_rule(
      selected_mtry_rule,
      length(selected_features)
    )
  } else {
    NA_integer_
  }

  if (algorithm == "ElasticNet") {
    selection_rule <-
      "Elastic Net lambda selected by grouped cv.glmnet"
  } else {
    selection_rule <-
      paste0(
        "RF mtry rule selected by repeated grouped inner CV: ",
        selected_mtry_rule
      )
  }

  final_fit <- if (algorithm == "ElasticNet") {
    fit_enet_processed(
      prepared$X_train[
        , selected_features, drop = FALSE
      ],
      y,
      prepared$X_test[
        , selected_features, drop = FALSE
      ],
      groups = groups,
      seed = seed + 3000L,
      alpha = ALPHA
    )
  } else {
    fit_rf_processed(
      prepared$X_train[
        , selected_features, drop = FALSE
      ],
      y,
      prepared$X_test[
        , selected_features, drop = FALSE
      ],
      mtry = selected_mtry,
      seed = seed + 3000L
    )
  }

  selection <- tibble(
    Stage = stage,
    OmicsType = omics_type,
    Algorithm = algorithm,
    Mode = mode,
    SelectedK = ifelse(
      mode == "FixedK",
      fixed_k,
      NA_integer_
    ),
    ActualK = length(selected_features),
    MtryRule = ifelse(
      algorithm == "RandomForest",
      selected_mtry_rule,
      NA_character_
    ),
    Mtry = ifelse(
      algorithm == "RandomForest",
      final_fit$mtry,
      NA_integer_
    ),
    Alpha = ifelse(
      algorithm == "ElasticNet",
      ALPHA,
      NA_real_
    ),
    Lambda1SE = ifelse(
      algorithm == "ElasticNet",
      final_fit$lambda_1se,
      NA_real_
    ),
    LambdaMin = ifelse(
      algorithm == "ElasticNet",
      final_fit$lambda_min,
      NA_real_
    ),
    SelectedInnerRMSE = selected_inner_rmse,
    SelectionRule = selection_rule
  )

  list(
    pred_train = final_fit$pred_train,
    pred_test = final_fit$pred_test,
    selected_features = selected_features,
    input_features = length(selected_features),
    effective_features = final_fit$n_features,
    selection = selection,
    inner_summary = inner_summary,
    inner_folds = inner_folds,
    prepared = prepared,
    final_fit = final_fit
  )
}

#  Fusion and stacking
fit_meta_nnls_matrix <- function(
  y,
  predictions
) {
  y <- as.numeric(y)
  P <- as.matrix(predictions)

  if (is.null(colnames(P))) {
    colnames(P) <- paste0(
      "Model",
      seq_len(ncol(P))
    )
  }

  ok <- is.finite(y) &
    apply(
      P,
      1,
      function(z) all(is.finite(z))
    )

  y <- y[ok]
  P <- P[ok, , drop = FALSE]

  y_bar <- mean(y)
  p_bar <- colMeans(P)

  yc <- y - y_bar
  Pc <- sweep(
    P,
    2,
    p_bar,
    "-"
  )

  objective <- function(w) {
    residuals <- yc -
      as.numeric(Pc %*% w)
    sum(residuals^2)
  }

  opt <- optim(
    par = rep(
      1 / ncol(P),
      ncol(P)
    ),
    fn = objective,
    method = "L-BFGS-B",
    lower = rep(0, ncol(P))
  )

  weights <- pmax(
    as.numeric(opt$par),
    0
  )
  names(weights) <- colnames(P)

  list(
    intercept =
      y_bar -
      sum(weights * p_bar),
    weights = weights
  )
}

predict_meta_nnls_matrix <- function(
  model,
  predictions
) {
  P <- as.matrix(predictions)

  if (is.null(dim(P))) {
    P <- matrix(P, ncol = 1L)
  }

  if (!is.null(colnames(P))) {
    P <- P[
      , names(model$weights), drop = FALSE
    ]
  }

  as.numeric(
    model$intercept +
      P %*% model$weights
  )
}

fit_early_fusion <- function(
  metab_fit,
  micro_fit,
  y_train,
  groups,
  seed
) {
  Xm_train <- metab_fit$prepared$X_train[
    , metab_fit$selected_features, drop = FALSE
  ]
  Xm_test <- metab_fit$prepared$X_test[
    , metab_fit$selected_features, drop = FALSE
  ]

  Xc_train <- micro_fit$prepared$X_train[
    , micro_fit$selected_features, drop = FALSE
  ]
  Xc_test <- micro_fit$prepared$X_test[
    , micro_fit$selected_features, drop = FALSE
  ]

  colnames(Xm_train) <- colnames(Xm_test) <-
    paste0("M_", colnames(Xm_train))

  colnames(Xc_train) <- colnames(Xc_test) <-
    paste0("C_", colnames(Xc_train))

  X_train <- cbind(Xm_train, Xc_train)
  X_test <- cbind(Xm_test, Xc_test)

  fit <- fit_enet_processed(
    X_train,
    y_train,
    X_test,
    groups = groups,
    seed = seed,
    alpha = ALPHA
  )

  list(
    pred_train = fit$pred_train,
    pred_test = fit$pred_test,
    input_features = ncol(X_train),
    effective_features = fit$n_features,
    lambda_1se = fit$lambda_1se,
    lambda_min = fit$lambda_min
  )
}

fit_stacking <- function(
  meta_train,
  X_metab_train_raw,
  X_micro_train_raw,
  X_metab_test_raw,
  X_micro_test_raw,
  mode,
  final_metab_fit,
  final_micro_fit,
  seed,
  inner_repeats = INNER_REPEATS
) {
  groups <- as.character(
    meta_train[[COL_ANIMAL]]
  )
  y <- as.numeric(
    meta_train[[COL_TARGET]]
  )

  foldid <- make_group_folds(
    groups,
    y,
    kfold = STACK_FOLDS,
    seed = seed
  )

  oof_metab <- rep(NA_real_, length(y))
  oof_micro <- rep(NA_real_, length(y))

  for (f in sort(unique(foldid))) {
    is_val <- foldid == f
    is_train <- !is_val

    sub_meta <- meta_train[
      is_train, , drop = FALSE
    ]

    sub_groups <- as.character(
      sub_meta[[COL_ANIMAL]]
    )
    sub_y <- as.numeric(
      sub_meta[[COL_TARGET]]
    )

    sub_inner_plan <- make_repeated_group_folds(
      sub_groups,
      sub_y,
      kfold = INNER_FOLDS,
      repeats = as.integer(inner_repeats),
      seed = seed + f * 100L
    )

    fit_m <- fit_omics_pipeline(
      meta_train = sub_meta,
      X_train_raw = X_metab_train_raw[
        is_train, , drop = FALSE
      ],
      X_test_raw = X_metab_train_raw[
        is_val, , drop = FALSE
      ],
      omics_type = "metab",
      mode = mode,
      inner_plan = sub_inner_plan,
      seed = seed + f * 10000L + 1L,
      stage = paste0(
        "StackOOF_Metab_",
        mode
      )
    )

    fit_c <- fit_omics_pipeline(
      meta_train = sub_meta,
      X_train_raw = X_micro_train_raw[
        is_train, , drop = FALSE
      ],
      X_test_raw = X_micro_train_raw[
        is_val, , drop = FALSE
      ],
      omics_type = "micro",
      mode = mode,
      inner_plan = sub_inner_plan,
      seed = seed + f * 10000L + 2L,
      stage = paste0(
        "StackOOF_Micro_",
        mode
      )
    )

    oof_metab[is_val] <- fit_m$pred_test
    oof_micro[is_val] <- fit_c$pred_test
  }

  if (
    any(!is.finite(oof_metab)) ||
      any(!is.finite(oof_micro))
  ) {
    stop(
      "Incomplete OOF predictions were generated ",
      "for within-organ stacking."
    )
  }

  oof_matrix <- cbind(
    Metabolomics = oof_metab,
    Microbiome = oof_micro
  )

  meta_model <- fit_meta_nnls_matrix(
    y,
    oof_matrix
  )

  final_train_matrix <- cbind(
    Metabolomics =
      final_metab_fit$pred_train,
    Microbiome =
      final_micro_fit$pred_train
  )

  final_test_matrix <- cbind(
    Metabolomics =
      final_metab_fit$pred_test,
    Microbiome =
      final_micro_fit$pred_test
  )

  list(
    pred_train = predict_meta_nnls_matrix(
      meta_model,
      final_train_matrix
    ),
    pred_test = predict_meta_nnls_matrix(
      meta_model,
      final_test_matrix
    ),
    input_features =
      final_metab_fit$input_features +
        final_micro_fit$input_features,
    effective_features =
      final_metab_fit$effective_features +
        final_micro_fit$effective_features,
    meta_model = meta_model,
    oof_metab = oof_metab,
    oof_micro = oof_micro
  )
}

fit_single_organ_fixedk_stack <- function(
  meta_train,
  X_metab_train_raw,
  X_micro_train_raw,
  X_metab_test_raw,
  X_micro_test_raw,
  tissue_label,
  seed,
  inner_repeats = INNER_REPEATS
) {
  groups <- as.character(
    meta_train[[COL_ANIMAL]]
  )
  y <- as.numeric(
    meta_train[[COL_TARGET]]
  )

  inner_plan <- make_repeated_group_folds(
    groups,
    y,
    kfold = INNER_FOLDS,
    repeats = as.integer(inner_repeats),
    seed = seed + 11L
  )

  metab_fit <- fit_omics_pipeline(
    meta_train = meta_train,
    X_train_raw = X_metab_train_raw,
    X_test_raw = X_metab_test_raw,
    omics_type = "metab",
    mode = "FixedK",
    inner_plan = inner_plan,
    seed = seed + 10000L,
    stage = paste0(
      tissue_label,
      "_Metab_FixedK"
    )
  )

  micro_fit <- fit_omics_pipeline(
    meta_train = meta_train,
    X_train_raw = X_micro_train_raw,
    X_test_raw = X_micro_test_raw,
    omics_type = "micro",
    mode = "FixedK",
    inner_plan = inner_plan,
    seed = seed + 20000L,
    stage = paste0(
      tissue_label,
      "_Micro_FixedK"
    )
  )

  stack_fit <- fit_stacking(
    meta_train = meta_train,
    X_metab_train_raw =
      X_metab_train_raw,
    X_micro_train_raw =
      X_micro_train_raw,
    X_metab_test_raw =
      X_metab_test_raw,
    X_micro_test_raw =
      X_micro_test_raw,
    mode = "FixedK",
    final_metab_fit = metab_fit,
    final_micro_fit = micro_fit,
    seed = seed + 40000L,
    inner_repeats = inner_repeats
  )

  list(
    pred_train = stack_fit$pred_train,
    pred_test = stack_fit$pred_test,
    stack_fit = stack_fit,
    metab_fit = metab_fit,
    micro_fit = micro_fit
  )
}

# Tissue-specific outer LOAO
prediction_rows <- function(
  meta_test,
  predictions,
  model,
  label,
  outer_fold,
  analysis = "Tissue-specific"
) {
  tibble(
    Analysis = analysis,
    Label = label,
    OuterFold = outer_fold,
    Model = model,
    AnimalID = as.character(
      meta_test[[COL_ANIMAL]]
    ),
    Observed = as.numeric(
      meta_test[[COL_TARGET]]
    ),
    Predicted = as.numeric(predictions)
  ) %>%
    group_by(
      Analysis,
      Label,
      OuterFold,
      Model,
      AnimalID
    ) %>%
    summarise(
      Observed = first(Observed),
      Predicted = mean(
        Predicted,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Residual = Observed - Predicted,
      AbsoluteError = abs(Residual),
      SquaredError = Residual^2
    )
}

feature_rows <- function(
  fit,
  label,
  outer_fold,
  omics_label
) {
  tibble(
    Label = label,
    OuterFold = outer_fold,
    Omics = omics_label,
    Mode = "FixedK",
    Algorithm =
      fit$selection$Algorithm[1],
    SelectedK =
      fit$selection$SelectedK[1],
    ActualK = fit$input_features,
    Feature = fit$selected_features,
    Rank = seq_along(
      fit$selected_features
    )
  )
}

run_one_outer_loao <- function(
  meta,
  X_metab,
  X_micro,
  split_info,
  label
) {
  outer_fold <- split_info$OuterFold
  seed <- model_partition_seed(
    "SingleOrganFixedKStack",
    label,
    split_info$TrainIDs,
    split_info$TestIDs
  )

  ids <- as.character(
    meta[[COL_ANIMAL]]
  )

  is_train <- ids %in%
    split_info$TrainIDs
  is_test <- ids %in%
    split_info$TestIDs

  if (!any(is_train) || !any(is_test)) {
    stop(
      "Empty training or test partition for ",
      label,
      ", outer fold ",
      outer_fold,
      "."
    )
  }

  meta_train <- meta[
    is_train, , drop = FALSE
  ]
  meta_test <- meta[
    is_test, , drop = FALSE
  ]

  y_train <- as.numeric(
    meta_train[[COL_TARGET]]
  )
  groups_train <- as.character(
    meta_train[[COL_ANIMAL]]
  )

  inner_plan <- make_repeated_group_folds(
    groups_train,
    y_train,
    kfold = INNER_FOLDS,
    repeats = INNER_REPEATS,
    seed = seed + 11L
  )

  metab_fixed <- fit_omics_pipeline(
    meta_train,
    X_metab[is_train, , drop = FALSE],
    X_metab[is_test, , drop = FALSE],
    "metab",
    "FixedK",
    inner_plan,
    seed + 10000L,
    "Metab_FixedK"
  )

  micro_fixed <- fit_omics_pipeline(
    meta_train,
    X_micro[is_train, , drop = FALSE],
    X_micro[is_test, , drop = FALSE],
    "micro",
    "FixedK",
    inner_plan,
    seed + 20000L,
    "Micro_FixedK"
  )

  ef_fixed <- fit_early_fusion(
    metab_fixed,
    micro_fixed,
    y_train,
    groups_train,
    seed + 30000L
  )

  stacking_fixed <- fit_stacking(
    meta_train,
    X_metab[is_train, , drop = FALSE],
    X_micro[is_train, , drop = FALSE],
    X_metab[is_test, , drop = FALSE],
    X_micro[is_test, , drop = FALSE],
    "FixedK",
    metab_fixed,
    micro_fixed,
    seed + 40000L
  )

  metab_all <- fit_omics_pipeline(
    meta_train,
    X_metab[is_train, , drop = FALSE],
    X_metab[is_test, , drop = FALSE],
    "metab",
    "ALL",
    inner_plan,
    seed + 50000L,
    "Metab_ALL"
  )

  micro_all <- fit_omics_pipeline(
    meta_train,
    X_micro[is_train, , drop = FALSE],
    X_micro[is_test, , drop = FALSE],
    "micro",
    "ALL",
    inner_plan,
    seed + 60000L,
    "Micro_ALL"
  )

  ef_all <- fit_early_fusion(
    metab_all,
    micro_all,
    y_train,
    groups_train,
    seed + 70000L
  )

  stacking_all <- fit_stacking(
    meta_train,
    X_metab[is_train, , drop = FALSE],
    X_micro[is_train, , drop = FALSE],
    X_metab[is_test, , drop = FALSE],
    X_micro[is_test, , drop = FALSE],
    "ALL",
    metab_all,
    micro_all,
    seed + 80000L
  )

  predictions <- bind_rows(
    prediction_rows(
      meta_test,
      metab_fixed$pred_test,
      "Metab_FixedK",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      micro_fixed$pred_test,
      "Micro_FixedK",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      ef_fixed$pred_test,
      "EF_FixedK",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      stacking_fixed$pred_test,
      "Stacking_FixedK",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      metab_all$pred_test,
      "Metab_ALL",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      micro_all$pred_test,
      "Micro_ALL",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      ef_all$pred_test,
      "EF_ALL",
      label,
      outer_fold
    ),
    prediction_rows(
      meta_test,
      stacking_all$pred_test,
      "Stacking_ALL",
      label,
      outer_fold
    )
  )

  selections <- bind_rows(
    metab_fixed$selection,
    micro_fixed$selection,
    metab_all$selection,
    micro_all$selection,
    tibble(
      Stage = c(
        "EF_FixedK",
        "EF_ALL",
        "Stacking_FixedK",
        "Stacking_ALL"
      ),
      OmicsType = c(
        "fusion",
        "fusion",
        "stacking",
        "stacking"
      ),
      Algorithm = c(
        EF_ALGORITHM,
        EF_ALGORITHM,
        "NNLS",
        "NNLS"
      ),
      Mode = c(
        "FixedK",
        "ALL",
        "FixedK",
        "ALL"
      ),
      SelectedK = c(
        K_METAB_PRIMARY +
          K_MICRO_PRIMARY,
        NA_integer_,
        K_METAB_PRIMARY +
          K_MICRO_PRIMARY,
        NA_integer_
      ),
      ActualK = c(
        ef_fixed$input_features,
        ef_all$input_features,
        stacking_fixed$input_features,
        stacking_all$input_features
      ),
      MtryRule = NA_character_,
      Mtry = NA_integer_,
      Alpha = c(
        ALPHA,
        ALPHA,
        NA_real_,
        NA_real_
      ),
      Lambda1SE = c(
        ef_fixed$lambda_1se,
        ef_all$lambda_1se,
        NA_real_,
        NA_real_
      ),
      LambdaMin = c(
        ef_fixed$lambda_min,
        ef_all$lambda_min,
        NA_real_,
        NA_real_
      ),
      SelectedInnerRMSE = NA_real_,
      SelectionRule = c(
        "Elastic Net lambda selected by grouped CV",
        "Elastic Net lambda selected by grouped CV",
        "Training-only OOF NNLS stacking",
        "Training-only OOF NNLS stacking"
      )
    )
  ) %>%
    mutate(
      Label = label,
      OuterFold = outer_fold,
      TestAnimal = as.character(
        meta_test[[COL_ANIMAL]][1]
      ),
      .before = 1
    )

  stack_weights <- bind_rows(
    tibble(
      Label = label,
      OuterFold = outer_fold,
      TestAnimal = as.character(
        meta_test[[COL_ANIMAL]][1]
      ),
      Mode = "FixedK",
      Intercept =
        stacking_fixed$meta_model$intercept,
      WeightMetab =
        stacking_fixed$meta_model$weights[
          "Metabolomics"
        ],
      WeightMicro =
        stacking_fixed$meta_model$weights[
          "Microbiome"
        ]
    ),
    tibble(
      Label = label,
      OuterFold = outer_fold,
      TestAnimal = as.character(
        meta_test[[COL_ANIMAL]][1]
      ),
      Mode = "ALL",
      Intercept =
        stacking_all$meta_model$intercept,
      WeightMetab =
        stacking_all$meta_model$weights[
          "Metabolomics"
        ],
      WeightMicro =
        stacking_all$meta_model$weights[
          "Microbiome"
        ]
    )
  )

  features <- bind_rows(
    feature_rows(
      metab_fixed,
      label,
      outer_fold,
      "Metabolomics"
    ),
    feature_rows(
      micro_fixed,
      label,
      outer_fold,
      "Microbiome"
    )
  )

  inner_summary <- bind_rows(
    metab_fixed$inner_summary,
    micro_fixed$inner_summary,
    metab_all$inner_summary,
    micro_all$inner_summary
  ) %>%
    mutate(
      Label = label,
      OuterFold = outer_fold,
      .before = 1
    )

  inner_folds <- bind_rows(
    metab_fixed$inner_folds,
    micro_fixed$inner_folds,
    metab_all$inner_folds,
    micro_all$inner_folds
  ) %>%
    mutate(
      Label = label,
      OuterFold = outer_fold,
      .before = 1
    )

  list(
    predictions = predictions,
    selections = selections,
    stack_weights = stack_weights,
    features = features,
    inner_summary = inner_summary,
    inner_folds = inner_folds
  )
}

run_loao_analysis <- function(
  meta,
  X_metab,
  X_micro,
  label
) {
  n_animals <- n_distinct(
    meta[[COL_ANIMAL]]
  )

  

  outer_splits <- make_outer_splits(
    meta
  )

  cat(
    sprintf(
      "  %-12s animals=%3d rows=%4d  ",
      label,
      n_animals,
      nrow(meta)
    )
  )

  outputs <- vector(
    "list",
    length(outer_splits)
  )

  for (i in seq_along(outer_splits)) {
    if (i %% 5L == 0L) cat(".")

    outputs[[i]] <- run_one_outer_loao(
      meta,
      X_metab,
      X_micro,
      outer_splits[[i]],
      label
    )
  }

  cat(" done\n")

  list(
    predictions = bind_rows(
      lapply(outputs, `[[`, "predictions")
    ),
    selections = bind_rows(
      lapply(outputs, `[[`, "selections")
    ),
    stack_weights = bind_rows(
      lapply(outputs, `[[`, "stack_weights")
    ),
    features = bind_rows(
      lapply(outputs, `[[`, "features")
    ),
    inner_summary = bind_rows(
      lapply(outputs, `[[`, "inner_summary")
    ),
    inner_folds = bind_rows(
      lapply(outputs, `[[`, "inner_folds")
    )
  )
}


# Performance summaries and comparisons
bootstrap_metric_ci <- function(
  observed,
  predicted,
  B = BOOT_B,
  seed = 42L
) {
  ok <- is.finite(observed) &
    is.finite(predicted)

  observed <- observed[ok]
  predicted <- predicted[ok]
  n <- length(observed)

  set.seed(seed)

  idx <- replicate(
    B,
    sample.int(
      n,
      n,
      replace = TRUE
    )
  )

  rmse_boot <- apply(
    idx,
    2,
    function(ii) {
      calc_rmse(
        observed[ii],
        predicted[ii]
      )
    }
  )

  mae_boot <- apply(
    idx,
    2,
    function(ii) {
      calc_mae(
        observed[ii],
        predicted[ii]
      )
    }
  )

  r2_boot <- apply(
    idx,
    2,
    function(ii) {
      calc_r2(
        observed[ii],
        predicted[ii]
      )
    }
  )

  tibble(
    RMSE_lower = quantile(
      rmse_boot,
      0.025,
      na.rm = TRUE
    ),
    RMSE_upper = quantile(
      rmse_boot,
      0.975,
      na.rm = TRUE
    ),
    MAE_lower = quantile(
      mae_boot,
      0.025,
      na.rm = TRUE
    ),
    MAE_upper = quantile(
      mae_boot,
      0.975,
      na.rm = TRUE
    ),
    R2_lower = quantile(
      r2_boot,
      0.025,
      na.rm = TRUE
    ),
    R2_upper = quantile(
      r2_boot,
      0.975,
      na.rm = TRUE
    )
  )
}

summarise_predictions <- function(
  predictions,
  group_columns
) {
  predictions %>%
    group_by(
      across(all_of(group_columns))
    ) %>%
    group_modify(function(df, key) {
      point <- tibble(
        NAnimals = n_distinct(
          df$AnimalID
        ),
        RMSE = calc_rmse(
          df$Observed,
          df$Predicted
        ),
        MAE = calc_mae(
          df$Observed,
          df$Predicted
        ),
        R2 = calc_r2(
          df$Observed,
          df$Predicted
        ),
        MeanError = mean(
          df$Residual,
          na.rm = TRUE
        ),
        MedianAbsoluteError = median(
          df$AbsoluteError,
          na.rm = TRUE
        )
      )

      key_text <- paste(
        unlist(key),
        collapse = "__"
      )

      ci <- bootstrap_metric_ci(
        df$Observed,
        df$Predicted,
        B = BOOT_B,
        seed = stable_seed(
          "MetricCI",
          key_text
        )
      )

      bind_cols(point, ci)
    }) %>%
    ungroup()
}

paired_prediction_comparison <- function(
  candidate_df,
  reference_df,
  comparison_label,
  B = BOOT_B,
  seed = 42L
) {
  paired <- candidate_df %>%
    select(
      AnimalID,
      Observed,
      PredCandidate = Predicted,
      AECandidate = AbsoluteError
    ) %>%
    inner_join(
      reference_df %>%
        select(
          AnimalID,
          PredReference = Predicted,
          AEReference = AbsoluteError
        ),
      by = "AnimalID"
    )

  n <- nrow(paired)

  point_delta <-
    calc_rmse(
      paired$Observed,
      paired$PredReference
    ) -
    calc_rmse(
      paired$Observed,
      paired$PredCandidate
    )

  set.seed(seed)

  boot_delta <- replicate(B, {
    ii <- sample.int(
      n,
      n,
      replace = TRUE
    )

    calc_rmse(
      paired$Observed[ii],
      paired$PredReference[ii]
    ) -
      calc_rmse(
        paired$Observed[ii],
        paired$PredCandidate[ii]
      )
  })

  perm_delta <- replicate(B, {
    swap <- sample(
      c(FALSE, TRUE),
      n,
      replace = TRUE
    )

    candidate_perm <- ifelse(
      swap,
      paired$PredReference,
      paired$PredCandidate
    )

    reference_perm <- ifelse(
      swap,
      paired$PredCandidate,
      paired$PredReference
    )

    calc_rmse(
      paired$Observed,
      reference_perm
    ) -
      calc_rmse(
        paired$Observed,
        candidate_perm
      )
  })

  p_perm <- (
    1 +
      sum(
        abs(perm_delta) >=
          abs(point_delta)
      )
  ) / (B + 1)

  tibble(
    Comparison = comparison_label,
    NAnimals = n,
    DeltaRMSE = point_delta,
    DeltaRMSELower = quantile(
      boot_delta,
      0.025,
      na.rm = TRUE
    ),
    DeltaRMSEUpper = quantile(
      boot_delta,
      0.975,
      na.rm = TRUE
    ),
    WinRate = mean(
      paired$AECandidate <
        paired$AEReference
    ),
    PPermutation = p_perm
  )
}

# Shared complete-cohort pairwise cross-organ analysis
build_complete_organ_dataset <- function(
  meta,
  X_metab,
  X_micro,
  tissues = TISSUES
) {
  id_sets <- lapply(
    tissues,
    function(tissue) {
      as.character(
        meta[[COL_ANIMAL]][
          as.character(meta[[COL_TISSUE]]) ==
            tissue
        ]
      )
    }
  )

  common_ids <- sort(
    Reduce(intersect, id_sets)
  )


  row_index <- list()
  metab_list <- list()
  micro_list <- list()

  for (tissue in tissues) {
    idx <- which(
      as.character(meta[[COL_TISSUE]]) ==
        tissue
    )

    tissue_ids <- as.character(
      meta[[COL_ANIMAL]][idx]
    )

    rows <- idx[
      match(common_ids, tissue_ids)
    ]

    if (anyNA(rows)) {
      stop(
        "Five-organ alignment failed for ",
        tissue,
        "."
      )
    }

    row_index[[tissue]] <- rows
    metab_list[[tissue]] <- X_metab[
      rows, , drop = FALSE
    ]
    micro_list[[tissue]] <- X_micro[
      rows, , drop = FALSE
    ]
  }

  y_matrix <- do.call(
    cbind,
    lapply(
      tissues,
      function(tissue) {
        as.numeric(
          meta[[COL_TARGET]][
            row_index[[tissue]]
          ]
        )
      }
    )
  )

  meta_complete <- meta[
    row_index[[tissues[1]]],
    , drop = FALSE
  ]

  meta_complete[[COL_ID]] <- paste0(
    common_ids,
    "__FiveOrgan"
  )
  meta_complete[[COL_TISSUE]] <-
    "FiveOrganComplete"

  list(
    meta = meta_complete,
    metab = metab_list,
    micro = micro_list,
    tissues = tissues,
    n_animals = length(common_ids)
  )
}

organ_set_label <- function(organs) {
  paste(organs, collapse = "+")
}

parse_organ_set <- function(label) {
  strsplit(
    as.character(label),
    "+",
    fixed = TRUE
  )[[1]]
}

make_pair_and_five_organ_subsets <- function(
  tissues = TISSUES
) {
  single_subsets <- lapply(
    tissues,
    function(tissue) tissue
  )

  pair_subsets <- combn(
    tissues,
    2L,
    simplify = FALSE
  )

  five_organ_subset <- list(tissues)

  c(
    single_subsets,
    pair_subsets,
    five_organ_subset
  )
}

run_one_shared_complete_outer <- function(
  complete_data,
  split_info
) {
  tissues <- complete_data$tissues
  subsets <- make_pair_and_five_organ_subsets(
    tissues
  )
  outer_fold <- split_info$OuterFold

  ids <- as.character(
    complete_data$meta[[COL_ANIMAL]]
  )

  is_train <- ids %in%
    split_info$TrainIDs
  is_test <- ids %in%
    split_info$TestIDs

  meta_train <- complete_data$meta[
    is_train, , drop = FALSE
  ]
  meta_test <- complete_data$meta[
    is_test, , drop = FALSE
  ]

  groups <- as.character(
    meta_train[[COL_ANIMAL]]
  )
  y <- as.numeric(
    meta_train[[COL_TARGET]]
  )

  train_rows <- which(is_train)

  final_fits <- setNames(
    vector("list", length(tissues)),
    tissues
  )

  test_organ_predictions <- setNames(
    rep(NA_real_, length(tissues)),
    tissues
  )

  for (tissue in tissues) {
    final_fits[[tissue]] <-
      fit_single_organ_fixedk_stack(
        meta_train,
        complete_data$metab[[tissue]][
          is_train, , drop = FALSE
        ],
        complete_data$micro[[tissue]][
          is_train, , drop = FALSE
        ],
        complete_data$metab[[tissue]][
          is_test, , drop = FALSE
        ],
        complete_data$micro[[tissue]][
          is_test, , drop = FALSE
        ],
        tissue,
        model_partition_seed(
          "SingleOrganFixedKStack",
          tissue,
          split_info$TrainIDs,
          split_info$TestIDs
        )
      )

    test_organ_predictions[tissue] <-
      mean(
        final_fits[[tissue]]$pred_test,
        na.rm = TRUE
      )
  }

   meta_foldid <- make_group_folds(
    groups,
    y,
    kfold = STACK_FOLDS,
    seed = stable_seed(
      "SharedOrganLevelOOF",
      id_set_key(split_info$TrainIDs),
      id_set_key(split_info$TestIDs)
    )
  )

  oof_matrix <- matrix(
    NA_real_,
    nrow = length(y),
    ncol = length(tissues),
    dimnames = list(NULL, tissues)
  )

  for (f in sort(unique(meta_foldid))) {
    is_meta_val <- meta_foldid == f
    is_meta_train <- !is_meta_val

    sub_meta <- meta_train[
      is_meta_train, , drop = FALSE
    ]

    for (tissue in tissues) {
      sub_fit <- fit_single_organ_fixedk_stack(
        sub_meta,
        complete_data$metab[[tissue]][
          train_rows[is_meta_train],
          , drop = FALSE
        ],
        complete_data$micro[[tissue]][
          train_rows[is_meta_train],
          , drop = FALSE
        ],
        complete_data$metab[[tissue]][
          train_rows[is_meta_val],
          , drop = FALSE
        ],
        complete_data$micro[[tissue]][
          train_rows[is_meta_val],
          , drop = FALSE
        ],
        tissue,
        model_partition_seed(
          "SingleOrganFixedKStack",
          tissue,
          groups[is_meta_train],
          groups[is_meta_val]
        ),
        inner_repeats =
          MULTI_ORGAN_OOF_INNER_REPEATS
      )

      oof_matrix[
        is_meta_val,
        tissue
      ] <- sub_fit$pred_test
    }
  }


  prediction_list <- list()
  weight_list <- list()
  pred_pos <- 1L
  weight_pos <- 1L

  for (organs in subsets) {
    n_organs <- length(organs)
    organ_set <- organ_set_label(organs)

    analysis_label <- if (
      n_organs == 1L
    ) {
      "Complete-cohort single-organ reference"
    } else if (
      n_organs == 2L
    ) {
      "Exploratory pairwise cross-organ integration"
    } else {
      "Supplementary five-organ sensitivity analysis"
    }

    if (n_organs == 1L) {
      tissue <- organs[1]
      one_prediction <- test_organ_predictions[tissue]

      prediction_list[[pred_pos]] <-
        prediction_rows(
          meta_test,
          one_prediction,
          paste0(
            tissue,
            "_Stacking_FixedK"
          ),
          organ_set,
          outer_fold,
          analysis = analysis_label
        ) %>%
        mutate(
          NOrgans = n_organs,
          OrganSet = organ_set,
          Method = "SingleOrgan",
          .before = 1
        )

      pred_pos <- pred_pos + 1L

      weight_list[[weight_pos]] <- tibble(
        OuterFold = outer_fold,
        TestAnimal = as.character(
          meta_test[[COL_ANIMAL]][1]
        ),
        NOrgans = n_organs,
        OrganSet = organ_set,
        Method = "SingleOrgan",
        Intercept = 0,
        Tissue = tissue,
        Weight = 1
      )

      weight_pos <- weight_pos + 1L
      next
    }

    equal_prediction <- mean(
      test_organ_predictions[organs]
    )

    prediction_list[[pred_pos]] <-
      prediction_rows(
        meta_test,
        equal_prediction,
        paste0(
          "EqualWeight_",
          n_organs,
          "Organs"
        ),
        organ_set,
        outer_fold,
        analysis = analysis_label
      ) %>%
      mutate(
        NOrgans = n_organs,
        OrganSet = organ_set,
        Method = "EqualWeight",
        .before = 1
      )

    pred_pos <- pred_pos + 1L

    weight_list[[weight_pos]] <- tibble(
      OuterFold = outer_fold,
      TestAnimal = as.character(
        meta_test[[COL_ANIMAL]][1]
      ),
      NOrgans = n_organs,
      OrganSet = organ_set,
      Method = "EqualWeight",
      Intercept = 0,
      Tissue = organs,
      Weight = rep(
        1 / n_organs,
        n_organs
      )
    )

    weight_pos <- weight_pos + 1L

    P_train <- oof_matrix[
      , organs, drop = FALSE
    ]

    meta_model <- fit_meta_nnls_matrix(
      y,
      P_train
    )

    P_test <- matrix(
      test_organ_predictions[organs],
      nrow = 1L,
      dimnames = list(NULL, organs)
    )

    hierarchical_prediction <-
      predict_meta_nnls_matrix(
        meta_model,
        P_test
      )

    prediction_list[[pred_pos]] <-
      prediction_rows(
        meta_test,
        hierarchical_prediction,
        paste0(
          "HierarchicalNNLS_",
          n_organs,
          "Organs"
        ),
        organ_set,
        outer_fold,
        analysis = analysis_label
      ) %>%
      mutate(
        NOrgans = n_organs,
        OrganSet = organ_set,
        Method = "HierarchicalNNLS",
        .before = 1
      )

    pred_pos <- pred_pos + 1L

    weight_list[[weight_pos]] <- tibble(
      OuterFold = outer_fold,
      TestAnimal = as.character(
        meta_test[[COL_ANIMAL]][1]
      ),
      NOrgans = n_organs,
      OrganSet = organ_set,
      Method = "HierarchicalNNLS",
      Intercept = meta_model$intercept,
      Tissue = names(
        meta_model$weights
      ),
      Weight = as.numeric(
        meta_model$weights
      )
    )

    weight_pos <- weight_pos + 1L
  }

  oof_predictions <- as_tibble(
    oof_matrix,
    .name_repair = "minimal"
  ) %>%
    mutate(
      MetaAnimal = groups,
      Observed = y,
      .before = 1
    ) %>%
    pivot_longer(
      cols = all_of(tissues),
      names_to = "Tissue",
      values_to = "OOFPredicted"
    ) %>%
    mutate(
      OuterFold = outer_fold,
      OuterTestAnimal = as.character(
        meta_test[[COL_ANIMAL]][1]
      ),
      .before = 1
    )

  list(
    predictions = bind_rows(prediction_list),
    weights = bind_rows(weight_list),
    oof_predictions = oof_predictions
  )
}

run_shared_complete_multi_organ_analysis <- function(
  complete_data
) {

  cat(
    sprintf(
      "  Shared pairwise cohort animals=%3d  ",
      complete_data$n_animals
    )
  )

  outer_splits <- make_outer_splits(
    complete_data$meta
  )

  outputs <- vector(
    "list",
    length(outer_splits)
  )

  for (i in seq_along(outer_splits)) {
    if (i %% 5L == 0L) cat(".")

    outputs[[i]] <-
      run_one_shared_complete_outer(
        complete_data,
        outer_splits[[i]]
      )
  }

  cat(" done\n")

  list(
    predictions = bind_rows(
      lapply(outputs, `[[`, "predictions")
    ),
    weights = bind_rows(
      lapply(outputs, `[[`, "weights")
    ),
    oof_predictions = bind_rows(
      lapply(outputs, `[[`, "oof_predictions")
    )
  )
}


# Shared multi-organ summaries and comparisons
summarise_shared_multi_organ <- function(
  predictions_multi
) {
  summary <- summarise_predictions(
    predictions_multi,
    c(
      "NOrgans",
      "OrganSet",
      "Method",
      "Model"
    )
  )

  candidate_rows <- summary %>%
    filter(NOrgans >= 2L)

  descriptive_list <- list()
  desc_pos <- 1L

  for (i in seq_len(nrow(candidate_rows))) {
    candidate <- candidate_rows[i, ]
    organs <- parse_organ_set(
      candidate$OrganSet
    )

    constituent_summary <- summary %>%
      filter(
        Method == "SingleOrgan",
        OrganSet %in% organs
      ) %>%
      arrange(RMSE, OrganSet)

    best_constituent <- constituent_summary %>%
      slice(1)

    descriptive_list[[desc_pos]] <- tibble(
      NOrgans = candidate$NOrgans,
      OrganSet = candidate$OrganSet,
      Method = candidate$Method,
      Model = candidate$Model,
      CandidateRMSE = candidate$RMSE,
      CandidateRMSELower =
        candidate$RMSE_lower,
      CandidateRMSEUpper =
        candidate$RMSE_upper,
      BestConstituent =
        best_constituent$OrganSet,
      BestConstituentRMSE =
        best_constituent$RMSE,
      DeltaRMSE_vs_BestConstituent =
        best_constituent$RMSE -
          candidate$RMSE,
      LowerRMSE_than_best_constituent_descriptive =
        candidate$RMSE <
          best_constituent$RMSE
    )

    desc_pos <- desc_pos + 1L
  }

  comparison_list <- list()
  comp_pos <- 1L

  candidate_keys <- predictions_multi %>%
    filter(NOrgans >= 2L) %>%
    distinct(
      NOrgans,
      OrganSet,
      Method,
      Model
    )

  for (i in seq_len(nrow(candidate_keys))) {
    key <- candidate_keys[i, ]
    organs <- parse_organ_set(
      key$OrganSet
    )

    candidate_df <- predictions_multi %>%
      filter(
        OrganSet == key$OrganSet,
        Method == key$Method
      )

    for (reference_tissue in organs) {
      reference_df <- predictions_multi %>%
        filter(
          NOrgans == 1L,
          Method == "SingleOrgan",
          OrganSet == reference_tissue
        )

      comparison_list[[comp_pos]] <-
        paired_prediction_comparison(
          candidate_df,
          reference_df,
          comparison_label = paste0(
            key$OrganSet,
            " ",
            key$Method,
            " vs ",
            reference_tissue
          ),
          B = MULTI_ORGAN_COMPARISON_B,
          seed = stable_seed(
            "SharedMultiOrganComparison",
            key$OrganSet,
            key$Method,
            reference_tissue
          )
        ) %>%
        mutate(
          NOrgans = key$NOrgans,
          OrganSet = key$OrganSet,
          CandidateMethod = key$Method,
          ReferenceType =
            "Constituent single organ",
          ReferenceModel =
            reference_tissue,
          .before = 1
        )

      comp_pos <- comp_pos + 1L
    }
  }

  learned_keys <- predictions_multi %>%
    filter(
      NOrgans >= 2L,
      Method == "HierarchicalNNLS"
    ) %>%
    distinct(NOrgans, OrganSet)

  for (i in seq_len(nrow(learned_keys))) {
    key <- learned_keys[i, ]

    comparison_list[[comp_pos]] <-
      paired_prediction_comparison(
        predictions_multi %>%
          filter(
            OrganSet == key$OrganSet,
            Method == "HierarchicalNNLS"
          ),
        predictions_multi %>%
          filter(
            OrganSet == key$OrganSet,
            Method == "EqualWeight"
          ),
        comparison_label = paste0(
          key$OrganSet,
          ": hierarchical NNLS vs equal weight"
        ),
        B = MULTI_ORGAN_COMPARISON_B,
        seed = stable_seed(
          "HierarchicalVsEqual",
          key$OrganSet
        )
      ) %>%
      mutate(
        NOrgans = key$NOrgans,
        OrganSet = key$OrganSet,
        CandidateMethod =
          "HierarchicalNNLS",
        ReferenceType =
          "Equal-weight benchmark",
        ReferenceModel =
          "EqualWeight",
        .before = 1
      )

    comp_pos <- comp_pos + 1L
  }

  comparisons <- bind_rows(
    comparison_list
  ) %>%
    group_by(
      CandidateMethod,
      ReferenceType,
      NOrgans
    ) %>%
    mutate(
      PHolmWithinMethodSize = p.adjust(
        PPermutation,
        method = "holm"
      )
    ) %>%
    ungroup()

  list(
    summary = summary,
    descriptive = bind_rows(
      descriptive_list
    ),
    comparisons = comparisons
  )
}

build_add_one_organ_deltas <- function(
  summary_multi
) {
  to_rows <- summary_multi %>%
    filter(NOrgans >= 2L)

  rows <- list()
  pos <- 1L

  for (i in seq_len(nrow(to_rows))) {
    to_row <- to_rows[i, ]
    to_organs <- parse_organ_set(
      to_row$OrganSet
    )

    for (added_organ in to_organs) {
      from_organs <- setdiff(
        to_organs,
        added_organ
      )
      from_set <- organ_set_label(
        from_organs
      )
      from_method <- if (
        length(from_organs) == 1L
      ) {
        "SingleOrgan"
      } else {
        to_row$Method
      }

      from_row <- summary_multi %>%
        filter(
          OrganSet == from_set,
          Method == from_method
        )

      if (nrow(from_row) != 1L) {
        next
      }

      rows[[pos]] <- tibble(
        FromNOrgans = length(from_organs),
        ToNOrgans = length(to_organs),
        FromOrganSet = from_set,
        AddedOrgan = added_organ,
        ToOrganSet = to_row$OrganSet,
        FromMethod = from_method,
        ToMethod = to_row$Method,
        FromRMSE = from_row$RMSE,
        ToRMSE = to_row$RMSE,
        DeltaRMSE =
          from_row$RMSE -
            to_row$RMSE,
        ImprovedAfterAddingOrgan =
          to_row$RMSE <
            from_row$RMSE
      )

      pos <- pos + 1L
    }
  }

  bind_rows(rows)
}


metab_raw <- read.csv("metabolite_matrix.csv",check.names = FALSE)
micro_raw <- read.csv("microbiome_PMI.csv",check.names = FALSE)
metab_map_raw <- read.csv("metabolit_mapping.csv",check.names = FALSE)

metab_raw[[COL_TARGET]] <- suppressWarnings(
  as.numeric(metab_raw[[COL_TARGET]])
)

micro_raw[[COL_TARGET]] <- suppressWarnings(
  as.numeric(micro_raw[[COL_TARGET]])
)

meta_metab_all <- metab_raw %>%
  select(all_of(META_REQUIRED))

meta_micro_all <- micro_raw %>%
  select(all_of(META_REQUIRED))

feat_metab_all <- metab_raw %>%
  select(-any_of(META_EXCLUDE)) %>%
  select(where(is.numeric))

feat_micro_all <- micro_raw %>%
  select(-any_of(META_EXCLUDE)) %>%
  select(where(is.numeric))

metab_map <- setNames(
  as.character(metab_map_raw[[2]]),
  as.character(metab_map_raw[[1]])
)

map_metab_names <- function(ids) {
  mapped <- metab_map[ids]

  ifelse(
    !is.na(mapped) &
      nchar(trimws(mapped)) > 0,
    mapped,
    ids
  )
}

joint_ids <- intersect(
  meta_metab_all[[COL_ID]],
  meta_micro_all[[COL_ID]]
)

meta_joint <- meta_metab_all %>%
  filter(
    .data[[COL_ID]] %in% joint_ids
  ) %>%
  arrange(.data[[COL_ID]])

match_metab <- match(
  meta_joint[[COL_ID]],
  meta_metab_all[[COL_ID]]
)

match_micro <- match(
  meta_joint[[COL_ID]],
  meta_micro_all[[COL_ID]]
)

feat_metab_joint <- feat_metab_all[
  match_metab, , drop = FALSE
]

feat_micro_joint <- feat_micro_all[
  match_micro, , drop = FALSE
]

metadata_check <- tibble(
  SampleID = meta_joint[[COL_ID]],
  Animal_metab = as.character(
    meta_metab_all[[COL_ANIMAL]][
      match_metab
    ]
  ),
  Animal_micro = as.character(
    meta_micro_all[[COL_ANIMAL]][
      match_micro
    ]
  ),
  Tissue_metab = as.character(
    meta_metab_all[[COL_TISSUE]][
      match_metab
    ]
  ),
  Tissue_micro = as.character(
    meta_micro_all[[COL_TISSUE]][
      match_micro
    ]
  ),
  PMI_metab = as.numeric(
    meta_metab_all[[COL_TARGET]][
      match_metab
    ]
  ),
  PMI_micro = as.numeric(
    meta_micro_all[[COL_TARGET]][
      match_micro
    ]
  )
) %>%
  mutate(
    PMIdiff = PMI_metab - PMI_micro
  ) %>%
  filter(
    Animal_metab != Animal_micro |
      Tissue_metab != Tissue_micro |
      !is.finite(PMI_metab) |
      !is.finite(PMI_micro) |
      abs(PMIdiff) > 1e-6
  )

duplicate_check <- meta_joint %>%
  count(
    .data[[COL_ANIMAL]],
    .data[[COL_TISSUE]],
    name = "NSamples"
  ) %>%
  filter(NSamples > 1L)



animal_pmi_check <- meta_joint %>%
  group_by(.data[[COL_ANIMAL]]) %>%
  summarise(
    NPMI = n_distinct(
      round(.data[[COL_TARGET]], 8)
    ),
    .groups = "drop"
  ) %>%
  filter(NPMI > 1L)



meta_joint <- meta_joint %>%
  filter(
    .data[[COL_TISSUE]] %in% TISSUES
  )

original_joint_ids <-
  meta_metab_all[[COL_ID]][match_metab]

row_index <- match(
  meta_joint[[COL_ID]],
  original_joint_ids
)

feat_metab_joint <- feat_metab_joint[
  row_index, , drop = FALSE
]

feat_micro_joint <- feat_micro_joint[
  row_index, , drop = FALSE
]

# Run tissue-specific LOAO
main_outputs <- lapply(
  setNames(TISSUES, TISSUES),
  function(tissue) {
    idx <- as.character(
      meta_joint[[COL_TISSUE]]
    ) == tissue

    run_loao_analysis(
      meta_joint[idx, , drop = FALSE],
      feat_metab_joint[
        idx, , drop = FALSE
      ],
      feat_micro_joint[
        idx, , drop = FALSE
      ],
      tissue
    )
  }
)

predictions_main <- bind_rows(
  lapply(
    main_outputs,
    `[[`,
    "predictions"
  )
)

selections_main <- bind_rows(
  lapply(
    main_outputs,
    `[[`,
    "selections"
  )
)

stack_weights_main <- bind_rows(
  lapply(
    main_outputs,
    `[[`,
    "stack_weights"
  )
)

features_main <- bind_rows(
  lapply(
    main_outputs,
    `[[`,
    "features"
  )
)

inner_summary_main <- bind_rows(
  lapply(
    main_outputs,
    `[[`,
    "inner_summary"
  )
)

inner_folds_main <- bind_rows(
  lapply(
    main_outputs,
    `[[`,
    "inner_folds"
  )
)

summary_main <- summarise_predictions(
  predictions_main,
  c("Label", "Model")
) %>%
  mutate(
    ModelLabel = unname(
      MODEL_LABELS[Model]
    ),
    AnalysisType = ifelse(
      Model %in% PRIMARY_MODELS,
      "Primary fixed-K",
      "Secondary all-retained-feature benchmark"
    )
  )


SINGLE_OMICS_FIXEDK <- c(
  "Metab_FixedK",
  "Micro_FixedK"
)

MULTI_OMICS_FIXEDK <- c(
  "EF_FixedK",
  "Stacking_FixedK"
)

fixedK_models_for_test <- c(
  SINGLE_OMICS_FIXEDK,
  MULTI_OMICS_FIXEDK
)

fixedK_duplicate_check <- predictions_main %>%
  filter(
    Model %in% fixedK_models_for_test
  ) %>%
  count(
    Label,
    Model,
    AnimalID,
    name = "NPredictions"
  ) %>%
  filter(
    NPredictions != 1L
  )

fixedK_animal_set_check <- predictions_main %>%
  filter(
    Model %in% fixedK_models_for_test
  ) %>%
  group_by(
    Label,
    Model
  ) %>%
  summarise(
    NAnimals = n_distinct(AnimalID),
    AnimalSet = paste(
      sort(
        unique(
          as.character(AnimalID)
        )
      ),
      collapse = "|"
    ),
    .groups = "drop"
  ) %>%
  group_by(Label) %>%
  summarise(
    NAnimalCounts = n_distinct(NAnimals),
    NAnimalSets = n_distinct(AnimalSet),
    .groups = "drop"
  ) %>%
  filter(
    NAnimalCounts != 1L |
      NAnimalSets != 1L
  )

fixedK_single_multi_comparisons <- tidyr::expand_grid(
  Label = TISSUES,
  SingleModel = SINGLE_OMICS_FIXEDK,
  MultiModel = MULTI_OMICS_FIXEDK
) %>%
  purrr::pmap_dfr(
    function(
      Label,
      SingleModel,
      MultiModel
    ) {
      tissue <- Label
      single_model <- SingleModel
      multi_model <- MultiModel

      candidate_df <- predictions_main %>%
        filter(
          .data$Label == .env$tissue,
          .data$Model == .env$multi_model
        )

      reference_df <- predictions_main %>%
        filter(
          .data$Label == .env$tissue,
          .data$Model == .env$single_model
        )

      paired_prediction_comparison(
        candidate_df = candidate_df,
        reference_df = reference_df,
        comparison_label = paste(
          single_model,
          "vs",
          multi_model
        ),
        B = SINGLE_MULTI_COMPARISON_B,
        seed = stable_seed(
          "FixedKSingleVsMulti",
          tissue,
          single_model,
          multi_model
        )
      ) %>%
        mutate(
          Label = tissue,
          SingleModel = single_model,
          MultiModel = multi_model,
          .before = 1
        )
    }
  ) %>%
  group_by(Label) %>%
  mutate(
    PHolmWithinTissue = p.adjust(
      PPermutation,
      method = SINGLE_MULTI_P_ADJUST
    ),
    Significant =
      is.finite(PHolmWithinTissue) &
      PHolmWithinTissue <
        SINGLE_MULTI_ALPHA,
    Significance = case_when(
      PHolmWithinTissue <= 0.0001 ~ "****",
      PHolmWithinTissue < 0.001 ~ "***",
      PHolmWithinTissue < 0.01 ~ "**",
      PHolmWithinTissue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Direction = case_when(
      DeltaRMSE > 0 ~
        "Multi-omics lower RMSE",
      DeltaRMSE < 0 ~
        "Single-omics lower RMSE",
      TRUE ~
        "Equal RMSE"
    )
  ) %>%
  ungroup()

feature_summary_main <- features_main %>%
  group_by(
    Label,
    Omics,
    Feature
  ) %>%
  summarise(
    SelectedFolds = n_distinct(
      OuterFold
    ),
    MeanRank = mean(Rank),
    MedianRank = median(Rank),
    .groups = "drop"
  ) %>%
  left_join(
    features_main %>%
      distinct(
        Label,
        Omics,
        OuterFold
      ) %>%
      count(
        Label,
        Omics,
        name = "NOuterFolds"
      ),
    by = c("Label", "Omics")
  ) %>%
  mutate(
    SelectionFrequency =
      SelectedFolds / NOuterFolds,
    DisplayName = ifelse(
      Omics == "Metabolomics",
      map_metab_names(Feature),
      Feature
    )
  ) %>%
  arrange(
    Label,
    Omics,
    desc(SelectionFrequency),
    MeanRank
  )

write.csv(
  predictions_main,
  file.path(
    output_dir,
    "LOAO_animal_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_main,
  file.path(
    output_dir,
    "LOAO_model_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  fixedK_single_multi_comparisons,
  file.path(
    output_dir,
    "FixedK_single_vs_multi_significance.csv"
  ),
  row.names = FALSE
)

write.csv(
  selections_main,
  file.path(
    output_dir,
    "LOAO_model_selections.csv"
  ),
  row.names = FALSE
)

write.csv(
  stack_weights_main,
  file.path(
    output_dir,
    "LOAO_stacking_weights.csv"
  ),
  row.names = FALSE
)

write.csv(
  feature_summary_main,
  file.path(
    output_dir,
    "FixedK_feature_selection_frequency.csv"
  ),
  row.names = FALSE
)

write.csv(
  inner_summary_main,
  file.path(
    output_dir,
    "RF_inner_CV_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  inner_folds_main,
  file.path(
    output_dir,
    "RF_inner_CV_fold_results.csv"
  ),
  row.names = FALSE
)


make_rmse_plot <- function(
  summary_df,
  models,
  title_text,
  subtitle_text,
  significance_df = NULL
) {
  model_position <- setNames(
    seq_along(models),
    models
  )

  plot_df <- summary_df %>%
    filter(
      Model %in% models
    ) %>%
    mutate(
      Label = factor(
        Label,
        levels = TISSUES
      ),
      Model = factor(
        Model,
        levels = models
      ),
      ModelPosition = unname(
        model_position[
          as.character(Model)
        ]
      )
    )

  p <- ggplot(
    plot_df,
    aes(
      x = ModelPosition,
      y = RMSE,
      ymin = RMSE_lower,
      ymax = RMSE_upper,
      fill = Model
    )
  ) +
    geom_col(
      width = 0.72,
      alpha = 0.84,
      color = "white"
    ) +
    geom_errorbar(
      width = 0.17,
      linewidth = 0.45
    ) +
    facet_wrap(
      ~Label,
      nrow = 1
    ) +
    scale_x_continuous(
      breaks = seq_along(models),
      labels = unname(
        MODEL_LABELS[models]
      ),
      expand = expansion(
        add = c(0.55, 0.55)
      )
    ) +
    scale_fill_manual(
      values = MODEL_COLORS[models],
      guide = "none"
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.02, 0.08)
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = "LOAO RMSE (days)"
    ) +
    theme_bw(
      base_size = 10.5,
      base_family = "Arial"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        size = 8.5
      ),
      strip.text = element_text(
        face = "bold"
      ),
      panel.grid.major.x =
        element_blank(),
      panel.grid.minor =
        element_blank(),
      plot.title = element_text(
        face = "bold"
      ),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10,
        unit = "pt"
      )
    )

  if (!is.null(significance_df)) {
    sig_df <- significance_df %>%
      filter(
        Significant,
        SingleModel %in% models,
        MultiModel %in% models
      ) %>%
      mutate(
        Label = as.character(Label),
        XStart = unname(
          model_position[SingleModel]
        ),
        XEnd = unname(
          model_position[MultiModel]
        ),
        Span = abs(
          XEnd - XStart
        )
      ) %>%
      arrange(
        Label,
        Span,
        XStart,
        XEnd
      ) %>%
      group_by(Label) %>%
      mutate(
        BracketLevel = row_number()
      ) %>%
      ungroup()

    if (nrow(sig_df) > 0L) {
      plot_top <- plot_df %>%
        mutate(
          Label = as.character(Label)
        ) %>%
        group_by(Label) %>%
        summarise(
          PlotTop = max(
            c(RMSE, RMSE_upper),
            na.rm = TRUE
          ),
          .groups = "drop"
        ) %>%
        mutate(
          PlotTop = if_else(
            is.finite(PlotTop),
            PlotTop,
            1
          ),
          Step = pmax(
            0.085 * PlotTop,
            0.12
          ),
          Tip = pmax(
            0.020 * PlotTop,
            0.035
          )
        )

      sig_df <- sig_df %>%
        left_join(
          plot_top,
          by = "Label"
        ) %>%
        mutate(
          Label = factor(
            Label,
            levels = TISSUES
          ),
          Y = PlotTop +
            Step * BracketLevel,
          YTip = Y - Tip,
          YText = Y +
            0.08 * Step,
          XText = (
            XStart + XEnd
          ) / 2
        )

      p <- p +
        geom_segment(
          data = sig_df,
          aes(
            x = XStart,
            xend = XEnd,
            y = Y,
            yend = Y
          ),
          inherit.aes = FALSE,
          linewidth = 0.45
        ) +
        geom_segment(
          data = sig_df,
          aes(
            x = XStart,
            xend = XStart,
            y = YTip,
            yend = Y
          ),
          inherit.aes = FALSE,
          linewidth = 0.45
        ) +
        geom_segment(
          data = sig_df,
          aes(
            x = XEnd,
            xend = XEnd,
            y = YTip,
            yend = Y
          ),
          inherit.aes = FALSE,
          linewidth = 0.45
        ) +
        geom_text(
          data = sig_df,
          aes(
            x = XText,
            y = YText,
            label = Significance
          ),
          inherit.aes = FALSE,
          size = 4.2,
          fontface = "bold",
          vjust = 0
        )
    }
  }

  p
}

p_fixed_rmse <- make_rmse_plot(
  summary_main,
  PRIMARY_MODELS,
  "Fixed-K PMI prediction performance",
  paste0(
    "Tissue-specific LOAO; metabolomics = Elastic Net (alpha = ",
    ALPHA,
    ", Top ",
    K_METAB_PRIMARY,
    "), microbiome = Random Forest (Top ",
    K_MICRO_PRIMARY,
    ", mtry selected by repeated grouped inner CV). ",
    "Error bars are animal-bootstrap 95% confidence intervals. ",
    "Stars indicate significant paired permutation tests after ",
    "within-tissue Holm correction."
  ),
  significance_df =
    fixedK_single_multi_comparisons
)


ggsave(
  file.path(
    output_dir,
    "Fig_5a.tiff"
  ),
  p_fixed_rmse,
  width = 17,
  height = 6,
  device = "tiff",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

p_all_rmse <- make_rmse_plot(
  summary_main,
  ALL_BENCHMARK_MODELS,
  "All-retained-feature PMI prediction performance",
  paste0(
    "Tissue-specific LOAO benchmark after training-only preprocessing; ",
    "metabolomics = Elastic Net, microbiome = Random Forest with ",
    "repeated grouped inner-CV mtry selection. ",
    "Error bars are animal-bootstrap 95% confidence intervals."
  )
)

ggsave(
  file.path(
    output_dir,
    "Fig_ALL_RMSE.tiff"
  ),
  p_all_rmse,
  width = 17,
  height = 6,
  device = "tiff",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)
COMPRESSION_MODELS <- c(
  "Metab_ALL",
  "Metab_FixedK",
  "Micro_ALL",
  "Micro_FixedK",
  "EF_ALL",
  "EF_FixedK",
  "Stacking_ALL",
  "Stacking_FixedK"
)

COMPRESSION_COLORS <- c(
  "Metabolomics" = unname(
    MODEL_COLORS["Metab_FixedK"]
  ),
  "Microbiome" = unname(
    MODEL_COLORS["Micro_FixedK"]
  ),
  "Early fusion" = unname(
    MODEL_COLORS["EF_FixedK"]
  ),
  "Stacking" = unname(
    MODEL_COLORS["Stacking_FixedK"]
  )
)

compression_plot_df <- summary_main %>%
  filter(
    Model %in% COMPRESSION_MODELS
  ) %>%
  mutate(
    Strategy = case_when(
      Model %in% c(
        "Metab_ALL",
        "Metab_FixedK"
      ) ~ "Metabolomics",
      
      Model %in% c(
        "Micro_ALL",
        "Micro_FixedK"
      ) ~ "Microbiome",
      
      Model %in% c(
        "EF_ALL",
        "EF_FixedK"
      ) ~ "Early fusion",
      
      Model %in% c(
        "Stacking_ALL",
        "Stacking_FixedK"
      ) ~ "Stacking",
      
      TRUE ~ NA_character_
    ),
    
    FeatureSet = case_when(
      stringr::str_ends(
        Model,
        "_ALL"
      ) ~ "ALL",
      
      stringr::str_ends(
        Model,
        "_FixedK"
      ) ~ "FixedK",
      
      TRUE ~ NA_character_
    ),
    
    Label = factor(
      Label,
      levels = TISSUES
    ),
    
    Strategy = factor(
      Strategy,
      levels = c(
        "Metabolomics",
        "Microbiome",
        "Early fusion",
        "Stacking"
      )
    ),
    
    FeatureSet = factor(
      FeatureSet,
      levels = c(
        "ALL",
        "FixedK"
      )
    )
  ) %>%
  filter(
    !is.na(Strategy),
    !is.na(FeatureSet)
  )

# Calculate numerical change after compression
compression_delta_df <- compression_plot_df %>%
  select(
    Label,
    Strategy,
    FeatureSet,
    RMSE
  ) %>%
  tidyr::pivot_wider(
    names_from = FeatureSet,
    values_from = RMSE
  ) %>%
  mutate(
    DeltaRMSE_FixedK_minus_ALL =
      FixedK - ALL,
    
    DeltaLabel = sprintf(
      "\u0394=%+.3f",
      DeltaRMSE_FixedK_minus_ALL
    )
  )

compression_panel_ranges <- compression_plot_df %>%
  group_by(Label) %>%
  summarise(
    PanelMin = min(
      RMSE_lower,
      na.rm = TRUE
    ),
    PanelMax = max(
      RMSE_upper,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    PanelSpan = pmax(
      PanelMax - PanelMin,
      0.50
    )
  )

spread_labels_up <- function(
    y,
    gap
) {
  if (length(y) <= 1L) {
    return(y)
  }
  
  ord <- order(y)
  adjusted <- y[ord]
  
  for (i in 2:length(adjusted)) {
    adjusted[i] <- max(
      adjusted[i],
      adjusted[i - 1L] + gap
    )
  }
  
  result <- numeric(length(y))
  result[ord] <- adjusted
  result
}

compression_delta_plot_df <- compression_delta_df %>%
  left_join(
    compression_panel_ranges,
    by = "Label"
  ) %>%
  mutate(
    LabelYInitial =
      pmax(ALL, FixedK) +
      0.035 * PanelSpan,
    
    LabelGap =
      0.065 * PanelSpan
  ) %>%
  group_by(Label) %>%
  mutate(
    LabelY = spread_labels_up(
      LabelYInitial,
      gap = first(LabelGap)
    )
  ) %>%
  ungroup() %>%
  mutate(
    FeatureSet = factor(
      "FixedK",
      levels = c(
        "ALL",
        "FixedK"
      )
    )
  )


compression_summary_table <- compression_delta_df %>%
  mutate(
    Interpretation = case_when(
      DeltaRMSE_FixedK_minus_ALL < 0 ~
        "Fixed-K lower RMSE",
      
      DeltaRMSE_FixedK_minus_ALL > 0 ~
        "ALL lower RMSE",
      
      TRUE ~
        "Equal RMSE"
    )
  ) %>%
  arrange(
    Label,
    Strategy
  )

write.csv(
  compression_summary_table,
  file.path(
    output_dir,
    "FixedK_vs_ALL_compression_summary.csv"
  ),
  row.names = FALSE
)


# Plot
plot_feature_compression <- ggplot(
  compression_plot_df,
  aes(
    x = FeatureSet,
    y = RMSE,
    group = Strategy,
    color = Strategy
  )
) +
  geom_line(
    linewidth = 0.90,
    alpha = 0.82
  ) +
  geom_errorbar(
    aes(
      ymin = RMSE_lower,
      ymax = RMSE_upper
    ),
    width = 0.07,
    linewidth = 0.55,
    alpha = 0.65
  ) +
  geom_point(
    aes(
      shape = FeatureSet
    ),
    size = 4.0,
    stroke = 1.0
  ) +
  
  # Delta-RMSE numerical labels.
  geom_text(
    data = compression_delta_plot_df,
    aes(
      x = FeatureSet,
      y = LabelY,
      label = DeltaLabel,
      color = Strategy
    ),
    inherit.aes = FALSE,
    nudge_x = 0.08,
    hjust = 0,
    size = 3.0,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~Label,
    nrow = 1
  ) +
  
  scale_color_manual(
    values = COMPRESSION_COLORS,
    drop = FALSE
  ) +
  
  # Circle = ALL; diamond = Fixed-K.
  scale_shape_manual(
    values = c(
      "ALL" = 16,
      "FixedK" = 18
    ),
    labels = c(
      "ALL" =
        "All retained features",
      "FixedK" =
        "Fixed-K features"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "ALL" =
        "All retained\nfeatures",
      "FixedK" =
        "Fixed-K\nfeatures"
    ),
    expand = expansion(
      add = c(0.35, 0.85)
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.04, 0.18)
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    title =
      "Feature compression effect: all retained features versus Fixed-K models",
    
    subtitle = paste0(
      "Tissue-specific LOAO comparison. Fixed-K models used the top ",
      K_METAB_PRIMARY,
      " metabolites and top ",
      K_MICRO_PRIMARY,
      " microbial features selected within each outer training partition."
    ),
    
    caption = paste0(
      "\u0394RMSE = Fixed-K RMSE - all-feature RMSE; ",
      "negative values favor the Fixed-K model. ",
      "Error bars are animal-bootstrap 95% confidence intervals."
    ),
    
    x = NULL,
    y = "LOAO RMSE (days)",
    color = "Strategy",
    shape = "Feature set"
  ) +
  
  theme_bw(
    base_size = 10.5,
    base_family = "Arial"
  ) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    
    panel.grid.major.x =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    axis.text.x = element_text(
      size = 9.5
    ),
    
    plot.title = element_text(
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      color = "grey30"
    ),
    
    plot.caption = element_text(
      color = "grey35",
      hjust = 0
    ),
    
    legend.position = "bottom",
    
    legend.box = "horizontal",
    
    plot.margin = ggplot2::margin(
      t = 8,
      r = 25,
      b = 8,
      l = 8
    )
  ) +
  
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(
        linewidth = 1.0,
        shape = 16
      )
    ),
    
    shape = guide_legend(
      order = 2
    )
  )


# Save figure
ggsave(
  file.path(
    output_dir,
    "Figure S5a.tiff"
  ),
  plot_feature_compression,
  width = 17,
  height = 6.6,
  device = "tiff",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

print(
  compression_summary_table,
  n = Inf,
  width = Inf
)

print(plot_feature_compression)


#  Fixed-K feature stability
fixed_feature_summary <- feature_summary_main

N_SHOW_METAB <- min(
  15L,
  K_METAB_PRIMARY
)
N_SHOW_MICRO <- min(
  10L,
  K_MICRO_PRIMARY
)
FREQ_CUT_SHOW <- 0.50
MIN_TISSUES_SHOW <- 2L


ORGAN_SPECIFIC_FREQ_CUT <- 0.50
N_SHOW_ORGAN_METAB <- min(
  15L,
  K_METAB_PRIMARY
)
N_SHOW_ORGAN_MICRO <- min(
  8L,
  K_MICRO_PRIMARY
)

N_METAB_TREND <- min(
  10L,
  K_METAB_PRIMARY
)
N_MICRO_TREND <- min(
  5L,
  K_MICRO_PRIMARY
)

OMICS_BAR_COLORS <- c(
  "Metabolomics" = "#2980B9",
  "Microbiome" = "#E67E22"
)
METAB_TREND_COLOR <- "#1A5276"
MICRO_TREND_COLOR <- "#784212"

save_supplementary_plot <- function(
  plot_object,
  filename_stem,
  width,
  height,
  limitsize = TRUE
) {
  ggsave(
    file.path(
      output_dir,
      paste0(filename_stem, ".png")
    ),
    plot_object,
    width = width,
    height = height,
    dpi = 220,
    limitsize = limitsize,
    bg = "white"
  )

  ggsave(
    file.path(
      output_dir,
      paste0(filename_stem, ".tiff")
    ),
    plot_object,
    width = width,
    height = height,
    device = "tiff",
    dpi = 300,
    compression = "lzw",
    limitsize = limitsize,
    bg = "white"
  )
}

make_feature_label_key <- function(
  feature_summary,
  omics_type,
  feature_pool = NULL
) {
  label_key <- feature_summary %>%
    filter(Omics == omics_type)

  if (!is.null(feature_pool)) {
    label_key <- label_key %>%
      filter(Feature %in% feature_pool)
  }

  label_key %>%
    distinct(
      Feature,
      DisplayName
    ) %>%
    group_by(DisplayName) %>%
    mutate(
      DisplayLabel = ifelse(
        n() > 1L,
        paste0(
          DisplayName,
          " [",
          Feature,
          "]"
        ),
        DisplayName
      )
    ) %>%
    ungroup()
}

get_feature_pool <- function(
  feature_summary,
  omics_type,
  n_show,
  frequency_cut,
  min_tissues
) {
  feature_summary %>%
    filter(Omics == omics_type) %>%
    group_by(Label) %>%
    slice_max(
      SelectionFrequency,
      n = n_show,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    group_by(Feature) %>%
    summarise(
      NQualifyingTissues = n_distinct(
        Label[
          SelectionFrequency >=
            frequency_cut
        ]
      ),
      .groups = "drop"
    ) %>%
    filter(
      NQualifyingTissues >=
        min_tissues
    ) %>%
    pull(Feature) %>%
    unique()
}

make_bubble_data <- function(
  feature_summary,
  feature_pool,
  omics_type
) {
  label_key <- make_feature_label_key(
    feature_summary,
    omics_type,
    feature_pool
  )

  feature_order <- feature_summary %>%
    filter(
      Omics == omics_type,
      Feature %in% feature_pool
    ) %>%
    left_join(
      label_key,
      by = c(
        "Feature",
        "DisplayName"
      )
    ) %>%
    group_by(
      Feature,
      DisplayLabel
    ) %>%
    summarise(
      AverageFrequency =
        mean(SelectionFrequency),
      NTissues = n_distinct(
        Label[
          SelectionFrequency >=
            FREQ_CUT_SHOW
        ]
      ),
      .groups = "drop"
    ) %>%
    arrange(
      desc(NTissues),
      desc(AverageFrequency),
      DisplayLabel
    )

  expand_grid(
    Feature = feature_pool,
    Label = TISSUES
  ) %>%
    left_join(
      feature_summary %>%
        filter(Omics == omics_type) %>%
        select(
          Feature,
          Label,
          SelectionFrequency
        ),
      by = c(
        "Feature",
        "Label"
      )
    ) %>%
    left_join(
      feature_order,
      by = "Feature"
    ) %>%
    mutate(
      SelectionFrequency =
        replace_na(
          SelectionFrequency,
          0
        ),
      Label = factor(
        Label,
        levels = TISSUES
      ),
      DisplayLabel = factor(
        DisplayLabel,
        levels = rev(
          feature_order$DisplayLabel
        )
      )
    )
}

make_bubble_panel <- function(
  df,
  panel_title,
  bubble_color
) {
  right_labels <- df %>%
    distinct(
      DisplayLabel,
      NTissues
    )

  ggplot(
    df,
    aes(
      x = Label,
      y = DisplayLabel
    )
  ) +
    geom_point(
      aes(
        size = ifelse(
          SelectionFrequency > 0,
          SelectionFrequency * 100,
          NA_real_
        ),
        alpha = SelectionFrequency
      ),
      fill = bubble_color,
      color = "white",
      shape = 21,
      stroke = 0.5,
      na.rm = TRUE
    ) +
    geom_text(
      data = right_labels,
      aes(
        x = length(TISSUES) + 0.62,
        y = DisplayLabel,
        label = sprintf(
          "%d/%d",
          NTissues,
          length(TISSUES)
        )
      ),
      inherit.aes = FALSE,
      size = 3.2,
      hjust = 0,
      color = bubble_color,
      family = "Arial"
    ) +
    scale_size_continuous(
      range = c(1.5, 9),
      breaks = c(25, 50, 75, 100),
      limits = c(0, 100)
    ) +
    scale_alpha_continuous(
      range = c(0.15, 0.92),
      limits = c(0, 1),
      guide = "none"
    ) +
    scale_x_discrete(
      expand = expansion(
        add = c(0.5, 1.1)
      )
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = panel_title,
      x = NULL,
      y = NULL,
      size = "Selection frequency (%)"
    ) +
    theme_bw(
      base_size = 11,
      base_family = "Arial"
    ) +
    theme(
      axis.text.x = element_text(
        face = "bold",
        size = 10
      ),
      axis.text.y =
        element_text(size = 9),
      panel.grid.major = element_line(
        color = "grey90",
        linewidth = 0.25
      ),
      panel.grid.minor =
        element_blank(),
      plot.title = element_text(
        size = 12,
        face = "bold",
        color = bubble_color,
        hjust = 0.5
      ),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(
        t = 8,
        r = 24,
        b = 8,
        l = 8
      )
    )
}

pool_metab <- get_feature_pool(
  fixed_feature_summary,
  "Metabolomics",
  N_SHOW_METAB,
  FREQ_CUT_SHOW,
  MIN_TISSUES_SHOW
)

pool_micro <- get_feature_pool(
  fixed_feature_summary,
  "Microbiome",
  N_SHOW_MICRO,
  FREQ_CUT_SHOW,
  MIN_TISSUES_SHOW
)

bubble_panels <- list()

if (length(pool_metab) > 0L) {
  df_metab_bubble <- make_bubble_data(
    fixed_feature_summary,
    pool_metab,
    "Metabolomics"
  )

  bubble_panels[["Metabolomics"]] <-
    make_bubble_panel(
      df_metab_bubble,
      "Metabolomics",
      OMICS_BAR_COLORS[
        "Metabolomics"
      ]
    )
}

if (length(pool_micro) > 0L) {
  df_micro_bubble <- make_bubble_data(
    fixed_feature_summary,
    pool_micro,
    "Microbiome"
  )

  bubble_panels[["Microbiome"]] <-
    make_bubble_panel(
      df_micro_bubble,
      "Microbiome",
      OMICS_BAR_COLORS[
        "Microbiome"
      ]
    )
}

if (length(bubble_panels) > 0L) {
  p_features <- if (
    length(bubble_panels) == 2L
  ) {
    bubble_panels[[1]] /
      bubble_panels[[2]]
  } else {
    bubble_panels[[1]]
  }

  p_features <- p_features +
    plot_annotation(
      title =
        "Cross-tissue stability of fixed-K selected features",
      subtitle = paste0(
        "Primary models use Top ",
        K_METAB_PRIMARY,
        " metabolites and Top ",
        K_MICRO_PRIMARY,
        " microbial features; identities are reselected ",
        "inside every outer LOAO training partition."
      ),
      theme = theme(
        plot.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5,
          family = "Arial"
        ),
        plot.subtitle = element_text(
          size = 9.5,
          color = "grey35",
          hjust = 0.5,
          family = "Arial"
        )
      )
    )

  n_bubble_rows <- sum(
    vapply(
      bubble_panels,
      function(p) {
        n_distinct(
          p$data$DisplayLabel
        )
      },
      integer(1)
    )
  )

  figure_height <- min(
    n_bubble_rows * 0.28 + 5,
    16
  )

  save_supplementary_plot(
    p_features,
    "Figure 5c",
    width = 10,
    height = figure_height,
    limitsize = FALSE
  )
} else {
  message(
    "No cross-tissue feature met the bubble-plot criteria: ",
    "selection frequency >= ",
    FREQ_CUT_SHOW,
    " in at least ",
    MIN_TISSUES_SHOW,
    " tissues."
  )
}


# Organ-specific stable-feature bubble visualization
get_organ_specific_feature_pool <- function(
    feature_summary,
    omics_type,
    n_show,
    frequency_cut
) {
  
 
  strictly_organ_specific_features <-
    feature_summary %>%
    filter(
      Omics == omics_type
    ) %>%
    group_by(Feature) %>%
    summarise(
      NStableTissues = n_distinct(
        Label[
          SelectionFrequency >= frequency_cut
        ]
      ),
      .groups = "drop"
    ) %>%

    filter(
      NStableTissues == 1L
    ) %>%
    pull(Feature)
  

  feature_summary %>%
    filter(
      Omics == omics_type,
      Feature %in%
        strictly_organ_specific_features,
      SelectionFrequency >=
        frequency_cut
    ) %>%
    group_by(Label) %>%
    arrange(
      desc(SelectionFrequency),
      MeanRank,
      .by_group = TRUE
    ) %>%
    slice_head(
      n = n_show
    ) %>%
    ungroup() %>%
    pull(Feature) %>%
    unique()
}

organ_pool_metab <-
  get_organ_specific_feature_pool(
    fixed_feature_summary,
    "Metabolomics",
    N_SHOW_ORGAN_METAB,
    ORGAN_SPECIFIC_FREQ_CUT
  )

organ_pool_micro <-
  get_organ_specific_feature_pool(
    fixed_feature_summary,
    "Microbiome",
    N_SHOW_ORGAN_MICRO,
    ORGAN_SPECIFIC_FREQ_CUT
  )

df_organ_metab_bubble <-
  make_bubble_data(
    fixed_feature_summary,
    organ_pool_metab,
    "Metabolomics"
  ) %>%
  mutate(
    SelectionFrequency = if_else(
      SelectionFrequency >=
        ORGAN_SPECIFIC_FREQ_CUT,
      SelectionFrequency,
      0
    )
  )

df_organ_micro_bubble <-
  make_bubble_data(
    fixed_feature_summary,
    organ_pool_micro,
    "Microbiome"
  ) %>%
  mutate(
    SelectionFrequency = if_else(
      SelectionFrequency >=
        ORGAN_SPECIFIC_FREQ_CUT,
      SelectionFrequency,
      0
    )
  )


p_organ_metab <-
  make_bubble_panel(
    df_organ_metab_bubble,
    "Metabolomics: tissue-specific stable features",
    OMICS_BAR_COLORS["Metabolomics"]
  )

p_organ_micro <-
  make_bubble_panel(
    df_organ_micro_bubble,
    "Microbiome: tissue-specific stable features",
    OMICS_BAR_COLORS["Microbiome"]
  )

p_organ_specific_features <-
  p_organ_metab /
  p_organ_micro +
  plot_annotation(
    title =
      "Organ-specific stability of fixed-K selected features",
    subtitle = paste0(
      "Only features reaching the stability threshold in exactly one tissue ",
      "are shown. Bubble size represents outer-LOAO selection frequency."
    ),
    theme = theme(
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        family = "Arial"
      ),
      plot.subtitle = element_text(
        size = 9.5,
        color = "grey35",
        hjust = 0.5,
        family = "Arial"
      )
    )
  )

n_organ_specific_rows <-
  n_distinct(
    df_organ_metab_bubble$DisplayLabel
  ) +
  n_distinct(
    df_organ_micro_bubble$DisplayLabel
  )

organ_specific_height <- min(
  max(
    8,
    n_organ_specific_rows * 0.27 + 5
  ),
  28
)

save_supplementary_plot(
  p_organ_specific_features,
  "Figure S5b",
  width = 11,
  height = organ_specific_height,
  limitsize = FALSE
)
make_frequency_bar <- function(
  df,
  omics_type
) {
  ggplot(
    df,
    aes(
      x = PlotLabel,
      y = SelectionFrequency * 100
    )
  ) +
    geom_col(
      fill =
        OMICS_BAR_COLORS[omics_type],
      width = 0.72,
      alpha = 0.88
    ) +
    geom_text(
      aes(
        label = sprintf(
          "%.0f%%",
          SelectionFrequency * 100
        )
      ),
      hjust = -0.12,
      size = 3
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 115),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) {
        paste0(x, "%")
      }
    ) +
    labs(
      title = omics_type,
      x = NULL,
      y = "Selection frequency"
    ) +
    theme_bw(
      base_size = 10.5,
      base_family = "Arial"
    ) +
    theme(
      panel.grid.major.y =
        element_blank(),
      panel.grid.minor =
        element_blank(),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

for (tissue in TISSUES) {
  bar_metab <- fixed_feature_summary %>%
    filter(
      Label == tissue,
      Omics == "Metabolomics"
    ) %>%
    slice_max(
      SelectionFrequency,
      n = K_METAB_PRIMARY,
      with_ties = FALSE
    ) %>%
    left_join(
      make_feature_label_key(
        fixed_feature_summary %>%
          filter(Label == tissue),
        "Metabolomics"
      ),
      by = c(
        "Feature",
        "DisplayName"
      )
    ) %>%
    arrange(
      SelectionFrequency,
      desc(MeanRank)
    ) %>%
    mutate(
      PlotLabel = factor(
        DisplayLabel,
        levels = DisplayLabel
      )
    )

  bar_micro <- fixed_feature_summary %>%
    filter(
      Label == tissue,
      Omics == "Microbiome"
    ) %>%
    slice_max(
      SelectionFrequency,
      n = K_MICRO_PRIMARY,
      with_ties = FALSE
    ) %>%
    left_join(
      make_feature_label_key(
        fixed_feature_summary %>%
          filter(Label == tissue),
        "Microbiome"
      ),
      by = c(
        "Feature",
        "DisplayName"
      )
    ) %>%
    arrange(
      SelectionFrequency,
      desc(MeanRank)
    ) %>%
    mutate(
      PlotLabel = factor(
        DisplayLabel,
        levels = DisplayLabel
      )
    )

  if (
    nrow(bar_metab) > 0L &&
      nrow(bar_micro) > 0L
  ) {
    p_bars <-
      make_frequency_bar(
        bar_metab,
        "Metabolomics"
      ) +
      make_frequency_bar(
        bar_micro,
        "Microbiome"
      ) +
      plot_layout(
        ncol = 2,
        widths = c(1.7, 1)
      ) +
      plot_annotation(
        title = paste(
          tissue,
          "fixed-K feature stability"
        )
      )

    save_supplementary_plot(
      p_bars,
      paste0(
        "SuppFig_fixedK_feature_frequency_",
        tissue
      ),
      width = 14,
      height = max(
        5.5,
        max(
          nrow(bar_metab),
          nrow(bar_micro)
        ) * 0.43 + 3.4
      )
    )
  }
}

trend_breaks <- sort(
  unique(
    as.numeric(
      meta_joint[[COL_TARGET]]
    )
  )
)
trend_breaks <- trend_breaks[
  is.finite(trend_breaks)
]

metab_trend_pool <- fixed_feature_summary %>%
  filter(Omics == "Metabolomics") %>%
  group_by(Label) %>%
  slice_max(
    SelectionFrequency,
    n = N_METAB_TREND,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  pull(Feature) %>%
  unique()

micro_trend_pool <- fixed_feature_summary %>%
  filter(Omics == "Microbiome") %>%
  group_by(Label) %>%
  slice_max(
    SelectionFrequency,
    n = N_MICRO_TREND,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  pull(Feature) %>%
  unique()

available_metab_trend <- intersect(
  metab_trend_pool,
  colnames(feat_metab_joint)
)
available_micro_trend <- intersect(
  micro_trend_pool,
  colnames(feat_micro_joint)
)

raw_metab_long <- tibble()
raw_micro_long <- tibble()

if (length(available_metab_trend) > 0L) {
  raw_metab_long <- feat_metab_joint %>%
    select(
      all_of(available_metab_trend)
    ) %>%
    mutate(
      Label =
        meta_joint[[COL_TISSUE]],
      TrendX =
        as.numeric(
          meta_joint[[COL_TARGET]]
        )
    ) %>%
    pivot_longer(
      cols = -c(
        Label,
        TrendX
      ),
      names_to = "Feature",
      values_to = "Value"
    ) %>%
    filter(
      is.finite(TrendX),
      is.finite(Value)
    )
}

if (length(available_micro_trend) > 0L) {
  micro_clr_all <- data.frame(
    clr_mat(feat_micro_joint),
    check.names = FALSE
  )

  raw_micro_long <- micro_clr_all %>%
    select(
      all_of(available_micro_trend)
    ) %>%
    mutate(
      Label =
        meta_joint[[COL_TISSUE]],
      TrendX =
        as.numeric(
          meta_joint[[COL_TARGET]]
        )
    ) %>%
    pivot_longer(
      cols = -c(
        Label,
        TrendX
      ),
      names_to = "Feature",
      values_to = "Value"
    ) %>%
    filter(
      is.finite(TrendX),
      is.finite(Value)
    )
}

make_feature_trend_plot <- function(
  trend_df,
  title_text,
  y_label,
  point_color,
  ncol_facets,
  markdown_strips = FALSE
) {
  p <- ggplot(
    trend_df,
    aes(
      x = TrendX,
      y = Value
    )
  ) +
    geom_point(
      color = point_color,
      size = 1.7,
      alpha = 0.62
    )

  if (
    nrow(trend_df) >= 4L &&
      n_distinct(
        trend_df$TrendX
      ) >= 4L
  ) {
    p <- p +
      geom_smooth(
        method = "loess",
        span = 0.85,
        se = TRUE,
        color = point_color,
        fill = point_color,
        alpha = 0.15,
        linewidth = 0.85
      )
  }

  p +
    scale_x_continuous(
      breaks = trend_breaks
    ) +
    facet_wrap(
      ~PanelLabel,
      scales = "free_y",
      ncol = ncol_facets
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Features are ranked by fixed-K LOAO selection frequency; ",
        "curves are descriptive only."
      ),
      x = "PMI (days)",
      y = y_label
    ) +
    theme_bw(
      base_size = 9.5,
      base_family = "Arial"
    ) +
    theme(
      strip.text = if (
        markdown_strips
      ) {
        ggtext::element_markdown(
          size = 8
        )
      } else {
        element_text(size = 7.8)
      },
      panel.grid.minor =
        element_blank(),
      plot.title =
        element_text(face = "bold")
    )
}

for (tissue in TISSUES) {
  trend_metab <- fixed_feature_summary %>%
    filter(
      Label == tissue,
      Omics == "Metabolomics"
    ) %>%
    slice_max(
      SelectionFrequency,
      n = N_METAB_TREND,
      with_ties = FALSE
    ) %>%
    arrange(
      desc(SelectionFrequency),
      MeanRank
    ) %>%
    left_join(
      make_feature_label_key(
        fixed_feature_summary %>%
          filter(Label == tissue),
        "Metabolomics"
      ),
      by = c(
        "Feature",
        "DisplayName"
      )
    ) %>%
    mutate(
      PanelLabel = sprintf(
        "%s\n(freq=%.0f%%)",
        DisplayLabel,
        100 * SelectionFrequency
      )
    )

  trend_micro <- fixed_feature_summary %>%
    filter(
      Label == tissue,
      Omics == "Microbiome"
    ) %>%
    slice_max(
      SelectionFrequency,
      n = N_MICRO_TREND,
      with_ties = FALSE
    ) %>%
    arrange(
      desc(SelectionFrequency),
      MeanRank
    ) %>%
    left_join(
      make_feature_label_key(
        fixed_feature_summary %>%
          filter(Label == tissue),
        "Microbiome"
      ),
      by = c(
        "Feature",
        "DisplayName"
      )
    ) %>%
    mutate(
      PanelLabel = sprintf(
        "<i>%s</i><br>(freq=%.0f%%)",
        DisplayLabel,
        100 * SelectionFrequency
      )
    )

  trend_metab_df <- raw_metab_long %>%
    filter(
      Label == tissue,
      Feature %in%
        trend_metab$Feature
    ) %>%
    left_join(
      trend_metab %>%
        select(
          Feature,
          PanelLabel
        ),
      by = "Feature"
    ) %>%
    mutate(
      PanelLabel = factor(
        PanelLabel,
        levels =
          trend_metab$PanelLabel
      )
    )

  trend_micro_df <- raw_micro_long %>%
    filter(
      Label == tissue,
      Feature %in%
        trend_micro$Feature
    ) %>%
    left_join(
      trend_micro %>%
        select(
          Feature,
          PanelLabel
        ),
      by = "Feature"
    ) %>%
    mutate(
      PanelLabel = factor(
        PanelLabel,
        levels =
          trend_micro$PanelLabel
      )
    )

  if (nrow(trend_metab_df) > 0L) {
    p_trend_metab <-
      make_feature_trend_plot(
        trend_metab_df,
        paste(
          tissue,
          "stable metabolite trends"
        ),
        "Metabolite abundance",
        METAB_TREND_COLOR,
        ncol_facets = 5L,
        markdown_strips = FALSE
      )

    save_supplementary_plot(
      p_trend_metab,
      paste0(
        "SuppFig_fixedK_trend_metab_",
        tissue
      ),
      width = 14,
      height = ceiling(
        max(
          1,
          nrow(trend_metab)
        ) / 5
      ) * 3.6 + 2.3
    )
  }

  if (nrow(trend_micro_df) > 0L) {
    p_trend_micro <-
      make_feature_trend_plot(
        trend_micro_df,
        paste(
          tissue,
          "stable microbial trends"
        ),
        "CLR abundance",
        MICRO_TREND_COLOR,
        ncol_facets = min(
          5L,
          max(
            1L,
            nrow(trend_micro)
          )
        ),
        markdown_strips = TRUE
      )

    save_supplementary_plot(
      p_trend_micro,
      paste0(
        "SuppFig_fixedK_trend_micro_",
        tissue
      ),
      width = max(
        10,
        max(
          1,
          nrow(trend_micro)
        ) * 2.7
      ),
      height = 5
    )
  }
}


#  Shared complete-cohort pairwise cross-organ analysis
complete_data <- build_complete_organ_dataset(
  meta_joint,
  feat_metab_joint,
  feat_micro_joint,
  TISSUES
)

multi_output <-
  run_shared_complete_multi_organ_analysis(
    complete_data
  )

predictions_multi <- multi_output$predictions
weights_multi <- multi_output$weights
oof_predictions_multi <-
  multi_output$oof_predictions

multi_results <- summarise_shared_multi_organ(
  predictions_multi
)


predictions_pairwise <- predictions_multi %>%
  filter(
    NOrgans %in% c(1L, 2L)
  )

summary_pairwise <- multi_results$summary %>%
  filter(
    NOrgans %in% c(1L, 2L)
  )

descriptive_pairwise <-
  multi_results$descriptive %>%
  filter(NOrgans == 2L)

comparisons_pairwise <-
  multi_results$comparisons %>%
  filter(NOrgans == 2L)

weights_pairwise <- weights_multi %>%
  filter(
    NOrgans %in% c(1L, 2L)
  )

pair_best_key <- descriptive_pairwise %>%
  transmute(
    NOrgans,
    OrganSet,
    CandidateMethod = Method,
    BestConstituent,
    BestConstituentRMSE,
    CandidateRMSE,
    DescriptiveDeltaRMSE =
      DeltaRMSE_vs_BestConstituent
  )

pair_vs_best_paired <-
  comparisons_pairwise %>%
  filter(
    ReferenceType ==
      "Constituent single organ"
  ) %>%
  inner_join(
    pair_best_key,
    by = c(
      "NOrgans",
      "OrganSet",
      "CandidateMethod"
    )
  ) %>%
  filter(
    ReferenceModel ==
      BestConstituent
  ) %>%
  mutate(
    ComparisonScope =
      "Post hoc lower-RMSE constituent comparison; descriptive",
    .after = ReferenceModel
  ) %>%
  arrange(
    OrganSet,
    CandidateMethod
  )

write.csv(
  predictions_pairwise,
  file.path(
    output_dir,
    "Pairwise_shared_LOAO_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_pairwise,
  file.path(
    output_dir,
    "Pairwise_shared_model_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  descriptive_pairwise,
  file.path(
    output_dir,
    "Pairwise_vs_best_constituent_descriptive.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_vs_best_paired,
  file.path(
    output_dir,
    "Pairwise_vs_best_constituent_paired_descriptive.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparisons_pairwise,
  file.path(
    output_dir,
    "Pairwise_paired_comparisons.csv"
  ),
  row.names = FALSE
)

write.csv(
  weights_pairwise,
  file.path(
    output_dir,
    "Pairwise_shared_weights.csv"
  ),
  row.names = FALSE
)

write.csv(
  oof_predictions_multi,
  file.path(
    output_dir,
    "Pairwise_shared_training_OOF_predictions.csv"
  ),
  row.names = FALSE
)


predictions_five_sensitivity <-
  predictions_multi %>%
  filter(
    NOrgans %in% c(
      1L,
      length(TISSUES)
    )
  )

summary_five_sensitivity <-
  multi_results$summary %>%
  filter(
    NOrgans %in% c(
      1L,
      length(TISSUES)
    )
  )

descriptive_five_sensitivity <-
  multi_results$descriptive %>%
  filter(
    NOrgans == length(TISSUES)
  )

comparisons_five_sensitivity <-
  multi_results$comparisons %>%
  filter(
    NOrgans == length(TISSUES)
  )

weights_five_sensitivity <-
  weights_multi %>%
  filter(
    NOrgans == length(TISSUES)
  )

write.csv(
  predictions_five_sensitivity,
  file.path(
    output_dir,
    "Supplementary_FiveOrgan_LOAO_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_five_sensitivity,
  file.path(
    output_dir,
    "Supplementary_FiveOrgan_model_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  descriptive_five_sensitivity,
  file.path(
    output_dir,
    "Supplementary_FiveOrgan_vs_best_constituent_descriptive.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparisons_five_sensitivity,
  file.path(
    output_dir,
    "Supplementary_FiveOrgan_paired_comparisons.csv"
  ),
  row.names = FALSE
)

write.csv(
  weights_five_sensitivity,
  file.path(
    output_dir,
    "Supplementary_FiveOrgan_weights.csv"
  ),
  row.names = FALSE
)

primary_complete_predictions <- predictions_main %>%
  filter(
    Model == "Stacking_FixedK",
    AnimalID %in%
      as.character(
        complete_data$meta[[COL_ANIMAL]]
      )
  ) %>%
  transmute(
    Tissue = Label,
    AnimalID,
    Observed,
    Predicted,
    Residual,
    AbsoluteError,
    SquaredError
  )

primary_complete_summary <-
  summarise_predictions(
    primary_complete_predictions,
    c("Tissue")
  ) %>%
  rename_with(
    ~paste0("PrimaryRestricted_", .x),
    -Tissue
  )

complete_refit_summary <- predictions_multi %>%
  filter(Method == "SingleOrgan") %>%
  mutate(Tissue = OrganSet) %>%
  summarise_predictions(
    c("Tissue")
  ) %>%
  rename_with(
    ~paste0("CompleteRefit_", .x),
    -Tissue
  )

primary_vs_complete_single <-
  primary_complete_summary %>%
  inner_join(
    complete_refit_summary,
    by = "Tissue"
  ) %>%
  mutate(
    DeltaRMSE_CompleteRefit_minus_PrimaryRestricted =
      CompleteRefit_RMSE -
        PrimaryRestricted_RMSE
  )

write.csv(
  primary_vs_complete_single,
  file.path(
    output_dir,
    "Pairwise_primary_vs_complete_single_bridge.csv"
  ),
  row.names = FALSE
)


# Pair-focused exploratory visualizations
method_labels_multi <- c(
  "SingleOrgan" = "Single-organ stacking",
  "EqualWeight" = "Equal-weight averaging",
  "HierarchicalNNLS" =
    "Hierarchical NNLS stacking"
)

# Color-blind-friendly palette.
PAIR_GAIN_COLOR <- "#009E73"
PAIR_LOSS_COLOR <- "#D55E00"
EQUAL_WEIGHT_COLOR <- "#E69F00"
HIERARCHICAL_COLOR <- "#0072B2"
SINGLE_ORGAN_COLOR <- "#68778D"
REFERENCE_LINE_COLOR <- "#6B7280"

theme_pmi_forest <- function(
  base_size = 11
) {
  theme_minimal(
    base_size = base_size,
    base_family = "Arial"
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.28),
        color = "#1F2937"
      ),
      plot.subtitle = element_text(
        color = "#4B5563",
        margin = ggplot2::margin(b = 10)
      ),
      plot.caption = element_text(
        size = rel(0.82),
        color = "#6B7280",
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      axis.title = element_text(
        color = "#273142"
      ),
      axis.text = element_text(
        color = "#4B5563"
      ),
      axis.text.y = element_text(
        face = "bold"
      ),
      panel.grid.major.y = element_line(
        color = "#E8ECF1",
        linewidth = 0.45
      ),
      panel.grid.major.x = element_line(
        color = "#EEF1F4",
        linewidth = 0.4
      ),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        color = "#374151"
      ),
      plot.margin = ggplot2::margin(
        t = 12,
        r = 20,
        b = 10,
        l = 12
      )
    )
}

pair_nnls_main_df <-
  pair_vs_best_paired %>%
  filter(
    CandidateMethod ==
      "HierarchicalNNLS"
  ) %>%
  mutate(
    Direction = ifelse(
      DeltaRMSE > 0,
      "Lower RMSE than best constituent",
      "Higher RMSE than best constituent"
    ),
    Direction = factor(
      Direction,
      levels = c(
        "Lower RMSE than best constituent",
        "Higher RMSE than best constituent"
      )
    ),
    ValueLabel = sprintf(
      "%+.2f",
      DeltaRMSE
    ),
    OrganSet = forcats::fct_reorder(
      OrganSet,
      DeltaRMSE
    )
  )

plot_pair_delta <- ggplot(
  pair_nnls_main_df,
  aes(
    x = DeltaRMSE,
    y = OrganSet
  )
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 0,
    ymin = -Inf,
    ymax = Inf,
    fill = PAIR_LOSS_COLOR,
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = PAIR_GAIN_COLOR,
    alpha = 0.055
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.75,
    color = REFERENCE_LINE_COLOR
  ) +
  geom_errorbarh(
    aes(
      xmin = DeltaRMSELower,
      xmax = DeltaRMSEUpper,
      color = Direction
    ),
    height = 0.17,
    linewidth = 0.9,
    alpha = 0.82
  ) +
  geom_point(
    aes(fill = Direction),
    shape = 21,
    size = 4.2,
    stroke = 0.9,
    color = "white"
  ) +
  geom_text(
    aes(
      label = ValueLabel,
      color = Direction,
      hjust = ifelse(
        DeltaRMSE >= 0,
        -0.38,
        1.38
      )
    ),
    size = 3.25,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Lower RMSE than best constituent" =
        PAIR_GAIN_COLOR,
      "Higher RMSE than best constituent" =
        PAIR_LOSS_COLOR
    ),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Lower RMSE than best constituent" =
        PAIR_GAIN_COLOR,
      "Higher RMSE than best constituent" =
        PAIR_LOSS_COLOR
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.10, 0.14)
    )
  ) +
  labs(
    title =
      "Incremental value of two-organ hierarchical stacking",
    subtitle = paste0(
      "Positive values indicate lower LOAO RMSE than the better constituent ",
      "single-organ model; all pairs were evaluated in the same ",
      complete_data$n_animals,
      " animals."
    ),
    caption = paste0(
      "Points show paired RMSE differences; horizontal bars are animal-bootstrap ",
      "95% confidence intervals. The better constituent was identified ",
      "descriptively within each pair."
    ),
    x =
      "Best constituent RMSE - hierarchical NNLS RMSE (days)",
    y = NULL,
    color = NULL,
    fill = NULL
  ) +
  theme_pmi_forest(
    base_size = 11
  )

save_supplementary_plot(
  plot_pair_delta,
  "Figure 5b",
  width = 10,
  height = 6.7
)

# Final compact output
final_summary <- summary_main %>%
  filter(Model %in% CORE_MODELS) %>%
  mutate(
    Label = factor(
      Label,
      levels = TISSUES
    ),
    ModelOrder = match(
      Model,
      CORE_MODELS
    )
  ) %>%
  arrange(
    Label,
    AnalysisType,
    ModelOrder
  ) %>%
  select(-ModelOrder)

write.csv(
  final_summary,
  file.path(
    output_dir,
    "FINAL_LOAO_summary_table.csv"
  ),
  row.names = FALSE
)


