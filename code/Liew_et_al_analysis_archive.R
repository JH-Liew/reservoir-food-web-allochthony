# R analyses accompanying Liew et al.
#
# 1. Packages and file locations ---------------------------------------------

library(readxl)
library(ggplot2)
library(vegan)
library(NetIndices)
library(igraph)
library(lme4)
library(rjags)
library(coda)
library(loo)

if (file.exists("reservoir_allochthony_data.xlsx")) {
  project_dir <- "."
} else if (file.exists("../reservoir_allochthony_data.xlsx")) {
  project_dir <- ".."
} else {
  stop("Cannot find reservoir_allochthony_data.xlsx")
}

project_dir <- normalizePath(project_dir, mustWork = TRUE)
workbook_file <- file.path(project_dir, "reservoir_allochthony_data.xlsx")
default_food_web_dir <- file.path(project_dir, "data", "food_webs")
configured_food_web_dir <- Sys.getenv("ALLOCHTHONY_MATRIX_DIR", unset = "")
food_web_dir <- if (nzchar(configured_food_web_dir))
  configured_food_web_dir else default_food_web_dir
output_dir <- file.path(project_dir, "outputs")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

reservoir_ids <- paste("Res", 1:12)
set.seed(20260904)

# 2. Analysis options --------------------------------------------------------

# Change a switch to TRUE to rerun that computationally intensive analysis.
run_primary_models <- FALSE
run_ndvi_models <- FALSE
run_area_ratio_model <- FALSE
run_res4_exclusion <- FALSE
run_biomass_model <- FALSE
run_feasibility_audit <- FALSE
run_exact_loo <- FALSE

# 3. Read reservoir and taxon data -------------------------------------------

read_table <- function(sheet) {
  as.data.frame(
    readxl::read_excel(workbook_file, sheet = sheet),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

reservoir_environment <- read_table("reservoir_environment")
taxon_allochthony_input <- read_table("taxon_allochthony")
taxon_metadata <- read_table("taxon_metadata")
mixing_inputs <- read_table("mixing_inputs")
feasibility_consumer <- read_table("feasibility_consumer")
taxon_replication <- read_table("taxon_replication")
stored_model_coefficients <- read_table("model_coefficients")
stored_model_comparison <- read_table("model_comparison")
stored_model_residuals <- read_table("model_residuals")
stored_predictor_correlations <- read_table("predictor_correlations")
ndvi_inputs <- read_table("ndvi_inputs")
stored_ndvi_coefficients <- read_table("ndvi_coefficients")
stored_ndvi_model_comparison <- read_table("ndvi_model_comparison")
stored_ndvi_correlations <- read_table("ndvi_correlations")
stored_exact_loo_comparison <- read_table("exact_loo_comparison")
stored_exact_loo_pointwise <- read_table("exact_loo_pointwise")

required_columns <- c(
  "reservoir", "mean_allochthony", "impervious_surfaces", "mean_rain",
  "max_rain", "tn", "tp", "biovolume", "cyano_prop", "res_area",
  "catch_area", "fish_biomass_top3_mean"
)
missing_columns <- setdiff(required_columns, names(reservoir_environment))
if (length(missing_columns)) {
  stop("reservoir_environment is missing: ",
       paste(missing_columns, collapse = ", "))
}
reservoir_data <- reservoir_environment[
  match(reservoir_ids, reservoir_environment$reservoir), , drop = FALSE
]
if (anyNA(reservoir_data$reservoir)) {
  stop("Reservoir identifiers are incomplete in reservoir_environment.")
}

# 4. Read food-web matrices ---------------------------------------------------

clean_taxon_name <- function(x) {
  x <- gsub("\\.", " ", trimws(x))
  x[x == "Cor"] <- "COr"
  x[x == "Omo"] <- "OMo"
  x
}

read_predation_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing prey–consumer matrix: ", path)
  x <- read.csv(path, row.names = 1, check.names = FALSE,
                na.strings = c("", "-", "NA"), stringsAsFactors = FALSE)
  row_names <- clean_taxon_name(rownames(x))
  col_names <- clean_taxon_name(colnames(x))
  x <- as.matrix(data.frame(lapply(x, as.numeric), check.names = FALSE))
  rownames(x) <- row_names
  colnames(x) <- col_names
  x[is.na(x)] <- 0
  if (!identical(rownames(x), colnames(x))) {
    stop("Prey and consumer names differ in ", basename(path))
  }
  if (any(x < 0)) stop("Negative diet fraction in ", basename(path))
  x
}

# The published release uses either spaces or underscores in the reservoir
# part of the filenames. Both forms are accepted without changing the inputs.
matrix_files <- character(length(reservoir_ids))
for (i in seq_along(reservoir_ids)) {
  possible_files <- file.path(
    food_web_dir,
    c(
      paste0("Bottom-up_predation_matrices_", reservoir_ids[i], ".csv"),
      paste0("Bottom-up_predation_matrices_Res_", i, ".csv")
    )
  )
  existing_file <- possible_files[file.exists(possible_files)]
  if (!length(existing_file)) {
    stop(
      "Posterior-median food-web matrix not found for ", reservoir_ids[i],
      ". Download the Wilkinson et al. matrices and place them in data/food_webs, ",
      "or set ALLOCHTHONY_MATRIX_DIR."
    )
  }
  matrix_files[i] <- existing_file[1]
}

propagate_allochthony <- function(web, reservoir) {
  taxa <- rownames(web)
  terrestrial <- grepl("riparian|terrestrial", tolower(taxa))
  if (!any(terrestrial)) {
    stop("No riparian or terrestrial producer rows found in ", reservoir)
  }
  direct <- as.numeric(colSums(web[terrestrial, , drop = FALSE]))
  names(direct) <- taxa
  current <- direct
  contributions <- list(current)
  for (step in seq_len(nrow(web) - 1L)) {
    next_step <- as.numeric(crossprod(web, current))
    names(next_step) <- taxa
    next_step[is.na(next_step)] <- 0
    if (max(abs(next_step), na.rm = TRUE) < 1e-12) break
    contributions[[length(contributions) + 1L]] <- next_step
    current <- next_step
  }
  total <- rowSums(do.call(cbind, contributions), na.rm = TRUE)
  names(total) <- taxa

  trophic_position <- tryCatch({
    tl <- NetIndices::TrophInd(web)[, "TL"]
    setNames(round(as.numeric(tl), 2), rownames(web))[taxa]
  }, error = function(e) rep(NA_real_, length(taxa)))

  data.frame(
    reservoir = reservoir,
    taxon = taxa,
    direct_allochthony = direct,
    allochthony = as.numeric(total),
    trophic_position = as.numeric(trophic_position),
    stringsAsFactors = FALSE
  )
}

predation_matrices <- vector("list", length(reservoir_ids))
names(predation_matrices) <- reservoir_ids
for (i in seq_along(reservoir_ids)) {
  predation_matrices[[i]] <- read_predation_matrix(matrix_files[i])
}

# 5. Calculate taxon-level allochthony ----------------------------------------

# Matrix rows are resources or prey and columns are consumers. Entries are
# posterior-median proportional diet fractions from the Wilkinson et al.
# food-web matrices. Terrestrial basal sources have value of 1, aquatic basal
# sources have value of 0. Macrophytes are classified as an aquatic resource.
# Links are propagated through the quantitative feeding network without
# abundance or biomass weighting.
taxon_allochthony <- do.call(
  rbind,
  Map(propagate_allochthony, predation_matrices, names(predation_matrices))
)
rownames(taxon_allochthony) <- NULL
write.csv(
  taxon_allochthony,
  file.path(output_dir, "taxon_allochthony_from_matrices.csv"),
  row.names = FALSE
)

# 6. Calculate reservoir-level indicators ------------------------------------

direct_taxa <- subset(taxon_allochthony, direct_allochthony > 0)
mean_direct_allochthony <- aggregate(
  direct_allochthony ~ reservoir, data = direct_taxa, FUN = mean
)
names(mean_direct_allochthony)[2] <- "mean_direct_allochthony"

consumer_taxa <- subset(taxon_allochthony, trophic_position > 1)
mean_overall_allochthony <- aggregate(
  allochthony ~ reservoir, data = consumer_taxa, FUN = mean
)
names(mean_overall_allochthony)[2] <- "matrix_mean_allochthony"

reservoir_indicator_check <- merge(
  reservoir_data[, c("reservoir", "mean_allochthony")],
  mean_overall_allochthony, by = "reservoir", all = TRUE
)
write.csv(
  mean_direct_allochthony,
  file.path(output_dir, "reservoir_direct_allochthony.csv"),
  row.names = FALSE
)
write.csv(
  mean_overall_allochthony,
  file.path(output_dir, "reservoir_allochthony_from_matrices.csv"),
  row.names = FALSE
)
write.csv(
  reservoir_indicator_check,
  file.path(output_dir, "matrix_to_workbook_check.csv"),
  row.names = FALSE
)

# 7. Food-web and ordination analyses ----------------------------------------

make_split_web <- function(web, taxon_values, output_file = NULL) {
  taxa <- colnames(web)
  consumers <- taxa[colSums(web) > 0]
  basal <- setdiff(taxa, consumers)
  node_a <- setNames(taxon_values$allochthony, taxon_values$taxon)[taxa]
  node_a[intersect(c("riparian plants", "riparian grasses"), taxa)] <- 1
  node_a[is.na(node_a)] <- 0
  edge_list <- do.call(rbind, lapply(consumers, function(consumer) {
    resources <- taxa[web[, consumer] > 0]
    if (!length(resources)) return(NULL)
    fractions <- web[resources, consumer]
    terrestrial <- fractions * node_a[resources]
    aquatic <- fractions * (1 - node_a[resources])
    rbind(
      data.frame(from = resources[terrestrial > 0], to = consumer,
                 diet_fraction = terrestrial[terrestrial > 0],
                 pathway = "allochthonous"),
      data.frame(from = resources[aquatic > 0], to = consumer,
                 diet_fraction = aquatic[aquatic > 0],
                 pathway = "autochthonous")
    )
  }))
  graph <- igraph::graph_from_data_frame(
    edge_list,
    vertices = data.frame(name = taxa, stringsAsFactors = FALSE),
    directed = TRUE
  )
  vertex_tl <- tryCatch(NetIndices::TrophInd(web)[, "TL"],
                        error = function(e) seq_along(taxa))
  vertex_tl <- setNames(as.numeric(vertex_tl), rownames(web))[taxa]
  level <- pmax(1, round(vertex_tl))
  x_coord <- numeric(length(taxa))
  for (lev in sort(unique(level))) {
    index <- which(level == lev)
    x_coord[index] <- if (length(index) == 1) 3.75 else
      seq(0, 7.5, length.out = length(index))
  }
  layout <- cbind(x_coord, vertex_tl - 1)
  pies <- lapply(taxa, function(taxon) {
    if (taxon %in% basal) return(c(1, 0))
    c(node_a[taxon], 1 - node_a[taxon])
  })
  V(graph)$shape <- ifelse(V(graph)$name %in% basal, "square", "pie")
  V(graph)$frame.color <- "black"
  V(graph)$color <- "white"
  V(graph)$size <- ifelse(V(graph)$name %in% basal, 18, 14)
  E(graph)$width <- E(graph)$diet_fraction * 10
  E(graph)$color <- ifelse(E(graph)$pathway == "allochthonous",
                           "seagreen", "lightblue")
  draw <- function() {
    plot(graph, layout = layout, vertex.pie = pies, vertex.label = NA,
         edge.arrow.size = 0.4, rescale = TRUE, asp = 0, margin = 0.05)
    legend("topleft", c("Primary producer", "Consumer"),
           pch = c(22, 21), pt.cex = 2, bty = "n")
    legend("topleft", c("Allochthonous", "Aquatic"),
           col = c("seagreen", "lightblue"), lwd = 5, bty = "n",
           inset = c(0, 0.14))
  }
  if (!is.null(output_file)) {
    tiff(output_file, width = 4000, height = 4000, res = 600, bg = "white")
    draw()
    dev.off()
  }
  invisible(graph)
}

make_split_web(
  predation_matrices[["Res 6"]],
  subset(taxon_allochthony, reservoir == "Res 6"),
  file.path(figure_dir, "Figure2_Res6_food_web.tif")
)

direct_long <- subset(
  taxon_allochthony, direct_allochthony > 0,
  select = c("reservoir", "taxon", "direct_allochthony")
)
community <- xtabs(direct_allochthony ~ reservoir + taxon, data = direct_long)
community <- unclass(community[reservoir_ids, , drop = FALSE])
ordination <- vegan::metaMDS(
  community, distance = "bray", k = 4, maxit = 10000,
  trymax = 5000, wascores = TRUE, trace = FALSE
)
site_scores <- as.data.frame(scores(ordination, display = "sites")[, 1:2])
names(site_scores) <- c("axis1", "axis2")
site_scores$reservoir <- rownames(site_scores)
taxon_scores <- as.data.frame(scores(ordination, display = "species")[, 1:2])
names(taxon_scores) <- c("axis1", "axis2")
taxon_scores$taxon <- rownames(taxon_scores)
taxon_scores$group <- ifelse(
  taxon_scores$taxon %in% taxon_metadata$taxon_code,
  "fish", "invertebrate"
)
nmds_plot <- ggplot() +
  geom_point(data = site_scores, aes(axis1, axis2), size = 3) +
  geom_point(data = taxon_scores, aes(axis1, axis2, shape = group), size = 2) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "nMDS1", y = "nMDS2", shape = NULL)
ggsave(
  file.path(figure_dir, "Figure3_nMDS_base.png"), nmds_plot,
  width = 7, height = 5, dpi = 600
)
write.csv(site_scores, file.path(output_dir, "nmds_site_scores.csv"), row.names = FALSE)
write.csv(taxon_scores, file.path(output_dir, "nmds_taxon_scores.csv"), row.names = FALSE)

# Taxon-level relationships and crossed random effects ----------------------

taxon_relationships <- NULL
mixed_coefficients <- NULL
if (!is.null(mean_direct_allochthony)) {
  taxon_relationships <- subset(taxon_allochthony, trophic_position > 1)
  taxon_relationships <- merge(
    taxon_relationships, mean_direct_allochthony, by = "reservoir", all.x = TRUE,
    sort = FALSE
  )
  taxon_relationships$relative_allochthony <-
    taxon_relationships$allochthony - taxon_relationships$mean_direct_allochthony
  taxon_relationships$absolute_value_relative_allochthony <-
    abs(taxon_relationships$relative_allochthony)
  taxon_relationships$z_trophic_position <-
    as.numeric(scale(taxon_relationships$trophic_position))
  signed_mixed <- lme4::lmer(
    relative_allochthony ~ z_trophic_position +
      (1 | reservoir) + (1 | taxon),
    data = taxon_relationships, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa")
  )
  absolute_mixed <- lme4::lmer(
    absolute_value_relative_allochthony ~ z_trophic_position +
      (1 | reservoir) + (1 | taxon),
    data = taxon_relationships, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa")
  )
  tidy_mixed <- function(fit, model_id) {
    co <- summary(fit)$coefficients
    data.frame(
      model_id = model_id, parameter = rownames(co), estimate = co[, "Estimate"],
      standard_error = co[, "Std. Error"], t_value = co[, "t value"],
      singular_fit = lme4::isSingular(fit, tol = 1e-4),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  mixed_coefficients <- rbind(
    tidy_mixed(signed_mixed, "signed_relative_allochthony"),
    tidy_mixed(absolute_mixed, "absolute_relative_allochthony")
  )
  write.csv(mixed_coefficients,
            file.path(output_dir, "mixed_effects_coefficients.csv"),
            row.names = FALSE)
  variance_components <- rbind(
    data.frame(model_id = "signed_relative_allochthony",
               as.data.frame(VarCorr(signed_mixed))),
    data.frame(model_id = "absolute_relative_allochthony",
               as.data.frame(VarCorr(absolute_mixed)))
  )
  write.csv(variance_components,
            file.path(output_dir, "mixed_effects_variance_components.csv"),
            row.names = FALSE)
}

# 8. Primary reservoir-level models -------------------------------------------

# Predictors are divided by their uncentred root-mean-square magnitudes and
# are not mean-centred. The denominator is n - 1 in the primary analysis.
n_reservoirs <- nrow(reservoir_data)
catchment_area_scale <- sqrt(
  sum(reservoir_data$catch_area^2) / (n_reservoirs - 1)
)
reservoir_area_scale <- sqrt(
  sum(reservoir_data$res_area^2) / (n_reservoirs - 1)
)
mean_rainfall_scale <- sqrt(
  sum(reservoir_data$mean_rain^2) / (n_reservoirs - 1)
)
peak_rainfall_scale <- sqrt(
  sum(reservoir_data$max_rain^2) / (n_reservoirs - 1)
)
phytoplankton_biovolume_scale <- sqrt(
  sum(reservoir_data$biovolume^2) / (n_reservoirs - 1)
)
total_nitrogen_scale <- sqrt(
  sum(reservoir_data$tn^2) / (n_reservoirs - 1)
)
total_phosphorus_scale <- sqrt(
  sum(reservoir_data$tp^2) / (n_reservoirs - 1)
)

catchment_area_scaled <- reservoir_data$catch_area / catchment_area_scale
reservoir_area_scaled <- reservoir_data$res_area / reservoir_area_scale
mean_rainfall_scaled <- reservoir_data$mean_rain / mean_rainfall_scale
peak_rainfall_scaled <- reservoir_data$max_rain / peak_rainfall_scale
phytoplankton_biovolume_scaled <-
  reservoir_data$biovolume / phytoplankton_biovolume_scale
total_nitrogen_scaled <- reservoir_data$tn / total_nitrogen_scale
total_phosphorus_scaled <- reservoir_data$tp / total_phosphorus_scale
impervious_cover <- reservoir_data$impervious_surfaces
cyanobacteria_proportion <- reservoir_data$cyano_prop

X_peak_runoff <- data.frame(
  catchment_size = catchment_area_scaled,
  rain_max = peak_rainfall_scaled,
  impervious = impervious_cover
)
X_average_runoff <- data.frame(
  catchment_size = catchment_area_scaled,
  rain_mean = mean_rainfall_scaled,
  impervious = impervious_cover
)
X_reservoir_null <- data.frame(reservoir_size = reservoir_area_scaled)
X_catchment_null <- data.frame(catchment_size = catchment_area_scaled)
X_potential_production <- data.frame(
  reservoir_size = reservoir_area_scaled,
  tn_scaled = total_nitrogen_scaled,
  tp_scaled = total_phosphorus_scaled
)
X_actual_production <- data.frame(
  reservoir_size = reservoir_area_scaled,
  phyto_biovolume = phytoplankton_biovolume_scaled,
  cyanobacteria = cyanobacteria_proportion
)

model_labels <- c(
  peak_catchment_runoff = "Peak catchment runoff",
  average_catchment_runoff = "Average catchment runoff",
  reservoir_size_null = "Reservoir-size null",
  catchment_size_null = "Catchment-size null",
  potential_aquatic_production = "Potential aquatic production",
  actual_aquatic_production = "Actual aquatic production",
  ndvi_only = "NDVI only",
  vegetation_runoff_alternative = "Vegetation-runoff alternative",
  peak_runoff_reference = "Peak-runoff reference",
  peak_runoff_plus_ndvi = "Peak runoff + NDVI",
  area_ratio_sensitivity = "Ratio sensitivity",
  overall_direct_allochthony = "Mean overall versus mean direct allochthony (Equation 5)"
)

make_gaussian_model <- function(predictor_names) {
  terms <- paste0(
    "b[", seq_along(predictor_names), "]*x",
    seq_along(predictor_names), "[i]"
  )
  paste0(
    "model {\n",
    "for (i in 1:N) {\n",
    "  y[i] ~ dnorm(mu[i], tau)\n",
    "  mu[i] <- b0 + ", paste(c(terms), collapse = " + "), "\n",
    "  loglik[i] <- logdensity.norm(y[i], mu[i], tau)\n",
    "}\n",
    "b0 ~ dnorm(0, 0.001)\n",
    "for (j in 1:K) { b[j] ~ dnorm(0, 0.001) }\n",
    "sigma ~ dunif(0, 1)\n",
    "tau <- pow(sigma, -2)\n",
    "}\n"
  )
}

fit_reservoir_model <- function(response, predictors, model_id, model_label,
                                reservoir, seed_offset, adaptation = 10000,
                                burn_in = 100000, retained = 10000) {
  predictor_names <- names(predictors)
  n_predictors <- length(predictors)
  model_data <- list(y = as.numeric(response), N = length(response), K = n_predictors)
  for (j in seq_len(n_predictors)) {
    model_data[[paste0("x", j)]] <- as.numeric(predictors[[j]])
  }
  inits <- lapply(seq_len(4), function(chain) {
    list(
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = 73100 + seed_offset * 100 + chain
    )
  })
  jags_fit <- rjags::jags.model(
    textConnection(make_gaussian_model(predictor_names)),
    data = model_data, inits = inits, n.chains = 4,
    n.adapt = adaptation, quiet = TRUE
  )
  update(jags_fit, burn_in, progress.bar = "none")
  samples <- rjags::coda.samples(
    jags_fit, variable.names = c("b0", "b", "sigma", "mu", "loglik"),
    n.iter = retained, thin = 1, progress.bar = "none"
  )
  sample_matrix <- as.matrix(samples)
  coefficient_columns <- c(
    "b0", if (n_predictors == 1) "b" else paste0("b[", seq_len(n_predictors), "]"),
    "sigma"
  )
  coefficient_names <- c("Intercept", predictor_names, "Residual SD")
  coefficient_summary <- summary(samples)
  rhat <- coda::gelman.diag(samples[, coefficient_columns, drop = FALSE],
                             autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."]
  ess <- coda::effectiveSize(samples[, coefficient_columns, drop = FALSE])
  coefficients <- data.frame(
    model_id = model_id,
    model = model_label,
    parameter = coefficient_names,
    mean = coefficient_summary$statistics[coefficient_columns, "Mean"],
    sd = coefficient_summary$statistics[coefficient_columns, "SD"],
    q2_5 = coefficient_summary$quantiles[coefficient_columns, "2.5%"],
    median = coefficient_summary$quantiles[coefficient_columns, "50%"],
    q97_5 = coefficient_summary$quantiles[coefficient_columns, "97.5%"],
    rhat = as.numeric(rhat), effective_sample_size = as.numeric(ess),
    row.names = NULL, stringsAsFactors = FALSE
  )
  loglik <- sample_matrix[, paste0("loglik[", seq_len(length(response)), "]"), drop = FALSE]
  chain_id <- rep(seq_along(samples), each = nrow(samples[[1]]))
  loo_fit <- loo::loo(loglik, r_eff = loo::relative_eff(exp(loglik), chain_id = chain_id))
  waic_fit <- loo::waic(loglik)
  comparison <- data.frame(
    model_id = model_id, model = model_label,
    predictors = paste(predictor_names, collapse = "; "),
    n_reservoirs = length(response), n_predictors = n_predictors,
    waic = waic_fit$estimates["waic", "Estimate"],
    waic_se = waic_fit$estimates["waic", "SE"],
    p_waic = waic_fit$estimates["p_waic", "Estimate"],
    looic = loo_fit$estimates["looic", "Estimate"],
    looic_se = loo_fit$estimates["looic", "SE"],
    p_loo = loo_fit$estimates["p_loo", "Estimate"],
    max_pareto_k = max(loo::pareto_k_values(loo_fit)),
    n_pareto_k_over_0_7 = sum(loo::pareto_k_values(loo_fit) > 0.7),
    max_rhat = max(rhat), min_effective_sample_size = min(ess),
    stringsAsFactors = FALSE
  )
  fitted <- colMeans(sample_matrix[, paste0("mu[", seq_len(length(response)), "]"), drop = FALSE])
  if (length(reservoir) != length(response)) {
    stop("Reservoir identifiers must match the fitted response length.")
  }
  residuals <- data.frame(
    model_id = model_id, model = model_label,
    reservoir = as.character(reservoir), observed_allochthony = response,
    posterior_mean_fitted = fitted, residual = response - fitted,
    stringsAsFactors = FALSE
  )
  list(jags_fit = jags_fit, samples = samples, comparison = comparison,
       coefficients = coefficients, residuals = residuals, loglik = loglik)
}

if (run_primary_models) {
  peak_runoff_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_peak_runoff,
    "peak_catchment_runoff", model_labels[["peak_catchment_runoff"]],
    reservoir_data$reservoir, seed_offset = 2
  )
  average_runoff_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_average_runoff,
    "average_catchment_runoff", model_labels[["average_catchment_runoff"]],
    reservoir_data$reservoir, seed_offset = 3
  )
  reservoir_null_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_reservoir_null,
    "reservoir_size_null", model_labels[["reservoir_size_null"]],
    reservoir_data$reservoir, seed_offset = 4
  )
  catchment_null_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_catchment_null,
    "catchment_size_null", model_labels[["catchment_size_null"]],
    reservoir_data$reservoir, seed_offset = 5
  )
  potential_production_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_potential_production,
    "potential_aquatic_production", model_labels[["potential_aquatic_production"]],
    reservoir_data$reservoir, seed_offset = 6
  )
  actual_production_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_actual_production,
    "actual_aquatic_production", model_labels[["actual_aquatic_production"]],
    reservoir_data$reservoir, seed_offset = 7
  )
  # Equation 5 uses rounded mean direct allochthony as the predictor.
  equation5_data <- merge(
    reservoir_data[, c("reservoir", "mean_allochthony")],
    mean_direct_allochthony,
    by = "reservoir", sort = FALSE
  )
  equation5_data <- equation5_data[
    match(reservoir_ids, equation5_data$reservoir), , drop = FALSE
  ]
  X_equation5 <- data.frame(
    mean_direct_allochthony = round(
      equation5_data$mean_direct_allochthony, digits = 2
    )
  )
  equation5_fit <- fit_reservoir_model(
    equation5_data$mean_allochthony, X_equation5,
    "overall_direct_allochthony", model_labels[["overall_direct_allochthony"]],
    equation5_data$reservoir, seed_offset = 8, retained = 10000
  )
  primary_fits <- list(
    peak_runoff_fit, average_runoff_fit, reservoir_null_fit,
    catchment_null_fit, potential_production_fit, actual_production_fit
  )
  model_comparison <- do.call(rbind, lapply(primary_fits, function(x) x$comparison))
  model_coefficients <- rbind(
    do.call(rbind, lapply(primary_fits, function(x) x$coefficients)),
    equation5_fit$coefficients
  )
  model_residuals <- do.call(rbind, lapply(primary_fits, function(x) x$residuals))
  model_comparison$delta_looic <- model_comparison$looic - min(model_comparison$looic)
  model_comparison$delta_waic <- model_comparison$waic - min(model_comparison$waic)
  model_comparison <- model_comparison[order(model_comparison$looic), , drop = FALSE]
  write.csv(model_comparison, file.path(output_dir, "model_comparison_results.csv"), row.names = FALSE)
  write.csv(model_coefficients, file.path(output_dir, "model_coefficients_results.csv"), row.names = FALSE)
  write.csv(model_residuals, file.path(output_dir, "model_residuals_results.csv"), row.names = FALSE)
} else {
  primary_fits <- NULL
  equation5_fit <- NULL
  model_comparison <- stored_model_comparison
  model_coefficients <- stored_model_coefficients
  model_residuals <- stored_model_residuals
}

# 9. Model comparison and diagnostics ----------------------------------------

correlation_variables <- c(
  "res_area", "catch_area", "mean_rain", "max_rain", "impervious_surfaces",
  "biovolume", "cyano_prop", "tn", "tp"
)
predictor_correlations <- do.call(rbind, combn(correlation_variables, 2,
  simplify = FALSE, FUN = function(pair) {
    data.frame(predictor_1 = pair[1], predictor_2 = pair[2],
               pearson_r = cor(reservoir_data[[pair[1]]], reservoir_data[[pair[2]]]),
               n = nrow(reservoir_data), stringsAsFactors = FALSE)
  }
))
write.csv(
  predictor_correlations,
  file.path(output_dir, "predictor_correlations_results.csv"),
  row.names = FALSE
)

model_rank_plot <- ggplot(model_comparison,
                          aes(x = reorder(model, looic), y = looic)) +
  geom_point(size = 2.5) + coord_flip() + theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "PSIS-LOOIC")
ggsave(file.path(figure_dir, "model_comparison.png"), model_rank_plot,
       width = 6.5, height = 4.5, dpi = 600)

peak_coefficient_data <- subset(
  model_coefficients,
  model_id == "peak_catchment_runoff" &
    !parameter %in% c("Intercept", "Residual SD")
)
peak_coefficient_plot <- ggplot(peak_coefficient_data,
                                aes(y = parameter, x = median)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = q2_5, xend = q97_5,
                   y = parameter, yend = parameter), linewidth = 0.6) +
  geom_point(size = 2.5) + theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Posterior median (95% credible interval)", y = NULL)
ggsave(file.path(figure_dir, "peak_model_coefficients.png"),
       peak_coefficient_plot, width = 6.5, height = 3.2, dpi = 600)

conditional_peak_plot <- function(fit, data, predictor) {
  sm <- as.matrix(fit$samples)
  peak_predictors <- c("catchment_size", "rain_max", "impervious")
  bcols <- paste0("b[", seq_along(peak_predictors), "]")
  grid <- seq(min(data[[predictor]]), max(data[[predictor]]), length.out = 100)
  fixed <- vapply(peak_predictors, function(name) mean(data[[name]]), numeric(1))
  x_grid <- matrix(rep(fixed, each = length(grid)), nrow = length(grid))
  x_grid[, match(predictor, peak_predictors)] <- grid
  mu <- vapply(seq_along(grid), function(i) {
    sm[, "b0"] + as.vector(sm[, bcols, drop = FALSE] %*% x_grid[i, ])
  }, numeric(nrow(sm)))
  data.frame(
    predictor_value = grid, median = apply(mu, 2, median),
    q2_5 = apply(mu, 2, quantile, 0.025),
    q97_5 = apply(mu, 2, quantile, 0.975),
    predictor = predictor, stringsAsFactors = FALSE
  )
}

# The following component plots are generated when the six models are refit.
if (!is.null(primary_fits)) {
  conditional_data <- do.call(rbind, lapply(
    c("catchment_size", "rain_max", "impervious"),
    function(predictor) conditional_peak_plot(
      peak_runoff_fit, X_peak_runoff, predictor
    )
  ))
  conditional_data$predictor <- factor(
    conditional_data$predictor,
    levels = c("catchment_size", "impervious", "rain_max"),
    labels = c("Catchment size", "Impervious cover", "Peak rainfall")
  )
  conditional_plot <- ggplot(conditional_data,
                             aes(x = predictor_value, y = median)) +
    geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.2) +
    geom_line() + facet_wrap(~ predictor, scales = "free_x") +
    theme_bw() + theme(panel.grid.minor = element_blank()) +
    labs(x = NULL, y = "Predicted reservoir allochthony")
  ggsave(file.path(figure_dir, "Figure4_conditional_relationships.png"),
         conditional_plot, width = 8, height = 3.5, dpi = 600)
}

# 10. NDVI sensitivity analysis ----------------------------------------------

ndvi_data <- reservoir_data
ndvi_data$catchment_size <- catchment_area_scaled
ndvi_data$reservoir_size <- reservoir_area_scaled
ndvi_data$rain_max <- peak_rainfall_scaled
ndvi_data$impervious <- impervious_cover
ndvi_data <- merge(
  ndvi_data,
  ndvi_inputs[, c("reservoir", "weighted_mean_ndvi")],
  by = "reservoir", sort = FALSE
)
ndvi_data <- ndvi_data[match(reservoir_ids, ndvi_data$reservoir), , drop = FALSE]
X_ndvi_only <- data.frame(weighted_mean_ndvi = ndvi_data$weighted_mean_ndvi)
X_vegetation_runoff <- data.frame(
  catchment_size = ndvi_data$catchment_size,
  rain_max = ndvi_data$rain_max,
  weighted_mean_ndvi = ndvi_data$weighted_mean_ndvi
)
X_peak_runoff_reference <- data.frame(
  catchment_size = ndvi_data$catchment_size,
  rain_max = ndvi_data$rain_max,
  impervious = ndvi_data$impervious
)
X_peak_runoff_plus_ndvi <- data.frame(
  catchment_size = ndvi_data$catchment_size,
  rain_max = ndvi_data$rain_max,
  impervious = ndvi_data$impervious,
  weighted_mean_ndvi = ndvi_data$weighted_mean_ndvi
)
if (run_ndvi_models) {
  ndvi_only_fit <- fit_reservoir_model(
    ndvi_data$mean_allochthony, X_ndvi_only, "ndvi_only",
    model_labels[["ndvi_only"]], ndvi_data$reservoir, seed_offset = 101
  )
  vegetation_runoff_fit <- fit_reservoir_model(
    ndvi_data$mean_allochthony, X_vegetation_runoff,
    "vegetation_runoff_alternative",
    model_labels[["vegetation_runoff_alternative"]],
    ndvi_data$reservoir, seed_offset = 102
  )
  peak_runoff_reference_fit <- fit_reservoir_model(
    ndvi_data$mean_allochthony, X_peak_runoff_reference,
    "peak_runoff_reference", model_labels[["peak_runoff_reference"]],
    ndvi_data$reservoir, seed_offset = 103
  )
  peak_runoff_plus_ndvi_fit <- fit_reservoir_model(
    ndvi_data$mean_allochthony, X_peak_runoff_plus_ndvi,
    "peak_runoff_plus_ndvi", model_labels[["peak_runoff_plus_ndvi"]],
    ndvi_data$reservoir, seed_offset = 104
  )
  ndvi_fits <- list(
    ndvi_only_fit, vegetation_runoff_fit, peak_runoff_reference_fit,
    peak_runoff_plus_ndvi_fit
  )
  ndvi_model_comparison <- do.call(
    rbind, lapply(ndvi_fits, function(x) x$comparison)
  )
  ndvi_coefficients <- do.call(
    rbind, lapply(ndvi_fits, function(x) x$coefficients)
  )
  ndvi_model_comparison$delta_looic <- ndvi_model_comparison$looic -
    min(ndvi_model_comparison$looic)
  ndvi_model_comparison$delta_waic <- ndvi_model_comparison$waic -
    min(ndvi_model_comparison$waic)
  ndvi_model_comparison <- ndvi_model_comparison[
    order(ndvi_model_comparison$looic), , drop = FALSE
  ]
  write.csv(
    ndvi_model_comparison,
    file.path(output_dir, "ndvi_model_comparison_results.csv"),
    row.names = FALSE
  )
  write.csv(
    ndvi_coefficients,
    file.path(output_dir, "ndvi_coefficients_results.csv"),
    row.names = FALSE
  )
} else {
  ndvi_fits <- NULL
  ndvi_model_comparison <- stored_ndvi_model_comparison
  ndvi_coefficients <- stored_ndvi_coefficients
}
ndvi_correlations <- stored_ndvi_correlations

# 11. Reservoir:catchment area-ratio sensitivity analysis --------------------

if (run_area_ratio_model) {
  log_area_ratio <- log(reservoir_data$res_area / reservoir_data$catch_area)
  log_area_ratio_scale <- sqrt(
    sum(log_area_ratio^2) / (n_reservoirs - 1)
  )
  scaled_log_area_ratio <- log_area_ratio / log_area_ratio_scale
  X_area_ratio <- data.frame(scaled_log_area_ratio = scaled_log_area_ratio)
  area_ratio_fit <- fit_reservoir_model(
    reservoir_data$mean_allochthony, X_area_ratio,
    "area_ratio_sensitivity", model_labels[["area_ratio_sensitivity"]],
    reservoir_data$reservoir, seed_offset = 200
  )
  ratio_coefficients <- area_ratio_fit$coefficients
  ratio_comparison <- area_ratio_fit$comparison
  write.csv(ratio_coefficients, file.path(output_dir, "area_ratio_coefficients_results.csv"), row.names = FALSE)
  write.csv(ratio_comparison, file.path(output_dir, "area_ratio_model_comparison_results.csv"), row.names = FALSE)
} else {
  area_ratio_fit <- NULL
}

# 12. Isotope-space feasibility assessment -----------------------------------

# The transformed parameters use the final source as a reference value. This
# produces non-negative proportions that sum to one.
source_proportions_from_parameters <- function(parameters) {
  unnormalised <- c(parameters, 0)
  unnormalised <- unnormalised - max(unnormalised)
  exponentiated <- exp(unnormalised)
  exponentiated / sum(exponentiated)
}
predict_consumer_isotopes <- function(source_proportions, source_d13c,
                                      source_d15n, concentration_d13c,
                                      concentration_d15n) {
  c(sum(source_proportions * concentration_d13c * source_d13c) /
      sum(source_proportions * concentration_d13c),
    sum(source_proportions * concentration_d15n * source_d15n) /
      sum(source_proportions * concentration_d15n))
}
calculate_mixture_variance <- function(source_proportions, source_sd_d13c,
                                       source_sd_d15n, concentration_d13c,
                                       concentration_d15n) {
  c(
    sum((source_proportions * concentration_d13c)^2 * source_sd_d13c^2) /
      sum(source_proportions * concentration_d13c)^2,
    sum((source_proportions * concentration_d15n)^2 * source_sd_d15n^2) /
      sum(source_proportions * concentration_d15n)^2
  )
}
calculate_squared_isotope_distance <- function(observed, predicted) {
  sum((observed - predicted)^2)
}

calculate_compatibility_score <- function(observed, predicted, variance) {
  variance[!is.finite(variance) | variance <= 0] <- 1e-8
  sum((observed - predicted)^2 / variance)
}

feasibility_objective <- function(z, obs, se, x, y, c13, c15, sx, sy,
                                  metric = c("euclidean", "chisq")) {
  metric <- match.arg(metric)
  p <- source_proportions_from_parameters(z)
  mu <- predict_consumer_isotopes(p, x, y, c13, c15)
  if (metric == "euclidean") {
    return(calculate_squared_isotope_distance(obs, mu))
  }
  variance <- calculate_mixture_variance(p, sx, sy, c13, c15) + se^2
  calculate_compatibility_score(obs, mu, variance)
}
make_feasibility_starts <- function(k, n_extra = 18L) {
  # A zero-valued transformed start corresponds to equal source proportions.
  starts <- list(rep(0, k - 1L))
  if (k > 1L) {
    starts[[2L]] <- rnorm(k - 1L, 0, 2)
    starts[[3L]] <- rnorm(k - 1L, 0, 4)
    starts <- c(starts, replicate(n_extra, rnorm(k - 1L, 0, 5), simplify = FALSE))
  }
  starts
}
make_preliminary_starts <- function(k) {
  # The draw-level starting parameter set comes from this three-start
  # ordinary-distance search: equal proportions and random scales 2 and 4.
  starts <- list(rep(0, k - 1L))
  if (k > 1L) {
    starts[[2L]] <- rnorm(k - 1L, 0, 2)
    starts[[3L]] <- rnorm(k - 1L, 0, 4)
  }
  starts
}
run_multistart_optimisation <- function(obs, se, x, y, c13, c15, sx, sy,
                                        metric, starts, maxit = 3000) {
  fits <- lapply(starts, function(start) {
    tryCatch(
      optim(start, feasibility_objective, obs = obs, se = se, x = x, y = y,
            c13 = c13, c15 = c15, sx = sx, sy = sy, metric = metric,
            method = "BFGS", control = list(maxit = maxit, reltol = 1e-11)),
      error = function(e) NULL
    )
  })
  values <- vapply(fits, function(fit) {
    if (is.null(fit) || !is.finite(fit$value)) NA_real_ else fit$value
  }, numeric(1))
  finite <- which(is.finite(values))
  if (!length(finite)) {
    return(list(value = NA_real_, par = rep(NA_real_, length(x) - 1L),
                convergence = 99L, p = rep(NA_real_, length(x)),
                pred13 = NA_real_, pred15 = NA_real_, n_finite = 0L,
                start_values = values,
                start_convergence = vapply(fits, function(fit) {
                  if (is.null(fit)) 99L else fit$convergence
                }, integer(1))))
  }
  selected <- finite[which.min(values[finite])]
  fit <- fits[[selected]]
  p <- source_proportions_from_parameters(fit$par)
  predicted <- predict_consumer_isotopes(p, x, y, c13, c15)
  list(value = fit$value, par = fit$par, convergence = fit$convergence,
       p = p, pred13 = predicted[1], pred15 = predicted[2],
       n_finite = length(finite), start_values = values,
       start_convergence = vapply(fits, function(fit) {
         if (is.null(fit)) 99L else fit$convergence
       }, integer(1)))
}
retry_compatibility <- function(old, obs, se, x, y, c13, c15, sx, sy,
                                cutoff) {
  k <- length(x) - 1L
  if (k <= 0L) return(old)
  starts <- c(list(old$par, rep(0, k)),
              lapply(seq_len(8), function(j) 4 * sin(seq_len(k) * j)))
  best <- old
  for (start in starts) {
    fit <- run_multistart_optimisation(obs, se, x, y, c13, c15, sx, sy,
                                  "chisq", list(start), maxit = 10000)
    if (is.finite(fit$value) &&
        (!is.finite(best$value) || fit$value < best$value)) best <- fit
    if (is.finite(fit$value) && fit$value <= cutoff && fit$convergence == 0L) break
  }
  best
}

run_feasibility_assessment <- function(inputs, consumer_records, n_draws = 300L,
                                       cutoff = qchisq(0.95, 2)) {
  set.seed(20260830)
  rows <- vector("list", nrow(consumer_records))
  for (i in seq_len(nrow(consumer_records))) {
    record <- consumer_records[i, , drop = FALSE]
    source <- inputs[inputs$reservoir == record$reservoir &
                       inputs$consumer == record$consumer, , drop = FALSE]
    if (!nrow(source)) stop("No mixing inputs for ", record$reservoir, " / ", record$consumer)
    x <- source$mean_d13c + source$trophic_discrimination_d13c_mean
    y <- source$mean_d15n + source$trophic_discrimination_d15n_mean
    sx <- sqrt(source$sd_d13c^2 + source$trophic_discrimination_d13c_sd^2)
    sy <- sqrt(source$sd_d15n^2 + source$trophic_discrimination_d15n_sd^2)
    c13 <- source$concentration_d13c
    c15 <- source$concentration_d15n
    obs <- c(record$consumer_mean_d13c, record$consumer_mean_d15n)
    se_obs <- c(record$consumer_se_d13c, record$consumer_se_d15n)
    has_consumer_se <- all(is.finite(se_obs))
    se <- if (has_consumer_se) se_obs else c(0, 0)
    starts21 <- make_feasibility_starts(length(x), 18L)
    # Both fixed-input objectives use the inclusive 21-start set. Ordinary
    # distance supplies geometric diagnostics; only compatibility determines
    # compatibility classification.
    euclid <- run_multistart_optimisation(obs, se, x, y, c13, c15, sx, sy,
                                     "euclidean", starts21)
    chisq <- run_multistart_optimisation(obs, se, x, y, c13, c15, sx, sy,
                                    "chisq", starts21)
    preliminary_euclid <- run_multistart_optimisation(
      obs, se, x, y, c13, c15, sx, sy, "euclidean",
      make_preliminary_starts(length(x)))
    draw_euclidean <- rep(NA_real_, n_draws)
    draw_chisq <- rep(NA_real_, n_draws)
    euclidean_failures <- logical(n_draws)
    compatibility_failures_before <- rep(FALSE, n_draws)
    compatibility_failures_after <- rep(FALSE, n_draws)
    compatibility_nonfinite_before <- rep(FALSE, n_draws)
    compatibility_nonfinite_after <- rep(FALSE, n_draws)
    draw_retries <- logical(n_draws)
    original_classification <- rep(NA, n_draws)
    draw_start <- if (all(is.finite(preliminary_euclid$par))) {
      list(preliminary_euclid$par)
    } else {
      list(rep(0, length(x) - 1L))
    }
    for (b in seq_len(n_draws)) {
      x_draw <- rnorm(length(x), source$mean_d13c, source$sd_d13c) +
        rnorm(length(x), source$trophic_discrimination_d13c_mean,
              source$trophic_discrimination_d13c_sd)
      y_draw <- rnorm(length(y), source$mean_d15n, source$sd_d15n) +
        rnorm(length(y), source$trophic_discrimination_d15n_mean,
              source$trophic_discrimination_d15n_sd)
      draw_e <- run_multistart_optimisation(obs, se, x_draw, y_draw, c13, c15,
                                       sx, sy, "euclidean", draw_start)
      euclidean_failures[b] <- draw_e$convergence != 0L || !is.finite(draw_e$value)
      if (is.finite(draw_e$value)) draw_euclidean[b] <- sqrt(draw_e$value)
      if (!has_consumer_se) next
      draw_c <- run_multistart_optimisation(obs, se, x_draw, y_draw, c13, c15,
                                       sx, sy, "chisq", draw_start)
      compatibility_failures_before[b] <- draw_c$convergence != 0L
      compatibility_nonfinite_before[b] <- !is.finite(draw_c$value)
      if (is.finite(draw_c$value)) {
        original_classification[b] <- draw_c$value <= cutoff
      }
      after_retry <- draw_c
      if (draw_c$convergence != 0L || !is.finite(draw_c$value)) {
        draw_retries[b] <- TRUE
        after_retry <- retry_compatibility(draw_c, obs, se, x_draw, y_draw,
                                           c13, c15, sx, sy, cutoff)
      }
      compatibility_failures_after[b] <- after_retry$convergence != 0L
      compatibility_nonfinite_after[b] <- !is.finite(after_retry$value)
      if (is.finite(after_retry$value)) draw_chisq[b] <- after_retry$value
    }
    finite_euclidean <- is.finite(draw_euclidean)
    finite_draws <- is.finite(draw_chisq)
    # This count records any issue in the initial ordinary-distance or
    # compatibility search. Finite compatibility scores are retained after
    # retry, so the final draw proportion uses the available finite scores.
    failure_union <- euclidean_failures | compatibility_failures_before |
      compatibility_nonfinite_before
    finite_original <- is.finite(original_classification)
    classification_changed <- finite_original & finite_draws &
      (original_classification != (draw_chisq <= cutoff))
    rows[[i]] <- data.frame(
      reservoir = record$reservoir, consumer = record$consumer,
      n_source_entries = nrow(source), n_consumer_obs = record$n_consumer_obs,
      n_unknown_producer_sources = record$n_unknown_producer_sources,
      consumer_mean_d13c = obs[1], consumer_mean_d15n = obs[2],
      consumer_sd_d13c = record$consumer_sd_d13c,
      consumer_sd_d15n = record$consumer_sd_d15n,
      consumer_se_d13c = record$consumer_se_d13c,
      consumer_se_d15n = record$consumer_se_d15n,
      point_min_euclidean = sqrt(euclid$value), point_pred_d13c = euclid$pred13,
      point_pred_d15n = euclid$pred15, point_min_chisq = chisq$value,
      point_chisq_reference_tail_area_df2 = pchisq(chisq$value, 2, lower.tail = FALSE),
      point_chisq_95pct = is.finite(chisq$value) && chisq$value <= cutoff,
      uncertainty_draws = n_draws,
      uncertainty_min_euclidean_median = if (any(finite_euclidean))
        median(draw_euclidean[finite_euclidean]) else NA_real_,
      uncertainty_min_euclidean_q025 = if (any(finite_euclidean))
        quantile(draw_euclidean[finite_euclidean], 0.025, names = FALSE) else NA_real_,
      uncertainty_min_euclidean_q975 = if (any(finite_euclidean))
        quantile(draw_euclidean[finite_euclidean], 0.975, names = FALSE) else NA_real_,
      uncertainty_p_euclidean_le_1 = if (any(finite_euclidean))
        mean(draw_euclidean[finite_euclidean] <= 1) else NA_real_,
      uncertainty_p_euclidean_le_2 = if (any(finite_euclidean))
        mean(draw_euclidean[finite_euclidean] <= 2) else NA_real_,
      uncertainty_p_chisq_95pct = if (any(finite_draws))
        mean(draw_chisq[finite_draws] <= cutoff) else NA_real_,
      uncertainty_n_optimizer_failures = sum(failure_union),
      point_optimisation_starts = length(starts21),
      euclidean_nonconvergence_count = sum(euclidean_failures),
      compatibility_nonconvergence_count_original = sum(compatibility_failures_before),
      compatibility_nonconvergence_count_final = sum(compatibility_failures_after),
      nonfinite_compatibility_count_original = sum(compatibility_nonfinite_before),
      nonfinite_compatibility_count_final = sum(compatibility_nonfinite_after),
      n_finite_compatibility_scores = sum(finite_draws),
      n_compatible_draws = if (any(finite_draws))
        sum(draw_chisq[finite_draws] <= cutoff) else 0L,
      n_compatibility_retries = sum(draw_retries),
      n_draw_classification_changes = sum(classification_changed),
      compatible_draw_proportion_original = if (any(finite_original))
        mean(original_classification[finite_original]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

if (run_feasibility_audit) {
  feasibility_results <- run_feasibility_assessment(
    mixing_inputs, feasibility_consumer
  )
  write.csv(feasibility_results,
            file.path(output_dir, "feasibility_consumer_results.csv"),
            row.names = FALSE)
} else {
  feasibility_results <- feasibility_consumer
}

feasibility_summary <- aggregate(
  uncertainty_p_chisq_95pct ~ reservoir, data = feasibility_results,
  FUN = median, na.action = na.pass
)
names(feasibility_summary)[2] <- "median_compatible_draw_proportion"

# 13. Exact leave-one-reservoir-out prediction --------------------------------

prepare_fold <- function(raw, heldout_index) {
  train_index <- setdiff(seq_len(nrow(raw)), heldout_index)
  train_data <- raw[train_index, , drop = FALSE]
  heldout_data <- raw[heldout_index, , drop = FALSE]
  n_training <- nrow(train_data)
  catchment_scale <- sqrt(sum(train_data$catch_area^2) / (n_training - 1))
  reservoir_scale <- sqrt(sum(train_data$res_area^2) / (n_training - 1))
  mean_rain_scale <- sqrt(sum(train_data$mean_rain^2) / (n_training - 1))
  peak_rain_scale <- sqrt(sum(train_data$max_rain^2) / (n_training - 1))
  phyto_scale <- sqrt(sum(train_data$biovolume^2) / (n_training - 1))
  tn_scale <- sqrt(sum(train_data$tn^2) / (n_training - 1))
  tp_scale <- sqrt(sum(train_data$tp^2) / (n_training - 1))
  train_data$catchment_size <- train_data$catch_area / catchment_scale
  heldout_data$catchment_size <- heldout_data$catch_area / catchment_scale
  train_data$reservoir_size <- train_data$res_area / reservoir_scale
  heldout_data$reservoir_size <- heldout_data$res_area / reservoir_scale
  train_data$rain_mean <- train_data$mean_rain / mean_rain_scale
  heldout_data$rain_mean <- heldout_data$mean_rain / mean_rain_scale
  train_data$rain_max <- train_data$max_rain / peak_rain_scale
  heldout_data$rain_max <- heldout_data$max_rain / peak_rain_scale
  train_data$phyto_biovolume <- train_data$biovolume / phyto_scale
  heldout_data$phyto_biovolume <- heldout_data$biovolume / phyto_scale
  train_data$tn_scaled <- train_data$tn / tn_scale
  heldout_data$tn_scaled <- heldout_data$tn / tn_scale
  train_data$tp_scaled <- train_data$tp / tp_scale
  heldout_data$tp_scaled <- heldout_data$tp / tp_scale
  train_data$impervious <- train_data$impervious_surfaces
  heldout_data$impervious <- heldout_data$impervious_surfaces
  train_data$cyanobacteria <- train_data$cyano_prop
  heldout_data$cyanobacteria <- heldout_data$cyano_prop
  list(train = train_data, heldout = heldout_data)
}
fit_exact_fold <- function(train, heldout, predictors, model_id,
                           model_label, seed_offset) {
  values <- lapply(train[predictors], as.numeric)
  names(values) <- predictors
  fit <- fit_reservoir_model(
    train$mean_allochthony, values, model_id, model_label,
    train$reservoir, seed_offset
  )
  x_holdout <- as.numeric(heldout[1, predictors, drop = TRUE])
  extension <- 0L
  initial_lpd <- NA_real_
  predictive_density_ess <- NA_real_
  repeat {
    sm <- as.matrix(fit$samples)
    bcols <- if (length(predictors) == 1) "b" else paste0("b[", seq_along(predictors), "]")
    mu <- sm[, "b0"] + as.vector(sm[, bcols, drop = FALSE] %*% x_holdout)
    sigma <- sm[, "sigma"]
    log_density <- stats::dnorm(heldout$mean_allochthony[[1]], mu, sigma, log = TRUE)
    max_lpd <- max(log_density)
    lpd <- max_lpd + log(mean(exp(log_density - max_lpd)))
    rhat <- coda::gelman.diag(fit$samples[, c("b0", bcols, "sigma"), drop = FALSE],
                              multivariate = FALSE, autoburnin = FALSE)$psrf[, "Point est."]
    ess <- coda::effectiveSize(fit$samples[, c("b0", bcols, "sigma"), drop = FALSE])
    z <- exp(log_density - max_lpd)
    chain_n <- nrow(fit$samples[[1]])
    z_ess <- coda::effectiveSize(coda::mcmc.list(lapply(seq_along(fit$samples), function(j)
      coda::mcmc(z[((j - 1) * chain_n + 1):(j * chain_n)]))))
    lpd_mcse <- stats::sd(z) / sqrt(z_ess) / mean(z)
    if (extension == 0L) initial_lpd <- lpd
    predictive_density_ess <- as.numeric(z_ess)
    if ((max(rhat) <= 1.01 && min(ess) >= 400 && lpd_mcse <= 0.0025) || extension >= 4L) break
    additional <- rjags::coda.samples(
      fit$jags_fit, c("b0", "b", "sigma"), n.iter = 50000,
      thin = 1, progress.bar = "none"
    )
    # Only coefficient and residual-dispersion draws are needed for the
    # held-out prediction. Keeping this common set of columns also allows
    # extension batches to be appended to the initial posterior draws.
    fit$samples <- coda::mcmc.list(lapply(seq_along(fit$samples), function(j) {
      additional_matrix <- as.matrix(additional[[j]])
      additional_matrix <- additional_matrix[, c("b0", bcols, "sigma"), drop = FALSE]
      coda::mcmc(rbind(
        as.matrix(fit$samples[[j]][, c("b0", bcols, "sigma"), drop = FALSE]),
        additional_matrix
      ))
    }))
    extension <- extension + 1L
  }
  predictive_draws <- stats::rnorm(length(mu), mu, sigma)
  predictive_q <- stats::quantile(predictive_draws, c(0.025, 0.975))
  list(lpd = lpd, predictive_mean = mean(mu),
       predictive_q2_5 = unname(predictive_q[1]),
       predictive_q97_5 = unname(predictive_q[2]),
       max_rhat = max(rhat), min_effective_sample_size = min(ess),
       n_draws = nrow(sm), lpd_mcse = lpd_mcse,
       predictive_density_ess = predictive_density_ess,
       initial_lpd = initial_lpd,
       extension_batches = extension)
}
run_exact_loo_refits <- function(raw) {
  primary_model_predictors <- list(
    peak_catchment_runoff = c("catchment_size", "rain_max", "impervious"),
    average_catchment_runoff = c("catchment_size", "rain_mean", "impervious"),
    reservoir_size_null = c("reservoir_size"),
    catchment_size_null = c("catchment_size"),
    potential_aquatic_production = c("reservoir_size", "tn_scaled", "tp_scaled"),
    actual_aquatic_production = c("reservoir_size", "phyto_biovolume", "cyanobacteria")
  )
  pointwise <- list()
  index <- 0L
  for (heldout_index in seq_len(nrow(raw))) {
    fold <- prepare_fold(raw, heldout_index)
    for (model_id in names(primary_model_predictors)) {
      index <- index + 1L
      fit <- fit_exact_fold(
        fold$train, fold$heldout, primary_model_predictors[[model_id]],
        model_id, model_labels[[model_id]], index
      )
      pointwise[[index]] <- data.frame(
        heldout_reservoir = fold$heldout$reservoir,
        model_id = model_id, model = model_labels[[model_id]],
        n_training_reservoirs = nrow(fold$train),
        observed_allochthony = fold$heldout$mean_allochthony,
        exact_log_predictive_density = fit$lpd,
        posterior_predictive_mean = fit$predictive_mean,
        posterior_predictive_q2_5 = fit$predictive_q2_5,
        posterior_predictive_q97_5 = fit$predictive_q97_5,
        max_rhat = fit$max_rhat,
        min_effective_sample_size = fit$min_effective_sample_size,
        posterior_draws = fit$n_draws,
        log_predictive_density_mcse = fit$lpd_mcse,
        predictive_density_ess = fit$predictive_density_ess,
        initial_log_predictive_density = fit$initial_lpd,
        extension_batches = fit$extension_batches,
        stringsAsFactors = FALSE
      )
    }
  }
  pointwise <- do.call(rbind, pointwise)
  comparison <- do.call(rbind, lapply(split(pointwise, pointwise$model_id), function(dat) {
    elpd <- sum(dat$exact_log_predictive_density)
    se <- sqrt(nrow(dat) * var(dat$exact_log_predictive_density))
    lpd_mcse <- sqrt(sum(dat$log_predictive_density_mcse^2, na.rm = TRUE))
    data.frame(model_id = dat$model_id[1], model = dat$model[1],
               n_reservoirs = nrow(dat), exact_elpd = elpd,
               exact_elpd_se = se, exact_looic = -2 * elpd,
               exact_looic_se = 2 * se,
               exact_looic_mcse = 2 * lpd_mcse,
               max_rhat_across_folds = max(dat$max_rhat),
               min_ess_across_folds = min(dat$min_effective_sample_size),
               stringsAsFactors = FALSE)
  }))
  peak_lpd <- pointwise$exact_log_predictive_density[
    pointwise$model_id == "peak_catchment_runoff"
  ]
  peak_ids <- pointwise$heldout_reservoir[
    pointwise$model_id == "peak_catchment_runoff"
  ]
  comparison$delta_exact_looic_se <- vapply(comparison$model_id, function(model_id) {
    if (model_id == "peak_catchment_runoff") return(0)
    other <- pointwise$exact_log_predictive_density[pointwise$model_id == model_id]
    other_ids <- pointwise$heldout_reservoir[pointwise$model_id == model_id]
    order_peak <- match(other_ids, peak_ids)
    2 * sqrt(length(other) * var(other - peak_lpd[order_peak]))
  }, numeric(1))
  comparison$delta_exact_looic <- comparison$exact_looic - min(comparison$exact_looic)
  comparison$rank_exact_looic <- rank(comparison$exact_looic, ties.method = "min")
  comparison <- comparison[, c(
    "model_id", "model", "n_reservoirs", "exact_elpd", "exact_elpd_se",
    "exact_looic", "exact_looic_se", "exact_looic_mcse",
    "max_rhat_across_folds", "min_ess_across_folds", "delta_exact_looic",
    "rank_exact_looic", "delta_exact_looic_se"
  )]
  list(pointwise = pointwise, comparison = comparison[order(comparison$rank_exact_looic), ])
}

if (run_exact_loo) {
  exact_loo <- run_exact_loo_refits(reservoir_data)
  write.csv(exact_loo$pointwise,
            file.path(output_dir, "exact_loo_pointwise_results.csv"),
            row.names = FALSE)
  write.csv(exact_loo$comparison,
            file.path(output_dir, "exact_loo_comparison_results.csv"),
            row.names = FALSE)
} else {
  exact_loo <- list(pointwise = stored_exact_loo_pointwise,
                    comparison = stored_exact_loo_comparison)
}

# 14. Res 4 exclusion sensitivity analysis -----------------------------------

reservoir_data_without_res4 <- reservoir_data[
  reservoir_data$reservoir != "Res 4", , drop = FALSE
]
n_without_res4 <- nrow(reservoir_data_without_res4)
catchment_area_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$catch_area^2) / (n_without_res4 - 1)
)
reservoir_area_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$res_area^2) / (n_without_res4 - 1)
)
mean_rainfall_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$mean_rain^2) / (n_without_res4 - 1)
)
peak_rainfall_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$max_rain^2) / (n_without_res4 - 1)
)
phytoplankton_biovolume_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$biovolume^2) / (n_without_res4 - 1)
)
total_nitrogen_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$tn^2) / (n_without_res4 - 1)
)
total_phosphorus_without_res4_scale <- sqrt(
  sum(reservoir_data_without_res4$tp^2) / (n_without_res4 - 1)
)
catchment_area_without_res4_scaled <-
  reservoir_data_without_res4$catch_area / catchment_area_without_res4_scale
reservoir_area_without_res4_scaled <-
  reservoir_data_without_res4$res_area / reservoir_area_without_res4_scale
mean_rainfall_without_res4_scaled <-
  reservoir_data_without_res4$mean_rain / mean_rainfall_without_res4_scale
peak_rainfall_without_res4_scaled <-
  reservoir_data_without_res4$max_rain / peak_rainfall_without_res4_scale
phytoplankton_biovolume_without_res4_scaled <-
  reservoir_data_without_res4$biovolume /
    phytoplankton_biovolume_without_res4_scale
total_nitrogen_without_res4_scaled <-
  reservoir_data_without_res4$tn / total_nitrogen_without_res4_scale
total_phosphorus_without_res4_scaled <-
  reservoir_data_without_res4$tp / total_phosphorus_without_res4_scale
impervious_without_res4 <- reservoir_data_without_res4$impervious_surfaces
cyanobacteria_without_res4 <- reservoir_data_without_res4$cyano_prop

X_peak_runoff_without_res4 <- data.frame(
  catchment_size = catchment_area_without_res4_scaled,
  rain_max = peak_rainfall_without_res4_scaled,
  impervious = impervious_without_res4
)
X_average_runoff_without_res4 <- data.frame(
  catchment_size = catchment_area_without_res4_scaled,
  rain_mean = mean_rainfall_without_res4_scaled,
  impervious = impervious_without_res4
)
X_reservoir_null_without_res4 <- data.frame(
  reservoir_size = reservoir_area_without_res4_scaled
)
X_catchment_null_without_res4 <- data.frame(
  catchment_size = catchment_area_without_res4_scaled
)
X_potential_production_without_res4 <- data.frame(
  reservoir_size = reservoir_area_without_res4_scaled,
  tn_scaled = total_nitrogen_without_res4_scaled,
  tp_scaled = total_phosphorus_without_res4_scaled
)
X_actual_production_without_res4 <- data.frame(
  reservoir_size = reservoir_area_without_res4_scaled,
  phyto_biovolume = phytoplankton_biovolume_without_res4_scaled,
  cyanobacteria = cyanobacteria_without_res4
)

if (run_res4_exclusion) {
  peak_runoff_without_res4_fit <- fit_reservoir_model(
    reservoir_data_without_res4$mean_allochthony,
    X_peak_runoff_without_res4, "peak_catchment_runoff",
    model_labels[["peak_catchment_runoff"]],
    reservoir_data_without_res4$reservoir, seed_offset = 301
  )
  average_runoff_without_res4_fit <- fit_reservoir_model(
    reservoir_data_without_res4$mean_allochthony,
    X_average_runoff_without_res4, "average_catchment_runoff",
    model_labels[["average_catchment_runoff"]],
    reservoir_data_without_res4$reservoir, seed_offset = 302
  )
  reservoir_null_without_res4_fit <- fit_reservoir_model(
    reservoir_data_without_res4$mean_allochthony,
    X_reservoir_null_without_res4, "reservoir_size_null",
    model_labels[["reservoir_size_null"]],
    reservoir_data_without_res4$reservoir, seed_offset = 303
  )
  catchment_null_without_res4_fit <- fit_reservoir_model(
    reservoir_data_without_res4$mean_allochthony,
    X_catchment_null_without_res4, "catchment_size_null",
    model_labels[["catchment_size_null"]],
    reservoir_data_without_res4$reservoir, seed_offset = 304
  )
  potential_production_without_res4_fit <- fit_reservoir_model(
    reservoir_data_without_res4$mean_allochthony,
    X_potential_production_without_res4, "potential_aquatic_production",
    model_labels[["potential_aquatic_production"]],
    reservoir_data_without_res4$reservoir, seed_offset = 305
  )
  actual_production_without_res4_fit <- fit_reservoir_model(
    reservoir_data_without_res4$mean_allochthony,
    X_actual_production_without_res4, "actual_aquatic_production",
    model_labels[["actual_aquatic_production"]],
    reservoir_data_without_res4$reservoir, seed_offset = 306
  )
  res4_exclusion_fits <- list(
    peak_runoff_without_res4_fit, average_runoff_without_res4_fit,
    reservoir_null_without_res4_fit, catchment_null_without_res4_fit,
    potential_production_without_res4_fit, actual_production_without_res4_fit
  )
  res4_exclusion_comparison <- do.call(
    rbind, lapply(res4_exclusion_fits, function(x) x$comparison)
  )
  res4_exclusion_coefficients <- do.call(
    rbind, lapply(res4_exclusion_fits, function(x) x$coefficients)
  )
  res4_exclusion_comparison$delta_looic <-
    res4_exclusion_comparison$looic - min(res4_exclusion_comparison$looic)
  res4_exclusion_comparison$delta_waic <-
    res4_exclusion_comparison$waic - min(res4_exclusion_comparison$waic)
  res4_exclusion_comparison <- res4_exclusion_comparison[
    order(res4_exclusion_comparison$looic), , drop = FALSE
  ]
  write.csv(
    res4_exclusion_comparison,
    file.path(output_dir, "res4_exclusion_model_comparison.csv"),
    row.names = FALSE
  )
  write.csv(
    res4_exclusion_coefficients,
    file.path(output_dir, "res4_exclusion_coefficients.csv"),
    row.names = FALSE
  )
} else {
  res4_exclusion_fits <- NULL
}

# 15. Fish-biomass analysis ---------------------------------------------------

fit_biomass_model <- function(reservoir_data) {
  dat <- reservoir_data[, c("reservoir", "mean_allochthony", "res_area",
                            "catch_area", "fish_biomass_top3_mean")]
  names(dat) <- c("reservoir", "allochthony", "reservoir_area",
                  "catchment_area", "fish_biomass_top3_mean")
  if (anyNA(dat)) stop("Missing values remain in the fish-biomass records.")
  n_biomass_reservoirs <- nrow(dat)
  fish_biomass_scale <- sqrt(
    sum(dat$fish_biomass_top3_mean^2) / (n_biomass_reservoirs - 1)
  )
  reservoir_area_biomass_scale <- sqrt(
    sum(dat$reservoir_area^2) / (n_biomass_reservoirs - 1)
  )
  catchment_area_biomass_scale <- sqrt(
    sum(dat$catchment_area^2) / (n_biomass_reservoirs - 1)
  )
  dat$fish_biomass <- dat$fish_biomass_top3_mean / fish_biomass_scale
  dat$reservoir_size <- dat$reservoir_area / reservoir_area_biomass_scale
  dat$catchment_size <- dat$catchment_area / catchment_area_biomass_scale
  model_code <- "
  model {
    for (i in 1:N) {
      fish_biomass[i] ~ dnorm(mu[i], tau)
      mu[i] <- b0 + b[1] * reservoir_size[i] +
        b[2] * catchment_size[i] + b[3] * allochthony[i]
    }
    b0 ~ dnorm(0, 0.001)
    for (m in 1:3) {
      b[m] ~ dnorm(0, 0.001)
    }
    sigma ~ dunif(0, 1)
    tau <- pow(sigma, -2)
  }
  "
  inits <- lapply(seq_len(4), function(chain) list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 731900 + chain))
  fit <- rjags::jags.model(textConnection(model_code), data = list(
    fish_biomass = dat$fish_biomass, reservoir_size = dat$reservoir_size,
    catchment_size = dat$catchment_size, allochthony = dat$allochthony,
    N = nrow(dat)), inits = inits, n.chains = 4, n.adapt = 10000, quiet = TRUE)
  update(fit, 100000, progress.bar = "none")
  samples <- rjags::coda.samples(fit, c("b0", "b", "sigma"), n.iter = 10000,
                                thin = 1, progress.bar = "none")
  stats <- summary(samples)
  sm <- as.matrix(samples)
  params <- c("b0", "b[1]", "b[2]", "b[3]", "sigma")
  out <- data.frame(parameter = c("Intercept", "Reservoir size", "Catchment size", "Allochthony", "Residual SD"),
                    mean = stats$statistics[params, "Mean"], sd = stats$statistics[params, "SD"],
                    q2_5 = stats$quantiles[params, "2.5%"], median = stats$quantiles[params, "50%"],
                    q97_5 = stats$quantiles[params, "97.5%"],
                    rhat = coda::gelman.diag(samples, multivariate = FALSE, autoburnin = FALSE)$psrf[params, "Point est."],
                    effective_sample_size = coda::effectiveSize(samples)[params])
  coefficient_plot <- ggplot(out[out$parameter != "Residual SD", , drop = FALSE],
                             aes(y = parameter, x = median)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_segment(aes(x = q2_5, xend = q97_5,
                     y = parameter, yend = parameter), linewidth = 0.6) +
    geom_point(size = 2.5) + theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Posterior median (95% credible interval)", y = NULL)
  ggsave(file.path(figure_dir, "SupplementaryFigureS2_biomass_coefficients.png"),
         coefficient_plot, width = 6.5, height = 3.5, dpi = 600)
  grid <- seq(min(dat$allochthony), max(dat$allochthony), length.out = 100)
  mean_reservoir_size <- mean(dat$reservoir_size)
  mean_catchment_size <- mean(dat$catchment_size)
  b <- sm[, c("b[1]", "b[2]", "b[3]")]
  conditional_mu <- vapply(grid, function(value) {
    sm[, "b0"] + b[, 1] * mean_reservoir_size +
      b[, 2] * mean_catchment_size + b[, 3] * value
  }, numeric(nrow(sm)))
  conditional <- data.frame(
    allochthony = grid, median = apply(conditional_mu, 2, median),
    q2_5 = apply(conditional_mu, 2, quantile, 0.025),
    q97_5 = apply(conditional_mu, 2, quantile, 0.975)
  )
  biomass_plot <- ggplot(dat, aes(x = allochthony, y = fish_biomass)) +
    geom_ribbon(data = conditional,
                aes(x = allochthony, ymin = q2_5, ymax = q97_5), alpha = 0.2,
                inherit.aes = FALSE) +
    geom_line(data = conditional, aes(x = allochthony, y = median), linewidth = 0.7) +
    geom_point(size = 2.5) + theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Mean food-web allochthony", y = "Scaled fish biomass")
  ggsave(file.path(figure_dir, "SupplementaryFigureS2_biomass_relationship.png"),
         biomass_plot, width = 7, height = 5, dpi = 600)
  write.csv(out, file.path(output_dir, "biomass_model_coefficients_results.csv"), row.names = FALSE)
  out
}

if (run_biomass_model) {
  biomass_coefficients <- fit_biomass_model(reservoir_data)
} else {
  biomass_coefficients <- NULL
}

# 16. Figures and exported results --------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
invisible(list(
  reservoir_data = reservoir_data,
  taxon_allochthony = taxon_allochthony,
  model_comparison = model_comparison,
  model_coefficients = model_coefficients,
  model_residuals = model_residuals,
  ndvi_model_comparison = ndvi_model_comparison,
  ndvi_coefficients = ndvi_coefficients,
  feasibility_results = feasibility_results,
  feasibility_summary = feasibility_summary,
  exact_loo = exact_loo
))
