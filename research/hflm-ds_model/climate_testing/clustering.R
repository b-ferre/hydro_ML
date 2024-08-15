### GLOBAL COMPUTE_CANADA PARAMS
num_cores <- 10 - 1

## note: throughout all of this, B(s, t) is assumed to be normalized such that all coeffs sum to 100                                                            # nolint
## the "shapes" of the hydrographs are what I aim to cluster on.
## FIXME: check w Dr. Ameli or Joe on whether this is a reasonable assumption

## plan:
###### 1. cluster all catchments based on full time-series inference on b(s, t) to get "true" universal hydrographs                                             # nolint
###### 2. cluster b(s, t) estimations based on first chunk using "true" universal hydrographs as "means"                                                        # nolint
###### 3. cluster b(s, t) estimations based on latter chunk using "true" universal hydrographs as "means"                                                       # nolint
###### 4. do sub-clustering within each cluster from (2), holding the cluster label fixed and then trying to find the average 4/5 ways that                     # nolint


##############################################################################
##                     CUSTOM K-MEANS USING KGE                             ##
##############################################################################
library(parallel)
library(hydroGOF)

gpu_kge <- function(observed, simulated) {
# Function to calculate Kling Gupta Efficiency (KGE)
  if (!is.numeric(simulated) | !is.numeric(observed)) {
    stop("Both simulated and observed data must be numeric vectors.")
  }
  
  # Convert vectors to GPU vectors
  sim_gpu <- gpuVector(simulated, type = "float")
  obs_gpu <- gpuVector(observed, type = "float")
  
  # Mean of observed data
  obs_mean <- mean(observed)
  
  # Calculate components
  r <- cor(sim_gpu, obs_gpu)
  alpha <- sd(sim_gpu) / sd(obs_gpu)
  beta <- mean(sim_gpu) / mean(obs_gpu)
  
  # Calculate KGE
  kge <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
  
  return(kge)
}

symmetric_kge <- function(vec1, vec2) {
    return((0.5 * (KGE(vec1, vec2) + KGE(vec2, vec1))))
}

custom_kmeans_parallel <- function(data, k, custom_loss_function, cl, max_iter = 100, centers = NA, tol = 1e-4) {                                          # nolint
  row_names <- rownames(data)

  # Initialize cluster centers (randomly select k points from data)
  if (length(dim(centers)) != 2) {
    set.seed(42)
    centers <- data[sample(seq_len(nrow(data)), k), ]
  }

# Function to assign points to the nearest center based on the custom loss function                                                                             # nolint
  assign_clusters <- function(data, centers) {                                                                                                        # nolint
    assignments <- parSapply(cl, seq_len(nrow(data)), function(i) {
      losses <- sapply(seq_len(nrow(centers)), function(j) symmetric_kge(data[i, ], centers[j, ]))                                                       # nolint
      which.min(losses)
    })
    assignments <- as.data.frame(assignments)
    rownames(assignments) <- row_names
    return(assignments)
  }

  # Function to update centers
  update_centers <- function(data, assignments, k) {
    new_centers <- sapply(1:k, function(j) {
      cluster_points <- data[assignments == j, , drop = FALSE]
      colMeans(cluster_points, na.rm = TRUE)
    })
    return(t(new_centers))
  }

  for (iter in 1:max_iter) {
    # Assignment step
    print(paste0("> calculating assignments for the ", iter, "th time..."))
    assignments <- assign_clusters(data, centers)

    # Update step
    print(paste0("> updating centers for the ", iter, "th time..."))
    new_centers <- update_centers(data, assignments, k)

    # Check for convergence
    if (sum((centers - new_centers)^2, na.rm = TRUE) < tol) {
      cat("Converged in", iter, "iterations.\n")
      break
    }

    centers <- new_centers
  }
  print("done!")

  return(list(centers = centers, assignments = assignments))
}

##############################################################################
##                                  1                                       ##
##############################################################################

#### PREP DATA
library(stringr)
extract_gridcode <- function(filename) {
    return(as.integer(str_remove_all(filename, "[^\\d]")))
}

to_test <-list.files("./results/normalized_raw_estimates/full_series", full.names = TRUE)                                                                               # nolint
n <- length(to_test)

normed_b_preds <- matrix(nrow = n, ncol = 365 * 150)        ## hard-coded matrix size based on M, max_lag used when estimating b(s, t)                          # nolint
gridcodes <- as.numeric(lapply(to_test, extract_gridcode))
rownames(normed_b_preds) <- gridcodes

print("building normalized beta(s, t) matrix")

for (i in 1:n) {
    pred_file <- to_test[i]
    b_pred <- as.numeric(read.csv(pred_file)$pred_beta)
    normed_b_pred <- b_pred / sum(b_pred, na.rm = TRUE)
    normed_b_preds[i, ] <- normed_b_pred
    if (i %% 100 == 0) {print(paste0(round(i / n * 100, digits = 3), "% done..."))}
}
####

#### INIT CLUSTER & DO (PARALLELIZED) CLUSTERING & SAVE RESULTS
print(paste(detectCores(), "cores detected..."))
print(paste0("initializing a ", num_cores, "-core cluster..."))
cl <- makeCluster(num_cores)   ##TODO: hardcode depending on clustering_run.sh number requested                            # nolint
clusterExport(cl, c("symmetric_kge", "KGE"))

library(clustlearn)
n_universal_hydrographs <- 8
n_iter_per_round <- 10
n_rounds <- 100

print(paste0("running ", n_rounds, " x ", n_iter_per_round, "-iteration epochs of k-means..."))                        # nolint
t0 <- Sys.time()

## load best (to date) centers from memory if they exist
if ("best_centers.rda" %in% list.files("./results/clustering_results")) {
  load("./results/clustering_results/best_centers.rda")
} else {
  best_centers <- NA
}

for (q in 1:n_rounds) {
    clustering <- custom_kmeans_parallel(data = normed_b_preds,
                                        k = n_universal_hydrographs,
                                        centers = best_centers,
                                        custom_loss_function = symmetric_kge,
                                        max_iter = n_iter_per_round,
                                        cl = cl)

    best_centers <- clustering$centers
    save(best_centers, file = "./results/clustering_results/best_centers.rda")

    print(paste0("finished with epoch ", q, " of ", n_rounds,
                ";  estimated time remaining : ", round((n_rounds - q) * ((Sys.time() - t0) / q))))                 # nolint
}
print("done!")

stopCluster(cl)

assignments <- as.data.frame(clustering$assignments)
save(assignments, file = "./results/clustering_results/grid_assignment.rda")
save(clustering, file = "./results/clustering_results/clustering.rda")
for (j in 1:n_universal_hydrographs) {
    b_pred <- best_centers[j, ]
    save(b_pred, file = paste0("./results/clustering_results/basis_hydrographs/", "b_pred_basis_", j, ".rda"))              # nolint
}
####

### PLOTTING (only do this locally)
do_plotting <- FALSE
if (do_plotting) {
    library(gridExtra)
    require(leaflet)
    require(leafpop)
    require(ggplot2)
    require(sf)
    require(here)
    require(stringr)
    require(jpeg)
    require(htmltools)
    require(base64enc)
    require(leaflet)
    require(leafpop)
    source("./cc_testing_funs.R")

    load("./results/clustering_results/clustering.rda")
    basis_hydrographs <- list.files("./results/clustering_results/basis_hydrographs", full.names = TRUE)
    basis_plots <- list()
    for (i in 1:length(basis_hydrographs)) {
        load(basis_hydrographs[i])
        b_plot <- plot_bvec_better(b_pred, M = 364, max_lag <- 150, title = paste0("basis function ", i))
        basis_plots[[i]] <- b_plot
    }

    do.call("grid.arrange", c(basis_plots, ncol = 2))
    ### MANUALLY SAVE OUTPUT FROM FUNCTION ABOVE

    map <- function() {
        load("./atlas.Rda")

        ## pre-process data for mapping
        locs <- data.frame(index = atlas$gridcode,
                            lng = atlas$longitude,
                            lat = atlas$latitude)

        ## trim locs such that only mappable 
        good_locs <- locs[locs$index %in% mappable, ]

        ## build map
        pb <- txtProgressBar(min = 1, max = nrow(good_locs), style = 3, width = 150, char = "~")                        # nolint
        map <- leaflet()    %>%
        addTiles()

        for (i in 1:nrow(good_locs)) {
            loc <- good_locs[i, ]
            suppressMessages({
            map <- addCircleMarkers(map, data = loc,
                        radius = 3,
                        group = paste(loc$index),
                        color = brewer.pal(n = 8, name = "Dark2")[assignments[paste0(loc$index), "cluster_id"]],
                        opacity = 1)
            map <- addPopupGraphs(map, graph = basis_plots[assignments[paste0(loc$index), "cluster_id"]],
                                group = paste(loc$index), width = 600, height = 600)})
            setTxtProgressBar(pb, i)
        }

        map <- addLegend(
            map = map,
            position = "topright",  # Choose the position
            colors = brewer.pal(n = n_universal_hydrographs, name = "Dark2"),
            labels = 1:n_universal_hydrographs,
            title = "Basis B(s, t) Functions")

        print("done!")

    return(map)
    }

    library(htmlwidgets)
    mappable <- rownames(assignments)
    map_ <- map()
    saveWidget(map_, "clustering_map.html")
}