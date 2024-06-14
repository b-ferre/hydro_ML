
## note: throughout all of this, B(s, t) is assumed to be normalized such that all coeffs sum to 100                                                            # nolint
## the "shapes" of the hydrographs are what I aim to cluster on.
## FIXME: check w Dr. Ameli or Joe on whether this is a reasonable assumption

## plan:
###### 1. cluster all catchments based on full time-series inference on b(s, t) to get "true" universal hydrographs                                             # nolint
###### 2. cluster b(s, t) estimations based on first chunk using "true" universal hydrographs as "means"                                                        # nolint
###### 3. cluster b(s, t) estimations based on latter chunk using "true" universal hydrographs as "means"                                                       # nolint
###### 4. do sub-clustering within each cluster from (2), holding the cluster label fixed and then trying to find the average 4/5 ways that                     # nolint


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
    if (i %% 100 == 0) {print(paste0(round(i / n, digits = 3), "% done..."))}
}
####

#### DO CLUSTERING & SAVE RESULTS
library(clustlearn)
n_universal_hydrographs <- 8
clustering <- kmeans(normed_b_preds,
                    centers = n_universal_hydrographs,
                    max_iterations = 1000, details = TRUE)

assignments <- as.data.frame(clustering$cluster)
save(assignments, file = "./results/clustering_results/grid_assignment.rda")
for (i in 1:n_universal_hydrographs) {
    save(as.numeric(clustering$centers[i, ]), file = paste0("./results/clustering_results/basis_hydrographs/", "basis_graph_", i, ".rda"))              # nolint
}
####
