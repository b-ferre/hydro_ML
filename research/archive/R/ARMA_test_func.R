## THIS SCRIPT IS MEANT TO BE USED TO EFFECTIVELY AND ACCURATELY TEST THE
## PERFORMANCE OF JIM KIRTSCHNER'S ARMA IRF MODEL USING BOTH FORWARD AND
## BACKWARD (currently N/A) VALIDATION TESTING ON A SINGLE CATCHMENT
validate_IRF <- function(catchment_no,
                        max_lag = 365,
                        complete = TRUE,
                        FD = FALSE,
                        train_size = 0.9,
                        dir_flag = 1,
                        random_subsample_size = NULL,
                        save_plot = TRUE) {
## load catchment data
catchment_data <- read_csv(
    paste(here(".", "data", "raw_data"), "/", catchment_no, ".csv", sep = ""),
    show_col_types = FALSE)
n <- nrow(catchment_data)

# if random_subsample_size != NULL, grab a semi-random, contiguous subsample of
# size m = (random_subsample_size * n).
# NOTE: clearly 0 < random_subasmple_size < 1.
# also important to note that using the word "random" here is a bit abusive -
# especially if m is "close" to 1, however simultaneously as m approaches 0,
# the robustness of the test on our selected random subsample goes down.
# CAUTION: it is obvious that m > 0, however more specifically...
# m * min(train_size, 1-train_size) > max_lag in order for the smaller of the
# two sets (train/validation) to have enough data in it to actually run testing.
if (!is.null(random_subsample_size)) {
m <- n * random_subsample_size
    if (m > n || (m * min(train_size, 1 - train_size) <= 3 * max_lag)) {
        print("ERROR: incompatible train_size and random_subsample_size;
        continuing testing with no random subsampling...")
    }
    else {
        rand <- ceiling(runif(1, 1, n - m))
        catchment_data <- catchment_data[rand:(rand+m), ]
        n <- m
    }
}

## split data into training/validation sets
if (dir_flag == 1) {
    train_set <- catchment_data[1:floor(train_size * n), ]
    validation_set <- catchment_data[ceiling(train_size * n):n, ]
}
if (dir_flag == -1) {
    train_set <- catchment_data[(n - floor(train_size * n)):n, ]
    validation_set <- catchment_data[1:ceiling(train_size * n), ]
}


## train IRF model on training set
# - NO FORWARD DIFFERENCING (ADDRESSES NON_STATIONARY BEHAVIOR) FOR NOW
# - complete = TRUE is forced because even I want the IRF model to report
#  ALL of the lags it is using for calculation
mod <- IRF(y = train_set$"Observed_Streamflow_mm/d",
        x = train_set$"Forcing_Precipitation_mm/d",
        nu = 0, FD = 0, m = max_lag, verbose = FALSE, h = NULL, complete = complete)

## following equation (20), rebuild b_i coeff's from IRF and AR coeffs
reg_coeffs <- mod$IRF
auto_coeffs <- mod$phi
h <- mod$h
if (h == 0) {
    auto_coeffs <- numeric(0)
}
conv_mtx <- toeplitz(c(1, -1 * auto_coeffs, rep(0, (max_lag))))
conv_mtx[upper.tri(conv_mtx)] <- 0
conv_vec <- as.matrix(c(reg_coeffs, rep(0, h)))
b <- conv_mtx %*% conv_vec

## predict validation set using b_i's and AR coeffs (eqn (15))
#   - ASSUMING alpha = 0 here, as I can't see how to extract a
#   nonzero alpha from the IRF model output. Just focusing on getting
#   the rough mechanics of prediction working for now, but this assumption
#   needs to be corrected or justified before I can move forward with rigor
pred_mtx <- matrix(nrow = nrow(validation_set) - (max_lag + h),
            ncol = (max_lag + (2 * h) + 1))

## FILL PRED_MTX COLUMN BY COLUMN (OPTIMIZE THIS LATER)
for (i in (1:nrow(pred_mtx))) {                                                                                         # nolint
    j <- i + h + max_lag
    if (h == 0) {
        pred_mtx[i, ] <- c(rev(validation_set$"Forcing_Precipitation_mm/d"[(j - (max_lag + h)):j]))                    # nolint
    } else {
        pred_mtx[i, ] <- c(rev(validation_set$"Forcing_Precipitation_mm/d"[(j - (max_lag + h)):j]),         # nolint
                        rev(validation_set$"Observed_Streamflow_mm/d"[(j - h):(j - 1)]))  ## FIXME: bruh I leaked test set data into my predictions. delete this or fix it but ignore all results from it until then.                   # nolint
    }
}
pred_vec <- c(b, auto_coeffs)     # SEE ASSUMPTION ABOVE
y_pred <- pred_mtx %*% pred_vec

## compare predicted y values to real y values
# (note truncation as y_i : i < max_lag + h is not predictable using this model)
y_true <- validation_set$"Observed_Streamflow_mm/d"
y_true <- as.matrix(y_true[(max_lag + h + 1):nrow(validation_set)])

## get kge for this validation test
kge <- KGE(y_pred, y_true)

## plot the model kernel (trained on train_set) if requested
if (save_plot) {
    pdf(here("data", "results", "ARMA_model", paste("kernel_", catchment_no, ".pdf", sep = "")))       # nolint
    plot(mod$IRF)
    dev.off()
}

neg_flag <- length(mod$IRF[mod$IRF < 0]) != 0

if (!neg_flag) {
    print(mod$IRF)
}

return(c(kge, h, neg_flag))
}

robust_validate_IRF <- function(catchment_no,
                                max_lag = 100,
                                complete = TRUE,
                                train_size = 0.9,
                                dir_flag = 1,
                                random_subsample_size = 0.6,
                                folds = 3) {
    results <- as.data.frame(matrix(data = NA, nrow = folds, ncol = 2))
    colnames(results) <- c("kge", "h")
    for (i in 1:folds) {
        results[i, ] <- validate_IRF(catchment_no = catchment_no,
                                    max_lag = max_lag,
                                    complete = complete,
                                    train_size = train_size,
                                    dir_flag = dir_flag,
                                    random_subsample_size = random_subsample_size,                      # nolint
                                    save_plot = FALSE)
    }

    summary <- data.frame(
        folds = folds,
        problems = length(results$kge[is.na(results$kge)]),
        median_kge = median(results$kge, na.rm = TRUE),
        mean_kge = mean(results$kge, na.rm = TRUE),
        median_AR_order = median(results$h, na.rm = TRUE))

    print(catchment_no)
    print(summary)
    return(summary)
}

plot_IRF <- function(catchment_no, max_lag = 100, complete = TRUE) {
    catchment_data <- read_csv(
    paste(here(".", "data", "raw_data"), "/", catchment_no, ".csv", sep = ""),
    show_col_types = FALSE)
    mod <- IRF(y = catchment_data$"Observed_Streamflow_mm/d",
        x = catchment_data$"Forcing_Precipitation_mm/d",
        nu = 0, FD = 0, m = max_lag, verbose = FALSE, h = NULL, complete = complete)
    
    print(mod$IRF)
    
    pdf(here("data", "results", "ARMA_model", paste("kernel_", catchment_no, ".pdf", sep = "")))       # nolint
    plot(mod$IRF)
    dev.off()
}