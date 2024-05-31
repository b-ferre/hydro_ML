## THIS SCRIPT IS MEANT TO BE USED TO EFFECTIVELY AND ACCURATELY TEST THE
## PERFORMANCE OF JIM KIRTSCHNER'S ARMA IRF MODEL USING BOTH FORWARD AND
## BACKWARD (currently N/A) VALIDATION TESTING ON A SINGLE CATCHMENT

## load packages
library(here)
library(data.table)             # data.table used inside model training
library(readr)                  # read_csv
library(hydroGOF)               # KGE

## load Jim's model code
source(here("kirschner_ar_stuff", "irfnnhs.r"))

## set global vars that control flow of script
catchment_no <- 2631
max_lag <- 365                    # maximum lag to consider when fitting model
train_size <- 0.90                # number >0, <1; determines training set size

## load catchment data
catchment_data <- read_csv(
    paste(here(".", "data", "raw_data"), "/", catchment_no, ".csv", sep = ""))
n <- nrow(catchment_data)

## split data into FORWARD training/validation sets
train_set <- catchment_data[1:floor(train_size * n), ]
validation_set <- catchment_data[ceiling(train_size * n):n, ]

## train IRF model on training set
# - NO FORWARD DIFFERENCING (ADDRESSES NON_STATIONARY BEHAVIOR) FOR NOW
# - complete = FALSE because I do not want to artificially force the model
#   to make the IRF coeffs converge to zero if they don't
mod <- IRF(y = train_set$"Observed_Streamflow_mm/d",
        x = train_set$"Forcing_Precipitation_mm/d",
        nu = 0, FD = 1, m = max_lag, verbose = TRUE, h = NULL, complete = FALSE)

## following equation (20), rebuild b_i coeff's from IRF and AR coeffs
reg_coeffs <- mod$IRF
auto_coeffs <- mod$phi
h <- mod$h
conv_mtx <- toeplitz(c(1, -1 * auto_coeffs, rep(0, (max_lag))))
conv_mtx[upper.tri(conv_mtx)] <- 0
conv_vec <- c(reg_coeffs, rep(0, h))

print(dim(conv_mtx))
print(max_lag)
print(h)
print(dim(conv_vec))

b <- conv_mtx %*% conv_vec

## predict validation set using b_i's and AR coeffs (eqn (15))
#   - ASSUMING alpha = 0 here, as I can't see how to extract a
#   nonzero alpha from the IRF model output. Just focusing on getting
#   the rough mechanics of prediction working for now, but this assumption
#   needs to be corrected or justified before I can move forward with rigor
pred_mtx <- matrix(nrow = nrow(validation_set) - (max_lag + h),
            ncol = (max_lag + (2 * h) + 1))
## FILL PRED_MTX COLUMN BY COLUMN (OPTIMIZE THIS LATER)
for (i in (1:nrow(pred_mtx))) {                                                                   # nolint
    j <- i + h + max_lag
    pred_mtx[i, ] <- c(rev(validation_set$"Forcing_Precipitation_mm/d"[(j - (max_lag + h)):j]),   # nolint
                        rev(validation_set$"Observed_Streamflow_mm/d"[(j - h):(j - 1)]))              # nolint
}
pred_vec <- c(b, auto_coeffs)     # SEE ASSUMPTION ABOVE
y_pred <- pred_mtx %*% pred_vec

## compare predicted y values to real y values
# (note truncation as y_i : i < max_lag + h is not predictable using this model)
y_true <- validation_set$"Observed_Streamflow_mm/d"
y_true <- as.matrix(y_true[(max_lag + h + 1):nrow(validation_set)])

kge <- KGE(y_pred, y_true)