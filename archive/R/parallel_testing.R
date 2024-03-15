library(parallel)
library(foreach)
library(doParallel)
library(testthat)
library(SlidingWindowReg)

print("runnning sequentially...")
t_0 <- Sys.time()
    s_out <- foreach(i = 1:2) %do% {
        Sys.sleep(1)
        data.frame(sqrts = sqrt(i * (1:100)))
    }
print(Sys.time() - t_0)

print("setting up parallelization...")
print(paste("detected cores:", detectCores()))
print("registering cluster...")
cl <- makeCluster(2)
registerDoParallel(cl, cores = 2)
expect_equal(2, getDoParWorkers())

print("runnning in parallel...")
t_0 <- Sys.time()
    p_out <- foreach(i = 1:2, .packages = "SlidingWindowReg", .combine = ) %dopar% {
        Sys.sleep(1)
        data.frame(sqrts = sqrt(i * (1:100)))
    }
print(Sys.time() - t_0)