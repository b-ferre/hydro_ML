import numpy as np

## restructuring Joe's raw code/functions into a callable model interface (in python)
class hflmds_model:

    ## model parameters (dictated by problem)
    T = None   ## cycle length; for our purposes it should be 365
    D = None   ## maximum considered lag - input variable from >D time-steps ago is assumed to have no impact on current target value

    ## hyperparameter-related model parameters
    D_h, D_v = None             ## matrices used for smoothness regularization
    wh_opt, wv_opt = None       ## hyper-parameters which control smoothing regularization on b(s, t)

    def format_data(self, X, Y):
        ## check for full cycles
        assert X.shape[0] % self.T == 0, "length(X) is not a multiple of cycle length."
        assert Y.shape[0] == X.shape[0], "length(Y) is not a multiple of cycle length."
        self.n = int(X.shape[0] / self.T)

        ## reformat both input and target time-series into matrix such that each row is one full cycle
        X_new, Y_new = X.reshape((self.n, self.T)), Y.reshape((self.n, self.T))

        ## remove cyclical signals such that mean for each time-step is zero across cycles
        X_new = X_new - np.mean(X_new, axis = 1)[:, np.newaxis]
        Y_new = Y_new - np.mean(Y_new, axis = 1)[:, np.newaxis]
        assert np.sum(np.mean(X_new, axis = 1)) == 0
        assert np.sum(np.mean(Y_new, axis = 1)) == 0

        return X_new, Y_new

    def build_norm_Mats(self):
        ## build D_h
        diff_vec = np.roll(np.eye(self.K)[0], 2)
        diff_vec[0] = -1
        D_h = np.tile(diff_vec, (self.K, 1))
        for i in range(self.K):
            D_h[i] = np.roll(D_h[i], i)
        self.D_h = D_h.astype(int)

        ## build D_v
        D_v = (-1 * (np.eye(self.K) - np.tri(self.K, k = 1) + np.tri(self.K, k = 0))).astype(int)
        D_v[(self.D - 1)::self.D, :] = -1 * D_v[(self.D - 1)::self.D, :]
        D_v[(self.D - 1)::self.D, :][D_v[(self.D - 1)::self.D, :] < 0] = 0
        self.D_v = D_v

    ## TODO: double check that this is working as needed
    def get_beta_gp_norms(self, beta_hat):
        assert len(beta_hat) == (self.T * (self.D - 1))
        beta_hat = beta_hat.reshape((self.T, self.D - 1))
        beta_gp_norms = np.cumsum(beta_hat ** 2, axis = 0).flatten()
        assert len(beta_gp_norms) == len(beta_hat)
        return beta_gp_norms

    def hyper_param_opt(self):
        NotImplementedError

    ## TODO: make sure wh/wv bounds are correct
    def __init__(self, D = 150, T = 365, wh_bounds = np.arange(10, 20), wv_bounds = np.arange(-5, 15)):
        self.T = T
        self.D = D
        self.K = D * T

        self.wh_bounds, self.wv_bounds = wh_bounds, wv_bounds

        self.build_norm_Mats()
        
    ## X, Y are time-series with length a multiple of self.T
    def train(self, X, Y):

        X, Y = self.format_data()

        NotImplementedError #TODO:


test_mod = hflmds_model(D = 2, T = 3)