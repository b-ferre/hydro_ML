import pandas as pd
import cupy as cp
import numpy as np

# Function to calculate KGE for all centers at once
def kge_vectorized(data, centers):
    n_samples, n_features = data.shape
    n_centers = centers.shape[0]

    # Calculate mean and std for data and centers
    mean_data = cp.mean(data, axis=1)
    std_data = cp.std(data, axis=1)
    mean_centers = cp.mean(centers, axis=1)
    std_centers = cp.std(centers, axis=1)

    # Calculate Pearson correlation coefficient
    centered_data = data - mean_data[:, cp.newaxis]
    centered_centers = centers - mean_centers[:, cp.newaxis]
    correlation = cp.dot(centered_data, centered_centers.T) / (n_features * std_data[:, cp.newaxis] * std_centers[cp.newaxis, :])
    
    # Calculate alpha and beta
    alpha = std_data[:, cp.newaxis] / std_centers[cp.newaxis, :]
    beta = mean_data[:, cp.newaxis] / mean_centers[cp.newaxis, :]

    # Calculate KGE
    kge = 1 - cp.sqrt((correlation - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)

    return kge

# Custom k-means-like algorithm using KGE
def kmeans_with_kge(data, k, max_iters=100, tol=1e-4):
    n_samples, n_features = data.shape
    centers = data[cp.random.choice(n_samples, k, replace=False)]
    prev_centers = cp.zeros_like(centers)
    labels = cp.zeros(n_samples, dtype=cp.int32)

    for iteration in range(max_iters):
        # Assign labels based on KGE
        kge_values = kge_vectorized(data, centers)
        labels = cp.argmax(kge_values, axis=1)

        # Update centers
        for j in range(k):
            members = data[labels == j]
            if len(members) > 0:
                centers[j] = cp.mean(members, axis=0)

        # Check for convergence
        if cp.allclose(centers, prev_centers, atol=tol):
            break
        prev_centers = cp.copy(centers)

    return labels, centers

# Read the data
x = pd.read_csv('data.csv')
data = x[['lat', 'lng']].values

# Convert data to GPU array
data_gpu = cp.asarray(data)

# Perform custom k-means clustering with KGE
k = 50
labels, centers = kmeans_with_kge(data_gpu, k)

# Output the result
x['cluster'] = cp.asnumpy(labels)
x.to_csv('output.csv', index=False)
