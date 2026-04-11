Official implementation of the paper "Dynamic spectral co-clustering of directed networks to unveil latent community paths in VAR-type models" by Younghoon Kim and Changryong Baek.
This repository provides a statistical framework to detect and track group structures within complex, directed networks over time.

Key Features:
Model Frameworks: Implements Periodic VAR (PVAR) models to capture seasonal network changes, and generalized Vector Heterogeneous Autoregressive (VHAR) models to analyze short-, medium-, and long-term network dependencies.

Network Estimation: Utilizes Lasso estimation to efficiently handle high-dimensional data and reliably estimate the initial network connections.

Community Tracking: Recovers hidden groups using spectral co-clustering and tracks their evolution, allowing users to observe exactly when communities split, merge, or remain stable.
Applications include analyzing time-series data, such as U.S. macroeconomic trends and global stock index movements.
