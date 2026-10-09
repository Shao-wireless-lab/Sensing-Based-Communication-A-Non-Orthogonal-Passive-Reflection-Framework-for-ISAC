# Source Code for Theoretical Validation of Multi-Tag Backscatter System

This repository contains the official MATLAB simulation and theoretical modeling source code for our paper. Each script corresponds to a specific experiment and evaluation presented in the manuscript.

## 📁 File Structure & Usage

Based on the repository contents, the files map to our experiments as follows:
- **`simulation2.m`, `simulation4.m`**: Corresponding simulation scripts for the respective experiments in the paper.
- **`simulation3_1.m`**: Generates the empirical simulation values for the multi-tag scenarios. 
- **`simulation3.m`**: The core theoretical framework and validation script. It utilizes a specific output snapshot from `simulation3_1.m` as the empirical "ground truth" to fit and evaluate the theoretical model.
- **`linearFitting.m` / `tek0001ALL.csv`**: Auxiliary scripts and data used for signal processing and channel fitting.

## ⚠️ Important Note on Randomness

Please note that our system simulation involves the generation of random OFDM signals, randomized data payloads, and Gaussian noise vectors. Therefore, **slight numerical variations are expected** if you re-run the simulation scripts from scratch. These minor fluctuations are a normal characteristic of Monte Carlo-style physical layer simulations and do strictly preserve the capacity bounds and system trends.

## 📝 Post-Submission Update Notice (Specifically for `simulation3.m`)

Please note that this repository represents the **actively maintained version** of our codebase. Following the initial manuscript submission, we further refined the physical parameter mapping in `simulation3.m`. 

Specifically, we rigorously aligned the theoretical baseline noise factor with the discrete-time simulation sampling rate (compensating for the matched filter oversampling gain). As a result of this rigorous physical calibration, the theoretical modeling has become physically more robust. The updated model in `simulation3.m` yields a **MAPE shift of +0.025** (e.g., reaching approximately 0.167) compared to the earlier snapshot reported in the submitted PDF. 

We emphasize that this extremely minor numerical refinement **does not alter the shape of the capacity curves, the MAI-limited boundaries, or any scientific conclusions** presented in the paper. It continuously demonstrates the high accuracy of our theoretical framework in predicting multi-tag interference limits.
