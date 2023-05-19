import numpy as np
import processing
import time
from sklearn.decomposition import FastICA
from sklearn import linear_model
from pymcr.mcr import McrAR
from pymcr.regressors import OLS, NNLS
from constants import SNR
from pymcr.constraints import ConstraintNonneg, ConstraintNorm
import matplotlib.pyplot as plt

def center(x):
    mean = np.mean(x, axis=1, keepdims=True)
    centered = x - mean
    return centered, mean


def covariance(x):
    mean = np.mean(x, axis=1, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean

    return (m.dot(m.T)) / n


def correlation(unknown, sources, nica):
    Cor = np.zeros((nica, sources.shape[0]))

    print(unknown.shape)
    print(sources.shape)
    # Matching
    for ii in range(nica):
        for jj in range(sources.shape[0]):
            y = np.array([unknown[ii, :] / unknown[ii, :].max(), sources[jj, :] / sources[jj, :].max()])
            Cor[ii, jj] = abs(np.corrcoef(y)[0, 1])
            if np.isnan(Cor[ii, jj]):
                Cor[ii, jj] = 0.
    return Cor


def whiten(X):
    # Calculate the covariance matrix
    print('Calculate the covariance matrix')
    coVarM = covariance(X)

    # Single value decomposition
    print('Single value decomposition')
    U, S, V = np.linalg.svd(coVarM)

    # Calculate diagonal matrix of eigenvalues
    print('Calculate diagonal matrix of eigenvalues')
    d = np.diag(1.0 / np.sqrt(S))

    # Calculate whitening matrix
    print('Calculate whitening matrix')
    whiteM = np.dot(U, np.dot(d, U.T))

    # Project onto whitening matrix
    print('Project onto whitening matrix')
    Xw = np.dot(whiteM, X)

    return Xw, whiteM


def fastIca(signals, alpha=1, thresh=1e-8, iterations=5000):
    m, n = signals.shape

    # Initialize random weights
    W = np.random.rand(m, m)

    for c in range(m):
        w = W[c, :].copy().reshape(m, 1)
        w = w / np.sqrt((w ** 2).sum())

        i = 0
        lim = 100
        while ((lim > thresh) & (i < iterations)):
            # Dot product of weight and signal
            ws = np.dot(w.T, signals)

            # Pass w*s into contrast function g
            wg = np.tanh(ws * alpha).T

            # Pass w*s into g prime
            wg_ = (1 - np.square(np.tanh(ws))) * alpha

            # Update weights
            wNew = (signals * wg.T).mean(axis=1) - wg_.mean() * w.squeeze()

            # Decorrelate weights
            wNew = wNew - np.dot(np.dot(wNew, W[:c].T), W[:c])
            wNew = wNew / np.sqrt((wNew ** 2).sum())

            # Calculate limit condition
            lim = np.abs(np.abs((wNew * w).sum()) - 1)

            # Update weights
            w = wNew

            # Update counter
            i += 1

        W[c, :] = w.T
    return W


def run_ica(source_spectra, sample_spectra):
    wn_space = np.linspace(source_spectra[0].WN_min,
                           source_spectra[0].WN_max,
                           source_spectra[0].num_points
                           )
    source_data, _ = processing.get_a(source_spectra)
    sample_data_preserve = np.array([sample.Absorbance for sample in sample_spectra])
    sample_data = sample_data_preserve * 1
    noisy_samples = np.array([processing.add_noise(sample, snr_db=SNR) for sample in sample_data_preserve])
    noisy_ICA = noisy_samples * 1
    noisy_data = [
        np.array([processing.add_noise(sample, snr_db=snr) for sample in sample_data_preserve])
        for snr in np.arange(20, 120, 10)
    ]
    sa_norm = np.array([sample / sample.max() for sample in sample_data])
    noisy_norm = np.array([sample / sample.max() for sample in noisy_samples])
    sf = np.trapz(sample_data[0, :], wn_space)
    means = noisy_samples.mean(1)
    stdsqr = np.sqrt(noisy_samples.std(1))
    for i in range(noisy_samples.shape[0]):
        print(f"Before scale max {noisy_samples[i].max()}")
        noisy_samples[i, :] /= noisy_samples[i, np.where(wn_space >= 2328)[0][0]]
        print(f"After scale max {noisy_samples[i].max()}")
        noisy_samples[i, :] = (noisy_samples[i, :] - means[i]) / stdsqr[i]
    for i in range(noisy_samples.shape[0]):
        print(noisy_samples[i].max())
    noisy_derivative = np.array([processing.derive(sample, order=2) for sample in noisy_samples])

    eps_fig, eps_ax = plt.subplots(1, 1)
    snr_range = np.arange(20, 120, 10)
    eps = []
    sk = []
    for idx, samples in enumerate(noisy_data):
        u, sv, v = np.linalg.svd(samples, full_matrices=True)
        s = np.zeros(samples.shape)
        np.fill_diagonal(s, sv, wrap=True)
        n = np.linalg.norm(samples, ord=1)
        e = np.zeros(samples.shape[0])
        sk_sum = []
        for nn in range(samples.shape[0]):
            Rec = u @ s[:, :nn + 1] @ v[:nn + 1, :]
            sk_sum.append(np.sum(s[:, :nn+1])/np.sum(s))
            e[nn] = np.linalg.norm(samples - Rec, ord=1) / n
        print(sk_sum)
        eps.append(e)
        sk.append(sk_sum)
    eps_k = np.array(eps).T
    for k_val in range(eps_k.shape[0]):
        color = (1, 0, 0) if k_val == 10 else (0, 0, k_val/eps_k.shape[0])
        eps_ax.plot(snr_range, np.log10(eps_k[k_val]), color=color)
    sk_fig, sk_ax = plt.subplots(1, 1)
    k_range = np.arange(1, noisy_data[0].shape[0]+1, 1)
    for idx, sk_val in enumerate(sk):
        sk_ax.plot(k_range, sk_val, color=(0, 0, idx/len(sk)))

    u, sv, v = np.linalg.svd(noisy_samples, full_matrices=True)
    s = np.zeros(noisy_samples.shape)
    np.fill_diagonal(s, sv, wrap=True)
    n = np.linalg.norm(noisy_samples, ord=1)
    e = np.zeros(noisy_samples.shape[0])
    for nn in range(noisy_samples.shape[0]):
        Rec = u @ s[:, :nn+1] @ v[:nn+1, :]
        e[nn] = np.linalg.norm(noisy_samples - Rec, ord=1) / n
    de = np.append(-e[0], np.diff(e))
    source_num = np.max([sum(e >= 1e-7) + 1, sum(de < -1e-7)])
    print(source_num)
    tic = time.time()
    ica = FastICA(fun='logcosh', n_components=10, tol=1e-10, max_iter=5000,
                  whiten='arbitrary-variance')
    ica.fit_transform(noisy_derivative.T)
    mixed_est = ica.mixing_
    toc = time.time()
    ica_runtime = toc - tic

    shift = np.linalg.pinv(mixed_est) @ (means.reshape(-1, 1) / stdsqr.reshape(-1, 1))
    source_ICA = (np.linalg.pinv(mixed_est) @ noisy_samples) + shift
    source_der_ICA = np.linalg.pinv(mixed_est) @ noisy_derivative
    mixed_est *= stdsqr.reshape(-1, 1)
    mixed_ICA = mixed_est
    # for i in range(mixed_est.shape[0]):
    #     fig = plt.figure("ICA test")
    #     plt.subplot(1, 1, 1)
    #     y_values = processing.y_vals(source_data, mixed_est[i])
    #     plt.plot(wn_space, y_values)
    # MCR
    tic = time.time()
    print(noisy_norm.max())
    mcrals = McrAR(c_regr=linear_model.ElasticNet(alpha=1e-5, l1_ratio=0.75), max_iter=700,
                   tol_err_change=noisy_norm.max() * 1e-8, st_regr='NNLS',
                   c_constraints=[ConstraintNonneg(), ConstraintNorm()])
    print(noisy_norm.shape)
    print(source_ICA.shape)
    mcrals.fit(noisy_norm, ST=source_ICA**2)
    toc = time.time()
    mcr_runtime = toc - tic
    source_MCR = mcrals.ST_opt_
    mixed_MCR = mcrals.C_opt_
    source_der_MCR = (np.linalg.pinv(mixed_MCR) @ noisy_derivative)

    print(f"source MCR shape {source_MCR.shape}")
    print(f"source data shape {source_data.shape}")

    Cor = correlation(source_MCR, source_data.T, source_num)
    I = Cor.argmax(1)
    Ivalues = Cor.max(1)
    I.sort()
    I = np.unique(I)
    L = source_data[I, :]
    reg = linear_model.Lasso(alpha=1e-1, max_iter=int(2e4), positive=True)
    reg.fit(L.T, noisy_norm.T)
    G = reg.coef_
    Gn = G / G.sum(1)[:, None] * 100
    print(Gn)
    # sample_centered, mean = center(sample_data)
    # sample_whitened, whiteCoeff = whiten(sample_centered)
    # print(np.round(covariance(sample_whitened)))
    # W = fastIca(sample_whitened, alpha=1)
    # unmixed = sample_whitened.T.dot(W.T)
    # unmixed = (unmixed.T - mean)
    # figure, axes = plt.subplots(nrows=2, ncols=1)
    # figure.tight_layout()
    # print(unmixed.shape)
    # for i in range(source_data.shape[0]):
    #     axes[0].plot(wn_space, source_data[i])
    # for i in range(unmixed.shape[0]):
    #     axes[1].plot(wn_space, unmixed[i])

