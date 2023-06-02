import scipy as sp
import numpy as np
import classes
import plotting
from .ica import run_ica, correlation
from .jade import jadeR
from constants import TRUE_SOLUTION, CONSTRAINT_TYPES, SNR, REQUIRED_SAMPLES
from typing import List, Dict
import matplotlib.pyplot as plt
from sklearn import linear_model
from functools import reduce
from pymcr.mcr import McrAR
from pymcr.constraints import ConstraintNonneg


def value_at_point(gasses, coefficients, data_point):
    value = 0
    for i in range(len(gasses)):
        print(f'\t{gasses[i].name}: {coefficients[i] * gasses[i].Absorbance[data_point]}')
        value += coefficients[i] * gasses[i].Absorbance[data_point]
    return value


def normalize(coefficients, nh3_pos):
    norm_coefficients = [coeff / (coefficients[nh3_pos]*1e8) for coeff in coefficients]
    return norm_coefficients


def get_a(spectra):
    array = []
    names = []

    for spectrum in spectra:
        names.append(spectrum.name)
        array.append(spectrum.Absorbance)

    a_horizontal = np.array(array)
    a = np.transpose(a_horizontal)

    return a, names


def y_vals(values, coefficients):
    return values.dot(coefficients)


def reorder_points(sample_name, names):
    reordered = []
    solution = TRUE_SOLUTION[sample_name]
    for name in names:
        reordered.append(solution[name])
    return reordered


def add_noise(spectrum, snr_db=20):
    mean = np.mean(spectrum)
    mean_db = 10*np.log10(mean)
    noise_db = mean_db - snr_db
    noise_std = 10 ** (noise_db/10)
    noise_mean = 0
    noise_spectrum = np.random.normal(noise_mean, noise_std, spectrum.shape[0])
    return spectrum+noise_spectrum


def lsql_batch_solve(a, b, constrain=True):
    if constrain:
        b = np.append(b, 1)
    result = sp.optimize.lsq_linear(
        a, b,
        bounds=[0, 1], lsq_solver='exact'
    )
    return result.x, result.cost, result.active_mask


def derive(spectrum, order):
    j = 1
    sg_window = 35
    sg_deg = 2
    while j <= order:
        spectrum = sp.signal.savgol_filter(spectrum, sg_window, sg_deg, deriv=1, mode='nearest')
        j += 1
    return spectrum


def get_spectra_array(spectra_classes: List[classes.Spectrum]):
    return np.array([spectrum.Absorbance for spectrum in spectra_classes])


def generate_coefficients() -> Dict[str, float]:
    ch4 = max(1e-12, np.random.normal(1.48e-6, 5e-8))
    co = max(1e-12, np.random.normal(4.7e-7, 4e-10))
    co2 = max(1e-12, np.random.normal(3.3e-4, 6e-7))
    h2o = max(1e-12, np.random.normal(1e-2, 5e-3))
    n2o = max(1e-12, np.random.normal(2.8e-7, 1e-9))
    nh3 = max(1e-12, np.random.normal(1e-8, 1e-12))
    o2 = max(1e-12, np.random.normal(2.07e-1, 5e-2))
    o3 = max(1e-12, np.random.normal(2e-8, 2e-11))
    so2 = max(1e-12, np.random.normal(8e-8, 8e-11))
    n2 = 1 - (ch4 + co + co2 + h2o + n2o + nh3 + o2 + o3 + so2)
    coeff_dict = {
        "Pure CH4": ch4,
        "Pure CO": co,
        "Pure CO2": co2,
        "Pure H2O": h2o,
        "Pure N2": n2,
        "Pure N2O": n2o,
        "Pure NH3": nh3,
        "Pure O2":  o2,
        "Pure O3": o3,
        "Pure SO2": so2
    }
    return coeff_dict


def generate_sample(sources, WN, i):
    coefficients = generate_coefficients()
    abs_spectrum = sources.dot(list(coefficients.values()))
    name = f"Randomly generated sample #{i}"
    return coefficients, \
        classes.Spectrum(name, np.array([WN, np.zeros(abs_spectrum.shape), abs_spectrum]))


def calculate(data: Dict[str, List[classes.Spectrum]]):
    a, names = get_a(data["sources"])
    a_constrained = np.append(a, [np.ones(a.shape[1])], axis=0)

    wn_space = np.linspace(data["samples"][0].WN_min, data["samples"][0].WN_max, data["samples"][0].num_points)

    known_solutions = TRUE_SOLUTION

    if len(data["samples"]) < REQUIRED_SAMPLES:
        for i in range(REQUIRED_SAMPLES - len(data["samples"])):
            coeffs, new_sample = generate_sample(a, data["sources"][0].WN, i+1)
            known_solutions.update({new_sample.name: coeffs})
            data["samples"].append(new_sample)

    sample_array = get_spectra_array(data["samples"])
    lowest_max = 1e9
    lowest_max_sample = -1
    for i in range(sample_array.shape[0]):
        print(f"Sample {i} max {sample_array[i].max()} at L={wn_space[np.argmax(sample_array[i])]}")
        if sample_array[i].max() < lowest_max:
            lowest_max = sample_array[i].max()
            lowest_max_sample = i
    lowest_max_wn = np.argmax(sample_array[lowest_max_sample])
    print(lowest_max_wn)
    rescaled_samples = np.array([sample_array[i] / sample_array[i, lowest_max_wn] for i in range(sample_array.shape[0])])
    for i in range(rescaled_samples.shape[0]):
        print(f"Sample {i} max {rescaled_samples[i].max()} at L={wn_space[np.argmax(rescaled_samples[i])]}")
    noisy_samples = np.array([add_noise(sample, snr_db=SNR) for sample in sample_array])

    batch_lsql = [lsql_batch_solve(a_constrained, sample) for sample in sample_array]
    solutions_lsql = []
    for index, (solution, residual, bound_type) in enumerate(batch_lsql):
        print(f"LSQ Linear solution for sample {data['samples'][index].name}")
        solutions_lsql.append(solution)
        for component_index, component in enumerate(data["sources"]):
            calculated_value = solution[component_index]
            component_name = names[component_index]
            known_value = known_solutions[data['samples'][index].name][component_name]
            print(f"\t{component.name:<8}: {calculated_value:<25}"
                  f"({CONSTRAINT_TYPES[bound_type[component_index]]:^24}) ",
                  f"known: {known_value:<23} "
                  f"(relative error {abs((known_value - calculated_value)/known_value)*100}%)")
        print(f"\tResidual: {residual}")

    batch_lsql_noisy = [lsql_batch_solve(a_constrained, sample) for sample in noisy_samples]
    solutions_lsql_noisy = []
    for index, (solution, residual, bound_type) in enumerate(batch_lsql_noisy):
        print(f"LSQ Linear solution for noisy sample {data['samples'][index].name} (SNR = {SNR})")
        solutions_lsql_noisy.append(solution)
        for component_index, component in enumerate(data["sources"]):
            calculated_value = solution[component_index]
            component_name = names[component_index]
            known_value = known_solutions[data['samples'][index].name][component_name]
            print(f"\t{component.name:<8}: {calculated_value:<25}"
                  f"({CONSTRAINT_TYPES[bound_type[component_index]]:^24}) ",
                  f"known: {known_value:<23} "
                  f"(relative error {abs((known_value - calculated_value) / known_value) * 100}%)")
        print(f"\tResidual: {residual}")

    U, Sv, V = np.linalg.svd(sample_array, full_matrices=True)
    S = np.zeros(sample_array.shape)
    np.fill_diagonal(S, Sv, wrap=True)
    N = np.linalg.norm(sample_array, ord=1)
    E = np.zeros(sample_array.shape[0])
    for nn in range(0, sample_array.shape[0]):
        Rec = U @ S[:, 0:nn + 1] @ V[0:nn + 1, :]
        E[nn] = np.linalg.norm(sample_array - Rec, ord=1) / N
    DE = np.append(-E[0], np.diff(E))
    nica = np.max([(sum(E >= 2e-8) + 1), sum(DE < -2e-8)])
    print(f"N of ICA components {nica}")

    # run_ica(data["sources"], data["samples"])
    jade_W = jadeR(sample_array, m=10)
    sources_JADE = jade_W.dot(sample_array)
    print(sources_JADE.shape)
    sources_JADE = abs(sources_JADE)

    print("Cors JADE")
    Cor = correlation(sources_JADE, a.T, 10)
    print(Cor)
    print("I")
    I = Cor.argmax(1)
    print(I)
    print("Identified species")
    confidence_JADE = []
    names_JADE = []
    print(a.T.shape)
    for spec, i in enumerate(I):
        print(f"{data['sources'][i].name}: {Cor[spec][i]*100:.3f}% confident")
        names_JADE.append(data['sources'][i].name)
        confidence_JADE.append(Cor[spec][i]*100)
        if Cor[spec][i]*100 > 95:
            sources_JADE[spec] = abs(sources_JADE[spec])*a.T[i].max()/abs(sources_JADE[spec]).max()
    Ivalues = Cor.max(1)
    I.sort()
    I = np.unique(I)
    print(I)

    L = a.T[I, :]
    print("Spectral matrix JADE")
    print(L)
    reg_JADE = linear_model.ElasticNet(alpha=1e-5, max_iter=int(2e6), positive=True)
    reg_JADE.fit(L.T, sample_array.T)
    G = reg_JADE.coef_
    print(G)
    print(G.sum(1))
    Gn = G / G.sum(1)[:, None]*100
    print("[JADE] The associated coefficients of the mixing matrix (proportional or equal to concentrations)")
    print(G / G.sum(1)[:, None] * 100)

    print(f"sources_JADE shape -- {sources_JADE.shape}")
    print(sources_JADE)
    print(np.asarray(sources_JADE)**2)
    mcrals = McrAR(c_regr=linear_model.ElasticNet(alpha=1e-5, l1_ratio=0.75), max_iter=700,
               tol_err_change=sample_array.max() * 1e-18, st_regr='NNLS', c_constraints=[ConstraintNonneg()])
    # mcrals.fit(sample_array, ST=np.asarray(sources_JADE)**2)
    mcrals.fit(sample_array, ST=a.T)
    sources_MCRALS = np.asmatrix(mcrals.ST_opt_)

    print("Cors MCRALS")
    Cor = correlation(sources_MCRALS, a.T, nica)
    print(Cor)
    print("I")
    I = Cor.argmax(1)
    print(I)
    print("Identified species")
    for spec, i in enumerate(I):
        print(f"{data['sources'][i].name}: {Cor[spec][i] * 100:.3f}% confident")
    Ivalues = Cor.max(1)
    I.sort()
    I = np.unique(I)
    print(I)

    L = a.T[I, :]
    print("Spectral matrix MCRALS")
    print(L)
    reg_MCRALS = linear_model.Lasso(alpha=1e-1, max_iter=int(2e4), positive=True)
    reg_MCRALS.fit(L.T, rescaled_samples.T)
    G = reg_MCRALS.coef_
    print(G)
    print(G.sum(1))
    Gn = G / G.sum(1)[:, None] * 100
    print("[MCRALS] The associated coefficients of the mixing matrix (proportional or equal to concentrations)")
    print(G / G.sum(1)[:, None] * 100)

    for idx, source in enumerate(sources_JADE):
        plt.figure(f"{names_JADE[idx]}: {confidence_JADE[idx]:.3f}% confident")
        plt.plot(wn_space, abs(source.T))

    for idx, source in enumerate(sources_MCRALS):
        plt.figure(f'MCRALS Source reconstruction {idx}')
        plt.plot(wn_space, abs(source.T))

    plt.show()
    solutions = [solutions_lsql, solutions_lsql_noisy]
    solution_titles = ["LSQL method", f"LSQL method with noisy samples (SNR = {SNR})"]

    # plotting.plot_data(data["samples"],
    #                    data["sources"],
    #                    solutions,
    #                    solution_titles)
