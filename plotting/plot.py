import matplotlib.pyplot as plt
import numpy as np
import constants
import processing


def plot_data_set(ax, x, y, label):
    ax.plot(x, y)
    ax.set_title(label)


def plot_data(sample_spectrum, source_spectra, solutions, labels):
    x = np.linspace(sample_spectrum.WN_min, sample_spectrum.WN_max, sample_spectrum.num_points)
    y_sample = np.array(sample_spectrum.KNa)
    data_fig = plt.figure()
    ax_source = data_fig.add_subplot(4, 1, 1)

    ax_source.plot(x, y_sample)
    ax_source.set_title("Original spectrum")

    ax_lsql = data_fig.add_subplot(4, 1, 2)
    ax_lsql_ub = data_fig.add_subplot(4, 1, 3)
    ax_nnls = data_fig.add_subplot(4, 1, 4)

    values = processing.get_a(source_spectra)
    y_lsql = processing.y_vals(values, solutions[0])
    y_lsql_ub = processing.y_vals(values, solutions[1])
    y_nnls = processing.y_vals(values, solutions[2])

    plot_data_set(ax_lsql, x, y_lsql, labels[0])
    plot_data_set(ax_lsql_ub, x, y_lsql_ub, labels[1])
    plot_data_set(ax_nnls, x, y_nnls, labels[2])

    error_fig = plt.figure()

    err_lsql = error_fig.add_subplot(3, 1, 1)
    err_lsql_ub = error_fig.add_subplot(3, 1, 2)
    err_nnls = error_fig.add_subplot(3, 1, 3)

    plot_data_set(err_lsql, x, y_sample-y_lsql, "LSQL error")
    plot_data_set(err_lsql_ub, x, y_sample-y_lsql_ub, "LSQL Unbound error")
    plot_data_set(err_nnls, x, y_sample-y_nnls, "NNLS error")

    meta_error_fig = plt.figure()

    meta_err_lsql = meta_error_fig.add_subplot(3, 1, 1)
    meta_err_lsql_ub = meta_error_fig.add_subplot(3, 1, 2)
    meta_err_nnls = meta_error_fig.add_subplot(3, 1, 3)

    plot_data_set(meta_err_lsql, x, y_lsql - y_lsql_ub, "LSQL-LSQL UB error")
    plot_data_set(meta_err_lsql_ub, x, y_lsql_ub - y_nnls, "LSQL UB-NNLS error")
    plot_data_set(meta_err_nnls, x, y_lsql - y_nnls, "LSQL-NNLS error")

    plt.show()
