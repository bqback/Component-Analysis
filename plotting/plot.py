import matplotlib.pyplot as plt
import numpy as np
import constants
import processing


def plot_data_set(ax, x, y, label):
    ax.plot(x, y)
    ax.set_title(label)


def plot_data(sample_spectra, source_spectra, solutions, solution_titles):
    values, _ = processing.get_a(source_spectra)
    for index, sample in enumerate(sample_spectra):
        x = np.linspace(sample.WN_min, sample.WN_max, sample.num_points)
        y_sample = np.array(sample.Absorbance)
        data_fig, axes = plt.subplots(nrows=1+2*len(solutions), ncols=1, sharex='all')
        data_fig.tight_layout()
        data_fig.supxlabel('$WN, cm^{-1}$')
        ax_source = axes[0]
        ax_source.set_ylabel('$A$', fontsize=12, rotation=0)

        ax_source.plot(x, y_sample)
        ax_source.set_title(sample.name)

        for idx, sol in enumerate(solutions):
            ax_sol = axes[2*idx+1]
            ax_sol.set_ylabel('$A$', fontsize=12, rotation=0)
            y_sol = processing.y_vals(values, sol[index])

            ax_sol.plot(x, y_sol)
            ax_sol.set_title(solution_titles[idx])

            ax_err = axes[2*idx+2]
            ax_err.set_ylabel('$A$', fontsize=12, rotation=0)
            y_err = y_sample - y_sol

            ax_err.plot(x, y_err)
            ax_err.set_title("$A_{model} - A_{calc}$")

    plt.show()
