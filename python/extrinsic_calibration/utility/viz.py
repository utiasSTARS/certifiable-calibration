import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import scipy.io as sio
import numpy as np
import utils

def plot_optimal():

    #Load the data
    prefixes = ["per10t0r", "per100t0r", "per500t0r", "per1000t0r",
                 "per0t100r", "per0t500r", "per0t1000r", "per0t5000r"]
    split_pre = [["per10t0r", "per100t0r", "per500t0r", "per1000t0r" ],
                 ["per0t100r", "per0t500r", "per0t1000r", "per0t5000r"]]
    dataset = "unscaled_pose_data"
    constraints = ["RCH", "RC", "RH", "R"]
    legend_labels = ["R+C+H", "R+C", "R+H", "R"]
    data = utils.results_loader(prefixes, dataset, constraints)

    #Initialize plot
    plt.rc('text', usetex=True)
    f = plt.figure()
    ax1 = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    axes_list = [ax1, ax2]
    x = np.arange(len(prefixes)/2)
    width = 0.20
    labels = ["1\%", "10\%", "50\%", "100\%"]
    x_labels = ["Translation Noise ", "Rotation Noise "]
    colours = ['#7DB8DE', '#F0A78E', '#ADC198', '#BF95C7']

    for i, list_pre in enumerate(split_pre):
        for j, constraint in enumerate(constraints):
            optimal = np.zeros((len(list_pre)))
            for k, prefix in enumerate(list_pre):
                results = data[prefix][constraint]["dual_time"].flatten()
                num_tests = results.shape
                optimal[k] = np.sum(np.not_equal(results, 1000.))/num_tests * 100
            axes_list[i].bar(x + (j - 3/2)*width, optimal, width, color=colours[j], edgecolor='#252525')
        axes_list[i].set_ylabel("Optimal \%")
        axes_list[i].grid(True, which='both', color='tab:grey', linestyle='--', alpha=0.5, linewidth=0.5)
    axes_list[0].set_xticks(x)
    axes_list[0].set_xticklabels(labels)
    axes_list[0].set_xlabel(x_labels[0])
    axes_list[1].set_xticks(x)
    axes_list[1].set_xticklabels(["10\%", "50\%", "100\%", "500\%"])
    axes_list[1].set_xlabel(x_labels[1])
    axes_list[0].legend(legend_labels, loc=1, prop={'size': 8})
    f.tight_layout()
    f.savefig("Optimal_plot.pdf")


    return

def plot_histogram_error(labels, data, axis_handles):
    if len(labels) == 0:
        labels = [''] * data.shape[0]
    for i in range(data.shape[0]):
        if i==0: #Translation Error
            axis_handles[0].hist(data[i, :], bins=np.arange(0, 0.5, 0.05), alpha = 0.5, label=labels[i])
            axis_handles[0].set_xlim((0, 0.5))
            #axis_handles[0].set_xticks((1.25/2 * np.arange(0, 3)).tolist())
            axis_handles[0].grid(True, which='both', color='tab:grey', linestyle='--', alpha=0.5, linewidth=0.5)
        elif i==1: #Rotation Error
            axis_handles[1].hist(data[i, :], bins=np.arange(0, 0.05, 0.005), alpha = 0.5, label=labels[i])
            axis_handles[1].set_xlim((0, 0.05))
            axis_handles[1].grid(True, which='both', color='tab:grey', linestyle='--', alpha=0.5, linewidth=0.5)
        else: #Alpha Error
            axis_handles[2].hist(data[i, :], bins=np.arange(0, 0.005, 0.0005), alpha = 0.5, label=labels[i])
            axis_handles[2].set_xlim((0, 0.005))
            axis_handles[2].grid(True, which='both', color='tab:grey', linestyle='--', alpha=0.5, linewidth=0.5)
    return

def create_error_fig():
    # Load the data
    # prefixes = ["per10t10r", "per25t10r", "per50t10r",
    #             "per10t50r", "per25t50r", "per50t50r",
    #             "per10t100r", "per25t100r", "per50t100r"]
    # split_pre = [["per10t10r", "per25t10r", "per50t10r"],
    #             ["per10t50r", "per25t50r", "per50t50r"],
    #             ["per10t100r", "per25t100r", "per50t100r"]]
    prefixes = ["matt_N_100_per5t10r", "matt_N_100_per5t20r",
                "matt_N_100_per10t10r", "matt_N_100_per10t20r"]
    split_pre = [["matt_N_100_per5t10r", "matt_N_100_per5t20r"],
                ["matt_N_100_per10t10r", "matt_N_100_per10t20r"]]
    sigma_t = [0.5, 0.5, 1, 1, 1.5, 2]
    sigma_r = [1, 1.5, 1, 1.5, 1, 2]
    dataset = "unscaled_pose_data"
    constraints = ["RCH"]
    data_of_interest = [["dual_trans_error", "dual_rot_error", "dual_alpha_error"], ["andreff_trans_error", "andreff_rot_error", "andreff_alpha_error"]]
    data_dict = utils.results_loader(prefixes, dataset, constraints)

    #Initialize parts of the figure
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
    f = plt.figure(figsize=(6.0, 9.0))
    axes_list = [f.add_subplot(4,3, 1 + m) for m in range(12)]

    #Loop through the data
    max = 0
    for i, data_list in enumerate(split_pre):
        for j, data_name in enumerate(data_list):
            plotting_axes = axes_list[3*j + 6 * i: 3*(j + 1) + 6 * i]
            data = 0
            for solver in data_of_interest:
                data = np.vstack([data_dict[data_name]['RCH'][error] for error in solver])
                labels= [ "Dual"  if 'dual' in error else "Linear" for error in solver ]
                data_mask = np.not_equal(data[0, :], 1000.)
                if np.sum(data_mask) > 0:
                    data = data[:, data_mask]
                    plot_histogram_error(labels, data, plotting_axes)
                    if np.max(data[2, :])>max:
                        max = np.max(data[2, :])
                else:
                    print("No Data")

    print(max)
    #Clean up plot a bit
    max_counts = 100
    [axis.set_ylim((0, max_counts)) for axis in axes_list]
    # Experiment Noise Values
    for i in range(4):
        # Experiment Labels
        axes_list[3*i+2].text(0.003, 0.65*max_counts, r' \['
                                         r' \begin{split} '
                                         r' \sigma_t &= ' + str(sigma_t[i]) + r' \% \\'
                                         r' \sigma_R &= ' + str(sigma_r[i]) + r' \% \\'
                                         r' \end{split}'
                                         r' \]', {'fontsize': 8},
                                 bbox=dict(boxstyle='round',fc="#FFFFFF"))
        # Legends
        axes_list[3*i].legend(loc=1, prop={'size': 8})
        # Y axis labels
        axes_list[3*i].set_ylabel("Counts")
    axes_list[-3].set_xlabel("Translation Error [m]")
    axes_list[-2].set_xlabel("Rotation Error [rad]")
    axes_list[-1].set_xlabel("Scale Error [unitless]")
    f.subplots_adjust(wspace=0.4, hspace=0.3)
    #f.suptitle(title)
    f.savefig("Error_Plot.pdf",  bbox_inches='tight')
    # plt.show()
    return

def create_histogram_error_fig(data_name, data, subslice, title):
    f = plt.figure()
    ax1 = f.add_subplot(331)
    ax2 = f.add_subplot(332)
    ax3 = f.add_subplot(333)
    ax4 = f.add_subplot(334)
    ax5 = f.add_subplot(335)
    ax6 = f.add_subplot(336)
    ax7 = f.add_subplot(337)
    ax8 = f.add_subplot(338)
    ax9 = f.add_subplot(339)
    plot_histogram_error(data_name[0], data[0][subslice, :], [ax1, ax2, ax3])
    plot_histogram_error(data_name[1], data[1][subslice, :], [ax4, ax5, ax6])
    plot_histogram_error(data_name[2], data[2][subslice, :], [ax7, ax8, ax9])
    ax7.set_xlabel("Translation Error [m]")
    ax8.set_xlabel("Rotation Error [rad]")
    ax9.set_xlabel("Scale Error [unitless]")
    ax1.set_ylabel("Counts")
    ax4.set_ylabel("Counts")
    ax7.set_ylabel("Counts")
    f.subplots_adjust(wspace=0.3, hspace=0.3)
    f.suptitle(title)
    f.savefig(title+".png",  bbox_inches='tight')
    #plt.show()
    return

def main():
    #plot_optimal()
    create_error_fig()
    data = sio.loadmat("plotting_data.mat")
    data_list = [("no_noise", data["no_noise"]), ("trans_noise", data["trans_noise"]), ("rot_noise", data["rot_noise"])]
    # create_histogram_error_fig([data_list[0][0], data_list[1][0], data_list[2][0]], [data_list[0][1], data_list[1][1], data_list[2][1]], [0, 2, 4], "Lagrangian Dual Solution Error")
    # create_histogram_error_fig([data_list[0][0], data_list[1][0], data_list[2][0]],
    #                            [data_list[0][1], data_list[1][1], data_list[2][1]], [1, 3, 5], "SDP Relaxation Solution Error")
    return

if __name__ == "__main__":
    main()
