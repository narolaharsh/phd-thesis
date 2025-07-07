import numpy as np
import matplotlib.pyplot as plt
import bilby
import corner

import scienceplots
plt.style.use(['science', 'bright'])

fontsize = 11

ticksize = 9
#plt.rc('xtick',labelsize=ticksize)
#plt.rc('ytick',labelsize=ticksize)

plt.rcParams.update({
    'font.size': fontsize,              # Base font size
    'axes.titlesize': fontsize,         # Axes title
    'axes.labelsize': fontsize,         # Axes labels
    'xtick.labelsize': ticksize,        # X-axis tick labels
    'ytick.labelsize': ticksize,        # Y-axis tick labels
    'legend.fontsize': fontsize,        # Legend text
    'figure.titlesize': fontsize        # Figure title
})
plt.rcParams['axes.labelsize'] = fontsize
data = np.genfromtxt('overlaps_0_posterior_samples.dat.txt', names = True)

signal_A = np.genfromtxt('single_injected_A_0_posterior_samples.dat.txt', names = True)
signal_B = np.genfromtxt('single_injected_B_0_posterior_samples.dat.txt', names = True)

data['luminosity_distance_A'] /= 1000
data['luminosity_distance_B'] /= 1000




#params = ['chirp_mass_A', 'chirp_mass_B', 'theta_jn_A', 'theta_jn_B', 'luminosity_distance_A', 'luminosity_distance_B']

params = ['chirp_mass_A', 'chirp_mass_B', 'mass_ratio_A', 'mass_ratio_B', 'luminosity_distance_A', 'luminosity_distance_B']



#latex_labels = [r'$\mathcal{M}^{\mathrm{A}}_c[M_{\odot}]$', r'$\mathcal{M}^{\mathrm{B}}_c[M_{\odot}]$', r'${\theta}^{\mathrm{A}}_{\mathrm{JN}}$', r'${\theta}^{\mathrm{B}}_{\mathrm{JN}}$', r'$d^{\mathrm{A}}_{\mathrm{L}}[\rm{Gpc}]$', r'$d^{\mathrm{B}}_{\mathrm{L}}[\rm{Gpc}]$']
latex_labels_A = [r'$\mathcal{M}^{\mathrm{A}}_c[M_{\odot}]$', r'$q^{\mathrm{A}}$',  r'$\theta^{\mathrm{A}}_{\mathrm{JN}}$', r'$d^{\mathrm{A}}_{\mathrm{L}}[\rm{Gpc}]$']
latex_labels_B = [r'$\mathcal{M}^{\mathrm{B}}_c[M_{\odot}]$', r'$q^\mathrm{B}$',  r'$\theta^{\mathrm{B}}_{\mathrm{JN}}$', r'$d^{\mathrm{B}}_{\mathrm{L}}[\rm{Gpc}]$']



injection_parameters = {
    "mass_ratio_A": 0.2231017059328386,
    "mass_ratio_B": 0.37172837551203913,
    "chirp_mass_A": 22.353653661893304,
    "chirp_mass_B": 47.64517481623601,
    "luminosity_distance_A": 7.268571745364448,
    "luminosity_distance_B": 50.09434936920892,
    "dec_A": 0.6557802485585996,
    "dec_B": -0.25296571570747095,
    "ra_A": 2.180684709384296,
    "ra_B": 4.80942360638223,
    "theta_jn_A": 1.8552485990854504,
    "theta_jn_B": 0.7722659493489856,
    "psi_A": 2.317315123447218,
    "psi_B": 0.3869169812264324,
    "phase_A": 0.736517637678661,
    "phase_B": 1.0168760933964565,
    "a_1_A": 0.5874657489189719,
    "a_1_B": 0.8251185871884739,
    "a_2_A": 0.7740290285427935,
    "a_2_B": 0.03742190139192977,
    "tilt_1_A": 1.8400317067401333,
    "tilt_1_B": 1.2114859809247314,
    "tilt_2_A": 2.4431533314845897,
    "tilt_2_B": 1.731705006294819,
    "phi_12_A": 1.4830050470970393,
    "phi_12_B": 1.2192989717535039,
    "phi_jl_A": 3.5537095497820763,
    "phi_jl_B": 4.538201491570135,
    "geocent_time_A": 0.0,
    "geocent_time_B": 0.013082053778854186
  }

truth_dictionary = [injection_parameters[p] for p in params]

joint_A = np.vstack([np.array(data['chirp_mass_A']), np.array(data['mass_ratio_A']), np.array(data['theta_jn_A']), np.array(data['luminosity_distance_A'])]).T
joint_B = np.vstack([np.array(data['chirp_mass_B']), np.array(data['mass_ratio_B']), np.array(data['theta_jn_B']), np.array(data['luminosity_distance_B'])]).T


signal_A_pe = np.vstack([np.array(signal_A['chirp_mass']), signal_A['mass_ratio'], signal_A['theta_jn'], signal_A['luminosity_distance']/1000]).T
signal_B_pe = np.vstack([np.array(signal_B['chirp_mass']), signal_B['mass_ratio'], signal_B['theta_jn'], signal_B['luminosity_distance']/1000]).T

fig = plt.figure(figsize=(5.5, 5.5))
total_bins = 40
linewidth = 1.5

params_A_tags = ['chirp_mass_A', 'mass_ratio_A', 'theta_jn_A', 'luminosity_distance_A']
truth_params = [injection_parameters[p] for p in params_A_tags]
fig1 = corner.corner(signal_A_pe, fig = fig,
                    bins = total_bins, truths = truth_params, truth_color = 'black',
                    color = 'rosybrown', 
                    density = 1,
                    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                   plot_density=True, 
                    plot_datapoints=False, 
                    max_n_ticks=3,
                    fill_contours=True, 
                    **{
                        "hist_kwargs": {"linewidth": linewidth, 'density': 1},         # 1D histograms
                        "contour_kwargs": {"linewidths": linewidth},     # KDE contours
                        "truth_kwargs": {"linewidth": linewidth},        # Truth lines
                    })


fig1 = corner.corner(joint_A, fig = fig,
                    labels = latex_labels_A, 
                    bins = total_bins, 
                    color = 'royalblue', 
                    density = 1,
                    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                   plot_density=True, 
                    plot_datapoints=False, 
                    max_n_ticks=3,
                    fill_contours=True, 
                    **{
                        "hist_kwargs": {"linewidth": linewidth, 'density': 1},         # 1D histograms
                        "contour_kwargs": {"linewidths": linewidth},     # KDE contours
                        "truth_kwargs": {"linewidth": linewidth},        # Truth lines
                    })



ndim = len(params_A_tags)
import matplotlib.lines as mpllines
lines = []

lines.append(mpllines.Line2D([0], [0], color='rosybrown'))
lines.append(mpllines.Line2D([0], [0], color='royalblue'))
labels = ['Signal A isolation', 'Signal A overlapping']
final_axes = fig1.get_axes()
final_axes[ndim - 1].legend(lines, labels, fontsize = fontsize)


fig1.savefig('overlapping_signals_A.pdf')
fig1.savefig('../../figures/overlapping_signals_A.pdf')



fig = plt.figure(figsize=(5.5, 5.5))
total_bins = 40
linewidth = 1.5
params_B_tags = ['chirp_mass_B', 'mass_ratio_B', 'theta_jn_B', 'luminosity_distance_B']
truth_params = [injection_parameters[p] for p in params_B_tags]
fig2 = corner.corner(signal_B_pe, fig = fig,
                    bins = total_bins, truths = truth_params, truth_color = 'black',
                    color = 'darkseagreen', 
                    density = 1,
                    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                   plot_density=True, 
                    plot_datapoints=False, 
                    max_n_ticks=3,
                    fill_contours=True, 
                    **{
                        "hist_kwargs": {"linewidth": linewidth, 'density': 1},         # 1D histograms
                        "contour_kwargs": {"linewidths": linewidth},     # KDE contours
                        "truth_kwargs": {"linewidth": linewidth},        # Truth lines
                    })



fig2 = corner.corner(joint_B, fig = fig,
                    labels = latex_labels_B, 
                    bins = total_bins, #truths = short_truth, truth_color = 'black',
                    color = 'royalblue', 
                    density = 1,
                    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                   plot_density=True, 
                    plot_datapoints=False, 
                    max_n_ticks=3,
                    fill_contours=True, 
                    **{
                        "hist_kwargs": {"linewidth": linewidth, 'density': 1},         # 1D histograms
                        "contour_kwargs": {"linewidths": linewidth},     # KDE contours
                        "truth_kwargs": {"linewidth": linewidth},        # Truth lines
                    })


ndim = len(params_B_tags)
import matplotlib.lines as mpllines
lines = []

lines.append(mpllines.Line2D([0], [0], color='darkseagreen'))
lines.append(mpllines.Line2D([0], [0], color='royalblue'))
labels = ['Signal B isolation', 'Signal B overlapping']

final_axes = fig2.get_axes()
final_axes[ndim - 1].legend(lines, labels, fontsize = fontsize)

fig2.savefig('overlapping_signals_B.pdf')
fig2.savefig('../../figures/overlapping_signals_B.pdf')



