import numpy as np
import matplotlib.pyplot as plt
import bilby
import corner

import scienceplots
plt.style.use(['science', 'bright'])

fontsize = 10

ticksize = 8
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
data['luminosity_distance_A'] /= 1000
data['luminosity_distance_B'] /= 1000




#params = ['chirp_mass_A', 'chirp_mass_B', 'theta_jn_A', 'theta_jn_B', 'luminosity_distance_A', 'luminosity_distance_B']

params = ['chirp_mass_A', 'chirp_mass_B', 'mass_ratio_A', 'mass_ratio_B', 'luminosity_distance_A', 'luminosity_distance_B']

#latex_labels = [r'$\mathcal{M}^{\mathrm{A}}_c[M_{\odot}]$', r'$\mathcal{M}^{\mathrm{B}}_c[M_{\odot}]$', r'${\theta}^{\mathrm{A}}_{\mathrm{JN}}$', r'${\theta}^{\mathrm{B}}_{\mathrm{JN}}$', r'$d^{\mathrm{A}}_{\mathrm{L}}[\rm{Gpc}]$', r'$d^{\mathrm{B}}_{\mathrm{L}}[\rm{Gpc}]$']
latex_labels = [r'$\mathcal{M}^{\mathrm{A}}_c[M_{\odot}]$', r'$\mathcal{M}^{\mathrm{B}}_c[M_{\odot}]$', r'$q^{\mathrm{A}}$', r'$q^\mathrm{B}$', r'$d^{\mathrm{A}}_{\mathrm{L}}[\rm{Gpc}]$', r'$d^{\mathrm{B}}_{\mathrm{L}}[\rm{Gpc}]$']


prod_pe = np.vstack([np.array(data[p]) for p in params]).T


fig = plt.figure(figsize=(6.5, 6.5))
total_bins = 40
linewidth = 1.25

fig = corner.corner(prod_pe, fig = fig,
                    labels = latex_labels, 
                    bins = total_bins, #truths = short_truth, truth_color = 'black',
                    color = 'royalblue', 
                    density = 1, 
                    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                   plot_density=True, 
                    plot_datapoints=False, 
                    max_n_ticks=3,
                    fill_contours=True, 
                    **{
                        "hist_kwargs": {"linewidth": linewidth},         # 1D histograms
                        "contour_kwargs": {"linewidths": linewidth},     # KDE contours
                        "truth_kwargs": {"linewidth": linewidth},        # Truth lines
                    })
fig.savefig('overlapping_signals.pdf')
fig.savefig('../../figures/overlapping_signals.pdf')


