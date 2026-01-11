import numpy as np
import bilby
import matplotlib.pyplot as plt


import scienceplots
plt.style.use(['science', 'bright'])

#import matplotlib as mpl

fontsize = 11
plt.rcParams['axes.labelsize'] = fontsize
plt.rc('xtick',labelsize=fontsize)
plt.rc('ytick',labelsize=fontsize)

plt.rcParams.update({
    'font.size': fontsize,              # Base font size
    'axes.titlesize': fontsize,         # Axes title
    'axes.labelsize': fontsize,         # Axes labels
    'xtick.labelsize': fontsize,        # X-axis tick labels
    'ytick.labelsize': fontsize,        # Y-axis tick labels
    'legend.fontsize': fontsize,        # Legend text
    'figure.titlesize': fontsize        # Figure title
})




injection_parameters = {
    "mass_ratio_A": 0.8231017059328386,
    "chirp_mass_A": 60.353653661893304,
    "luminosity_distance_A": 7.268571745364448 * 500 * 100,
    "dec_A": 0.6557802485585996,
    "ra_A": 2.180684709384296,
    "theta_jn_A": 1.8552485990854504,
    "psi_A": 2.317315123447218,
    "phase_A": 0.736517637678661,
    "a_1_A": 0.5874657489189719,
    "a_2_A": 0.7740290285427935,
    "tilt_1_A": 0.,
    "tilt_2_A": 0.,
    "phi_12_A": 1.4830050470970393,
    "phi_jl_A": 3.5537095497820763,
    "geocent_time_A": 3600.0,
  }


signal_A = {key.replace('_A', ''): injection_parameters[key] for key in injection_parameters.keys() if '_A' in key}

duration = 32
sampling_frequency = 4096
minimum_frequency = 20.
approximant = 'IMRPhenomXPHM'
reference_frequency = 50. 


waveform_arguments = dict(
    waveform_approximant=approximant,
    reference_frequency=reference_frequency,
    minimum_frequency=minimum_frequency,
)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,
)


ifos = bilby.gw.detector.InterferometerList(["L1"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=signal_A["geocent_time"] - duration + 2,
)

ifos[0].inject_signal(
    waveform_generator=waveform_generator, parameters=signal_A
)







ifos_A = bilby.gw.detector.InterferometerList(["L1"])
ifos_A.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=signal_A["geocent_time"] - duration + 2,
)

ifos_A[0].inject_signal(
    waveform_generator=waveform_generator, parameters=signal_A
)

fig, ax = plt.subplots(1, 1, figsize = (8.5, 4))
trigger_time = injection_parameters['geocent_time_A'] + 34e-3

ax.plot(1e3*(ifos_A[0].time_array-trigger_time), ifos_A[0].whitened_time_domain_strain, color = 'royalblue', lw = 2, alpha = 1)
ax.set_xlim(1e3*-.30, 1e3*(1 * .01))

# Make axes and tick labels invisible
ax.axis('off')

fig.savefig('signal.pdf', facecolor='none', transparent=True)
fig.savefig('signal.png', facecolor='none', transparent=True, dpi = 250)

