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
    "mass_ratio_A": 0.2231017059328386,
    "mass_ratio_B": 0.37172837551203913,
    "chirp_mass_A": 22.353653661893304,
    "chirp_mass_B": 47.64517481623601,
    "luminosity_distance_A": 7.268571745364448 * 500,
    "luminosity_distance_B": 50.09434936920892 * 100,
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
    "geocent_time_A": 3600.0,
    "geocent_time_B": 3600.013082053778854186
  }


signal_A = {key.replace('_A', ''): injection_parameters[key] for key in injection_parameters.keys() if '_A' in key}
signal_B = {key.replace('_B', ''): injection_parameters[key] for key in injection_parameters.keys() if '_B' in key}

duration = 32
sampling_frequency = 4096
minimum_frequency = 5.
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


ifos = bilby.gw.detector.InterferometerList(["ET"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=signal_A["geocent_time"] - duration + 2,
)

ifos[0].inject_signal(
    waveform_generator=waveform_generator, parameters=signal_A
)

ifos[0].inject_signal(waveform_generator = waveform_generator, 
                        parameters = signal_B)






ifos_A = bilby.gw.detector.InterferometerList(["ET"])
ifos_A.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=signal_A["geocent_time"] - duration + 2,
)

ifos_A[0].inject_signal(
    waveform_generator=waveform_generator, parameters=signal_A
)


ifos_B= bilby.gw.detector.InterferometerList(["ET"])
ifos_B.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=signal_A["geocent_time"] - duration + 2,
)

ifos_B[0].inject_signal(
    waveform_generator=waveform_generator, parameters=signal_B
)



fig, ax = plt.subplots(1, 1, figsize = (8, 4))
trigger_time = injection_parameters['geocent_time_A'] + 34e-3

ax.plot(1e3*(ifos_A[0].time_array-trigger_time), ifos_A[0].whitened_time_domain_strain, color = 'royalblue', lw = 2, alpha = 1, label = 'Signal 1')
ax.plot(1e3*(ifos_B[0].time_array-trigger_time), ifos_B[0].whitened_time_domain_strain, color = 'darkseagreen', lw = 2, alpha = 1, label = 'Signal 2')
ax.plot(1e3*(ifos[0].time_array-trigger_time), ifos[0].whitened_time_domain_strain, color = 'black', lw = 2, ls = '--', label = 'Signal 1 + Signal 2')
ax.legend(frameon = True, fancybox = True)
ax.set_xlim(1e3*-0.25, 1e3*.03)
ax.set_ylabel('Whitened strain [$\sigma$]')
ax.set_xlabel('Time[ms]')
fig.savefig('vis_overlap_signals.pdf')
fig.savefig('../../figures/vis_overlap_signals.pdf')