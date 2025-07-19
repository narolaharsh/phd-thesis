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


injection_parameters = {'chirp_mass': 15.121562495409771, 
                        'mass_ratio': 0.24993916746465974, 
                        'a_1': 0.4122192199971295, 'a_2': 0.9602728646975162, 
                        'phi_12': 5.449220181039508, 'phi_jl': 0.11262825760945509, 'tilt_1': 0.18058635171182916, 'tilt_2': 2.2123761042070944, 
                        'ra': 3.819356832760839, 'dec': 0.6216145465611382, 'luminosity_distance': 724.4973450195514, 
                        'theta_jn': 0.7908280708598625, 'phase': 2.0890792533810787, 'psi': 2.5259323400438882, 'geocent_time': 1239082262.1819956}


"""

injection_parameters = {
    "mass_ratio": 8,
    'chirp_mass': 36,
    "luminosity_distance": 7.268571745364448 * 500,
    "dec": 0.6557802485585996,
    "ra": 2.180684709384296,
    "theta_jn": 1.8552485990854504,
    "psi": 2.317315123447218,
    "phase": 0.736517637678661,
    "a_1": 0.3874657489189719,
    "a_2": 0.8740290285427935,
    "tilt_1": 1.8400317067401333,
    "tilt_2": 1.4431533314845897,
    "phi_12": 1.4830050470970393,
    "phi_jl": 3.5537095497820763,
    "geocent_time": 3600.0,
  }
"""

phenomd_parameters = injection_parameters.copy()
pop_params = ['tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'a_1', 'a_2']
for p in pop_params:
    phenomd_parameters.pop(p)
phenomd_parameters['chi_1'] = injection_parameters['a_1']
phenomd_parameters['chi_2'] = injection_parameters['a_2']

duration = 8
sampling_frequency = 4096
minimum_frequency = 20
approximant = 'IMRPhenomD'
reference_frequency = 100. 


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


ifos = bilby.gw.detector.InterferometerList(["V1"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=phenomd_parameters["geocent_time"] - duration + 2,
)

ifos[0].inject_signal(waveform_generator=waveform_generator, parameters=phenomd_parameters)
phenomd_strain = ifos[0].whitened_time_domain_strain





waveform_arguments = dict(
    waveform_approximant='IMRPhenomXPHM',
    reference_frequency=reference_frequency,
    minimum_frequency=minimum_frequency,
)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,)


ifos = bilby.gw.detector.InterferometerList(["V1"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - duration + 2,
)

ifos[0].inject_signal(
    waveform_generator=waveform_generator, parameters=injection_parameters
)
xphm_strain = ifos[0].whitened_time_domain_strain


fig, ax = plt.subplots(1, 1, figsize = (8.5, 2.5))
trigger_time = injection_parameters['geocent_time'] - 10e-3

ax.plot(1e3*(ifos[0].time_array-trigger_time), phenomd_strain, color = 'royalblue', lw = 2, alpha = 1, label = 'IMRPhenomD', zorder = 100)
ax.plot(5 + 1e3*(ifos[0].time_array-trigger_time), xphm_strain, color = 'salmon', lw = 2, label = 'IMRPhenomXPHM', alpha = 0.7, zorder = 100)

times = 1e3*(ifos[0].time_array-trigger_time)
ax.fill_between(x = [-300, -55], y1 = -2, y2 = 2, color = 'grey', alpha = 0.2, zorder = 0)
ax.fill_between(x = [-55, 0], y1 = -2, y2 = 2, color = 'grey', alpha = 0.4, zorder = 0)
ax.fill_between(x = [0, 100], y1 = -2, y2 = 2, color = 'grey', alpha = 0.6, zorder = 0)


ax.legend(frameon = True, fancybox = True)
ax.set_xlim(1e3*-0.30, 1e3*.07)
ax.set_ylim(-1.75, 1.75)
ax.set_ylabel('Whitened strain [$\sigma$]')
ax.set_xlabel('Time [ms]')
ax.text(x = -150, y = -1.6, s = 'Inspiral', color='black')
ax.text(x = -52, y = -1.6, s = 'Intermediate', color='black')
ax.text(x = 10, y = -1.6, s = 'Merger-\nringdown', color='black')


fig.savefig('phenomd_xphm.pdf')
#fig.savefig('../../figures/phenomd_xphm.pdf')