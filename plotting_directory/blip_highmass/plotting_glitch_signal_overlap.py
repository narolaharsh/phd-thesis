import numpy as np
import gengli
import matplotlib.pyplot as plt
import bilby
from argparse import ArgumentParser
from scipy.signal.windows import tukey
from scipy.interpolate import UnivariateSpline
from astropy.cosmology import Planck18
from bilby.core.utils import logger
import scienceplots
import sys
import os
import json
from gwpy.timeseries import TimeSeries
sys.path.append('../../../null_stream_glitch_mitigation/')
import utils
plt.style.use(['science'])


parser = ArgumentParser()

parser.add_argument('--outdir', type = str, default = 'delete_me')
parser.add_argument('--label', type = str, default='delete_me')
parser.add_argument('--minimum-frequency', type = float, default = 20.)
parser.add_argument('--seed', type = int, default = 2025)
parser.add_argument('--duration', type =float, default = 256.)
parser.add_argument('--signal-delta-t', type = float, default = 0.)
parser.add_argument('--glitch-trigger-time', type = float, default=36462)
parser.add_argument('--injection-label', type = str, help = 'Name of the GW signal we want to inject', default='gw150914')
parser.add_argument('--injection-file', type = str, default = '../../../null_stream_glitch_mitigation/parameter_estimation/phenomd_parameters.json')
parser.add_argument('--glitch-snr', type = float, default = 12)
parser.add_argument('--redshift', type = float, help = 'Redshift of the event', default=2.)
parser.add_argument('--gengli-seed', type = int, default=2025)
args = parser.parse_args()

bilby.core.utils.setup_logger(outdir=args.outdir, label=args.label)

###### setting up some usual arguments for bilby #########
sampling_frequency = 2048.
waveform_approximant = 'IMRPhenomD'
reference_frequency = 20.

injection_parameters = json.load(open(args.injection_file, 'r'))[args.injection_label]
z = args.redshift

injection_parameters['mass_1'] = injection_parameters['mass_1'] * (1 + z)
injection_parameters['mass_2'] = injection_parameters['mass_2'] * (1 + z)
injection_parameters['luminosity_distance'] = 15924 * 3 #Planck18.luminosity_distance(z).to('Mpc').to_value()

injection_parameters['geocent_time'] = args.glitch_trigger_time + args.signal_delta_t

print('Injection parameters', injection_parameters)

approx_duration = bilby.gw.utils.calculate_time_to_merger(args.minimum_frequency, injection_parameters['mass_1'], injection_parameters['mass_2'], safety = 1.2)
print('Approximate duration of the signal:', approx_duration)

start_time= args.glitch_trigger_time - args.duration + 3
print('Start time of the frame', start_time)

waveform_arguments = dict(
    waveform_approximant=waveform_approximant,
    reference_frequency=reference_frequency,
    minimum_frequency=args.minimum_frequency,
)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=args.duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,
)

np.random.seed(args.seed)
bilby.core.utils.random.seed(args.seed)

ifos = bilby.gw.detector.InterferometerList(["ET"])

for ifo in ifos:
    ifo.minimum_frequency = args.minimum_frequency
    ifo.maximum_frequency = .5 * sampling_frequency

ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=args.duration,
    start_time = start_time)

coloured_noise_only_timeseries = ifos[0].time_domain_strain.copy()
ifos.inject_signal(waveform_generator=waveform_generator, parameters=injection_parameters)

################# set up a gengli glitch and injecting it into the data ############################
gengli_glitch = True

print('Glitch SNR input', args.glitch_snr)
if gengli_glitch:
    generator = gengli.glitch_generator('L1')
    l1_to_ET_factor = 1
    glitch_snr = args.glitch_snr / l1_to_ET_factor #args.glitch_snr_ratio*np.abs(ifos[0].meta_data['matched_filter_SNR'])
    print('Glitch snr for ET', glitch_snr)

    _glitch = generator.get_glitch(seed = args.gengli_seed, snr = glitch_snr, srate=sampling_frequency)

    glitch = utils.center_maxima(_glitch)
    

glitch_delta_t_array = np.arange(0, len(glitch), 1)*1/sampling_frequency

fig, ax = plt.subplots(1, 1)
ax.plot(glitch_delta_t_array, glitch)
ax.set_ylim(-10, 10)
fig.savefig(f'{args.outdir}/{args.label}_gengli_glitch.pdf', dpi = 250)

store_glitch_time_series = open(f'{args.outdir}/{args.label}_gengli_glitch.dat', 'w')

for ii in range(len(glitch)):
    store_glitch_time_series.write(f"{glitch_delta_t_array[ii]}\t{glitch[ii]}\n")
store_glitch_time_series.close()
######################################################################

############# injecting the glitch ###############################

asd = np.array([ifos[0].frequency_array, ifos[0].amplitude_spectral_density_array]).T
white_noise_only_strain = utils.coloured_timeseries_to_white_timeseries(coloured_noise_only_timeseries, asd, sampling_frequency, args.minimum_frequency, duration=args.duration)
white_glitch_only_timeseries, glitch_start_time = utils.add_glitch(glitch = glitch, noise = white_noise_only_strain, t_inj = args.glitch_trigger_time, time_array = ifos[0].time_array)
        
fig.savefig(f"{args.outdir}/{args.label}_white_noise_glitch_only.pdf", dpi = 250)

coloured_glitch_only_time_array = utils.white_timeseries_to_coloured_timeseries(white_glitch_only_timeseries, asd, sampling_frequency, args.minimum_frequency, duration=args.duration)

######################################################################


############### saving glitch only data ######################
glitch_only_time_series = TimeSeries(coloured_glitch_only_time_array, times = ifos[0].time_array, channel = 'ET1:STRAIN')
signal_only_timeseries = TimeSeries(ifos[0].time_domain_strain, times = ifos[0].time_array, channel = 'ET1:STRAIN')


plt.rcParams['axes.labelsize'] = 11.
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)

###### making the q scans for the ET1 and for null stream ################
fig, axes = plt.subplots(1, 2, figsize = (8, 4), sharey = True, sharex = True, gridspec_kw = {'wspace': 0.05})

maximum_frequency = 512
qrange = (16, 32)
tres = 0.01
fres = 0.01



ax = axes[0]
_, ax, imshow = utils.plot_q_transform(signal_only_timeseries.crop(args.glitch_trigger_time-1.5, args.glitch_trigger_time+1), 
                    ax = ax, whiten = True, minimum_frequency = args.minimum_frequency, 
                    trigger_time = args.glitch_trigger_time, 
                    maximum_frequency = maximum_frequency,
                    qrange = qrange, tres = tres, fres = fres)

ax.set_xlabel('Time [s]')
ax.set_ylabel('Frequency [Hz]')
ax.grid()



ax = axes[1]
_, ax, imshow = utils.plot_q_transform(glitch_only_time_series.crop(args.glitch_trigger_time-1.5, args.glitch_trigger_time+1), 
                    ax = ax, whiten = True, minimum_frequency = args.minimum_frequency, 
                    trigger_time = args.glitch_trigger_time, 
                    maximum_frequency = maximum_frequency,
                    qrange = qrange, tres = tres, fres = fres)

ax.set_xlim(-0.8, 0.8)


ax.grid()
ax.set_ylim(20, 300)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.15, 0.015, 0.7])
ax.set_xlabel('Time [s]')

fig.colorbar(imshow, cax=cbar_ax, label = "Normalized energy")
fig.savefig(f'../../figures/{args.label}_highmass_blip_qscan.png', bbox_inches = 'tight', dpi = 2500)
#fig.savefig(f'./{args.label}_highmass_blip_qscan.png', bbox_inches = 'tight', dpi = 250)

##########################################################################

