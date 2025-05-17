import numpy as np
import gengli
import matplotlib.pyplot as plt
import bilby
from argparse import ArgumentParser
from astropy.cosmology import Planck18
from bilby.core.utils import logger
import scienceplots
import sys
import os
import json
from tqdm import tqdm
from gwpy.timeseries import TimeSeries
sys.path.append('../../../null_stream_glitch_mitigation')
import utils
plt.style.use(['science'])

plt.rcParams['axes.labelsize'] = 11.
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)

red = (0.80, 0.36, 0.36)
green = (0.62, 0.79, 0.62)
blue = (0.27, 0.51, 0.71)

parser = ArgumentParser()

parser.add_argument('--outdir', type = str, default = 'delete_me')
parser.add_argument('--label', type = str, default='delete_me')
parser.add_argument('--minimum-frequency', type = float, default = 20.)
parser.add_argument('--seed', type = int, default = 2025)
parser.add_argument('--duration', type =float, default = 256.)
parser.add_argument('--signal-delta-t', type = float, default = 0.)
parser.add_argument('--glitch-trigger-time', type = float, default=36462)
parser.add_argument('--injection-label', type = str, help = 'Name of the GW signal we want to inject', default='gw150914')
parser.add_argument('--injection-file', type = str, default = '/data/gravwav/hnarola/null_stream_glitch_mitigation/parameter_estimation/phenomd_parameters.json')
parser.add_argument('--glitch-snr', type = float, default = 12)
parser.add_argument('--redshift', type = float, help = 'Redshift of the event', default=2.)
parser.add_argument('--gengli-seed', type = int, default=2025)
args = parser.parse_args()

bilby.core.utils.setup_logger(outdir=args.outdir, label=args.label)

plot_colorbar = False
plot_tcp = True
plot_null_stream = True
plot_timeseries = True
gengli_glitch = True
plot_signal_only_qscan = False
make_likelihood_curve = False



###### setting up some usual arguments for bilby #########
sampling_frequency = 2048.
waveform_approximant = 'IMRPhenomD'
reference_frequency = 20.

injection_parameters = json.load(open(args.injection_file, 'r'))[args.injection_label]
z = args.redshift

injection_parameters['mass_1'] = 15
injection_parameters['mass_2'] = 14

injection_parameters['mass_1'] = injection_parameters['mass_1'] * (1 + z)
injection_parameters['mass_2'] = injection_parameters['mass_2'] * (1 + z)
injection_parameters['luminosity_distance'] = 15924.566651659155*1 #Planck18.luminosity_distance(z).to('Mpc').to_value()

injection_parameters['geocent_time'] = args.glitch_trigger_time

approx_duration = bilby.gw.utils.calculate_time_to_merger(args.minimum_frequency, injection_parameters['mass_1'], injection_parameters['mass_2'], safety = 1.2)


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

ifos.inject_signal(waveform_generator=waveform_generator, parameters=injection_parameters)

if make_likelihood_curve:
    priors = bilby.gw.prior.BBHPriorDict()
    likelihood = bilby.gw.GravitationalWaveTransient(interferometers=ifos,
    waveform_generator=waveform_generator,
    priors=priors)

    truth = injection_parameters['luminosity_distance']

    xvals = np.linspace(40, 50, int(500))


    collect_likelihoods = []

    for xx in tqdm(xvals):
        p = injection_parameters.copy()
        p['mass_1'] = xx
        likelihood.parameters = p
        collect_likelihoods.append(likelihood.log_likelihood_ratio())
    

    fig, ax = plt.subplots(1, 1, figsize = (4,3))
    ax.plot(xvals, collect_likelihoods, color = 'black', lw = 1.5)
    ax.set_xlabel('GW Parameter')
    ax.set_title('Measurement', fontsize = 11)

    ax.set_ylabel('PDF')
    #ax.set_yscale('log')
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.grid(alpha = 0.5)
    fig.savefig(f'{args.outdir}/{args.label}_pe_vals.pdf', bbox_inches = 'tight')
    fig.savefig(f'../../figures/{args.label}_pe_vals.pdf', bbox_inches = 'tight')
    exit()



coloured_signal_only_data = ifos[0].time_domain_strain.copy()

################# set up a gengli glitch and injecting it into the data ############################

if gengli_glitch:
    generator = gengli.glitch_generator('L1')
    l1_to_ET_factor = 1
    glitch_snr = args.glitch_snr / l1_to_ET_factor #args.glitch_snr_ratio*np.abs(ifos[0].meta_data['matched_filter_SNR'])
    glitch_snr_timeseries = 3.5

    args.gengli_seed = 108
    _glitch = generator.get_glitch(seed = args.gengli_seed, snr = 12, srate=sampling_frequency)

    glitch = utils.center_maxima(_glitch)
    


plot_gengli_glitch = False
if plot_gengli_glitch:
    glitch_delta_t_array = np.arange(0, len(glitch), 1)*1/sampling_frequency

    fig, ax = plt.subplots(1, 1)
    ax.plot(glitch_delta_t_array, glitch)
    ax.set_ylim(-10, 10)
    fig.savefig(f'{args.outdir}/{args.label}_gengli_glitch.pdf', dpi = 250)

######################################################################

############# injecting the glitch ###############################

asd = np.array([ifos[0].frequency_array, ifos[0].amplitude_spectral_density_array]).T
white_signal_noise_timeseries = utils.coloured_timeseries_to_white_timeseries(ifos[0].time_domain_strain.copy(), asd, sampling_frequency, args.minimum_frequency, duration=args.duration)
white_signal_noise_glitch_timeseries, glitch_start_time = utils.add_glitch(glitch = glitch, noise = white_signal_noise_timeseries, t_inj = args.glitch_trigger_time - args.signal_delta_t, time_array = ifos[0].time_array)
        

coloured_signal_noise_glitch_timeseries = utils.white_timeseries_to_coloured_timeseries(white_signal_noise_glitch_timeseries, asd, sampling_frequency, args.minimum_frequency, duration=args.duration)


############### saving glitch only data ######################
if plot_signal_only_qscan:

    maximum_frequency = 512
    qrange = (16, 32)
    tres = 0.01
    fres = 0.01
    signal_only_timeseries = TimeSeries(coloured_signal_only_data, times = ifos[0].time_array, channel = 'ET1:STRAIN')
    fig, ax = plt.subplots(1, 1, figsize = (5, 3))

    _, ax, imshow = utils.plot_q_transform(signal_only_timeseries.crop(args.glitch_trigger_time-1.5, args.glitch_trigger_time+1), 
                        ax = ax, whiten = True, minimum_frequency = args.minimum_frequency, 
                        trigger_time = args.glitch_trigger_time, 
                        maximum_frequency = maximum_frequency,
                        qrange = qrange, tres = tres, fres = fres)


    ax.set_xlabel('Time [s]')
    ax.grid()

    ax.set_title('Cleaned $\\mathrm{ET}_1$', fontsize = 11)

    ax.set_ylabel('Frequency [Hz]')
    ax.set_ylim(10, 300)
    ax.set_xlim(-1.1, 0.3)

    fig.savefig(f'{args.outdir}/{args.label}_cleaned_ET1.png', bbox_inches = 'tight', dpi = 500)
    fig.savefig(f'../../figures/{args.label}_cleaned_ET1.png', bbox_inches = 'tight', dpi = 500)







if plot_tcp:
    fig, axes = plt.subplots(1, 3, figsize = (15, 3), sharey = True, sharex = True, gridspec_kw = {'wspace': 0.05})

    maximum_frequency = 512
    qrange = (16, 32)
    tres = 0.01
    fres = 0.01


    for ii, ifo in enumerate(ifos):
        if ii==0:
            timeseries_to_plot = TimeSeries(coloured_signal_noise_glitch_timeseries, times = ifos[0].time_array, channel = 'ET1:STRAIN')
        else:
            timeseries_to_plot = TimeSeries(ifo.time_domain_strain.copy(), times = ifos[0].time_array, channel = f'{ifo.name}:STRAIN')

        ax = axes[ii]
        _, ax, imshow = utils.plot_q_transform(timeseries_to_plot.crop(args.glitch_trigger_time-1.5, args.glitch_trigger_time+1), 
                        ax = ax, whiten = True, minimum_frequency = args.minimum_frequency, 
                        trigger_time = args.glitch_trigger_time, 
                        maximum_frequency = maximum_frequency,
                        qrange = qrange, tres = tres, fres = fres)


        ax.set_xlabel('Time [s]')
        ax.grid()


    axes[0].set_title('$\\mathrm{ET}_1$', fontsize = 11)
    axes[1].set_title('$\\mathrm{ET}_2$', fontsize = 11)
    axes[2].set_title('$\\mathrm{ET}_3$', fontsize = 11)
    axes[0].set_ylabel('Frequency [Hz]')
    axes[0].set_ylim(10, 300)
    axes[0].set_xlim(-1.1, 0.3)

    if plot_colorbar:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.81, 0.15, 0.01, 0.7])

        fig.colorbar(imshow, cax=cbar_ax, label = "Normalized energy")

    fig.savefig(f'{args.outdir}/{args.label}_tcp.png', bbox_inches = 'tight', dpi = 250)
    fig.savefig(f'../../figures/{args.label}_tcp.png', bbox_inches = 'tight', dpi = 250)



if plot_null_stream:
    null_stream_timeseries = coloured_signal_noise_glitch_timeseries + ifos[1].time_domain_strain + ifos[2].time_domain_strain
    timeseries_to_plot = TimeSeries(null_stream_timeseries, times = ifos[0].time_array, channel = f'{ifo.name}:STRAIN')


    maximum_frequency = 512
    qrange = (16, 32)
    tres = 0.01
    fres = 0.01



    fig, ax = plt.subplots(1, 1, figsize = (5, 3))

    _, ax, imshow = utils.plot_q_transform(timeseries_to_plot.crop(args.glitch_trigger_time-1.5, args.glitch_trigger_time+1), 
                    ax = ax, whiten = True, minimum_frequency = args.minimum_frequency, 
                    trigger_time = args.glitch_trigger_time, 
                    maximum_frequency = maximum_frequency,
                    qrange = qrange, tres = tres, fres = fres)


    ax.grid()
    ax.set_title('Null Stream', fontsize = 11)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Frequency [Hz]')
    ax.set_ylim(10, 300)
    ax.set_xlim(-1.1, 0.3)
    if plot_colorbar:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.81, 0.15, 0.02, 0.7])

        fig.colorbar(imshow, cax=cbar_ax, label = "Normalized energy")

    fig.savefig(f'{args.outdir}/{args.label}_null_stream.png', bbox_inches = 'tight', dpi = 250)
    fig.savefig(f'../../figures/{args.label}_null_stream.png', bbox_inches = 'tight', dpi = 250)    


if plot_timeseries:
    background_data = white_signal_noise_glitch_timeseries

    glitch_delta_t_array = np.arange(0, len(glitch), 1)*1/sampling_frequency

    roll = glitch_start_time - injection_parameters['geocent_time']

    fig, ax  = plt.subplots(1, 1, figsize = (4, 3))
    lw = 1.5

    ax.plot(ifos[0].time_array - injection_parameters['geocent_time'], background_data, color = 'grey', lw = 1.5, alpha = 0.4)

    
    green = (0.24, 0.9, 0.44)	#'green' #(0.62, 0.5, 0.5, 1.0)

    ax.plot(glitch_delta_t_array+roll, glitch/4, color = 'black', ls = '--', lw = 1.5)

    total_width = 0.3 + 0.1

    old_right = (0.5) / (1.1 + 0.3)
    old_left = (1.1-args.signal_delta_t) / (1.1 + 0.3)
  
    #ax.set_xlim(-0.2 - (total_width*old_left),  -0.2 + (total_width * old_right))
    ax.set_xlim(-0.7, 0.05)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Whitened Strain [$\\sigma$]')
    ax.set_ylim(-30, 30)
    ax.set_title('Reconstructed Glitch', fontsize = 11)
    #ax.set_xticks(-1*np.array([0.35, 0.25, 0.15]))
    ax.grid(alpha = 0.2)
    fig.savefig(f'{args.outdir}/{args.label}_timeseries.pdf', bbox_inches = 'tight')
    fig.savefig(f'../../figures/{args.label}_timeseries.pdf', bbox_inches = 'tight') 
    



