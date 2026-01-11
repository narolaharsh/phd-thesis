import numpy as np
import gengli
import scienceplots
import matplotlib.pyplot as plt
plt.style.use(['science'])





def gaussian(x, mu, sigma):
    fac1 = (x-mu)**2
    fac2 = sigma**2

    exp_factor = np.exp(-0.5 * fac1 /fac2)

    exp_factor *= 1/(np.sqrt(2*np.pi) * sigma)
    return exp_factor



x = np.arange(-3, 3, 0.01)
y = gaussian(x, 0, 1)

fig, ax = plt.subplots(1, 1)
ax.plot(x, y, color = 'lightgreen')
ax.axis("off")
fig.savefig('gaussian.png', transparent=True, dpi=500)



exit()



generator = gengli.glitch_generator('L1')
seed = 2026
sampling_frequency = 4096



# Generate and save five different glitches
for i in range(5):
    glitch = generator.get_glitch(seed=seed + i, snr=12, srate=sampling_frequency)

    fig, ax = plt.subplots(1, 1)
    ax.plot(np.arange(0, len(glitch), 1)/sampling_frequency, glitch, color='salmon')
    ax.axis("off")

    fig.savefig(f"glitch_{i+1}.png", transparent=True, dpi=500)
    plt.close(fig)

