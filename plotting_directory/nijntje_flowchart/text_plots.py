import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])



plt.rcParams['axes.labelsize'] = 11.
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)



fig, ax = plt.subplots(1, 2, figsize = (8, 3))
ax[0].set_title('2. RJMCMC', fontsize = 11)
ax[0].set_xlabel('1.Construct null stream', fontsize = 11)
ax[0].set_ylabel('3. Subtract glitch', fontsize = 11)
ax[1].set_title('4. Measure GW parameters +', fontsize = 11)

ax[1].set_xlabel('2. Unmodelled reconstruction RJMCMC', fontsize = 11)
ax[1].set_title('4. Modelled reconstruction of GW Nested sampling', fontsize = 11)

fig.savefig('./inkscape_text.pdf')
fig.savefig('././../../figures/inkscape_text.pdf')