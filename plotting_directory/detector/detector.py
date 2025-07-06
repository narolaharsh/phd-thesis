import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scienceplots
import matplotlib.transforms as transforms

plt.rcParams['axes.labelsize'] = 11.
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)

plt.style.use(['science'])


fig, ax = plt.subplots(1, 1, figsize = (4, 4))

block_width = 0.1/2
block_height = 0.05/2

mirror_loc = [.2, 1.]
input_mirror_loc = 0.4

mirror1 = patches.Rectangle((0.2 - block_width, mirror_loc[1] - block_height), 0.1, 0.05, linewidth=1, edgecolor='black', facecolor='oldlace')
input_mirror1 = patches.Rectangle((0.2 - block_width, input_mirror_loc - block_height), 0.1, 0.05, linewidth=1, edgecolor='black', facecolor='oldlace')


mirror2 = patches.Rectangle((mirror_loc[1] - block_height, 0.2 - block_width), 0.05, 0.1, linewidth=1, edgecolor='black', facecolor='oldlace')
input_mirror2 = patches.Rectangle((input_mirror_loc - block_height, 0.2 - block_width), 0.05, 0.1, linewidth=1, edgecolor='black', facecolor='oldlace')

laser = patches.Rectangle((-0.1, 0.2 - block_height), 0.1, 0.05, linewidth=1, edgecolor='black', facecolor = 'salmon')


################## beam splitter ########################
loc = [0.2 - block_width, 0.2 - block_height]
beam_splitter = patches.Rectangle((loc[0], loc[1]), 0.1, 0.05, linewidth=1, edgecolor = 'black', facecolor = 'mintcream')
t = transforms.Affine2D().rotate_deg_around(loc[0] + 0.1/2, loc[1] + 0.05/2, 45) + ax.transData
beam_splitter.set_transform(t)
ax.add_patch(beam_splitter)
########################################################################################################################


############# add lasers ####################
N = 20
x = np.full(N, .2)
y = np.linspace(0.2, 1., N)

ax.plot(x, y, color = 'red', alpha = 0.5, zorder = 0)
ax.plot(y, x, color = 'red', alpha = 0.5, zorder = 0)

y_shrink = np.linspace(0.4, 1., N)
ax.plot(x, y_shrink, color = 'red', alpha = 0.8, zorder = 1)
ax.plot(y_shrink, x, color = 'red', alpha = 0.8, zorder = 1)

#### source laser 

ax.plot(np.linspace(-0.1, .2, N), np.full(N, .2), color = 'red', alpha = 0.5, zorder = 0)
ax.plot(np.full(N, .2), np.linspace(-0.1, .2, N),  color = 'red', alpha = 0.5, zorder = 0)

dx = 0.02
ax.annotate('', xy=(0.4+dx, 0.2-0.08), xytext=(1-dx, 0.2-.08),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1))



ax.annotate('', xy=(0.2-0.08, 0.4+dx), xytext=(0.2-0.08, 1-dx),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1))




############## Adding text ###############################
fontsize  = 7.5
#ax.text(1., .3, 'End mirror', color='black', fontsize=fontsize, ha='center', va='center')  
#ax.text(.2, 1.1, 'End mirror', color='black', fontsize=fontsize, ha='center', va='center')  
ax.text(-0.045, 0.195, 'LS', color='black', fontsize=fontsize, ha='center', va='center')  
ax.text(0.202, 0.195, 'B', color='black', fontsize=fontsize, ha='center', va='center')  
ax.text(0.4, 0.195, 'M', color='black', fontsize=fontsize, ha='center', va='center')  
ax.text(1., 0.195, 'M', color='black', fontsize=fontsize, ha='center', va='center')  

ax.text(.2, 0.395, 'M', color='black', fontsize=fontsize, ha='center', va='center')  
ax.text(0.2, 1-0.005, 'M', color='black', fontsize=fontsize, ha='center', va='center') 

ax.text(0.205, -0.1+0.005, 'P', color='black', fontsize=fontsize, ha='center', va='center')   

ax.text(0.7, 0.07, 'L', color='black', fontsize=fontsize, ha='center', va='center')
ax.text(0.07, 0.7, 'L', color='black', fontsize=fontsize, ha='center', va='center')   




##############################################################



# Add the rectangle to the Axes
ax.add_patch(mirror1)
ax.add_patch(input_mirror1)
ax.scatter(0.2, -0.1, marker = 'v', s = 150, color = 'lemonchiffon', edgecolor = 'black')
ax.add_patch(mirror2)
ax.add_patch(input_mirror2)

ax.add_patch(laser)



# Set limits and aspect ratio
ax.set_xlim(-0.2, mirror_loc[1]+0.2)
ax.set_ylim(-0.2, mirror_loc[1]+0.2)
ax.set_aspect('equal')
ax.axis('off')

fig.savefig('simple_detector.pdf')
fig.savefig('../../figures/simple_detector.pdf')
fig.savefig('../../figures/simple_detector.svg')