import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import matplotlib.patheffects as path_effects
from matplotlib import transforms
from PIL import Image

plt.style.use(['science'])


np.random.seed(22)
def add_equilateral_triangle(ax, x, y, size, aspect_ratio, rotation=0, color='lightgreen', linewidth=2):
    """
    Add an equilateral triangle to the plot.

    Parameters:
    -----------
    ax : matplotlib axes
        The axes to draw on
    x, y : float
        Center position of the triangle (in axis coordinates 0-1)
    size : float
        Side length of the triangle (in y-direction)
    aspect_ratio : float
        Width/height ratio of the page (width_mm/length_mm)
    rotation : float
        Rotation angle in degrees (default: 0)
    color : str
        Border color (default: 'white')
    linewidth : float
        Width of the border line (default: 2)
    """
    # Height of equilateral triangle in y-direction
    height = size * np.sqrt(3) / 2

    # Vertices of equilateral triangle (pointing up) relative to center
    # Define in "square" coordinates first
    vertices = np.array([
        [0, 2*height/3],           # Top vertex
        [-size/2, -height/3],      # Bottom left
        [size/2, -height/3],       # Bottom right
        [0, 2*height/3]            # Close the triangle
    ])

    # Apply rotation in square coordinates
    if rotation != 0:
        angle_rad = np.radians(rotation)
        cos_a = np.cos(angle_rad)
        sin_a = np.sin(angle_rad)
        rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        vertices = vertices @ rotation_matrix.T

    # Scale x-coordinates for aspect ratio
    vertices[:, 0] /= aspect_ratio

    # Translate to position
    vertices[:, 0] += x
    vertices[:, 1] += y

    ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=linewidth,
            transform=ax.transAxes, alpha = 0.2, zorder = 0)

def add_v_shape(ax, x, y, size, aspect_ratio, rotation=0, color='lightgreen', linewidth=2):
    """
    Add a V-shape with 90-degree angle to the plot.

    Parameters:
    -----------
    ax : matplotlib axes
        The axes to draw on
    x, y : float
        Center position of the V (in axis coordinates 0-1)
    size : float
        Length of each arm of the V
    aspect_ratio : float
        Width/height ratio of the page (width_mm/length_mm)
    rotation : float
        Rotation angle in degrees (default: 0)
    color : str
        Border color (default: 'white')
    linewidth : float
        Width of the border line (default: 2)
    """
    # V-shape with 90-degree angle (45 degrees on each side from vertical)
    # Define vertices from bottom-left, to center point, to bottom-right
    vertices = np.array([
        [-size/np.sqrt(2), -size/np.sqrt(2)],  # Left arm end
        [0, 0],                                  # Center point (apex)
        [size/np.sqrt(2), -size/np.sqrt(2)]     # Right arm end
    ])

    # Apply rotation in square coordinates
    if rotation != 0:
        angle_rad = np.radians(rotation)
        cos_a = np.cos(angle_rad)
        sin_a = np.sin(angle_rad)
        rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        vertices = vertices @ rotation_matrix.T

    # Scale x-coordinates for aspect ratio
    vertices[:, 0] /= aspect_ratio

    # Translate to position
    vertices[:, 0] += x
    vertices[:, 1] += y

    ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=linewidth,
            transform=ax.transAxes, alpha = 0.2, zorder = 0)

def add_image(ax, image_path, x, y, width, height, rotation=0, zorder=10, alpha=1.0):
    """
    Add an image (PNG/PDF) to the plot at a specified location with rotation.

    Parameters:
    -----------
    ax : matplotlib axes
        The axes to draw on
    image_path : str
        Path to the image file (e.g., './signal.png')
    x, y : float
        Center position in axis coordinates (0-1)
    width, height : float
        Width and height in axis coordinates (0-1)
    rotation : float
        Rotation angle in degrees (default: 0)
    zorder : int
        Drawing order (default: 10, above background shapes)
    alpha : float
        Transparency (default: 1.0)
    """
    # Load the image
    img = plt.imread(image_path)

    # Create extent for the image centered at (x, y)
    extent = [x - width/2, x + width/2, y - height/2, y + height/2]

    # Display the image
    im = ax.imshow(img, extent=extent, aspect='auto', zorder=zorder,
                   transform=ax.transAxes, interpolation='bilinear', alpha=alpha)

    # Apply rotation around the center point
    if rotation != 0:
        # Create a rotation transform around the center point (x, y) in axes coordinates
        trans_data = transforms.Affine2D().rotate_deg_around(x, y, rotation) + ax.transAxes
        im.set_transform(trans_data)

length_mm = 240
width_mm = 63.5#170

# Convert mm to inches (1 inch = 25.4 mm)
width_inch = width_mm / 25.4
length_inch = length_mm / 25.4

fig, ax = plt.subplots(1, 1, figsize = (width_inch, length_inch))
backfig, backax = plt.subplots(1, 1, figsize = (width_inch, length_inch))

# Dark blue background (RGB normalized to 0-1)
fig.patch.set_facecolor('black')
ax.set_facecolor('black')

# Setup background figure
backfig.patch.set_facecolor('black')
backax.set_facecolor('black')
backax.axis('off')

invite = "Invitation \n to \n attend \n my thesis \n defense "
title = "\n Precision \n science \n with \n Einstein \n Telescope"
author = "Harsh Narola"
address = "12th February, 2026 \n 16.15 University Hall \n Domplein 29 \n 3512 JE, Utrecht"

# Remove axes
ax.axis('off')

# Add text elements in white with highlighting effect
title_text = ax.text(0.5, 0.8, title, fontsize=28, color='white', ha='center', va='center',
                     weight='bold', transform=ax.transAxes, zorder = 20)
title_text.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground='black', alpha=0.9)])

author_text = ax.text(0.5, 0.3, author, fontsize=17, color='white', ha='center', va='center',
                      transform=ax.transAxes, zorder = 20)
author_text.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground='black', alpha=0.9)])

address_text = ax.text(0.5, 0.1, address, fontsize=17, color='white', ha='center', va='center',
                      transform=ax.transAxes, zorder = 20)
address_text.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground='black', alpha=0.9)])

# Calculate aspect ratio
aspect_ratio = width_mm / length_mm

# Configuration for background shapes
N = 100
shape_size = 0.05

# Add triangles with beta distribution (concentrated at edges)
triangle_locations = np.random.beta(0.2, 0.2, size=(2, N))
triangle_rotations = np.random.uniform(0, 360, N)
for ii in range(N):
    add_equilateral_triangle(ax, triangle_locations[0, ii], triangle_locations[1, ii],
                            size=shape_size, aspect_ratio=aspect_ratio,
                            rotation=triangle_rotations[ii])

# Add V-shapes with beta distribution
v_locations = np.random.beta(0.2, 0.2, size=(2, N))
v_rotations = np.random.uniform(0, 360, N)
for ii in range(N):
    add_v_shape(ax, v_locations[0, ii], v_locations[1, ii],
                size=shape_size, aspect_ratio=aspect_ratio,
                rotation=v_rotations[ii])

# Add signal waveform images
signal_locations = np.random.beta(0.2, 0.2, size=(2, N))
signal_rotations = np.random.uniform(-90, 90, N)
for ii in range(N):
    add_image(ax, "./signal.png", x=signal_locations[0, ii], y=signal_locations[1, ii],
              width=0.8, height=0.03, rotation=signal_rotations[ii], zorder = 20)
    

# Add glitch waveform images
glitch_locations = np.random.beta(0.2, 0.2, size=(2, N))
glitch_rotations = np.random.uniform(0, 360, N)
for ii in range(N):
    rand_int = np.random.randint(1, 6)  # Random integer from 1 to 5
    glitch_name = f"./glitch_{rand_int}.png"
    add_image(ax, glitch_name, x=glitch_locations[0, ii], y=glitch_locations[1, ii],
              width=0.8, height=0.1, rotation=glitch_rotations[ii], alpha = 0.5, zorder=20)


# Add gaussian images
gaussian_locations = np.random.beta(0.2, 0.2, size=(2, N))
gaussian_rotations = np.random.uniform(-45, 45, N)
for ii in range(N):
    add_image(ax, "./gaussian.png", x=gaussian_locations[0, ii], y=gaussian_locations[1, ii],
              width=0.2, height=0.05, rotation=gaussian_rotations[ii])

fig.savefig('bookmark.pdf', facecolor=fig.get_facecolor())
