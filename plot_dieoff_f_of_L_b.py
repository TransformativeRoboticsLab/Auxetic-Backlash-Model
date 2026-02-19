import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Generate data
backlash = np.linspace(0.001, 0.20, 100)
cell_sizes = np.linspace(10, 50, 30)

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Calculate die-off values for each cell size and create curves
die_off_values = []
for cell_size in cell_sizes:
    # Model: die-off decreases exponentially with backlash
    # and increases with cell size
    die_off = cell_size * 3.5 * np.exp(-20 * backlash)
    die_off_values.append(die_off[50])  # Use middle value for color
    
# Normalize die-off values for colormap
die_off_array = np.array(die_off_values)
norm = plt.Normalize(vmin=20, vmax=140)
cmap = cm.viridis

# Plot curves
for i, cell_size in enumerate(cell_sizes):
    die_off = cell_size * 3.5 * np.exp(-20 * backlash)
    color = cmap(norm(die_off_values[i]))
    ax.plot(backlash, die_off, color=color, linewidth=1.5)

# Labels and formatting
ax.set_xlabel('Normalized Backlash [n.d.]', fontsize=11)
ax.set_ylabel('Cell Size [mm]', fontsize=11)
ax.set_title('Die-off Distance given Cell Size and Backlash', fontsize=12)
ax.set_xlim(0, 0.20)
ax.set_ylim(10, 50)
ax.grid(False)

# Add colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Die-off [pix]', fontsize=11)

plt.tight_layout()
plt.show()