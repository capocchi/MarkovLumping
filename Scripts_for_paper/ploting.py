import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the data for the bar graph
matrix = np.random.rand(11, 4, 4)  # Replace with your own data
reduction_methods = ['method1', 'method2', 'method3', 'method4'] * 11
value_criteria = np.tile(np.arange(4), 11 * len(reduction_methods))[:176].reshape(matrix.shape)

# Create the 3D bar graph
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xpos, ypos = np.meshgrid(np.arange(matrix.shape[1]), np.arange(matrix.shape[2]))
xpos = xpos.flatten()   # Convert to 1D array
ypos = ypos.flatten()   # Convert to 1D array
zpos = np.zeros_like(xpos)
dx = 0.5 * np.ones_like(zpos)
dy = dx.copy()
dz = matrix.flatten()

colors = ['r', 'g', 'b', 'c']  # Replace with your own color scheme
for i, method in enumerate(reduction_methods):
    ax.bar3d(xpos + i * 0.25, ypos, zpos, dx, dy, dz[i::len(reduction_methods)],
             color=colors[i % len(colors)], alpha=0.8)

# Set the labels and title for the bar graph
ax.set_xlabel('Reduction Method')
ax.set_ylabel('Value Criteria')
ax.set_zlabel('Matrix')
ax.set_xticks(np.arange(matrix.shape[1]) + len(reduction_methods) / 2)
ax.set_xticklabels(reduction_methods[:matrix.shape[1]])
ax.set_yticks(np.arange(4) + 0.5)
ax.set_yticklabels(np.arange(4))
ax.set_title('3D Bar Graph Example')

# Display the bar graph
plt.show()
