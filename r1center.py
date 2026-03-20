#Centered on the Milky Way 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

G = 6.67430e-11

# Masses (kg)
m1 = 1.5e42  # Milky Way
m2 = 1.0e42  # Andromeda

# Initial positions (m)
r1 = np.array([0.0, 0.0])
r2 = np.array([2.0e21, 0.0])  # Current distance between galaxies

# Initial velocities (m/s)
v1 = np.array([0.0, 0.0])
v2 = np.array([0.0, -4e5])  # 400000 m/s towards MW typical value for galaxy that is gravitationally bound

# Time settings
t_max = 2.5e16  # Total time (792 million years) needed smaller time scale for simulation to work
dt = 5e12       # Time step (160000 years) balance for accuracy

# RK4 integration
def rk4_step(r1, v1, r2, v2, dt):
    def force(r1, r2):
        r = r2 - r1
        dist = np.linalg.norm(r)
        if dist == 0:
            return np.zeros(2)
        return G * m1 * m2 * r / dist**3

    def derivatives(r1, v1, r2, v2):
        f = force(r1, r2)
        a1 = f / m1
        a2 = -f / m2
        return np.hstack((v1, a1, v2, a2))

    k1 = dt * derivatives(r1, v1, r2, v2)
    k2 = dt * derivatives(r1 + 0.5 * k1[0:2], v1 + 0.5 * k1[2:4],
                          r2 + 0.5 * k1[4:6], v2 + 0.5 * k1[6:8])
    k3 = dt * derivatives(r1 + 0.5 * k2[0:2], v1 + 0.5 * k2[2:4],
                          r2 + 0.5 * k2[4:6], v2 + 0.5 * k2[6:8])
    k4 = dt * derivatives(r1 + k3[0:2], v1 + k3[2:4],
                          r2 + k3[4:6], v2 + k3[6:8])

    r1_new = r1 + (k1[0:2] + 2 * k2[0:2] + 2 * k3[0:2] + k4[0:2]) / 6
    v1_new = v1 + (k1[2:4] + 2 * k2[2:4] + 2 * k3[2:4] + k4[2:4]) / 6
    r2_new = r2 + (k1[4:6] + 2 * k2[4:6] + 2 * k3[4:6] + k4[4:6]) / 6
    v2_new = v2 + (k1[6:8] + 2 * k2[6:8] + 2 * k3[6:8] + k4[6:8]) / 6

    return r1_new, v1_new, r2_new, v2_new

# Simulation loop
positions_r1 = [r1.copy()]
positions_r2 = [r2.copy()]
t = 0.0

while t < t_max:
    r1, v1, r2, v2 = rk4_step(r1, v1, r2, v2, dt)
    positions_r1.append(r1.copy())
    positions_r2.append(r2.copy())
    t += dt

positions_r1 = np.array(positions_r1)
positions_r2 = np.array(positions_r2)

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlabel("x-position (m)")
ax.set_ylabel("y-position (m)")
ax.set_title("Milky Way–Andromeda Collision (Accelerated Timescale)")
ax.grid(True)
ax.set_aspect('equal')

# Trajectories and current positions
mw_path, = ax.plot([], [], 'b-', lw=1, label='Milky Way')
and_path, = ax.plot([], [], 'r-', lw=1, label='Andromeda')
mw_dot, = ax.plot([], [], 'bo')
and_dot, = ax.plot([], [], 'ro')

plt.legend()

# Animation
def update(frame):
    # Update trajectories
    mw_path.set_data(positions_r1[:frame, 0], positions_r1[:frame, 1])
    and_path.set_data(positions_r2[:frame, 0], positions_r2[:frame, 1])
    mw_dot.set_data([positions_r1[frame, 0]], [positions_r1[frame, 1]])
    and_dot.set_data([positions_r2[frame, 0]], [positions_r2[frame, 1]])

    # Center plot on Milky Way
    center_x, center_y = positions_r1[frame]
    width = 3e21
    ax.set_xlim(center_x - width / 1, center_x + width / 1)
    ax.set_ylim(center_y - width / 1, center_y + width / 1)

    return mw_path, and_path, mw_dot, and_dot

# Animation
ani = FuncAnimation(fig, update, frames=len(positions_r1),
                    blit=True, interval=1, repeat=False)

plt.show()