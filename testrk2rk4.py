import numpy as np
import matplotlib.pyplot as plt

def rk4(f, y0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(1, len(t)):
        k1 = dt * f(t[i-1], y[i-1])
        k2 = dt * f(t[i-1] + 0.5*dt, y[i-1] + 0.5*k1)
        k3 = dt * f(t[i-1] + 0.5*dt, y[i-1] + 0.5*k2)
        k4 = dt * f(t[i-1] + dt, y[i-1] + k3)
        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
    return t, y
def rk2(f, y0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(1, len(t)):
        k1 = dt * f(t[i-1], y[i-1])
        k2 = dt * f(t[i-1] + 0.5*dt, y[i-1] + 0.5*k1)
        y[i] = y[i-1] + k2
    return t, y
def calculate_energy(y, m1, m2, G):
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = y.T
    # Calculate distances and velocities
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    v1 = np.sqrt(vx1**2 + vy1**2)
    v2 = np.sqrt(vx2**2 + vy2**2)
    # Kinetic energy
    ke1 = 0.5 * m1 * v1**2
    ke2 = 0.5 * m2 * v2**2
    # Potential energy
    pe = -G * m1 * m2 / r
    # Total energy
    total_energy = ke1 + ke2 + pe
    return total_energy

def two_body_motion(t, y):
    # Constants
    G = 6.67430e-11  # gravitational constant in m^3 kg^-1 s^-2
    m1 = 5.972e24    # mass of the Earth in kg
    m2 = 7.342e22    # mass of the Moon in kg
    # Unpack positions and velocities
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = y
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    # Gravitational forces
    f1x = G * m2 * (x2 - x1) / r**3
    f1y = G * m2 * (y2 - y1) / r**3
    f2x = G * m1 * (x1 - x2) / r**3
    f2y = G * m1 * (y1 - y2) / r**3
    # Derivatives
    return np.array([vx1, vy1, vx2, vy2, f1x, f1y, f2x, f2y])

# Initial conditions (positions (m), velocities (m/s))
y0 = np.array([0, 0, 384400000, 0, 0, 0, 0, 1022])  # Earth at origin, Moon 384,400 km away
t0, tf, dt = 0, 10*24*3600, 100  # simulate for 10 days with a step of 100 seconds


# Constants
G = 6.67430e-11  # gravitational constant in m^3 kg^-1 s^-2
m1 = 5.972e24    # mass of the Earth in kg
m2 = 7.342e22    # mass of the Moon in kg

# Run simulations
t_rk2, y_rk2 = rk2(two_body_motion, y0, t0, tf, dt)
t_rk4, y_rk4 = rk4(two_body_motion, y0, t0, tf, dt)

# Calculate energies
energy_rk2 = calculate_energy(y_rk2, m1, m2, G)
energy_rk4 = calculate_energy(y_rk4, m1, m2, G)

# Plotting trajectories
plt.figure(figsize=(10, 8))
plt.plot(y_rk2[:, 0], y_rk2[:, 1], 'b-', label='Earth RK2')
plt.plot(y_rk2[:, 2], y_rk2[:, 3], 'b--', label='Moon RK2')
plt.plot(y_rk4[:, 0], y_rk4[:, 1], 'r-', label='Earth RK4')
plt.plot(y_rk4[:, 2], y_rk4[:, 3], 'r--', label='Moon RK4')
plt.title('Trajectory Comparison between RK2 and RK4')
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.grid(True)
plt.axis('equal')

# Plotting energy
plt.figure(figsize=(10, 5))
plt.plot(t_rk2, energy_rk2, 'b-', label='Energy RK2')
plt.plot(t_rk4, energy_rk4, 'r-', label='Energy RK4')
plt.title('Energy Conservation Comparison')
plt.xlabel('Time (s)')
plt.ylabel('Total Energy (J)')
plt.legend()
plt.grid(True)
plt.show()
