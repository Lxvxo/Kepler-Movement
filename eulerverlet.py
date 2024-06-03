import numpy as np
import matplotlib.pyplot as plt

def gravitational_force(m1, m2, r):
    """Calculate the gravitational force exerted on m2 by m1 at distance r."""
    G = 6.67430e-11  # gravitational constant in m^3 kg^-1 s^-2
    return -G * m1 * m2 / np.linalg.norm(r)**3 * r

def energy(m1, m2, r, v):
    """Calculate the total mechanical energy of the system (kinetic + potential)."""
    kinetic = 0.5 * m2 * np.linalg.norm(v)**2
    potential = -np.linalg.norm(gravitational_force(m1, m2, r)) * np.linalg.norm(r)
    return kinetic + potential

def verlet(r, v, dt, m1, m2, total_time):
    """Integrate using the Verlet method."""
    positions = [r]
    energies = []
    a = gravitational_force(m1, m2, r) / m2
    r_prev = r - v * dt + 0.5 * a * dt**2
    for _ in range(int(total_time / dt)):
        r_new = 2 * r - r_prev + gravitational_force(m1, m2, r) / m2 * dt**2
        v = (r_new - r_prev) / (2 * dt)  # Estimate velocity for energy calculation
        r_prev, r = r, r_new
        positions.append(r)
        energies.append(energy(m1, m2, r, v))
    return positions, energies

def euler(r, v, dt, m1, m2, total_time):
    """Integrate using the Euler method."""
    positions = [r]
    energies = []
    for _ in range(int(total_time / dt)):
        a = gravitational_force(m1, m2, r) / m2
        r = r + v * dt
        v = v + a * dt
        positions.append(r)
        energies.append(energy(m1, m2, r, v))
    return positions, energies

# Simulation parameters
m1 = 5.972e24  # Earth mass in kg
m2 = 1000  # Satellite mass in kg
r_initial = np.array([6.371e6 + 400e3, 0], dtype=float)  # initial position in meters
v_initial = np.array([0, 7660], dtype=float)  # initial velocity for a low Earth orbit
dt = 10  # time step in seconds
total_time = 3600  # total simulation time in seconds

# Run simulations
positions_verlet, energies_verlet = verlet(r_initial, v_initial, dt, m1, m2, total_time)
positions_euler, energies_euler = euler(r_initial, v_initial, dt, m1, m2, total_time)

# Plotting results
plt.figure(figsize=(12, 6))

# Plot energy conservation
plt.subplot(1, 2, 1)
plt.plot(energies_verlet, label='Verlet', linestyle='-')
plt.plot(energies_euler, label='Euler', linestyle='--')
plt.title('Energy Conservation: Verlet vs Euler')
plt.xlabel('Time Steps')
plt.ylabel('Total Energy (Joules)')
plt.legend()
plt.grid(True)

# Plot orbits
plt.subplot(1, 2, 2)
positions_verlet = np.array(positions_verlet)
positions_euler = np.array(positions_euler)
plt.plot(positions_verlet[:, 0], positions_verlet[:, 1], label='Verlet')
plt.plot(positions_euler[:, 0], positions_euler[:, 1], label='Euler')
plt.title('Orbital Paths: Verlet vs Euler')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.legend()
plt.axis('equal')

plt.tight_layout()
plt.show()
