import numpy as np
import matplotlib.pyplot as plt

def gravitational_acceleration(m, r):
    G = 6.67430e-11  # m^3 kg^-1 s^-2, constante gravitationnelle
    return -G * m / np.linalg.norm(r)**3 * r

def verlet(r, v, dt, total_time):
    positions = [r]
    r_prev = r - v * dt + gravitational_acceleration(5.972e24, r) * dt**2 / 2
    for _ in range(int(total_time / dt)):
        r_new = 2 * r - r_prev + gravitational_acceleration(5.972e24, r) * dt**2
        r_prev, r = r, r_new
        positions.append(r)
    return np.array(positions)

def rk4(r, v, dt, total_time):
    positions = [r]
    for _ in range(int(total_time / dt)):
        k1_v = gravitational_acceleration(5.972e24, r)
        k1_r = v
        k2_v = gravitational_acceleration(5.972e24, r + 0.5 * dt * k1_r)
        k2_r = v + 0.5 * dt * k1_v
        k3_v = gravitational_acceleration(5.972e24, r + 0.5 * dt * k2_r)
        k3_r = v + 0.5 * dt * k2_v
        k4_v = gravitational_acceleration(5.972e24, r + dt * k3_r)
        k4_r = v + dt * k3_v
        v += (dt / 6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        r += (dt / 6) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
        positions.append(r)
    return np.array(positions)
import numpy as np

def initial_previous_position(r, v, dt, acceleration):
    """ Calculer la position précédente basée sur la position et la vitesse actuelles. """
    return r - v * dt + 0.5 * acceleration * dt**2

# Exemple d'utilisation
r_initial = np.array([6.371e6 + 400e3, 0], dtype=float)  # Position initiale
v_initial = np.array([0, 7660], dtype=float)  # Vitesse initiale
dt = 10  # Pas de temps
acceleration = gravitational_acceleration(1, r_initial)  # Utiliser votre fonction d'accélération

r_prev = initial_previous_position(r_initial, v_initial, dt, acceleration)

# Conditions initiales, explicitement flottantes
r_initial = np.array([6.371e6 + 400e3, 0.0], dtype=float)  # 400 km au-dessus de la surface de la Terre
v_initial = np.array([0.0, 7660.0], dtype=float)  # vitesse initiale nécessaire pour une orbite basse circulaire
dt = 10  # pas de temps en secondes
total_time = 3600  # temps total de simulation en secondes

positions_verlet = verlet(r_initial, v_initial, dt, total_time)
positions_rk4 = rk4(r_initial, v_initial, dt, total_time)

plt.figure(figsize=(10, 5))
plt.plot(positions_verlet[:, 0], positions_verlet[:, 1], label='Verlet')
plt.plot(positions_rk4[:, 0], positions_rk4[:, 1], label='RK4')
plt.title('Comparison of Verlet and RK4 Integration Methods')
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()
