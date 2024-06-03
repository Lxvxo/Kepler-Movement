import numpy as np
import matplotlib.pyplot as plt

def euler(f, x0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    x = np.zeros((len(t), len(x0)))
    x[0] = x0
    for i in range(1, len(t)):
        x[i] = x[i-1] + dt * f(t[i-1], x[i-1])
    return t, x

def runge_kutta2(f, x0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    x = np.zeros((len(t), len(x0)))
    x[0] = x0
    for i in range(1, len(t)):
        k1 = dt * f(t[i-1], x[i-1])
        k2 = dt * f(t[i-1] + dt/2, x[i-1] + k1/2)
        x[i] = x[i-1] + k2
    return t, x

# Fonction qui représente les équations du mouvement
def orbital_motion(t, y):
    G, M = 4 * np.pi**2, 1  # Constantes
    r = np.sqrt(y[0]**2 + y[1]**2)
    return np.array([y[2], y[3], -G*M*y[0]/r**3, -G*M*y[1]/r**3])

# Conditions initiales
x0 = [1, 0, 0, 2*np.pi]
t0, tf, dt = 0, 1, 0.001  # De 0 à 2 périodes avec un dt de 0.01

# Simulation
t_euler, x_euler = euler(orbital_motion, x0, t0, tf, dt)
t_rk2, x_rk2 = runge_kutta2(orbital_motion, x0, t0, tf, dt)

# Tracé des résultats
plt.figure(figsize=(10, 5))
plt.plot(x_euler[:,0], x_euler[:,1], label='Euler')
plt.plot(x_rk2[:,0], x_rk2[:,1], label='Runge-Kutta 2')
plt.legend()
plt.title('Comparaison des trajectoires Orbitales')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.grid(True)
plt.show()
