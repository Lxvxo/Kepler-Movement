# Kepler Movement Simulation

## Introduction
This project simulates the motion of a planet around a star, following Kepler's laws of planetary motion. The simulation uses various numerical methods to solve the differential equations governing the system. The primary methods implemented include the Euler method, an implicit method, and the Runge-Kutta 4th order method (RK4).

### Example of trajectory for distinct initial conditions: 

![Capture](https://github.com/user-attachments/assets/6c8a1a7b-f7d2-4ee7-8d52-8796d9984729)

### Example of comparison of methods:

![Capture2](https://github.com/user-attachments/assets/474a7c2f-5f07-480e-a80b-759069d2f965)

## Methods

### 1. Implicit Method (`Kepler_implicite.py`)
This method uses an implicit approach to solve the differential equations. It iteratively computes the positions and velocities of a planet based on gravitational forces. This method is generally more stable but can be computationally intensive.

### 2. Euler Method (`euler.py`)
The Euler method is a straightforward numerical approach to solve ordinary differential equations. It updates the position and velocity at each time step using the current values, making it simple but less accurate for larger time steps.

### 3. Runge-Kutta 4th Order Method (`rk4.py`)
The RK4 method is a more accurate method for solving differential equations. It calculates intermediate steps to provide a more accurate update for position and velocity, making it more reliable than the Euler method, especially for larger time steps.

## Implementation

- **Gravitational Constant:** The simulation assumes a normalized gravitational constant \(G = 4\pi^2\) in astronomical units.
- **Time Steps:** The simulation runs over a fixed number of iterations, adjusting positions and velocities at each step.
- **Initial Conditions:** The planet starts at a specific position and velocity, typically set to simulate an orbit.

## Visualization
The simulation results can be visualized using the Matplotlib library. The scripts generate plots of the planetary orbits, showcasing how different numerical methods affect the accuracy and stability of the simulation.

## How to Run
1. Clone the repository.
2. Install the required Python libraries using `pip install -r requirements.txt`.
3. Run any of the scripts (`Kepler_implicite.py`, `euler.py`, or `rk4.py`) to simulate the orbit.
4. View the generated plots to analyze the results.

## Disclaimer
> [!IMPORTANT]
> In this project, other methods such as "Verlet" or other files used only to compare the methods are also used. To view this and learn more about the mathematical method do not hesitate to take a look at the file "KEPLER.pdf" which is basically a support for an oral.
