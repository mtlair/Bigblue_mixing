# Pseudocode logic for 1D Momentum Decay along 15 deg Cone
import numpy as np

# Inputs
U0 = 4.0          # Eductor nozzle exit velocity (m/s)
d0 = 0.025        # Nozzle diameter or equivalent slot height (m)
theta = np.radians(15)  # Cone angle
R_tank = 1.143    # Tank radius (90 in diameter = 45 in radius = 1.143 m)
x = np.linspace(0, 1.2, 100) # Distance down the cone slope (m)

# 1. Axisymmetric wall jet growth rate (Launder & Rodi, 1983)
dy12_dx = 0.073   # Jet spreading rate

# 2. Loop down cone path length (x)
for xi in x:
    r_local = R_tank - xi * np.sin(theta)  # Shrinking radius of cone
    y12 = d0 + dy12_dx * xi                # Layer growth
    
    # Velocity decay combining jet expansion and shrinking cone perimeter
    U_max = U0 * np.sqrt(d0 / (y12 * (r_local / R_tank)))
    
    if U_max < 2.0:
        print(f"CRITICAL: Scouring velocity lost at x = {xi:.2f} meters down cone!")
        break