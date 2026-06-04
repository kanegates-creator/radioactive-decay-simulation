# --- Libraries ---
import os
import math
import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
mass_sample = 40.0
molar_mass_ca47 = 46.9545465
n_avogadro = 6.0221415e23
n0 = mass_sample * n_avogadro / molar_mass_ca47
DaytoSec = 24 * 60 * 3600

# Decay constants in s^-1
lam = np.zeros(2)
lam[0] = math.log(2) / (4.536 * DaytoSec)
lam[1] = math.log(2) / (3.349 * DaytoSec)
lam2_modified = lam[1] / 10 # I made the decay constant for Scandium much smaller

# --- Simulation settings ---
t_half_parent = np.log(2) / lam[0]
t_max = 10 * t_half_parent
steps = 1000
t = np.linspace(0, t_max, steps)

# --- Population & Activity Equations ---

# Parent: Ca-47
n1 = n0 * np.exp(-lam[0] * t)

# Daughter: Sc-47
n2 = (lam[0] / (lam[1] - lam[0])) * n0 * (np.exp(-lam[0] * t) - np.exp(-lam[1] * t))
n2_modified = (lam[0] / (lam2_modified - lam[0])) * n0 * (np.exp(-lam[0] * t) - np.exp(-lam2_modified * t))

# Granddaughter: Ti-47 (Stable)
n3 = n0 * (1 + (lam[0] * np.exp(-lam[1] * t) - lam[1] * np.exp(-lam[0] * t)) / (lam[1] - lam[0]))
n3_modified = n0 * (1 + (lam[0] * np.exp(-lam2_modified * t) - lam[1] * np.exp(-lam[0] * t)) / (lam2_modified - lam[0]))

# --- Graph Configuration ---

# Population Graph
plt.figure(1)
plt.plot(t, n1, label='$^{47}$Ca (Parent)', color = 'blue')
plt.plot(t, n2, label='$^{47}$Sc (Daughter)', color = 'red')
plt.plot(t, n3, label='$^{47}$Ti (Granddaughter)', color = 'green')
plt.xlim(0, t_max)
plt.title('Nuclide Population')
plt.xlabel('Time (s)')
plt.ylabel('Nuclide Population (atoms)')
plt.legend()
plt.grid(True)
plt.show()

# Sensitivity Test
plt.figure(2)
plt.plot(t, n1, label='$^{47}$Ca (Parent)', color = 'blue')
plt.plot(t, n2_modified, label='$^{47}$Sc (Daughter Modified)', color = 'red')
plt.plot(t, n3_modified, label='$^{47}$Ti (Granddaughter Modified)', color = 'green')
plt.xlim(0, t_max)
plt.title('Sensitivity Test')
plt.xlabel('Time (s)')
plt.ylabel('Nuclide Population (atoms)')
plt.legend()
plt.grid(True)
plt.show()
