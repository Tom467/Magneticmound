import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft

# ==========================================================
# PARAMETERS
# ==========================================================

L = 80
N_grains = 600
g = 0.01
kappa = 0.5
steps = 4000
lambda_values = np.linspace(0, 3, 15)

# ==========================================================
# ENERGY FUNCTION
# ==========================================================

def energy(h, Lambda):
    E = 0.0
    
    # Gravity
    E += g * np.sum(h**2)
    
    # Surface smoothing
    E += kappa * np.sum((np.roll(h, -1) - h)**2)
    
    # Dipolar long-range repulsion
    for i in range(L):
        for j in range(i+1, L):
            r = abs(i-j)
            E += Lambda * h[i]*h[j] / r**3
            
    return E

# ==========================================================
# RELAXATION
# ==========================================================

def relax_system(Lambda):
    h = np.zeros(L, dtype=int)
    
    # initial mound
    for _ in range(N_grains):
        i = L//2 + np.random.randint(-5, 6)
        h[i] += 1

    for _ in range(steps):
        i = np.random.randint(1, L-1)
        if h[i] > 0:
            direction = np.random.choice([-1,1])
            new_h = h.copy()
            new_h[i] -= 1
            new_h[i+direction] += 1
            
            if energy(new_h, Lambda) < energy(h, Lambda):
                h = new_h

    return h

# ==========================================================
# ORDER PARAMETER
# ==========================================================

def column_order_parameter(h):
    spectrum = np.abs(fft(h))
    return np.max(spectrum[1:len(h)//2])

# ==========================================================
# SWEEP LAMBDA
# ==========================================================

order_parameters = []
profiles = []

print("Running simulations...")

for Lambda in lambda_values:
    print(f"Lambda = {Lambda:.2f}")
    h = relax_system(Lambda)
    profiles.append(h.copy())
    M = column_order_parameter(h)
    order_parameters.append(M)

order_parameters = np.array(order_parameters)

# Estimate critical Lambda from gradient
gradient = np.gradient(order_parameters)
crit_index = np.argmax(gradient)
Lambda_c_sim = lambda_values[crit_index]

print("\nEstimated Lambda_c (simulation) ≈", Lambda_c_sim)

# ==========================================================
# COMPUTE A FROM DISCRETE DIPOLE KERNEL
# ==========================================================

print("\nComputing A from discrete kernel...")

L_kernel = 300
r = np.arange(1, L_kernel)

def D(k):
    return np.sum(np.cos(k*r)/r**3)

k_vals = np.linspace(0, 0.4, 200)
D_vals = np.array([D(k) for k in k_vals])

# Quadratic fit near k=0
fit = np.polyfit(k_vals[:20], D_vals[:20], 2)
A = -fit[0]   # D ≈ D0 - A k^2

print("Estimated A =", A)

# Theoretical critical Lambda
Lambda_c_theory = 4*kappa/A

print("Predicted Lambda_c (theory) ≈", Lambda_c_theory)

# ==========================================================
# PLOTS
# ==========================================================
a
plt.figure()
plt.plot(lambda_values, order_parameters, marker='o')
plt.axvline(Lambda_c_sim, linestyle='--', label='Simulated Λc')
plt.axvline(Lambda_c_theory, linestyle=':', label='Theory Λc')
plt.xlabel("Dipolar strength")
plt.ylabel("Order parameter")
#plt.title("Field-driven instability")
plt.legend()
plt.show()

# Example profiles
plt.figure(figsize=(12,4))
for idx in [0, len(lambda_values)//2, -1]:
    plt.plot(profiles[idx], label=f"Lambda={lambda_values[idx]:.2f}")
    plt.xlabel("Horizontal position")
plt.ylabel("Height of mound")
plt.legend()
#plt.title("Mound profiles")
plt.show()
