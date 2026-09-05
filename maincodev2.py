import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# PARAMETERS
# ==========================================================

L = 80
N_grains = 600

g = 0.01
kappa = 0.5

steps = 4000

lambda_values = np.linspace(0, 3, 15)

rng = np.random.default_rng(12345)


# ==========================================================
# FINITE-LATTICE OPERATORS
# ==========================================================

def laplacian_matrix(L):
    """
    Open/finite boundary discrete Laplacian.

    No periodic wrapping is used.
    """

    A = np.zeros((L, L))

    for i in range(L):
        if i > 0:
            A[i, i-1] = -1
        if i < L-1:
            A[i, i+1] = -1

        A[i, i] = 2

    # Open boundaries:
    A[0, 0] = 1
    A[-1, -1] = 1

    return A


def dipolar_matrix(L):
    """
    Finite long-range dipolar-like kernel

        K_ij = 1 / |i-j|^3,  i != j
        K_ii = 0

    The finite lattice explicitly represents the lateral
    confinement imposed by the experimental chamber.
    """

    K = np.zeros((L, L))

    for i in range(L):
        for j in range(L):
            if i != j:
                r = abs(i - j)
                K[i, j] = 1.0 / r**3

    return K


LAP = laplacian_matrix(L)
KERNEL = dipolar_matrix(L)


# ==========================================================
# ENERGY FUNCTION
# ==========================================================

def energy(h, Lambda):
    """
    Total energy:

        E = g sum_i h_i^2
          + kappa sum_i (h_{i+1} - h_i)^2
          + Lambda sum_{i<j} h_i h_j / |i-j|^3

    The finite lattice provides the lateral boundary.
    """

    E_gravity = g * np.sum(h**2)

    E_surface = kappa * np.sum(np.diff(h)**2)

    E_magnetic = 0.0

    for i in range(L):
        for j in range(i + 1, L):
            r = j - i
            E_magnetic += Lambda * h[i] * h[j] / r**3

    return E_gravity + E_surface + E_magnetic


# ==========================================================
# INITIAL MOUND
# ==========================================================

def initial_mound():
    """
    Construct the initial finite mound by depositing grains
    near the centre of the chamber.
    """

    h = np.zeros(L, dtype=int)

    for _ in range(N_grains):
        i = L // 2 + rng.integers(-5, 6)
        h[i] += 1

    return h


# ==========================================================
# RELAXATION
# ==========================================================

def relax_system(Lambda):
    """
    Zero-temperature Monte-Carlo-like energy-lowering dynamics.

    A grain is moved only if the move lowers the total energy.

    The total mass is conserved.
    """

    h = initial_mound()

    current_energy = energy(h, Lambda)

    for _ in range(steps):

        # Do not select boundary sites because grains cannot
        # leave the finite chamber.
        i = rng.integers(1, L - 1)

        if h[i] == 0:
            continue

        direction = rng.choice([-1, 1])
        j = i + direction

        new_h = h.copy()

        new_h[i] -= 1
        new_h[j] += 1

        new_energy = energy(new_h, Lambda)

        if new_energy < current_energy:
            h = new_h
            current_energy = new_energy

    return h


# ==========================================================
# COLUMNAR / ROUGHNESS MEASURES
# ==========================================================

def roughness(h):
    """
    Nearest-neighbour roughness:

        R = 1/(L-1) sum_i (h_{i+1} - h_i)^2

    Small R  -> smooth mound
    Large R  -> strongly modulated / columnar profile
    """

    return np.mean(np.diff(h)**2)


def normalized_roughness(h):
    """
    Roughness normalized by the mean height squared.
    This helps compare profiles with slightly different
    overall height scales.
    """

    mean_h = np.mean(h)

    if mean_h == 0:
        return 0.0

    return roughness(h) / mean_h**2


def spectral_columnarity(h, cutoff_fraction=0.25):
    """
    Fraction of nonzero spectral power at relatively
    short wavelengths.

    This is a secondary measure; roughness is used as the
    primary columnarity measure.
    """

    h_fluct = h - np.mean(h)

    spectrum = np.abs(np.fft.rfft(h_fluct))**2

    if np.sum(spectrum[1:]) == 0:
        return 0.0

    k = np.arange(len(spectrum))

    cutoff = int(cutoff_fraction * len(k))

    high_k = np.sum(spectrum[cutoff:])

    total = np.sum(spectrum[1:])

    return high_k / total


# ==========================================================
# EXACT FINITE-LATTICE ENERGETIC STABILITY
# ==========================================================

def theoretical_lambda_c():
    """
    Calculate the energetic instability threshold for the
    finite quadratic model under fixed total mass.

    The Hessian is

        H = 2g I + 2 kappa L + Lambda K.

    Because total mass is fixed, the uniform vector is excluded.
    The threshold is obtained from the first eigenvalue that
    becomes zero in the mass-conserving subspace.
    """

    I = np.eye(L)

    # Non-magnetic part of Hessian
    B = 2 * g * I + 2 * kappa * LAP

    # Construct an orthonormal basis perpendicular to uniform mode
    u = np.ones(L) / np.sqrt(L)

    Q, _ = np.linalg.qr(
        np.column_stack((u, np.eye(L)[:, 1:]))
    )

    # Numerical alternative: use SVD to obtain null space
    U, S, Vt = np.linalg.svd(u.reshape(1, -1))

    # Null space of the mass constraint
    Z = Vt[1:].T

    B_red = Z.T @ B @ Z
    K_red = Z.T @ KERNEL @ Z

    # Solve generalized eigenvalue problem
    # B_red + Lambda K_red = 0
    #
    # Transform using B^(-1/2)

    eig_B, vec_B = np.linalg.eigh(B_red)

    B_inv_sqrt = (
        vec_B
        @ np.diag(1.0 / np.sqrt(eig_B))
        @ vec_B.T
    )

    C = B_inv_sqrt @ K_red @ B_inv_sqrt

    eig_C = np.linalg.eigvalsh(C)

    negative_modes = eig_C[eig_C < 0]

    if len(negative_modes) == 0:
        return np.inf

    # B + Lambda K becomes unstable when
    #
    # 1 + Lambda * mu = 0
    #
    # for the most negative eigenvalue mu.

    Lambda_c = -1.0 / np.min(negative_modes)

    return Lambda_c


# ==========================================================
# RUN SIMULATIONS
# ==========================================================

roughness_values = []
spectral_values = []
profiles = []
energies = []

print("Running simulations...")

for Lambda in lambda_values:

    print(f"Lambda = {Lambda:.2f}")

    h = relax_system(Lambda)

    profiles.append(h.copy())

    roughness_values.append(normalized_roughness(h))

    spectral_values.append(
        spectral_columnarity(h)
    )

    energies.append(
        energy(h, Lambda)
    )


roughness_values = np.array(roughness_values)
spectral_values = np.array(spectral_values)
energies = np.array(energies)


# ==========================================================
# THEORETICAL INSTABILITY
# ==========================================================

Lambda_c_theory = theoretical_lambda_c()

print("\nFinite-lattice energetic instability:")
print(f"Lambda_c ≈ {Lambda_c_theory:.3f}")


# ==========================================================
# FIGURE 1: ROUGHNESS / COLUMNARITY
# ==========================================================

plt.figure(figsize=(7, 5))

plt.plot(
    lambda_values,
    roughness_values,
    marker='o',
    label='Normalized roughness'
)

plt.axvline(
    Lambda_c_theory,
    linestyle=':',
    label=r'$\Lambda_c$ (energetic stability)'
)

plt.xlabel(r'Dipolar strength $\Lambda$')
plt.ylabel(r'Normalized roughness $R/\langle h\rangle^2$')

plt.legend()
plt.tight_layout()
plt.show()


# ==========================================================
# FIGURE 2: SPECTRAL COLUMNARITY
# ==========================================================

plt.figure(figsize=(7, 5))

plt.plot(
    lambda_values,
    spectral_values,
    marker='o'
)

plt.axvline(
    Lambda_c_theory,
    linestyle=':',
    label=r'$\Lambda_c$ (energetic stability)'
)

plt.xlabel(r'Dipolar strength $\Lambda$')
plt.ylabel('High-wavenumber spectral fraction')

plt.legend()
plt.tight_layout()
plt.show()


# ==========================================================
# FIGURE 3: MOUND PROFILES
# ==========================================================

plt.figure(figsize=(12, 4))

indices = [
    0,
    len(lambda_values) // 2,
    -1
]

for idx in indices:

    plt.plot(
        profiles[idx],
        label=f'Λ = {lambda_values[idx]:.2f}'
    )

plt.xlabel('Horizontal position')
plt.ylabel('Height')

plt.legend()
plt.tight_layout()
plt.show()
