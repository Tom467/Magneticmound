# 1D Dipolar Instability Simulation

Monte Carlo simulation of a one-dimensional granular height field with competing interactions:
- Gravity (local confinement)
- Surface tension (nearest-neighbor smoothing)
- Long-range dipolar repulsion (∝ 1/r³)

The system is relaxed via energy-lowering stochastic moves. By sweeping the dipolar strength (Λ), the code detects the onset of a pattern-forming instability and estimates the critical coupling numerically. The result is compared with a theoretical prediction derived from the discrete dipolar kernel.

---

## Physical Model

The total energy consists of:

- **Gravity term**: penalizes large heights  
- **Surface smoothing term**: penalizes local gradients  
- **Dipolar repulsion**: long-range interaction decaying as 1/r³  

An order parameter is computed from the dominant Fourier mode of the height profile.

---

## What the Code Does

1. Initializes a localized mound of grains
2. Relaxes the system using Monte Carlo energy minimization
3. Sweeps dipolar strength Λ
4. Computes:
   - Order parameter via FFT
   - Numerical estimate of critical Λ
   - Theoretical prediction from small-k expansion
5. Plots:
   - Order parameter vs Λ
   - Example height profiles

---

## How to Run

```bash
python simulation.py
