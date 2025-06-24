#  Heat Equation Solver Using Neural Approximation and Fourier Series

This repository contains a MATLAB implementation for solving the 1D heat equation with homogeneous boundary conditions. It includes both an exact solution using a Fourier sine series and a neural-inspired approximation method based on parameter optimization.

---

## 📁 Files

- `hh1_H2.m` — Main function that computes:
  - The exact solution via Fourier expansion
  - The approximate solution via a neural-like iterative method
  - Mesh plots comparing both
- `Cn_nh.m` — Computes the Fourier sine coefficients for a function defined in the interval [0, 1]
- Additional helper files may be added for different experiments

---

## 📐 Problem Description

The heat equation:


is solved in the interval [0, 1] with homogeneous boundary conditions \( u(0,t) = u(1,t) = 0 \), and initial condition \( u(x,0) = f(x) \).  

Two approaches are presented:
- 📈 **Exact solution**: Using truncated Fourier sine series  
- 🧠 **Neural approximation**: Parameters are adjusted via gradient descent to minimize the difference with the exact solution

---

##  How to Use

1. Open MATLAB and navigate to the folder containing the `.m` files
2. Run the function in the Command Window:

```matlab
[exa, Y2] = hh1_H2(0, 1, 10, 1e-6);
