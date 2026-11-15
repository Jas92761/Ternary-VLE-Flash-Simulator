# Multi-Component Vapor-Liquid Equilibrium (VLE) & Flash Simulator

A rigorous scientific computing and numerical optimization tool developed in Python to solve multi-component thermodynamic phase equilibria using the Non-Random Two-Liquid (NRTL) activity coefficient model.

---

## 🧪 Problem Statement & System Overview
This project solves **Problem 12.20** for a non-ideal ternary system consisting of:
1. **Acetone (1)**
2. **Methanol (2)**
3. **Water (3)**

Evaluated at a fixed temperature of **$T = 338.15\text{ K}$ ($65^\circ\text{C}$)** under three core phase-behavior regimes:
- **Part (a) - BUBL P:** Calculates total bubble pressure given liquid mole fractions ($x_1=0.3, x_2=0.4, x_3=0.3$).
- **Part (b) - DEW P:** Solves for dew pressure and liquid composition given vapor mole fractions ($y_1=0.3, y_2=0.4, y_3=0.3$).
- **Part (c) - P,T-Flash:** Evaluates phase split fractions and equilibrium compositions at $P = (P_{\text{bubble}} + P_{\text{dew}})/2$ using Rachford-Rice mass balance constraints.

---

## 📐 Thermodynamic Equations & Mathematical Formulations

### 1. Saturation Pressure (Antoine Equation)
Pure component vapor pressures ($P_i^{\text{sat}}$) are modeled as:
$$P_i^{\text{sat}} = \exp\left(A_i - \frac{B_i}{T + C_i}\right)$$

### 2. NRTL Activity Coefficient Model ($\gamma_i$)
To account for liquid-phase non-ideality, excess Gibbs energy is evaluated using the multicomponent NRTL model:
- **Binary Interaction Energy Parameter:** 
  $$\tau_{ij} = \frac{b_{ij}}{RT}$$
- **Non-Randomness Factor:** 
  $$G_{ij} = \exp(-\alpha_{ij}\tau_{ij})$$
- **Activity Coefficient Expression ($\ln \gamma_i$):**
  $$\ln \gamma_i = \sum_j \frac{x_j \tau_{ji} G_{ji}}{\sum_k x_k G_{ki}} + \sum_j \frac{x_j G_{ij}}{\sum_k x_k G_{ik}} \left( \tau_{ij} - \frac{\sum_m x_m \tau_{mj} G_{mj}}{\sum_k x_k G_{kj}} \right)$$

### 3. Phase Equilibrium & Flash Constraints
- **Modified Raoult's Law (Bubble/Dew Pressure):**
  $$P = \sum_{i} x_i \gamma_i P_i^{\text{sat}} \quad \text{and} \quad \frac{1}{P} = \sum_{i} \frac{y_i}{\gamma_i P_i^{\text{sat}}}$$
- **Rachford-Rice Mass Balance (Flash Routine):**
  Solved using K-values ($K_i = \frac{\gamma_i P_i^{\text{sat}}}{P}$) and bounded vapor fractions ($\beta$).

---

## 💻 Numerical Implementation & Tech Stack
- **Language:** Python
- **Libraries:** `NumPy` (matrix-based array operations), `SciPy` (`scipy.optimize.root` for non-linear root finding, `scipy.optimize.least_squares` for bounded flash optimization).
