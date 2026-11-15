import numpy as np
from scipy.optimize import root, least_squares

print("=" * 60)
print(" TERNARY VLE & FLASH SIMULATOR (NRTL MODEL)")
print(" System: Acetone(1) / Methanol(2) / Water(3)")
print("=" * 60)

# Given data
T = 338.15                # Temperature in K (65°C)
R = 8.314                 # Universal gas constant (J/mol·K)

# Antoine constants for Acetone(1), Methanol(2), Water(3)
A = np.array([14.3145, 16.5785, 16.3872])
B = np.array([2756.22, 3638.27, 3885.70])
C = np.array([-45.090, -33.650, -42.980])

# NRTL binary parameters (b in J/mol, alpha dimensionless)
b12, b21 = 773.3, 932.1
b13, b31 = 2642.1, 5013.3
b23, b32 = -1062.9, 3538.7

a12, a13, a23 = 0.3084, 0.5343, 0.2994

# Saturation pressure function (Antoine equation)
def Psat(T_val):
    return np.exp(A - B / (T_val + C))

# NRTL activity coefficient calculator
def gamma_from_x(T_val, x):
    X = np.array([0.0, x[0], x[1], x[2]])
    
    t = np.zeros((4, 4))
    t[1][2] = b12 / (R * T_val)
    t[2][1] = b21 / (R * T_val)
    t[1][3] = b13 / (R * T_val)
    t[3][1] = b31 / (R * T_val)
    t[2][3] = b23 / (R * T_val)
    t[3][2] = b32 / (R * T_val)
    t[1][1] = t[2][2] = t[3][3] = 0.0

    G = np.zeros((4, 4))
    G[1][2] = np.exp(-a12 * t[1][2])
    G[2][1] = np.exp(-a12 * t[2][1])
    G[1][3] = np.exp(-a13 * t[1][3])
    G[3][1] = np.exp(-a13 * t[3][1])
    G[2][3] = np.exp(-a23 * t[2][3])
    G[3][2] = np.exp(-a23 * t[3][2])
    G[1][1] = G[2][2] = G[3][3] = 1.0

    ln_gamma = np.zeros(4)
    for i in range(1, 4):
        Di = sum(X[k] * G[i][k] for k in range(1, 4))
        Num_i = sum(X[k] * t[i][k] * G[i][k] for k in range(1, 4))
        
        sum1 = 0.0
        sum2 = 0.0
        for j in range(1, 4):
            Dj = sum(X[k] * G[j][k] for k in range(1, 4))
            sum1 += X[j] * t[j][i] * G[j][i] / Dj
            sum2 += (X[j] * G[i][j] / Di) * (t[i][j] - Num_i / Di)
            
        ln_gamma[i] = sum1 + sum2

    return np.exp(ln_gamma[1:4])

Ps = Psat(T)

# ==========================================
# PART (a): BUBL P Calculation
# ==========================================
x_given = np.array([0.3, 0.4, 0.3])  # x1=0.3, x2=0.4, x3=0.3
gamma_bubble = gamma_from_x(T, x_given)
P_bubble = np.sum(x_given * gamma_bubble * Ps)

print("\n--- (a) BUBL P ---")
print(f"Liquid Composition (x): {x_given}")
print(f"Bubble Pressure (P_bubble): {P_bubble:.4f} kPa")


# ==========================================
# PART (b): DEW P Calculation
# ==========================================
y_given = np.array([0.3, 0.4, 0.3])  # y1=0.3, y2=0.4, y3=0.3

def dew_equations(v):
    P = v[0]
    x = v[1:4]
    gamma = gamma_from_x(T, x)
    return np.append(x - y_given * P / (gamma * Ps), np.sum(x) - 1.0)

dew_guess = [70.0, 0.1, 0.2, 0.7]
dew_solution = root(dew_equations, dew_guess)

if not dew_solution.success:
    raise ValueError("Dew pressure calculation did not converge.")

P_dew = dew_solution.x[0]
x_dew = dew_solution.x[1:4]

print("\n--- (b) DEW P ---")
print(f"Vapor Composition (y): {y_given}")
print(f"Dew Pressure (P_dew): {P_dew:.4f} kPa")
print(f"Equilibrium Liquid Comp (x): {x_dew}")


# ==========================================
# PART (c): P,T-Flash Calculation
# ==========================================
P_flash = (P_bubble + P_dew) / 2.0
z_given = np.array([0.3, 0.4, 0.3])  # z1=0.3, z2=0.4, z3=0.3

def flash_equations(v):
    x1, x2, beta = v
    x = np.array([x1, x2, 1.0 - x1 - x2])
    gamma = gamma_from_x(T, x)
    K = gamma * Ps / P_flash
    y = K * x

    eq1 = (1.0 - beta) * x[0] + beta * y[0] - z_given[0]
    eq2 = (1.0 - beta) * x[1] + beta * y[1] - z_given[1]
    eq3 = np.sum(y) - 1.0
    return [eq1, eq2, eq3]

flash_guess = [0.2, 0.3, 0.6]
flash_solution = least_squares(
    flash_equations,
    flash_guess,
    bounds=([1e-8, 1e-8, 1e-8], [0.999, 0.999, 0.999])
)

if not flash_solution.success:
    raise ValueError("Flash calculation did not converge.")

x1_f, x2_f, beta = flash_solution.x
x_flash = np.array([x1_f, x2_f, 1.0 - x1_f - x2_f])
gamma_flash = gamma_from_x(T, x_flash)
K_flash = gamma_flash * Ps / P_flash
y_flash = K_flash * x_flash

print("\n--- (c) P,T-FLASH ---")
print(f"Feed Composition (z): {z_given}")
print(f"Flash Pressure (P): {P_flash:.4f} kPa  [(P_bubble + P_dew) / 2]")
print(f"Vapor Fraction (beta): {beta:.4f}")
print(f"Liquid Fraction (1 - beta): {1.0 - beta:.4f}")
print(f"Resulting Liquid Composition (x): {x_flash}")
print(f"Resulting Vapor Composition (y): {y_flash}")
print("=" * 60)