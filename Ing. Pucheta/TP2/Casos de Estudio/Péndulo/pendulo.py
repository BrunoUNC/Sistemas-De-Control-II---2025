import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are

# ============================
# Parámetros físicos del sistema
# ============================
m = 0.1
F = 0.1
l = 1.6
g = 9.8
M = 1.5

# ============================
# Matrices del sistema continuo linealizado
# ============================
Ac = np.array([
    [0, 1, 0, 0],
    [0, -F/M, -m*g/M, 0],
    [0, 0, 0, 1],
    [0, -F/(l*M), -g*(m+M)/(l*M), 0]
])
Bc = np.array([[0], [1/M], [0], [1/(l*M)]])
Cc = np.array([[1, 0, 0, 0], [0, 0, 1, 0]])
Dc = np.array([[0]])

# ============================
# Discretización
# ============================
Ts = 1e-2
Ad = np.eye(4) + Ac * Ts
Bd = Bc * Ts

# ============================
# Sistema ampliado para referencia en delta (no en phi)
# ============================
Cref = Cc[0:1, :]
Aamp = np.block([
    [Ad, np.zeros((4, 1))],
    [-Cref @ Ad, np.eye(1)]
])
Bamp = np.vstack((Bd, -Cref @ Bd))

Q = np.diag([0.1, 1e-2, 1, 0.1, 0.0000093])
R = np.array([[1.9e-4]])
P = solve_discrete_are(Aamp, Bamp, Q, R)
K = np.linalg.inv(Bamp.T @ P @ Bamp + R) @ (Bamp.T @ P @ Aamp)
Kp = K[:, :4]
Ki = -K[:, 4]

# ============================
# Simulación
# ============================
T = 25
h = Ts / 20
t = np.arange(0, T, h)
ref = np.ones_like(t) * 10

x = np.array([0.0, 0.0, np.pi, 0.0])
x_hat = np.array([0.0, 0.0, np.pi, 0.0])
phi = [x[2]]
omega = [x[3]]
delta = [x[0]]
delta_dot = [x[1]]
v = [0]
u_ctrl = []

phi_pp = 0
bool_switch = False
mass_multiplied = False

for i in range(len(t)):
    y = Cc @ x
    y_obs = Cc @ (x_hat - np.array([0, 0, np.pi, 0]))
    v.append(v[-1] + (ref[i] - y[0]))
    u = -Kp @ (x - np.array([0, 0, np.pi, 0])) + Ki * v[-1]


    u_ctrl.append(u)

    # dinámica no lineal
    deltaP = float(x[1])
    phi_val = float(x[2])
    omega_val = float(x[3])

    for _ in range(int(Ts / h)):
        p_pp = (1 / (M + m)) * (u - m * l * phi_pp * np.cos(phi_val) + m * l * omega_val**2 * np.sin(phi_val) - F * deltaP)
        phi_pp = (1 / l) * (g * np.sin(phi_val) - p_pp * np.cos(phi_val))

        deltaP += h * p_pp
        delta_val = delta[-1] + h * deltaP
        omega_val += h * phi_pp
        phi_val += h * omega_val

        delta.append(delta_val)
        delta_dot.append(deltaP)
        omega.append(omega_val)
        phi.append(phi_val)

        x = np.array([delta_val, deltaP, phi_val, omega_val])
        if not mass_multiplied and delta_val >= 9.99:
            m *= 10
            mass_multiplied = True

    x_hat = Ad @ x_hat + Bd.flatten() * u + np.array([0, 0, 0, 0])  # sin observador aún

# ============================
# Gráficos
# ============================
plt.figure(figsize=(10, 8))
plt.subplot(3, 1, 1)
plt.plot(t, delta[:len(t)], label='Desplazamiento')
plt.plot(t, ref[:len(t)], '--', label='Referencia')
plt.legend()
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, phi[:len(t)], label='Ángulo φ')
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t, u_ctrl[:len(t)], label='u(t)')
plt.grid()

plt.tight_layout()
plt.show()
