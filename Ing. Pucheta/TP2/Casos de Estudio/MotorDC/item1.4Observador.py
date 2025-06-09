#Comparación del sistema con y sin observador de Luenberger

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
from scipy.linalg import solve_continuous_are

plt.rcParams.update({'font.size': 18})

# Parámetros físicos del motor
Ra = 2.27
La = 0.0047
Ki = 0.25
Km = 0.25
J = 0.00233
B = 0.00131

# Matrices del sistema
A = np.array([
    [-Ra / La, -Km / La, 0],
    [Ki / J, -B / J, 0],
    [0, 1, 0]
])
Bmat = np.array([
    [1 / La],
    [0],
    [0]
])
C = np.array([[0, 0, 1]])

# Sistema ampliado
Aamp = np.block([
    [A, np.zeros((3, 1))],
    [-C, np.zeros((1, 1))]
])
Bamp = np.vstack([Bmat, [[0]]])
Camp = np.hstack([C, [[0]]])

# Controlador LQR sobre sistema ampliado
Q = np.diag([0.1, 0.1, 100, 10000])  # Penalización en ángulo + integrador
R = np.array([[5]])
P = solve_continuous_are(Aamp, Bamp, Q, R)
K = np.linalg.inv(R) @ (Bamp.T @ P)
Kx = K[:, :3]
Ki_int = K[:, 3]

# Tiempo de simulación y paso de integración
h = 1e-04
tF = 20
t = np.arange(0, tF, h)
N = len(t)

# Señal de referencia (alternancia cada 2s)
ref = (np.pi / 2) * square(2 * np.pi * (1 / 10) * t)

# === SIMULACIÓN SIN OBSERVADOR === #
X = np.zeros((3, N))
e_int = np.zeros(N)
u = np.zeros(N)
ia_real = np.zeros(N)

for k in range(N - 1):
    ia, omega, theta = X[:, k]
    err = ref[k] - theta
    e_int[k + 1] = e_int[k] + h * err
    u_k = (-Kx @ X[:, k] - Ki_int * e_int[k + 1]).item()
    u[k] = u_k

    ia_dot = (-Ra * ia - Km * omega + u_k) / La
    omega_dot = (Ki * ia - B * omega) / J
    theta_dot = omega

    X[0, k + 1] = ia + h * ia_dot
    X[1, k + 1] = omega + h * omega_dot
    X[2, k + 1] = theta + h * theta_dot
    ia_real[k] = ia

# === SIMULACIÓN CON OBSERVADOR === #
# Observador de Luenberger
Ao = A.T
Bo = C.T
Qo = np.diag([1, 0.1, 0.1])
Ro = np.array([[5]])
Po = solve_continuous_are(Ao, Bo, Qo, Ro)
Ko = np.linalg.inv(Ro) @ (Bo.T @ Po)  # ganancia del observador

X_obs = np.zeros((3, N))
X_hat = np.zeros((3, N))
e_int_obs = np.zeros(N)
u_obs = np.zeros(N)
ia_hat = np.zeros(N)

X_hat[:, 0] = np.array([0.2, -0.1, 0.1])  # Condición inicial estimada

for k in range(N - 1):
    theta_hat = X_hat[2, k]
    err2 = ref[k] - theta_hat
    e_int_obs[k + 1] = e_int_obs[k] + h * err2
    u_k_obs = (-Kx @ X_hat[:, k] - Ki_int * e_int_obs[k + 1]).item()
    u_obs[k] = u_k_obs

    # Sistema real
    ia, omega, theta = X_obs[:, k]
    ia_dot = (-Ra * ia - Km * omega + u_k_obs) / La
    omega_dot = (Ki * ia - B * omega) / J
    theta_dot = omega
    X_obs[0, k + 1] = ia + h * ia_dot
    X_obs[1, k + 1] = omega + h * omega_dot
    X_obs[2, k + 1] = theta + h * theta_dot

    # Observador
    y = C @ X_obs[:, k]
    y_hat = C @ X_hat[:, k]
    x_hat_dot = A @ X_hat[:, k] + Bmat.flatten() * u_k_obs + Ko.flatten() * float(y - y_hat)
    X_hat[:, k + 1] = X_hat[:, k] + h * x_hat_dot
    ia_hat[k] = X_hat[0, k]

# === GRÁFICOS === #
plt.figure(figsize=(12, 10))

plt.subplot(4, 1, 1)
plt.plot(t, X[2], label='Ángulo sin observador')
plt.plot(t, X_obs[2], '--', label='Ángulo con observador')
plt.plot(t, ref, ':', label='Referencia')
plt.ylabel('Ángulo [rad]', fontsize=12)
plt.title('Ángulo')
plt.legend(fontsize=12)
plt.grid()

plt.subplot(4, 1, 2)
plt.plot(t, ia_real, label='Corriente real')
plt.plot(t, ia_hat, '--', label='Corriente estimada')
plt.ylabel('Corriente [A]', fontsize=12)
plt.title('Corriente de armadura')
plt.legend(fontsize=12)
plt.grid()

plt.subplot(4, 1, 3)
plt.plot(t, u, label='u sin observador')
plt.plot(t, u_obs, '--', label='u con observador')
plt.ylabel('Acción de control [V]', fontsize=12)
plt.title('Acción de control (Va)')
plt.legend(fontsize=12)
plt.grid()

plt.subplot(4, 1, 4)
plt.plot(t, ia_real - ia_hat, label='Error de observación de $i_a$')
plt.ylabel('Error [A]', fontsize=12)
plt.xlabel('Tiempo [s]', fontsize=12)
plt.legend(fontsize=12)
plt.grid()

plt.tight_layout()
plt.show()


plt.figure(figsize=(12, 10))

plt.subplot(2, 1, 1)
plt.plot(t, ia_real, label='Corriente real')
plt.plot(t, ia_hat, '--', label='Corriente estimada')
plt.ylabel('Corriente [A]', fontsize=12)
plt.title('Corriente de armadura')
plt.legend(fontsize=12)
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(t, u, label='u sin observador')
plt.plot(t, u_obs, '--', label='u con observador')
plt.ylabel('Acción de control [V]', fontsize=12)
plt.title('Acción de control (Va)')
plt.legend(fontsize=12)
plt.grid()

plt.tight_layout()
plt.show()
