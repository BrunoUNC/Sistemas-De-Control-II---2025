# Controlador LQR para un motor de corriente continua
# Evolución al modificar parámetro R del controlador LQR

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
from scipy.linalg import solve_continuous_are

# ============================ #
# Parámetros físicos del motor
# ============================ #
Ra = 2.27
La = 0.0047
Ki = 0.25
Km = 0.25
J = 0.00233
B = 0.00131

# ============================ #
# Modelo lineal en espacio de estados
# ============================ #
# x = [ia, omega, theta]
A = np.array([
    [-Ra / La,    -Km / La,    0],
    [Ki / J,     -B / J,       0],
    [0,           1,           0]
])
Bmat = np.array([
    [1 / La],
    [0],
    [0]
])
C = np.array([[0, 0, 1]])
D = np.array([[0]])

# ============================ #
# Definición de torque cuadrado periódico
# ============================ #
def torque_periodico(t):
    t_mod = t % 5
    return 0.12 if 2 <= t_mod < 3 else 0.0

# ============================ #
# Tiempo y paso de integración
# ============================ #
# Fijar h como el paso más restrictivo (R=10)
autovalores = np.linalg.eigvals(A)
polo_mas_rapido = np.min(np.abs(np.real(autovalores[autovalores != 0])))
Tmin = 1 / polo_mas_rapido
h = Tmin / 20  # 20 muestras por constante de tiempo
print(f"Paso de integración fijo: h = {h:.2e} s")

tF = 20
N = int(tF / h)
t = np.linspace(0, tF, N)
ref = (np.pi / 2) * square(2 * np.pi * (1 / 10) * t)

# ============================ #
# Simulaciones para varios R
# ============================ #
R_values = [2, 5, 10, 20, 30]
resultados_angulo = {}
resultados_control = {}

for Rval in R_values:
    Q = np.diag([1, 1, 40])
    R = np.array([[Rval]])
    P = solve_continuous_are(A, Bmat, Q, R)
    K = np.linalg.inv(R) @ (Bmat.T @ P)

    X = np.zeros((3, N))
    u = np.zeros(N)
    for k in range(N - 1):
        ia, omega, theta = X[:, k]
        torque = torque_periodico(t[k])
        x_ref = np.array([0, 0, ref[k]])
        u_k = (-K @ (X[:, k] - x_ref)).item()   
        u_k = np.clip(u_k, -5, 5)
        u[k] = u_k

        # Dinámica del sistema
        ia_dot = (-Ra * ia - Km * omega + u_k) / La
        omega_dot = (Ki * ia - B * omega - torque) / J
        theta_dot = omega

        X[0, k+1] = ia + h * ia_dot
        X[1, k+1] = omega + h * omega_dot
        X[2, k+1] = theta + h * theta_dot

    resultados_angulo[Rval] = X[2, :]
    resultados_control[Rval] = u

# ============================ #
# Gráficos
# ============================ #
plt.figure(figsize=(10, 6))
for Rval in R_values:
    plt.plot(t, resultados_angulo[Rval], label=f'R={Rval}')
plt.plot(t, ref, '--k', label='Referencia')
plt.title('Ángulo del motor para distintos R en LQR')
plt.xlabel('Tiempo [s]')
plt.ylabel('Ángulo θ [rad]')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
for Rval in R_values:
    plt.plot(t, resultados_control[Rval], label=f'R={Rval}')
plt.title('Acción de control u(t) para distintos R')
plt.xlabel('Tiempo [s]')
plt.ylabel('Tensión [V]')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
