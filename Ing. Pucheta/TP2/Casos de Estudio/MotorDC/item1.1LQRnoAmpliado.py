#Control de un motor de corriente continua con LQR
# Referencia alternante cada 5s (π/2, -π/2)


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
C = np.array([[0, 0, 1]])  # salida: ángulo
D = np.array([[0]])       

# ============================ #
# Diseño del controlador LQR
# ============================ #
Q = np.diag([1, 1, 40])     # penaliza el error de ángulo
R = np.array([[10]])        # penaliza el esfuerzo de control
                            # a mayot R, menor esfuerzo de control
                              
P = solve_continuous_are(A, Bmat, Q, R) # Resolución de la ecuación de Riccati
K = np.linalg.inv(R) @ (Bmat.T @ P)     # Ganancia de control LQR 

# ============================ #
# Definición de torque cuadrado periódico
# ============================ #
def torque_periodico(t):
    t_mod = t % 5
    return 0.12 if 2 <= t_mod < 3 else 0.0

# ============================ #
# Simulación Euler 
# ============================ #

# Cálculo del paso de integración h en función de los polos
autovalores = np.linalg.eigvals(A)
polo_mas_rapido = np.min(np.abs(np.real(autovalores[autovalores != 0])))
Tmin = 1 / polo_mas_rapido
h = Tmin / 20  # Resolución mínima de 20 muestras por constante de tiempo
print(f"Paso de integración h = {h:.2e} s")
tF = 20
N = int(tF / h)
t = np.linspace(0, tF, N)

# Referencia alternante cada 5s (π/2, -π/2)
ref = (np.pi / 2) * square(2 * np.pi * (1 / 10) * t)  # periodo 10s

# Inicialización
X = np.zeros((3, N))  # x = [ia, omega, theta]
u = np.zeros(N)
torque_ext = np.zeros(N)

for k in range(N - 1):
    ia, omega, theta = X[:, k]
    torque_ext[k] = torque_periodico(t[k])
    x_ref = np.array([0, 0, ref[k]])

    # Control (escalar puro)
    u_k = (-K @ (X[:, k] - x_ref)).item()
    u_k = np.clip(u_k, -5, 5) # Limitación de la acción de control
    u[k] = u_k

    # Derivadas físicas (estilo modmotor)
    ia_dot = (-Ra * ia - Km * omega + u_k) / La
    omega_dot = (Ki * ia - B * omega - torque_ext[k]) / J
    theta_dot = omega

    # Integración por Euler
    X[0, k+1] = ia + h * ia_dot
    X[1, k+1] = omega + h * omega_dot
    X[2, k+1] = theta + h * theta_dot

# ============================ #
# Gráficos
# ============================ #
plt.figure(figsize=(10, 10))

plt.subplot(4, 1, 1)
plt.plot(t, X[2, :], label='Ángulo θ(t)')
plt.plot(t, ref, '--', label='Referencia')
plt.ylabel('[rad]')
plt.title('Respuesta del sistema')
plt.legend()
plt.grid()

plt.subplot(4, 1, 2)
plt.plot(t, u, label='u(t) [V]')
plt.ylabel('[V]')
plt.title('Acción de control (Va)')
plt.grid()

plt.subplot(4, 1, 3)
plt.plot(t, X[0, :], label='Corriente i_a(t)')
plt.ylabel('[A]')
plt.xlabel('Tiempo [s]')
plt.title('Corriente de armadura')
plt.grid()

plt.subplot(4, 1, 4)
plt.plot(t, torque_ext, label='Torque TL')
plt.ylabel('[Nm]')
plt.title('Perturbación externa')
plt.grid()

plt.tight_layout()
plt.show()


# === Plano de fase: velocidad angular vs ángulo ===
plt.figure(figsize=(6, 5))
plt.plot(X[2, :], X[1, :], linewidth=1.5)
plt.xlabel('Ángulo [rad]')
plt.ylabel('Velocidad angular [rad/s]')
plt.title('Plano de fase')
plt.grid(True)
plt.tight_layout()
plt.show()