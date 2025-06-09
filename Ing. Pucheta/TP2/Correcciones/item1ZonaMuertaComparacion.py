import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.signal import square

# Parámetros físicos del motor
Ra = 2.27
La = 0.0047
Ki = 0.25
Km = 0.25
J = 0.00233
B = 0.00131

# Matrices del sistema original
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
C = np.array([[0, 0, 1]])  # salida: ángulo

# Sistema ampliado con integrador del error
Ae = np.block([
    [A, np.zeros((3, 1))],
    [-C, np.zeros((1, 1))]
])
Be = np.vstack([Bmat, [[0]]])

# LQR sobre el sistema ampliado
Q = np.diag([0.1, 0.1, 100, 1000])  # mayor peso al ángulo y al integrador
R = np.array([[5]])
P = solve_continuous_are(Ae, Be, Q, R)
K = np.linalg.inv(R) @ Be.T @ P
Kx = K[:, :3]
Ki_int = K[:, 3]

# Configuración temporal
h = 1e-4
tF = 5
t = np.arange(0, tF, h)
N = len(t)
ref = (np.pi / 2) * square(2 * np.pi * (1 / 10) * t)  # referencia alternante

def simular_con_zona_muerta(zm_threshold):
    X = np.zeros((3, N))  # ia, omega, theta
    e_int = np.zeros(N)
    u = np.zeros(N)

    for k in range(N - 1):
        ia, omega, theta = X[:, k]
        err = ref[k] - theta
        e_int[k+1] = e_int[k] + h * err

        # Control LQR con zona muerta
        ue = (-Kx @ X[:, k] - Ki_int * e_int[k+1]).item()
        if abs(ue) > zm_threshold:
            u_k = ue - zm_threshold * np.sign(ue)
        else:
            u_k = 0.0
        u[k] = u_k

        # Dinámica física del motor
        ia_dot = (-Ra * ia - Km * omega + u_k) / La
        omega_dot = (Ki * ia - B * omega) / J
        theta_dot = omega

        # Integración por Euler
        X[0, k+1] = ia + h * ia_dot
        X[1, k+1] = omega + h * omega_dot
        X[2, k+1] = theta + h * theta_dot

    return X, u

# Simulación para dos zonas muertas distintas
zm1 = 0.5
zm2 = 5
X1, u1 = simular_con_zona_muerta(zm1)
X2, u2 = simular_con_zona_muerta(zm2)

# Gráficos
plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(t, X1[2], label=f'θ(t), ZM={zm1}')
plt.plot(t, X2[2], '--', label=f'θ(t), ZM={zm2}')
plt.plot(t, ref, ':', label='Referencia')
plt.ylabel('Ángulo [rad]')
plt.title('Respuesta angular con distintas zonas muertas')
plt.legend()
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(t, u1, label=f'u(t), ZM={zm1}')
plt.plot(t, u2, '--', label=f'u(t), ZM={zm2}')
plt.ylabel('Tensión [V]')
plt.xlabel('Tiempo [s]')
plt.title('Acción de control')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
