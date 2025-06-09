import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
from scipy.linalg import solve_continuous_are

plt.rcParams.update({'font.size': 16})

# Tiempo de simulación y paso
h = 1e-4
tF = 20
t = np.arange(0, tF, h)
N = len(t)

# Referencia cuadrada alternante
ref = (np.pi / 2) * square(2 * np.pi * (1 / 10) * t)

# --- Definición de simulación con parámetros ajustables --- #
def simulate_motor(Ra, La, Ki, Km, J, B, label):
    A = np.array([
        [-Ra / La, -Km / La, 0],
        [Ki / J, -B / J, 0],
        [0, 1, 0]
    ])
    Bmat = np.array([[1 / La], [0], [0]])
    C = np.array([[0, 0, 1]])

    Aamp = np.block([[A, np.zeros((3, 1))], [-C, np.zeros((1, 1))]])
    Bamp = np.vstack([Bmat, [[0]]])
    Q = np.diag([0.1, 0.1, 100, 10000])
    R = np.array([[5]])
    P = solve_continuous_are(Aamp, Bamp, Q, R)
    K = np.linalg.inv(R) @ (Bamp.T @ P)
    Kx = K[:, :3]
    Ki_int = K[:, 3]

    Ao = A.T
    Bo = C.T
    Qo = np.diag([1, 0.1, 0.1])
    Ro = np.array([[5]])
    Po = solve_continuous_are(Ao, Bo, Qo, Ro)
    Ko = np.linalg.inv(Ro) @ (Bo.T @ Po)

    X = np.zeros((3, N))
    X_hat = np.zeros((3, N))
    e_int = np.zeros(N)
    u = np.zeros(N)

    X_hat[:, 0] = np.array([0.2, -0.1, 0.1])

    for k in range(N - 1):
        err = ref[k] - X_hat[2, k]
        e_int[k + 1] = e_int[k] + h * err
        u_k = (-Kx @ X_hat[:, k] - Ki_int * e_int[k + 1]).item()
        u[k] = u_k

        ia, omega, theta = X[:, k]
        ia_dot = (-Ra * ia - Km * omega + u_k) / La
        omega_dot = (Ki * ia - B * omega) / J
        theta_dot = omega

        X[0, k + 1] = ia + h * ia_dot
        X[1, k + 1] = omega + h * omega_dot
        X[2, k + 1] = theta + h * theta_dot

        y = C @ X[:, k]
        y_hat = C @ X_hat[:, k]
        x_hat_dot = A @ X_hat[:, k] + Bmat.flatten() * u_k + Ko.flatten() * float(y - y_hat)
        X_hat[:, k + 1] = X_hat[:, k] + h * x_hat_dot

    return X[2, :], label

# Parámetros físicos de dos motores distintos
params_motor_1 = (2.27, 0.0047, 0.25, 0.25, 0.00233, 0.00131, "Motor 1")
params_motor_2 = (2.25, 0.005, 0.259, 0.25, 0.00284, 0.0014, "Motor 2")

theta1, label1 = simulate_motor(*params_motor_1)
theta2, label2 = simulate_motor(*params_motor_2)

# Gráfico comparativo
plt.figure(figsize=(10, 6))
plt.plot(t, theta1, label=label1)
plt.plot(t, theta2, '--', label=label2)
plt.plot(t, ref, 'k:', label='Referencia')
plt.title("Respuesta angular de los dos motores (Real vs identificado)")
plt.xlabel("Tiempo [s]")
plt.ylabel("Ángulo [rad]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

