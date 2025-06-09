# Controlador LQR para ángulo del motor DC con referecia distinta de cero y perturbación periódica
# Se ha agregado un integrador de error para generar el sistema ampliado
# Se presenta ajuste de zona muerta a la acción de control


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.signal import square

plt.rcParams.update({'font.size': 18})

# ============================
# Parámetros físicos del motor
# ============================
Ra = 2.27
La = 0.0047
Ki = 0.25
Km = 0.25
J = 0.00233
B = 0.00131

# ============================
# Modelo en espacio de estados
# ============================
A = np.array([
    [-Ra / La,    -Km / La,    0],
    [Ki / J,     -B / J,       0],
    [0,           1,           0]
])

Bmat = np.array([
    [1 / La], 
    [0], [0]
])

C = np.array([[0, 0, 1]])

# Matriz de salida del ángulo (para integrador de error)
Cy = C

#Generacion de sistema ampliado
# Matrices ampliadas para agregar integrador del error
Ae = np.block([
    [A, np.zeros((3, 1))],
    [-Cy, np.zeros((1, 1))]
])
Be = np.vstack([Bmat, np.zeros((1, 1))])


# ============================
# Controlador LQR
# ============================
# Matrices de penalización para LQR
Q = np.diag([0.1, 0.1, 100, 10000])  # penaliza error en ángulo + integrador (sist ampliado)
R = np.array([[5]])          # Contol sobre la accion de control
                             # A mayor R, menor accion de control
# Resolución de Riccati
P = solve_continuous_are(Ae, Be, Q, R)
K = np.linalg.inv(R) @ Be.T @ P

# Separación del controlador
Kx = K[:, :3]
Ki_int = K[:, 3]


################################################
# Zona muerta como función continua por partes #
################################################
def aplicar_zona_muerta(u, umbral):
    if abs(u) > umbral:
        return u - umbral * np.sign(u)
    else:
        return 0.0


# ============================
# Simulación con Euler
# ============================
autovalores = np.linalg.eigvals(A)
polo_mas_rapido = np.min(np.abs(np.real(autovalores[autovalores != 0])))
Tmin = 1 / polo_mas_rapido
h = Tmin / 20
print(f"Paso de integración h = {h:.2e} s")

# Simulación
tF = 20
N = int(tF / h)
t = np.linspace(0, tF, N)
ref = (np.pi / 2) * square(2 * np.pi * (1 / 10) * t)

X = np.zeros((3, N))        # Estados: ia, omega, theta
e_int = np.zeros(N)         # Estado adicional: integrador del error
u = np.zeros(N)             # Acción de control
corriente = np.zeros(N)     # ia

for k in range(N - 1):
    ia, omega, theta = X[:, k]
 
    y = theta
    error = ref[k] - y
    e_int[k + 1] = e_int[k] + h * error

    # Control ampliado
    u_k = (-Kx @ X[:, k] - Ki_int * e_int[k + 1]).item()
    u_k = aplicar_zona_muerta(u_k, 0.5)      # Aplicación de zona muerta corregida
    u[k] = u_k

    corriente[k] = ia

    # Dinámica del sistema (modelo físico)
    ia_dot = (-Ra * ia - Km * omega + u_k) / La
    omega_dot = (Ki * ia - B * omega) / J
    theta_dot = omega

    # Euler
    X[0, k + 1] = ia + h * ia_dot
    X[1, k + 1] = omega + h * omega_dot
    X[2, k + 1] = theta + h * theta_dot

# Plots
plt.figure(figsize=(10, 10))

plt.subplot(2, 1, 1)
plt.plot(t, X[2, :], label='θ(t)')
plt.plot(t, ref, '--', label='Referencia')
plt.ylabel('Ángulo [rad]')
plt.title('Ángulo')
plt.legend()
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(t, u, label='u(t)')
plt.ylabel('Tensión [V]')
plt.title('Acción de control (Va)')
plt.grid()


plt.tight_layout()
plt.show()
