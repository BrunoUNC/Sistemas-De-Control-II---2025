# Controlador LQR para ángulo del motor DC con referecia distinta de cero y perturbación periódica
# Se ha agregado un integrador de error para generar el sistema ampliado
# Se agrega una zona muerta a la acción de control


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.signal import square

plt.rcParams.update({'font.size': 18})

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
Q = np.diag([0.1, 0.1, 500])     # penaliza el error de ángulo
R = np.array([[10]])        # penaliza el esfuerzo de control
                            # a mayot R, menor esfuerzo de control
                              
P = solve_continuous_are(A, Bmat, Q, R) # Resolución de la ecuación de Riccati
K = np.linalg.inv(R) @ (Bmat.T @ P)     # Ganancia de control LQR 

# ============================
# Perturbación: torque periódico
# ============================
#def torque_periodico(t):
#    t_mod = t % 5
#    return 0.12 if 2 <= t_mod < 3 else 0.0

# ============================
# Simulación con Euler
# ============================
autovalores = np.linalg.eigvals(A)
polo_mas_rapido = np.min(np.abs(np.real(autovalores[autovalores != 0])))
Tmin = 1 / polo_mas_rapido
h = 1e-04
print(f"Paso de integración h = {h:.2e} s")

# Simulación
tF = 10
N = int(tF / h)
t = np.linspace(0, tF, N)

# Referencia alternante cada 5s (π/2, -π/2)
ref = (np.pi / 2) * square(2 * np.pi * (1 / 1000) * t)  # periodo 10s

X = np.zeros((3, N))        # Estados: ia, omega, theta
e_int = np.zeros(N)         # Estado adicional: integrador del error
u = np.zeros(N)             # Acción de control
torque_ext = np.zeros(N)    # Perturbación
corriente = np.zeros(N)     # ia


# Inicialización
X = np.zeros((3, N))  # x = [ia, omega, theta]

torque_ext = np.zeros(N)


################################################
# Zona muerta como función continua por partes #
################################################
def aplicar_zona_muerta(u, umbral):
    if abs(u) > umbral:
        return u - umbral * np.sign(u)
    else:
        return 0.0



for k in range(N - 1):
    ia, omega, theta = X[:, k]
    #torque_ext[k] = torque_periodico(t[k])
    x_ref = np.array([0, 0, ref[k]])


  
  # Control
    u_k = (-K @ (X[:, k] - x_ref)).item()
    #u_k = np.clip(u_k, -5, 5)               # Limitación de la acción de control
    u_k = aplicar_zona_muerta(u_k, 0.5)      # Aplicación de zona muerta corregida
    u[k] = u_k




    # Derivadas físicas (estilo modmotor)
    ia_dot = (-Ra * ia - Km * omega + u_k) / La
    omega_dot = (Ki * ia - B * omega ) / J
    theta_dot = omega

    # Integración por Euler
    X[0, k+1] = ia + h * ia_dot
    X[1, k+1] = omega + h * omega_dot
    X[2, k+1] = theta + h * theta_dot

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
plt.xlabel('Tiempo [s]')
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

