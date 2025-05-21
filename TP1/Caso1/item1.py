import numpy as np
from scipy import signal
import matplotlib.pyplot as plt



'''
Se simula un circuito RLC en serie utilizando variables de estado. 

1-Importación de bibliotecas:
    numpy para operaciones matemáticas.
    scipy.signal para modelado y simulación de sistemas dinámicos.
    matplotlib.pyplot para graficar resultados.

2-Definición de parámetros del circuito:
    R: Resistencia en ohmios.
    L: Inductancia en henrios.
    C: Capacitancia en faradios.

3-Modelado en variables de estado:
    Define las variables de estado:
    x1 = i: Corriente en la malla.
    x2 = Vc: Tensión en el capacitor.
    Definición de matrices del sistema:
        Matriz A: Representa la dinámica del sistema.
        Matriz B: Representa la influencia de la entrada en el sistema.
        Matriz C: Selecciona la salida
        Matriz D: Matriz de transmisión directa (en este caso es cero).

5- Generación de la señal de entrada:
    Se genera una señal cuadrada de -12V con un ciclo de trabajo del 50% y una frecuencia de 500Hz.

6- Simulación del sistema:
    Se utiliza el método de Euler para simular el sistema en el tiempo definido.
    Se inicializan las variables y se almacenan los estados y salidas en cada paso de tiempo.


Finalmente, se realiza la visualización de sus respuestas de los sistemas.
'''


# Parámetros del circuito
R = 220        # Ohm
L = 0.5        # Henry
c = 2.2e-6     # Farad



#######################################################
###         Modelado en variables de estado         ###
#######################################################

# Varaibles de estado:
# x1 = i = Corriente de la malla
# x2 = Vc = Tensión en el capacitor


# Matrices del sistema
A = np.array([[-R/L, -1/L],
              [1/c, 0]])

B = np.array([[1/L],
                [0]])

D= np.array([[0]])       # Matriz D = 0, ya que no hay transmisión directa en el sistema


# Matriz C para la salida (Selector de salida) 

C1 = np.array([[R, 0]])  # Matriz C para la salida =  Tensión sobre el resistor
system1 = signal.StateSpace(A, B, C1, D)

C2 = np.array([[1, 0]])  # Selector de salida = Corriente
system2 = signal.StateSpace(A, B, C2, D)

C3 = np.array([[0, 1]])  # Selector de salida = Tensión sobre el capacitor
system3 = signal.StateSpace(A, B, C3, D)



#######################################################
####            Simulación del sistema              ###
#######################################################

# Parámetros del tiempo de simulación

# Duración total de la señal en s
total_duration = 0.100  # 100ms

dt = 10e-6  # Tiempo de muestreo (10 µs)
num_steps = int(total_duration / dt)  # Número de pasos de simulación
time = np.linspace(0, total_duration, num_steps)  # Recalcular el tiempo para que coincida con num_steps


#######################################################
###                 Señal de entrada                ###
#######################################################

# Frecuencia y Duty Cycle de la señal en Hz 
frequency_hz = 50
duty = 0.5
u_square = -12 * signal.square(2 * np.pi * frequency_hz * time, duty)
u = np.copy(u_square)
u[time < 0.01] = 0  # Los primeros 10ms (0.01s) se establecen en 0

# Inicialización de variables
x = np.array([0, 0])  # Estado inicial [i(0), Vc(0)]
x_history = np.zeros((num_steps, 2))  # Historial de estados
y1 = np.zeros(num_steps)  # Salida del sistema 1 (VR)
y2 = np.zeros(num_steps)  # Salida del sistema 2 (i)
y3 = np.zeros(num_steps)  # Salida del sistema 3 (VC)



#######################################################
###        Simulacióin con Método de Euler          ###
#######################################################

for k in range(num_steps):
    # Entrada en el instante actual
    u_k = u[k]

    # Derivada del estado: dx/dt = A*x + B*u
    dx_dt = A @ x + B.flatten() * u_k    #B.flatten() convierte la matriz B (que es 2D) en un vector 1D para que pueda
                                         #multiplicarse correctamente con el escalar u_k. 
                                         #Esto asegura que la operación sea compatible en términos de dimensiones.
 
    # Actualización del estado: x[k+1] = x[k] + dx/dt * dt
    x = x + dx_dt * dt

    # Guardar el estado y las salidas
    x_history[k, :] = x
    y1[k] = (C1 @ x).item()  # Extraer el valor escalar
    y2[k] = (C2 @ x).item()  # Extraer el valor escalar
    y3[k] = (C3 @ x).item()  # Extraer el valor escalar

# Crear la figura con cuatro subplots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 12))

# Graficar la señal de entrada u
ax1.plot(time, u, 'g')
ax1.set_title('Señal de entrada y salidas de los sistemas')
ax1.set_ylabel('Entrada u (V)')
ax1.grid(True)

# Graficar la respuesta de system1
ax2.plot(time, y1, 'b')
ax2.set_ylabel('Salida system1 = VR [V]')
ax2.grid(True)

# Graficar la respuesta de system2
ax3.plot(time, y2, 'r')
ax3.set_ylabel('Salida system2 = i [A]')
ax3.grid(True)

# Graficar la respuesta de system3
ax4.plot(time, y3, 'b')
ax4.set_xlabel('Tiempo [s]')
ax4.set_ylabel('Salida system3 = VC [V]')
ax4.grid(True)

plt.tight_layout()
plt.show()