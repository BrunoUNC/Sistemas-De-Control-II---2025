import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sympy as sym
from scipy.interpolate import interp1d
import os


#######################################################
###                Lectura de Datos                 ###
#######################################################
# Archivo y hojas de datos
archivo_excel = 'Curvas_Medidas_RLC_2025.xlsx'
hoja_datos_crudos = 'Hoja1'         # Datos crudos
hoja_nombres_columnas = 'Hoja2'     # Nombres de las columnas

if not os.path.isfile(archivo_excel):
    raise FileNotFoundError(f"The file '{archivo_excel}' does not exist.")

# Leer los nombres de las columnas desde la hoja 2
nombres_columnas = pd.read_excel(archivo_excel, sheet_name=hoja_nombres_columnas, engine='openpyxl').columns
print("Nombres de las columnas leídas:", nombres_columnas)

# Leer los datos crudos desde la hoja 1 y asignar los nombres de las columnas
data = pd.read_excel(archivo_excel, sheet_name=hoja_datos_crudos, engine='openpyxl', header=None)
data.columns = nombres_columnas
print("Columnas asignadas al DataFrame:", data.columns)

# Extraer columnas según las etiquetas de la hoja 2
try:
    tiempo = data['Tiempo [Seg.]'].to_numpy()  # Nota el punto al final
    corriente = data['Corriente [A]'].to_numpy()
    tension_capacitor = data['Tensión en el capacitor [V]'].to_numpy()
    tension_entrada = data['Tensión de entrada [V]'].to_numpy()
    tension_salida = data['Tensión de salida [V]'].to_numpy()
except KeyError as e:
    raise KeyError(f"Error: No se encontró la columna {e} en el DataFrame. Verifica los nombres de las columnas.")

#######################################################
###                Gráfico de señales               ###
#######################################################
plt.figure(figsize=(12, 10))

# Gráfica de la corriente
plt.subplot(4, 1, 1)
plt.plot(tiempo, corriente, label="Corriente [A]", color='blue')
plt.title("Señales del circuito RLC")
plt.ylabel("Corriente [A]")
plt.grid()
plt.legend()

# Gráfica de la tensión en el capacitor
plt.subplot(4, 1, 2)
plt.plot(tiempo, tension_capacitor, label="Tensión en el capacitor [V]", color='red')
plt.ylabel("Tensión [V]")
plt.grid()
plt.legend()

# Gráfica de la tensión de entrada
plt.subplot(4, 1, 3)
plt.plot(tiempo, tension_entrada, label="Tensión de entrada [V]", color='green')
plt.ylabel("Tensión [V]")
plt.grid()
plt.legend()

# Gráfica de la tensión de salida
plt.subplot(4, 1, 4)
plt.plot(tiempo, tension_salida, label="Tensión de salida [V]", color='purple')
plt.xlabel("Tiempo [s]")
plt.ylabel("Tensión [V]")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()


#######################################################
###       Cálculo de Función de Transferencia       ###
#######################################################

# Definir variables simbólicas
s, I, R, L, C, Ve, Vc = sym.symbols('s I R L C Ve Vc')

# Definir las ecuaciones diferenciales en el dominio de Laplace
eq1 = sym.Eq(s * I, (1 / L) * (Ve - Vc - R * I))  # Primera ecuación diferencial
eq2 = sym.Eq(s * Vc, (1 / C) * I)                 # Segunda ecuación diferencial

# Resolver el sistema de ecuaciones para Vc en función de Ve
soluciones = sym.solve([eq1, eq2], (I, Vc))
Vc_s = soluciones[Vc]  # Solución para Vc

# Calcular la función de transferencia H(s) = Vc(s) / Ve(s)
H_s = sym.simplify(Vc_s / Ve)

# Mostrar la función de transferencia
print("Función de transferencia H(s):", H_s)

#######################################################
###      Calculo de dVc/dt y I promedio             ###
###        (entre 10ms y 10.8ms)                    ###
#######################################################
t1 = 0.0105
t2 = 0.0108

# Interpolación para obtener Vc(t1), Vc(t2)
interp_vc = interp1d(tiempo, tension_capacitor, kind='linear')
vc_t1 = interp_vc(t1)
vc_t2 = interp_vc(t2)

# Regla del paralelogramo: pendiente
dvc_dt = (vc_t2 - vc_t1) / (t2 - t1)

# Corriente promedio entre 10ms y 11.5ms
mask_c = (tiempo >= t1) & (tiempo <= t2)
i_avg = np.mean(corriente[mask_c])

# Cálculo de C
C_est = i_avg / dvc_dt

#######################################################
###           Método de Chen (10ms a 20ms)          ###
#######################################################
# Definir los instantes t1 y t2 
t1_chen = 0.0103  # En segundos 
t2_chen = 0.0106
t3_chen = 0.0109


# Interpolación para obtener y_t1 y y_t2 directamente
interp_vc = interp1d(tiempo, tension_capacitor, kind='linear', fill_value="extrapolate")
interp_vin = interp1d(tiempo, tension_entrada, kind='linear', fill_value="extrapolate")

vc_t1 = interp_vc(t1_chen)
vc_t2 = interp_vc(t2_chen)
vc_t3 = interp_vc(t3_chen)
vin_t1 = interp_vin(t1_chen)
vin_t2 = interp_vin(t2_chen)
vin_t3 = interp_vin(t3_chen)


K = 12     # Ganancia del circuito

vc_t1 = interp_vc(t1_chen)/K
vc_t2 = interp_vc(t2_chen)/K
vc_t3 = interp_vc(t3_chen)/K

Kp=K/12


# Método de Chen
k1=(vc_t1/Kp)-1
k2=(vc_t2/Kp)-1
k3=(vc_t3/Kp)-1

b=4*k1**3*k3-3*k1**2*k2**2-4*k2**3+k3**2+6*k1*k2*k3

alpha_1=(k1*k2+k3-np.sqrt(b))/(2*(k1**2+k2))
alpha_2=(k1*k2+k3+np.sqrt(b))/(2*(k1**2+k2))

beta=(2*k1**3+3*k1*k2+k3-np.sqrt(b))/(np.sqrt(b))

t1_chen = t1_chen - 0.010                       # Origen del escalón en 10ms; normalizacion de la escala de tiempo.

T_1=-t1_chen/np.log(alpha_1)
T_2=-t1_chen/np.log(alpha_2)
T_3 = beta*(T_1-T_2)+T_1


# Calcular R y L con C estimada
L_est = (T_1 * T_2)/C_est
R_est = (T_1 + T_2)/C_est


#######################################################
###               Resultados Finales                ###
#######################################################
print("\n===== MÉTODO DE CHEN =====")
print(f"t1 = {t1_chen*1e3:.2f} ms, t2 = {t2_chen*1e3:.2f} ms")
print(f"alpha_1 = {alpha_1}, alpha_2 = {alpha_2}")
print(f"beta = {beta}")
print(f"y(t1) = {vc_t1:.4f}, y(t2) = {vc_t2:.4f}, y(t3) = {vc_t3:.4f}")
print(f"T1 = {T_1*1e3:.4f}, T2 = {T_2*1e3:.4f}, T3 = {T_3*1e3:.4f}")
print(f"K estimada = {K:.4f}")

print("\n===== COMPONENTES ESTIMADOS =====")
print(f"R estimada = {R_est:.2f} Ω")
print(f"L estimada = {L_est*1e3:.4f} mH")
print(f"C estimada = {C_est*1e6:.4f} µF")



#######################################################
###   Comparación: Tensión Calculada vs Experimental ###
#######################################################
from scipy.signal import lti, lsim

# Crear el sistema LTI a partir de los valores estimados
numerador = [1]  # Numerador de H(s) = 1 / (L C s^2 + R C s + 1)
denominador = [L_est * C_est, R_est * C_est, 1]  # Denominador de H(s)
sistema = lti(numerador, denominador)

# Simular la respuesta del sistema con la tensión de entrada experimental
t_sim, vc_calculada, _ = lsim(sistema, tension_entrada, tiempo)

# Crear la figura con dos subplots
plt.figure(figsize=(12, 8))

# Subplot 1: Tensión de entrada
plt.subplot(2, 1, 1)
plt.plot(tiempo, tension_entrada, label="Tensión de entrada", color='green', linewidth=2)
plt.title("Tensión de Entrada")
plt.xlabel("Tiempo [s]")
plt.ylabel("Tensión [V]")
plt.grid()
plt.legend()

# Subplot 2: Comparación de la tensión sobre el capacitor
plt.subplot(2, 1, 2)
plt.plot(tiempo, tension_capacitor, label="Tensión experimental (Excel)", color='red', linestyle='--', linewidth=3)
plt.plot(t_sim, vc_calculada, label="Tensión calculada (Método de Chen)", color='blue')
plt.title("Comparación: Tensión sobre el Capacitor")
plt.xlabel("Tiempo [s]")
plt.ylabel("Tensión [V]")
plt.grid()
plt.legend()

# Ajustar el diseño y mostrar el gráfico
plt.tight_layout()
plt.show()