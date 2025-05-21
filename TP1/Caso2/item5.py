import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.signal
from numpy.linalg import eig
from scipy.linalg import expm
from scipy.interpolate import interp1d
from scipy.signal import lti, lsim, step
from sympy import *
import control as ctl
import cmath as cm

print(plt.style.available)

#######################################################
###                Lectura de Datos                 ###
#######################################################
# Archivo y hojas de datos
archivo_excel = 'Curvas_Medidas_Motor_2025.xlsx'  
hoja_datos_crudos = 'Hoja1'         # Datos crudos
hoja_nombres_columnas = 'Hoja2'     # Nombres de las columnas

if not os.path.isfile(archivo_excel):
    raise FileNotFoundError(f"The file '{archivo_excel}' does not exist.")

# Leer los datos crudos y los nombres de las columnas
df_nombres_columnas = pd.read_excel(archivo_excel, sheet_name=hoja_nombres_columnas, engine='openpyxl')
df_datos = pd.read_excel(archivo_excel, sheet_name=hoja_datos_crudos, engine='openpyxl', header=None)

# Asignar nombres de columnas al DataFrame de datos
df_datos.columns = df_nombres_columnas.columns

# Extraer columnas directamente del DataFrame
tiempo = df_datos['Tiempo [Seg.]']
Velocidad_angular = df_datos['Velocidad angular [rad /seg]']
Corriente_en_armadura = df_datos['Corriente en armadura [A]']
Tensión = df_datos['Tensión [V]']
torque = df_datos['Torque']

#######################################################
###                Gráfico de señales               ###
#######################################################
plt.figure(figsize=(12, 10))

# Gráfica del tension
plt.subplot(4, 1, 1)
plt.plot(tiempo, Tensión, label="Tensión [V]", color='red')
plt.xlabel("Tiempo [s]")
plt.ylabel("Tensión [V]")
plt.title("Curvas medidas del motor")
plt.grid()
plt.legend()


# Gráfica de la corriente
plt.subplot(4, 1, 2)
plt.plot(tiempo, Corriente_en_armadura, label="Corriente en armadura [A]", color='blue')
plt.ylabel("Corriente [A]")
plt.grid()
plt.legend()

# Gráfica de la velocidad
plt.subplot(4, 1, 3)
plt.plot(tiempo, Velocidad_angular, label="Velocidad angular [rad/s]", color='orange')
plt.ylabel("Velocidad [rad/s]")
plt.grid()
plt.legend()

# Gráfica del torque
plt.subplot(4, 1, 4)
plt.plot(tiempo, torque, label="Torque TL [Nm]", color='green')
plt.xlabel("Tiempo [s]")
plt.ylabel("Torque [Nm]")
plt.grid()
plt.legend()



plt.tight_layout()
plt.show()

#######################################################
###            Funciones de Transferencia           ###
#######################################################

# Definición de variables simbólicas

Va, Ra, La, Ia, Ki, Wr, Jm, Km, Bm, TL, s=symbols('Va Ra La Ia Ki Wr Jm Km Bm TL s')

# Pasivando TL, se tiene la Wr/VA como
# FdT: Wr/Va
Eqn1=Eq(s*Ia,-(Ra/La)*Ia-(Km/La)*Wr+(1/La)*Va)
res01=solve(Eqn1,Ia)
Ia=res01[0]  #Despejo Ia
Eqn2=Eq(s*Wr, Ki/Jm*res01[0]-Bm/Jm*Wr) #Sustituyo Ia en la ecuacion de Wr
Wr_Va=solve(Eqn2,Wr)   #Cuidado, Va queda en el numerador
print('Wr/Va=')
pretty_print(Wr_Va)

# Pasivando Va, se tiene la Wr/TL como
# FdT: Wr/TL

Eqn3=Eq(s*Ia,-(Ra/La)*Ia - (Km/La)*Wr)      #Primera ecuacion diferencial pasivando Va
res03=solve(Eqn3,Ia)                        #Despejo Ia 
Eqn4=Eq(s*Wr, (Ki/Jm)*res03[0] - (Bm/Jm)*Wr-(1/Jm)*TL ) #Sustituyo Ia en la ecuacion de Wr
Wr_TL=solve(Eqn4,Wr)   #Cuidado, TL queda en el numerador
print('Wr/TL=')   
pretty_print(Wr_TL)


# Normalización de la ecuacion de Wr/Va y Wr/TL para tener a0=1
# En este paso ya elimino Va y TL de los numeradores
Wr_Va=Ki/(Jm*La)/(s**2+(Bm*La+Jm*Ra)/(Jm*La)*s+(Bm*Ra  + Ki*Km)/(Jm*La))

Wr_TL=(La/Ra*s+1)/(Ra*Jm*La)/(s**2+(Bm*La+Jm*Ra)/(Jm*La)*s+(Bm*Ra  + Ki*Km)/(Jm*La))

pretty_print(Wr_Va)
pretty_print(Wr_TL)


#######################################################
###                  Método de Chen                 ###
#######################################################

#Cálculo de Ra:
pico_corriente = Corriente_en_armadura.max()
max_tension = Tensión.max()
Ra_est = max_tension / pico_corriente

print(f"Pico de Corriente: {pico_corriente} [A]")
print(f"Tensión máxima: {max_tension} [V]")
print(f"Ra_est calculado: {Ra_est} [Ohm]")

# Definición de variables simbólicas para el método de Chen
#Variables simbólicas p se corresponden con el Chen de wr/TL

T1,T2,T3,T3p,K1_C,K2_C, s=symbols('T1 T2 T3 T3p K1_C K2_C s')

Wr_Va_Chen=K1_C/(T1*T2*s**2+(T1+T2)*s+1)

Wr_TL_Chen=K2_C*(T3*s+1)/(T1*T2*s**2+(T1+T2)*s+1)

print('Wr/Va Chen=')
pretty_print(Wr_Va_Chen)
print('Wr/TL Chen=')
pretty_print(Wr_TL_Chen)

# Normalizaciión de la ecuacion de Wr/Va y Wr/TL para tener a0=1
Wr_Va_Chen=K1_C/(T1*T2)/(s**2+(T1+T2)/(T1*T2)*s+1/(T1*T2))
Wr_TL_Chen=K2_C*(T3*s+1)/(T1*T2)/( s**2+(T1+T2)/(T1*T2)*s+1/(T1*T2))
print('Wr/Va Chen=')
pretty_print(Wr_Va_Chen)
print('Wr/TL Chen=')
pretty_print(Wr_TL_Chen)

#Igualando iguales potencias de s entre FdT del motor y FdT del modelo de Chen:
Eqn5=Eq(K1_C/(T1*T2), Ki/(Jm*La))
Eqn6=Eq(K2_C/(T1*T2), 1/(Ra*Jm*La))
Eqn7=Eq(T3p, La/Ra)
Eqn8=Eq((Bm*La + Jm*Ra)/(Jm*La),(T1+T2)/(T1*T2))
Eqn9=Eq((Bm*Ra + Ki*Km)/(Jm*La),1/(T1*T2))

#Solución de los parámetros del motor en función de los parámetros del método de Chen:
sol = solve((Eqn5, Eqn6, Eqn7, Eqn8, Eqn9),(Ra, La, Ki, Jm, Km, Bm))
print('Ra, La, Ki, Jm, Km, Bm')
pretty_print(sol[0])

##########################################
# Aplicación práctica del Método de Chen #
##########################################
#########       wr/VA:      ##############
# Definir los instantes t1 y t2 
t1_chen = 0.102  # En segundos 
t2_chen = 0.103
t3_chen = 0.104

# Interpolación para obtener y_t1 y y_t2 directamente
interp_wr = interp1d(tiempo, Velocidad_angular, kind='linear', fill_value="extrapolate")
interp_va = interp1d(tiempo, Tensión, kind='linear', fill_value="extrapolate")

wr_t1 = interp_wr(t1_chen)
wr_t2 = interp_wr(t2_chen)
wr_t3 = interp_wr(t3_chen)

va_t1 = interp_va(t1_chen)
va_t2 = interp_va(t2_chen)
va_t3 = interp_va(t3_chen)

K = 7.6246
K1C = K / 2  # Ganancia del sistema

#Parámetros del Método de Chen:
k1=(wr_t1/K1C)-1
k2=(wr_t2/K1C)-1
k3=(wr_t3/K1C)-1

be =  4 * (k1**3) * k3 - 3 * (k1**2) * (k2**2) - 4 * (k2**3) + (k3**2) + 6 * k1 * k2 * k3

if be < 0:
    alfa1 = (k1 * k2 + k3 - np.sqrt(be)) / (2 * (k1**2 + k2))
    alfa2 = (k1 * k2 + k3 + np.sqrt(be)) / (2 * (k1**2 + k2))
else:
    alfa1 = (k1 * k2 + k3 - np.sqrt(be)) / (2 * (k1**2 + k2))
    alfa2 = (k1 * k2 + k3 + np.sqrt(be)) / (2 * (k1**2 + k2))
beta = (2 * k1**3 + 3 * k1 * k2 + k3 - np.sqrt(be)) / np.sqrt(be)

t1_chen_n = t1_chen - 0.100             # Origen del escalón en 10ms; normalizacion de la escala de tiempo.

T1 = np.real(-t1_chen_n / np.log(alfa1))
T2 = np.real(-t1_chen_n / np.log(alfa2))
T3 = np.real( beta * (T1 - T2) + T1)


print("\n===== MÉTODO DE CHEN Va=====")
print(f"t1 = {t1_chen*1e3:.2f} ms, t2 = {t2_chen*1e3:.2f} ms, t2 = {t3_chen*1e3:.2f} ms")
print(f"y(t1) = {wr_t1:.4f}, y(t2) = {wr_t2:.4f}, y(t3) = {wr_t3:.4f}")
print(f"alpha_1 = {alfa1}, alpha_2 = {alfa2}")
print(f"beta = {beta}")
print(f"T1 = {T1*1e3:.4f}, T2 = {T2*1e3:.4f}")
print(f"K estimada = {K1C:.4f}")

FT = K1C  / ((T1 * s + 1) * (T2 * s + 1))
print("Función de transferencia FT:")
pretty_print(FT)

num_FT = [K1C]  # Numerador de FT
den_FT = [T1 * T2, (T1 + T2), 1]  # Denominador de FT
FT_lti = lti(num_FT, den_FT)  # Crear el sistema LTI

# Simular la respuesta de FT con lsim
t_sim, y_sim1, _ = lsim(FT_lti, U=Tensión, T=tiempo)


# Graficar la respuesta al escalón de FT escalón 2V
plt.figure(figsize=(10, 6))
plt.plot(tiempo, Tensión, 'g--', label="Entrada escalón (2 V)")
plt.plot(t_sim, y_sim1, label="Respuesta al escalón wr/Va(FT, simulada)", linestyle='-', color='blue')
plt.plot(tiempo, Velocidad_angular, label="Datos medidos (Excel)", linestyle='-', color='red')
plt.title('Respuesta al escalón de la función de transferencia')
plt.xlabel('Tiempo (s)')
plt.ylabel('Salida')
plt.grid()
plt.legend()
plt.show()


##########################################
# Aplicación práctica del Método de Chen #
##########################################
#########       wr/TL:      ##############

t1_chen = 0.702  # En segundos 
t2_chen = 0.712
t3_chen = 0.722

# Interpolación para obtener y_t1 y y_t2 directamente
interp_wr = interp1d(tiempo, Velocidad_angular, kind='linear', fill_value="extrapolate")
interp_TL = interp1d(tiempo, torque, kind='linear', fill_value="extrapolate")

wr_t1 = -(interp_wr(t1_chen)-7.6245)
wr_t2 = -(interp_wr(t2_chen)-7.6245)
wr_t3 = -(interp_wr(t3_chen)-7.6245)
TL_t1 = interp_TL(t1_chen)
TL_t2 = interp_TL(t2_chen)
TL_t3 = interp_TL(t3_chen)

K2 = 0.12
K2C = (7.6245 - 3.6379)/K2 # Ganancia del circuito para wr/TL
    
# Parámetros del método de Chen:
k1=(1/K2)*(wr_t1/K2C)-1
k2=(1/K2)*(wr_t2/K2C)-1
k3=(1/K2)*(wr_t3/K2C)-1

be =  4 * k1**3 * k3 - 3 * k1**2 * k2**2 - 4 * k2**3 + k3**2 + 6 * k1 * k2 * k3

if be > 0:
    alfa1 = (k1 * k2 + k3 - np.sqrt(be)) / (2 * (k1**2 + k2))
    alfa2 = (k1 * k2 + k3 + np.sqrt(be)) / (2 * (k1**2 + k2))
else:
    alfa1 = (k1 * k2 + k3 - np.sqrt(be)) / (2 * (k1**2 + k2))
    alfa2 = (k1 * k2 + k3 + np.sqrt(be)) / (2 * (k1**2 + k2))
beta = (2 * k1**3 + 3 * k1 * k2 + k3 - np.sqrt(be)) / np.sqrt(be)

t1_chen_n = t1_chen - 0.701                      # Normalizacion de la escala de tiempo.

T_1p=np.real(-t1_chen_n/np.log(alfa1))
T_2p=np.real(-t1_chen_n/np.log(alfa2))
T_3p =np.real( beta*(T_1p-T_2p)+T_1p)

print("\n===== MÉTODO DE CHEN TL=====")
print(f"t1 = {t1_chen*1e3:.2f} ms, t2 = {t2_chen*1e3:.2f} ms, t2 = {t3_chen*1e3:.2f} ms")
print(f"y(t1) = {wr_t1:.4f}, y(t2) = {wr_t2:.4f}, y(t3) = {wr_t3:.4f}")
print(f"alpha_1 = {alfa1}, alpha_2 = {alfa2}")
print(f"beta = {beta}")
print(f"T1_p = {T_1p*1e3:.4f}, T2_p = {T_2p*1e3:.4f}, T3_p = {T_3p*1e3:.4f}")
print(f"K estimada = {K2C:.4f}")



FT_1 = (K2C * (T_3p*s + 1) )/ ((T1 * s + 1) * (T2 * s + 1))
print("Función de transferencia wr/TL:")
pretty_print(FT_1)

num_FT_1 = [K2C * T_3p, K2C]  # Numerador de FT_1
den_FT_1 = [T1 * T2, (T1 + T2), 1]  # Denominador de FT_1
FT_1_lti = lti(num_FT_1, den_FT_1)  # Crear el sistema LTI


# Simular la respuesta de FT_1 con lsim
t_sim, y_sim2, _ = lsim(FT_1_lti, U=torque, T=tiempo)

plt.figure(figsize=(10, 8))

# Subplot 1: Torque = entrada
plt.subplot(2, 1, 1)
plt.plot(tiempo, torque, 'g--', label="Torque (0.12)")
plt.title('Entrada: Torque')
plt.xlabel('Tiempo (s)')
plt.ylabel('Torque [Nm]')
plt.grid()
plt.legend()

# Subplot 2: Respuesta al escalón y datos medidos de wr/TL
plt.subplot(2, 1, 2)
plt.plot(tiempo, y_sim2, label="Respuesta al escalón wr/TL (FdT calculada)", linestyle='-', color='blue')
plt.plot(tiempo, Velocidad_angular, label="Datos medidos (Excel)", linestyle='-', color='red')
plt.title('Respuesta al escalón de la función de transferencia')
plt.xlabel('Tiempo (s)')
plt.ylabel('Velocidad angular [rad/s]')
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()


######################################################################################


# # Cálculo de los parámetros del motor a partir de los parámetros del método de Chen y las ecuaciones sol.

Ra_est = Ra_est                                     
La_est = Ra_est*T_3p
Ki_est = (K1C)/(K2C*Ra_est)
Jm_est = (T1*T2)/(K2C*Ra_est**2*T_3p)
Km_est = ((T1-T_3p)*(T2-T_3p))/(K1C*T_3p**2)
Bm_est = ((-T1*T2+T1*T_3p+T2*T_3p))/(K2C*Ra_est**2*T_3p**2)


# #######################################################
# ###               Resultados Finales                ###
# #######################################################


print("\n===== COMPONENTES ESTIMADOS =====")


print(f"Ra estimada = {Ra_est:.2f}")
print(f"La estimada = {La_est:.4f}")
print(f"Ki estimada = {Ki_est:.4f}")
print(f"Km estimada = {Km_est:.2f}")
print(f"Jm estimada = {Jm_est:.4f}")
print(f"Bm estimada = {Bm_est:.4f}")


Wr_Va_est= (Ki_est/(Jm_est*La_est)) / (s**2  +  ((Bm_est*La_est+Jm_est*Ra_est)/(Jm_est*La_est))*s  +  ((Bm_est*Ra_est+Ki_est*Km_est)/(Jm_est*La_est)))

Wr_TL_est=(((La_est/Ra_est)*s+1)/(Ra_est*Jm_est*La_est)) / (s**2  +  ((Bm_est*La_est+Jm_est*Ra_est)/(Jm_est*La_est))*s  +  ((Bm_est*Ra_est + Ki_est*Km_est)/(Jm*La)))

# Definir los coeficientes de las funciones de transferencia Wr_Va_est y Wr_TL_est
# Wr_Va_est = (Ki_est / (Jm_est * La_est)) / (s^2 + ((Bm_est * La_est + Jm_est * Ra_est) / (Jm_est * La_est)) * s + ((Bm_est * Ra_est + Ki_est * Km_est) / (Jm_est * La_est)))
num_Wr_Va = [Ki_est / (Jm_est * La_est)]
den_Wr_Va = [1, (Bm_est * La_est + Jm_est * Ra_est) / (Jm_est * La_est), (Bm_est * Ra_est + Ki_est * Km_est) / (Jm_est * La_est)]

# Wr_TL_est = ((La_est / Ra_est) * s + 1) / (Ra_est * Jm_est * La_est) / (s^2 + ((Bm_est * La_est + Jm_est * Ra_est) / (Jm_est * La_est)) * s + ((Bm_est * Ra_est + Ki_est * Km_est) / (Jm_est * La_est)))
num_Wr_TL = [(La_est / Ra_est)/(Jm_est*La_est*Ra_est), 1/(Jm_est*La_est*Ra_est)]
den_Wr_TL = [1, (Bm_est * La_est + Jm_est * Ra_est) / (Jm_est * La_est), (Bm_est * Ra_est + Ki_est * Km_est) / (Jm_est * La_est)]

# Crear sistemas LTI para ambas funciones de transferencia
system_Wr_Va = lti(num_Wr_Va, den_Wr_Va)
system_Wr_TL = lti(num_Wr_TL, den_Wr_TL)
print("\nFunción de transferencia teórica Wr_Va con parámetros estimados:")
print(num_Wr_Va)
print(den_Wr_Va)
print("\nFunción de transferencia teórica Wr_TL con parámetros estimados:")
print(num_Wr_TL)
print(den_Wr_TL)

# Simular la respuesta de Wr_Va_est con la entrada Tensión
_, y_Wr_Va, _ = lsim(system_Wr_Va, Tensión, tiempo)  # Usar directamente las señales del archivo

# Simular la respuesta de Wr_TL_est con la entrada torque
_, y_Wr_TL, _ = lsim(system_Wr_TL, torque, tiempo)  # Usar directamente las señales del archivo

# Graficar las respuestas en subplots separados
plt.figure(figsize=(12, 12))

# Entrada y salida de Wr_Va_kest
plt.subplot(5, 1, 1)
plt.plot(tiempo, Tensión, label="Entrada: Tensión [V]", color="red", linestyle="--")
plt.xlabel("Tiempo [s]")
plt.ylabel("Tensión [V]")
plt.title("Entrada de Wr_Va_est")
plt.grid()
plt.legend()

plt.subplot(5, 1, 2)
plt.plot(tiempo, y_Wr_Va, label="Salida: Velocidad angular [rad/s] (Wr_Va_est)", color="blue")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad angular [rad/s]")
plt.title("Salida de Wr_Va_est")
plt.grid()
plt.legend()

# Entrada y salida de Wr_TL_est
plt.subplot(5, 1, 3)
plt.plot(tiempo, torque, label="Entrada: Torque [Nm]", color="green", linestyle="--")
plt.xlabel("Tiempo [s]")
plt.ylabel("Torque [Nm]")
plt.title("Entrada de Wr_TL_est")
plt.grid()
plt.legend()

plt.subplot(5, 1, 4)
plt.plot(tiempo, y_Wr_TL, label="Salida: Velocidad angular [rad/s] (Wr_TL_est)", color="orange")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad angular [rad/s]")
plt.title("Salida de Wr_TL_est")
plt.grid()
plt.legend()

plt.subplot(5, 1, 5)
plt.plot(tiempo, (y_Wr_Va-y_Wr_TL), label="Modelo teórico con parámetros estimados", color="blue")
plt.plot(tiempo, Velocidad_angular, label="velocidad angular medida", color="red")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad angular [rad/s]")
plt.title("Salida de Wr_TL_est")
plt.grid()
plt.legend()


plt.tight_layout()
plt.show()

# Calcular la respuesta al escalón para Wr_Va_est (escalón de amplitud 2)
t_Wr_Va_step, y_Wr_Va_step = step(system_Wr_Va)
y_Wr_Va_step = 2 * y_Wr_Va_step  # Escalón de amplitud 2

# Calcular la respuesta al escalón para Wr_TL_est (escalón de amplitud 0.12)
t_Wr_TL_step, y_Wr_TL_step = step(system_Wr_TL)
y_Wr_TL_step = 0.12 * y_Wr_TL_step  # Escalón de amplitud 0.12

# Graficar las respuestas
plt.figure(figsize=(12, 8))

# Respuesta de Wr_Va_est
plt.subplot(2, 1, 1)
plt.plot(t_Wr_Va_step, y_Wr_Va_step, label="Salida: Velocidad angular [rad/s] (Wr_Va_est)", color="blue")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad angular [rad/s]")
plt.title("Respuesta al escalón de Wr_Va_est (amplitud 2)")
plt.grid()
plt.legend()

# Respuesta de Wr_TL_est
plt.subplot(2, 1, 2)
plt.plot(t_Wr_TL_step, y_Wr_TL_step, label="Salida: Velocidad angular [rad/s] (Wr_TL_est)", color="orange")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad angular [rad/s]")
plt.title("Respuesta al escalón de Wr_TL_est (amplitud 0.12)")
plt.grid()
plt.legend()


plt.tight_layout()
plt.show()

#################################################################

# Simulación de la respuesta de FT y FT_1 usando torque como entrada

# Simulación de FT con entrada de tensión
t_FT, y_FT, _ = lsim(FT_lti, U=Tensión, T=tiempo, X0=[0, 0])

# Simulación de FT_1 con entrada de torque
t_FT1, y_FT1, _ = lsim(FT_1_lti, U=torque, T=tiempo, X0=[0, 0])

# Modelo combinado: suma de ambas respuestas
y_modelo = y_FT - y_FT1

plt.figure(figsize=(12, 10))


# Subplot 1: FT wr/Va
plt.subplot(3, 1, 1)
plt.plot(tiempo, y_FT, '--', color='green', label='FT wr/Va')
plt.ylabel('Velocidad [rad/s]')
plt.title('FT wr/Va - Respuesta al escalón para escalón de 2V')
plt.grid()
plt.legend()

# Subplot 2: FT_1 wr/TL
plt.subplot(3, 1, 2)
plt.plot(tiempo, -y_FT1, '--', color='orange', label='FT wr/TL')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad [rad/s]')
plt.title('FT wr/TL - Respuesta al escalón para escalón de 0.12Nm')	
plt.grid()
plt.legend()

# Subplot 3: Velocidad angular medida (Excel) y Modelo combinado superpuestos
plt.subplot(3, 1, 3)
plt.plot(tiempo, Velocidad_angular, color='red', label='Velocidad angular (Datos Excel)')
plt.plot(tiempo, y_modelo, color='blue', label='Modelo (FT wr/Va + FT wr/TL)')
plt.ylabel('Velocidad [rad/s]')
plt.title('Comparación: Curva dato vs Modelo')
plt.grid()
plt.legend()


plt.tight_layout()
plt.show()