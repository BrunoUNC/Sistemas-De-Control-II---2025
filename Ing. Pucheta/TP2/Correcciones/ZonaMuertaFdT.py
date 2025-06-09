import numpy as np
import matplotlib.pyplot as plt

# Vector de entrada (tensión de entrada)
ue = np.arange(-2.5, 2.5 + 0.1, 0.1)
N = len(ue)
uo = np.zeros(N)

# Zona muerta simétrica centrada en cero
ZM = 0.5  # Umbral de zona muerta

for i in range(N):
    if abs(ue[i]) > ZM:
        uo[i] = ue[i] - ZM * np.sign(ue[i])

# Gráfico
plt.plot(ue, uo, 'k')
plt.xlabel('Tensión de entrada')
plt.ylabel('$V_o$')
plt.title('Zona muerta aplicada')
plt.grid(True)
plt.show()
