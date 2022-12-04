import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


#constantes
r_t = 6371000 #metros
k = 10**-6 #metros/segundos
t_0 = 2000 + 273 #kelvin
t_s = 273 #kelvin
grad = 0.025 #grados/metro

def suma(n, t):
    """
    sumatoria de la expresión del gradiente
    :param n: cantidad de n para la sumatoria
    :param t: tiempo en millones de años
    :return: sumatoria
    """
    t = t * 60**2 * 24 * 365 * 10**6 # en segundos
    i = np.arange(1, n+1)
    cte = k * np.pi**2 / r_t**2
    potencia = np.exp(-cte * i**2 * t)
    sumatoria = np.sum(potencia)
    return sumatoria


def gradiente_sup(t, n=10**4):
    """
    gradiente en superficie para el modelo esférico
    :param n: numero de n
    :param t: tiempo en millones de años
    :return: gradiente superficial
    """
    grad = -2 * (t_0 - t_s) * suma(n, t) / r_t
    return grad

def funcion_cero(t, solucion=0.025):
    """
    funcion a la que hay que encontrarle el cero para encontrar el tiempo en el que el gradiente
    es el pedido (0.025 grados/metro)
    :param t: tiempo
    :param solucion: gradiente en el que queremos el tiempo
    :return: tiempo
    """
    return gradiente_sup(t) + solucion

def gradiente_semi_sup(t):
    """
    gradiente en superficie para el modelo del semiespacio
    :param t: tiempo
    :return: gradiente en superficie
    """
    t = t * 60**2 * 24 * 365 * 10**6 # en segundos
    grad = (t_0 - t_s) / np.sqrt(k * np.pi * t)

    return -grad

def funcion_cero_semi(t, solucion=0.025):
    """
    funcion a la que hay que encontrarle el cero para encontrar el tiempo en el que el gradiente
    es el pedido (0.025 grados/metro)
    :param t: tiempo
    :param solucion: gradiente en el que queremos el tiempo
    :return: tiempo
    """
    return gradiente_semi_sup(t) + solucion

sol_grad = fsolve(funcion_cero, 50) # tiempo en el que el gradiente es 25 grados por km
sol_grad_semi = fsolve(funcion_cero_semi, 50) # tiempo en el que el gradiente es 25 grados por km

print('tiempo para un gradiente de 25 grados por km(modelo esferico): ' + str(sol_grad) + 'Ma')
print('tiempo para un gradiente de 25 grados por km(modelo de semiespacio): ' + str(sol_grad_semi) + 'Ma')


n = 10**4
gradientes = []
gradientes_semi = []
tiempo = np. arange(30, 100)
for time in range(30, 100):
    gradientes.append(gradiente_sup(time))
    gradientes_semi.append(gradiente_semi_sup(time))


plt.figure(figsize=(10, 6))
plt.plot(tiempo, gradientes, label='modelo esférico', color='g')
plt.plot(tiempo, gradientes_semi, label='modelo del semiespacio', color='r')
plt.axhline(-0.025, linewidth=0.8, color='k')
plt.axvline(sol_grad, color='g', linewidth=0.8)
plt.axvline(sol_grad_semi, color='r', linewidth=0.8)
plt.title('Gráfico tiempo v/s gradiente superficial')
plt.xlabel('tiempo [millones de años]')
plt.ylabel('gradiente de temperatura superficial [$\dfrac{grados}{metro}$]')
plt.savefig('imagenes\grafico_gradiente_vs_tiempo')

plt.legend()
plt.show()
