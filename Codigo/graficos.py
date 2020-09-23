from prob import *
import matplotlib.pyplot as plt

#codigo para graficar el integrando original de Gamma
x = np.linspace(0, 10, num=500)
y = np.zeros(len(x))
counter = 0
for xi in x:
    y[counter] = xi**(K/2-1)*np.exp(-xi)
    counter += 1
plt.figure(1)
plt.clf()
plt.plot(x, y)
plt.fill_between(x, y)
plt.xlabel('x', fontsize=13)
plt.ylabel('$y = x^{k/2 - 1}e^{-x}$', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

#codigo para graficar el integrando modificado de Gamma
x = np.linspace(0, 1)
y = np.zeros(len(x))
counter = 0
for xi in x:
    y[counter] = Gamma_du(K/2, xi)
    counter += 1
plt.figure(2)
plt.clf()
plt.plot(x, y)
plt.xlabel('x', fontsize=13)
plt.ylabel(r'$y = (-\ln (x))^{K/2 - 1}$', fontsize=15)
plt.fill_between(x, y)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

#codigo para graficar Gamma
x = np.linspace(1, 3)
y = np.zeros(len(x))
counter = 0
for xi in x:
    y[counter] = Gamma(xi, dx = 1e-3)
    counter += 1
plt.figure(3)
plt.clf()
plt.plot(x, y)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

#graficos de tolerancias comentados porque demoran mucho

#codigo para graficar distintas tolerancias relativas y determinar el valor a utilizar 
""" tolerancias_de_prueba = np.logspace(-6, -1)
error = np.zeros(len(tolerancias_de_prueba))

plt.clf()

counter = 0
for tol in tolerancias_de_prueba:
    dif1 = np.abs(Gamma(1, dx = tol) - 1) #Gamma(1) = 1
    dif2 = np.abs(Gamma(2, dx = tol) - 1) #Gamma(2) = 1
    dif3 = np.abs(Gamma(3, dx = tol) - 2) #Gamma(3) = 2
    dif4 = np.abs(Gamma(4, dx = tol) - 6) #Gamma(4) = 6
    error[counter] = np.mean(np.array([dif1, dif2, dif3, dif4]))
    counter += 1
plt.figure(4)
plt.plot(tolerancias_de_prueba, error)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Tolerancia relativa para Simpson', fontsize=15)
plt.ylabel('Promedio de errores absolutos cometido', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show() """

""" tolerancias_de_prueba_2 = np.logspace(-12, -2)
error_2 = np.zeros(len(tolerancias_de_prueba_2))

plt.clf()

counter = 0
for tol_2 in tolerancias_de_prueba_2:
    a = biseccion(prob_menos_95, K, 5, 15, tol_abs=tol_2)
    dif = np.abs(prob_menos_95(K, a3))
    error_2[counter] = dif
    counter += 1
plt.figure(5)
plt.plot(tolerancias_de_prueba_2, error_2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Tolerancia absoluta para bisección', fontsize=15)
plt.ylabel('Error absoluto cometido', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show() """


#finalmente se graficará la forma de la función chi2, así como el valor a obtenido
x = np.linspace(0, 1.1*a3)
y = chi2(K, x)

x1 = np.linspace(0, a3)
y1 = chi2(K, x1)

plt.figure(6)
plt.clf()
plt.title(r'$\chi ^2 (x)$ para $k=4.495$')
plt.xlabel('x', fontsize=15)
plt.ylabel(r'y = $\chi ^2 (x)$', fontsize=15)
plt.axvline(a3, label = 'a = {}'.format(a3), color = 'r')
plt.fill_between(x1, 0, y1, label = r'$P(x < a) = \int_0^a \chi^2 (x) dx$')
plt.plot(x, y, label=r'$\chi^2(x)$')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.show()

x=np.linspace(0,20)
y=np.zeros(len(x))
i=0
for xi in x:
    y[i]=prob_chi2(K, xi)
    i+=1
plt.figure(7)
plt.plot(x,y, color='green')
plt.xlabel('x', fontsize=15)
plt.ylabel(r'y=$\int_0^x \chi^2(t)dt$', fontsize=15)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.axhline(y=0.95, color='r', label='y=0.95')
plt.axvline(x=a3, color='b', label='x={}'.format(a3))
plt.legend()
plt.show()