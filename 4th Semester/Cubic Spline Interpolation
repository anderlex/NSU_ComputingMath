import math
import matplotlib.pyplot as plt
import numpy as np
a,b =  10,15
n = 10
A = -20
B = 6

def function(x):
    return 100

x = np.linspace(a, b, n)
y = np.zeros(n)
for i in range(n):
    y[i] = (function(x[i]))

#Создаем массив h длин подынтервалов:
h = np.zeros(n - 1)
for i in range(n-1):
    h[i] = x[i + 1] - x[i]

#Метод прогонки:
def Progonka_method(h):
    C = np.zeros(n)
    Alphai = np.zeros(n)
    Betai = np.zeros(n)
    Alphai[1] = -0.5
    Betai[1] = ((y[1]-y[0])/h[0]**2 - A/h[0])*3
    for i in range(1, n - 1):
        Ai = h[i - 1]/6
        Bi = (h[i - 1]+h[i])/3
        Ci = h[i]/6
        Fi = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        Alphai[i + 1] = -Ci/(Ai*Alphai[i]+Bi)
        Betai[i + 1] = (Fi - Ai*Betai[i])/(Ai*Alphai[i]+Bi)
    Ai = h[n - 2]/6
    Bi = h[n - 2]/3
    C[n - 1] = (B - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]) - Ai*Betai[n-1])/(Ai*Alphai[n-1]+Bi)
    for i in reversed(range(n - 1)):
        C[i] = Alphai[i + 1]*C[i + 1] + Betai[i + 1]
    return C

def Spline_interpolation(x0):
    if x0 <= x[0]:
        i = 0
    elif x0 >= x[n - 1]:
        i = n - 2
    else:
        i = 0
        j = n - 1
        while i + 1 < j:
            k = i + (j - i)//2
            if x0 <= x[k]:
                j = k
            else:
                i = k
    h = x[i+1]-x[i]
    u = -C[i]/(6*h)
    v = C[i+1]/(6*h)
    C1 = ((y[i+1]-y[i])/h+h/6*(C[i]-C[i+1]))
    C2 = y[i] - C[i]*h*h/6
    Si = u*(x0 - x[i+1])**3+v*(x0-x[i])**3+C1*(x0-x[i])+C2
    return Si

#С помощью метода прогонки находим неизвестные коэффициенты С:
C = Progonka_method(h)
dots = 3000

#Создаем массивы для рисования интерполируемой функции:
x_ris = np.linspace(a, b, dots)
y_ris = np.zeros(dots)
S = np.zeros(dots)
for i in range(dots):
    y_ris[i] = function(x_ris[i])
    S[i] = Spline_interpolation(x_ris[i])

#Собственно, рисуем:
plt.style.use('ggplot')
fig, ax = plt.subplots()
ax.plot(x_ris, y_ris, x_ris, S)
ax.scatter(x,y)
ax.set_title("Интерполяция сплайном", fontsize = 24)
ax.set_xlabel("Узлы", fontsize = 14)
ax.set_ylabel("Значения", fontsize = 14)
ax.tick_params(axis = 'both', labelsize = 10)
plt.show()

#Задаем эпсилон > 0 и приближенно вычисляем производные на концах отрезка
eps = 1.0E-5
dfda = 1/eps*(Spline_interpolation(a+eps)-Spline_interpolation(a))
dfdb = 1/eps*(Spline_interpolation(b)-Spline_interpolation(b-eps))
print(dfda)
print(dfdb)
