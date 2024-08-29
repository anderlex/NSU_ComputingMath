import math
import matplotlib.pyplot as plt
import copy

def f(x):
	if x > 0:
		return 100
	else: return 101

def Newton(x, y, z):
    sum = y[0]
    n = len(x)
    for i in range(1, n):
        F = 0
        for j in range(i + 1):
            den = 1
            for k in range(i + 1):
                if k != j:
                    den *= x[j] - x[k]
            F += y[j]/den
        for k in range(i):
            F *= z - x[k]
        sum += F
    return sum


# Задаём начальные значения
N = 6					
Xi = [0] * (N + 1)
Yi = [0]*(N + 1)
a = -15
b = 20

Xi[0] = a
Yi[0] = f(a)
Xi[N] = b
Yi[N] = f(b)


for i in range(1, N):
    Xi[i] = a + ((b - a)/(N))*(i + 0.45*math.cos(5*i))
    Yi[i] = f(Xi[i])

num = 500
x = [0] * (num + 1)
y = [0] * (num + 1)

for i in range(num):
	x[i] = a + ((b - a)/num) * i
	y[i] = f(x[i])

x[num] = b
y[num] = f(b)

Yn = [0] * (num + 1)
axs = copy.deepcopy(Xi)
for i in range(num + 1):
    Yn[i] = Newton(Xi, Yi, x[i])

plt.style.use('ggplot')
fig, ax = plt.subplots()
ax.plot(x, y,c='red', linewidth = 2)
ax.plot(x, Yn,c='blue', linewidth = 1)
for i in range(N + 1):
	ax.scatter(axs[i],Yi[i],c='green',linewidth = 2)
ax.set_title("Интерполяция", fontsize = 24)
ax.set_xlabel("Узлы", fontsize = 14)
ax.set_ylabel("Значения", fontsize = 14)
ax.tick_params(axis = 'both', labelsize = 10)

plt.show()
