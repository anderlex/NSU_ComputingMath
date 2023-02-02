import math
import numpy as np
import matplotlib.pyplot as plt

a,b = 0, math.pi
N = 21
h = (b - a)/(N - 1)

def K(x,s):
	return s + math.sin(x)
	#return x*math.exp(x*s)
	#return -x*s/2

def f(x):
	return x + math.pi**2*(1/2*math.sin(x)+1/3*math.pi)
	#return math.exp(x)
	#return 5/6*x

def U(x):
	return x
	#return 1
	#return x

x = np.linspace(a, b, N)
F = np.zeros(N)
for i in range(N):
	F[i] = f(x[i])

U_exact = np.zeros(N)
for i in range(N):
	U_exact[i] = U(x[i])

def reflection_method(A, b):
	F = b.copy()
	U = np.zeros(A.shape[1])
	R = A.shape[1]
	N = A.shape[0]
	#Метод отражений
	for i in range(N - 1):
		E = np.eye(N - i)
		x = np.around(A[i:, i], decimals = 4)
		e = np.zeros(x.size)
		e[0] = 1
		if np.dot(x,x) == 0:
			continue
		else:
			if np.array_equal(x/math.sqrt(np.dot(x, x)),e):
				u = x/math.sqrt(np.dot(x, x))
			else:
				if (x[0] > 0):
					ui = x + math.sqrt(np.dot(x, x))*e
				else:
					ui = x - math.sqrt(np.dot(x, x))*e
				u = ui/math.sqrt(np.dot(ui,ui))
			A[i:N,i:N] = A[i:N,i:N] - np.dot(2*np.outer(u,u), A[i:N,i:N])
			F[i:N] = F[i:N] - np.dot(2*np.outer(u,u),F[i:N])
	#Обратный ход Гаусса
	U[R - 1] = F[R - 1]/A[R-1,R-1]
	for k in reversed(range(R - 1)):
		summa = 0
		for j in range(k + 1, R):
			summa += A[k][j]*U[j]
		U[k] = (F[k] - summa)/A[k][k]
	return U

#Матрица для метода правых прямоугольников

A_right_rectangle = np.eye((N))
for i in range(N):
	for j in range(1, N):
		A_right_rectangle[i][j] += h*K(x[i],x[j])
U_right_rectangle = reflection_method(A_right_rectangle, F)

#Матрица для метода Симпсона

A_simpson = np.eye((N))
for i in range(N):
	for j in range(N):
		if j == 0 or j == N - 1:
			A_simpson[i][j] += (h/3)*K(x[i],x[j])
		elif j % 2 != 0:
			A_simpson[i][j] += (4*h/3)*K(x[i],x[j])
		elif j % 2 == 0:
			A_simpson[i][j] += (2*h/3)*K(x[i],x[j])
U_simpson = reflection_method(A_simpson, F)
print(math.pi)
plt.style.use('ggplot')
fig, ax = plt.subplots()
ax.scatter(x, U_right_rectangle, c = 'red')
ax.scatter(x, U_simpson, c = 'blue')
ax.plot(x, U_exact, c = 'black')
ax.set_title("Красный - метод правых прям-ов, синий - метод Симпсона", fontsize = 14)
ax.set_xlabel("Узлы", fontsize = 14)
ax.set_ylabel("Значения", fontsize = 14)
ax.tick_params(axis = 'both', labelsize = 10)
#plt.ylim([0.9,1.02])
plt.show()