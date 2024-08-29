import math
import numpy as np
import matplotlib.pyplot as plt

a, b = 1, 1.2

eps = 10E-6

n = 20

h = 1/n

p = int(3*n*n/4+3*n/2)

def f(x,y):
	return 0.2*np.exp(x)*np.cos(y)

def phi(x,y):
	return np.exp(x)*np.cos(y)

def w(i,j):
	if j <= int(n/2):
		return int((n + 1 + j)*j/2 + i - n/2 + j)
	
	else:
		return int((3*n*n+10*n)/8 + (12*j*n - 14*n - 4*j*j + 12*j - 5*n*n)/8 + i)

def max_sob(A, n):
	x = np.ones(n)
	mlam_prev = 0
	x = x/np.linalg.norm(x)
	y = np.dot(A,x)
	mlam = np.dot(y,x)
	x = y/np.linalg.norm(y)
	while(abs(mlam - mlam_prev) > eps):
		mlam_prev = mlam
		x = x/np.linalg.norm(x)
		y = np.dot(A,x)
		mlam = np.dot(y,x)
		x = y/np.linalg.norm(y)
	return mlam

def min_sob(A, n, max_lam):
	x = np.ones(n)
	mlam_prev = 0
	B = max_lam*np.eye(n) - A
	x = x/np.linalg.norm(x)
	y = np.dot(B,x)
	mlam = np.linalg.norm(np.dot(A,y))/np.linalg.norm(y)
	x = y/np.linalg.norm(y)
	while(abs(mlam - mlam_prev) > eps):
		mlam_prev = mlam
		x = x/np.linalg.norm(x)
		y = np.dot(B,x)
		mlam = np.linalg.norm(np.dot(A,y))/np.linalg.norm(y)
		x = y/np.linalg.norm(y)
	return mlam

def upper_relaxation(A, b, N, cond):
	x = np.zeros(N)
	k = 0
	while True:
		x_cur = np.copy(x)
		k += 1
		for i in range(N):
			suma = 0
			for j in range(N):
				if j < i:
					suma += A[i][j]*x_cur[j]
				elif j > i:
					suma += A[i][j]*x[j]
			x_cur[i] = (1-cond)*x[i]+cond*(b[i] - suma)/A[i][i]
			"""
		if np.linalg.norm(np.dot(A,x_cur) - b) < eps:
			print(np.linalg.norm(np.dot(A,x_cur) - b))
			break
			"""
		if np.linalg.norm(x_cur - x) < eps:
			print(np.linalg.norm(x_cur - x))
			break
		x = x_cur
	return x_cur 

n_i_j = np.zeros((p+1, 3))
boundary = []
area = []

m = 0
for j in range(n+1):
	for i in range(n+1):
		if i + j >= int(n/2) and i + j <= int(3*n/2):
			n_i_j[m][0] = m
			n_i_j[m][1] = i
			n_i_j[m][2] = j
			if i + j == int(n/2) or i + j == int(3*n/2) or i == 0 or j == 0 or i == n or j == n:
				boundary.append(m)
			else: area.append(m)
			m += 1

k = len(area)

u_exact = np.zeros(k)
for i in range(k):
	u_exact[i] = phi(n_i_j[area[i]][1]*h,n_i_j[area[i]][2]*h)

matrix = np.zeros((k,k))
F = np.zeros(k)

for i in range(k):
	l = area[i]
	matrix[i][i] = 2*(a+b)
	F[i] = h*h*f(n_i_j[l][1]*h,n_i_j[l][2]*h)

	if l - 1 in boundary:

		F[i] += a*phi(n_i_j[l - 1][1]*h,n_i_j[l][2]*h)

	else:
		matrix[i][i - 1] = -a

	if l + 1 in boundary:

		F[i] += a*phi(n_i_j[l + 1][1]*h,n_i_j[l][2]*h)

	else:
		matrix[i][i + 1] = -a

	if w(n_i_j[l][1],n_i_j[l][2]-1) in boundary:

		F[i] += b*phi(n_i_j[l][1]*h,n_i_j[l - 1][2]*h)

	else:
		matrix[i][area.index(w(n_i_j[l][1],n_i_j[l][2]-1))] = -b

	if w(n_i_j[l][1],n_i_j[l][2]+1) in boundary:

		F[i] += b*phi(n_i_j[l][1]*h,n_i_j[l + 1][2]*h)

	else:
		matrix[i][area.index(w(n_i_j[l][1],n_i_j[l][2]+1))] = -b


max_lambda = max_sob(matrix, k)

min_lambda = min_sob(matrix, k, max_lambda)

cond_A = max_lambda/min_lambda

tau_opt = 2 - 4/math.sqrt(cond_A)

u_numerical = np.zeros(p+1)

u_numerical = upper_relaxation(matrix, F, k, tau_opt)

print(np.linalg.norm(u_numerical - u_exact, ord = np.inf))
