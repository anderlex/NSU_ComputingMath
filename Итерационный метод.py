from math import *
import numpy as np

eps = 10e-8
def eigenvalues_rotate(M):
	N = M.shape[0]
	A = M.copy()
	k = 0
	while True:
		k += 1
		m = 0
		T = np.eye((N))
		m_i = 0
		m_j = 1
		for i in range(N):
			for j in range(i + 1):
				if i == j:
					continue
				else:
					if abs(A[i][j]) > m:
						m = abs(A[i][j])
						m_i = i
						m_j = j
		if A[m_i][m_i] == A[m_j][m_j]:
			c = 1/sqrt(2)
			s = 1/sqrt(2)
		else:
			phi = atan(2*A[m_i][m_j]/(A[m_i][m_i] - A[m_j][m_j]))/2
			c = cos(phi)
			s = sin(phi)
		T[m_i][m_i], T[m_i][m_j], T[m_j][m_i], T[m_j][m_j] = c, -s, s, c
		A = np.dot(T.transpose(),np.dot(A, T))

		non_diag_sum = 0
		diag_min = abs(A[0][0])
		for i in range(1, N):
			if abs(A[i][i]) < diag_min:
				diag_min = abs(A[i][i])

		for i in range(1, N):
			for j in range(i):
				non_diag_sum += A[i][j]**2
		non_diag_sum = sqrt(non_diag_sum)

		#print(diag_min)
		#print(non_diag_sum)
		if non_diag_sum < eps*diag_min:
			break
	#print(k)
	eigen = np.zeros(N)
	for i in range(N):
		eigen[i] = A[i][i]
	return eigen

def richardson_method(A, b, N, eigenvalues):
	x_cur = np.ones(N)
	E = np.eye(N)
	m = min(eigenvalues)
	M = max(eigenvalues)
	tau_optimal = 2/(M + m)
	k = 0
	while True:
		k += 1
		x = x_cur
		x_cur = np.dot((E - tau_optimal*A), x) + tau_optimal*b
		if np.linalg.norm(np.dot(A,x_cur) - b, ord = 2) < eps * m:
			print(np.linalg.norm(np.dot(A,x_cur) - b, ord = 2))
			break
		#print(np.linalg.norm(np.dot(A,x_cur) - b, ord = 2))
	print(k)
	return x_cur

def ziedel_method(A, b, N, eigenvalues):
	m = min(eigenvalues)
	print(m)
	x = np.ones(N)
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
			x_cur[i] = (b[i] - suma)/A[i][i]
		if np.linalg.norm(np.dot(A,x_cur) - b, ord = 2) < eps * m:
			print(np.linalg.norm(np.dot(A,x_cur) - b, ord = 2))
			break
		#print(np.linalg.norm(np.dot(A,x_cur) - b, ord = 2))
		x = x_cur
	print(k)
	return x_cur 

INPUT = 'test3.dat'
Values = np.loadtxt(INPUT)
N = int(Values[0])
A = np.eye((N))
b = np.zeros(N)
for i in range(N):
	for j in range(N):
		A[i][j] = Values[i*N + j + 1]
for i in range(N):
	b[i] = Values[N*N + i + 1]


eigenvalues = eigenvalues_rotate(A)
#print(eigenvalues)

x_richardson = richardson_method(A, b, N, eigenvalues)
print(x_richardson)

x_ziedel = ziedel_method(A, b, N, eigenvalues)
print(x_ziedel)
