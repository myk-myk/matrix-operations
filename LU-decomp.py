import numpy as np
import scipy.linalg


def LU_decomposition(A):
    n = len(A[0])
    L = np.zeros([n, n])
    U = np.zeros([n, n])
    for i in range(n):
        L[i][i] = 1
        if i == 0:
            U[0][0] = A[0][0]
            for j in range(1, n):
                U[0][j] = A[0][j]
                L[j][0] = A[j][0] / U[0][0]
        else:
            for j in range(i, n):
                temp = 0
                for k in range(0, i):
                    temp = temp + L[i][k] * U[k][j]
                U[i][j] = A[i][j] - temp
            for j in range(i + 1, n):
                temp = 0
                for k in range(0, i):
                    temp = temp + L[j][k] * U[k][i]
                L[j][i] = (A[j][i] - temp) / U[i][i]
    print("L = ", L)
    print()
    print("U = ", U)
    gauss(U, gauss(L, BB))


def gauss(AA, BB):
    n = len(BB)
    A = []
    for aa in AA:
        A.append(aa[:])
    B = BB[:]
    # прямий хід алгоритму
    for k in range(n - 1):
        # визначення найбільшого елемента
        value = 0
        index = 0
        for i in range(k, n):
            if abs(A[i][k]) > value:
                value = abs(A[i][k])
                index = i
        # перестановка рядків
        A[k], A[index] = A[index], A[k]
        B[k], B[index] = B[index], B[k]
        # перевірка умови
        if A[k][k] == 0:
            print("Error: A[k][k] == 0")
            return
        # побудова трикутної матриці
        for i in range(k + 1, n):
            m = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - A[k][j] * m
            B[i] = B[i] - B[k] * m
    # зворотній хід алгоритму
    X = B[:]
    X[n - 1] = B[n - 1] / A[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, n):
            X[i] = X[i] - A[i][j] * X[j]
        X[i] = X[i] / A[i][i]
    print("X = ", X)
    return X


if __name__ == '__main__':
    B = np.array([[-11, 1, 6, 3],
                  [5, 14, 3, 5],
                  [1, -1, 9, 4],
                  [-4, 4, 2, 12]])
    BB = [0, -8, 1, 1]
    LU_decomposition(B)

    P, L, U = scipy.linalg.lu(B)
    print("L = ", L)
    print("U = ", U)
