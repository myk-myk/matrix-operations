import numpy as np


def seidel(A, B, delta):
    n = len(A)
    x = np.zeros(n)  # zero vector

    converge, it = False, 0
    while not converge:
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (B[i] - s1 - s2) / A[i][i]
        print(it, ": ", x)
        converge = np.sqrt(sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= delta
        x = x_new
        it += 1


if __name__ == '__main__':
    A = [[12, 1, -3, -1],
         [4, -12, 1, -2],
         [-1, -1, 14, 1],
         [-3, 5, -4, -11]]
    B = [-3, -10, 8, -5]
    delta = 1e-11
    seidel(A, B, delta)
