import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def table(x_, y):
    quotient = [[0] * len(x_) for _ in range(len(x_))]
    for n_ in range(len(x_)):
        quotient[n_][0] = y[n_]
    for i in range(1, len(x_)):
        for j in range(i, len(x_)):
            quotient[j][i] = (quotient[j][i - 1] - quotient[j - 1][i - 1]) / (x_[j] - x_[j - i])
    return quotient


def get_corner(result):
    link = []
    for i in range(len(result)):
        link.append(result[i][i])
    return link


def newton(data_set, x_p, x_7):
    result = data_set[0]
    for i in range(1, len(data_set)):
        p = data_set[i]
        for j in range(i):
            p *= (x_p - x_7[j])
        result += p
    return result


def draw_picture(x_list, y_list, node):
    plt.title("newton")
    plt.xlabel("x")
    plt.ylabel("y")
    for i in range(len(x_list)):
        plt.scatter(x_list[i], y_list[i], color="purple", linewidths=2)
    plt.scatter(node[0], node[1], color="blue", linewidth=2)
    plt.show()


def fun(x):
    return pow(math.e, x) * math.sin(x)


if __name__ == '__main__':
    x = 11.63
    x_1 = [9, 9.75, 10.5, 11.25, 12]
    y_1 = [3339.43, -5481.1, -31946.59, -74405, -87329.8]
    middle = table(x_1, y_1)
    n = get_corner(middle)
    newton = newton(n, x, x_1)
    print("Real value: ", fun(x))
    print("Newton's interpolation: {}".format(newton))
    draw_picture(x_1, y_1, (x, newton))
