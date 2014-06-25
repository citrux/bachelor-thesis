from matplotlib import pyplot as plt
from numpy import linspace, e, array

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=10)
plt.rc("text.latex", unicode=True)
plt.rc("text.latex", preamble="\\usepackage[utf8x]{inputenc}\
                               \\usepackage[T2A]{fontenc}\
                               \\usepackage[russian]{babel}\
                               \\usepackage{amsmath}\
                               \\usepackage{pscyr}")
plt.rc("figure", figsize=(4, 4))

# внешние параметры и постоянные
q = 1.6e-19
k = 1.38e-23
T = 298
F = 96485
eps0 = 8.85e-12
eps = 1

# параметры мембраны
c_out = 430
c_in = 50
thickness = 1e-8
phi = 0.1
# начальное распределение концентраций
start = lambda x: 0.0 * x

# параметры ионов
u = 4.4e-8
z = 1
D = u * k * T / q / abs(z)

# параметры разностной схемы
nodes = 101
r = .25  # <.5 для сходимости схемы
xi = linspace(0, 1, nodes)
dxi = 1 / (nodes - 1)
tau = 0
dtau = r * dxi * dxi

t = lambda tau: tau * thickness ** 2 / D
x = thickness * xi
abscissa = 1e9 * x


# теперь очередь разностного оператора
#
# j+1 |           i
#  j  | …   i-1   i   i+1   …
#
# он рассчитывает значение
def linear_field_operator(c, i):
     w = q * abs(z) * phi / k / T
     s = .5 * w * dxi
     return c[i] * (1 - 2 * r) + c[i + 1] * r * (1 - s) + c[i - 1] * r * (1 + s)


def linear_conc_operator(c, i):
     w = q * abs(z) * phi / k / T
     s = .5 * w * dxi
     b = w * dtau
     g = lambda i: 1 / (xi[i] - c_out / (c_out - c_in))
     return (c[i] * (1 - 2 * r + b * g(i)**2) +
             c[i + 1] * r * (1 - s * g(i)) +
             c[i - 1] * r * (1 + s * g(i)))


def linear_field_stat_conc(x):
     a = z * q * phi / k / T
     num = c_out * e**a - c_in - (c_out - c_in) * e**(a * x / thickness)
     return num / (e**a - 1)


def linear_conc_stat_conc(x):
     return c_out - (c_out - c_in) * x / thickness


def main(name):
     if name == "linear_conc":
         stat_conc = linear_conc_stat_conc
         operator = linear_conc_operator

     if name == "linear_field":
         stat_conc = linear_field_stat_conc
         operator = linear_field_operator

     # устанавливаем начальное условие
     c = start(linspace(0, thickness, nodes))
     c[0] = c_out
     c[-1] = c_in
     tau = 0

     # теперь считаем:
     step = 5e-9 * D / thickness ** 2
     new_c = c

     for j in range(1, 5):
         while tau < step * j:
             for i in range(1, nodes-1):
                 new_c[i] = operator(c, i)
             c = new_c
             tau += dtau

         plt.plot(abscissa, c, color="black", linewidth=.5, linestyle="-")
         if j == 4:
             plt.annotate("\\(t = %.0f~\\text{нс}\\)" % (t(tau) * 1e9),
                          xy=(abscissa[50], c[50]),
                          xytext=(abscissa[50], c[50] + 30),
                          arrowprops = {"arrowstyle": "->",
                                        "color": "black",
                                        "shrinkB": .01})

     s_conc = stat_conc(x)
     plt.plot(abscissa, s_conc, color="red")
     plt.annotate("\\(t \\to \\infty \\)", xy=(abscissa[70], s_conc[70]),
                  xytext=(abscissa[70], s_conc[70] + 30),
                  arrowprops={"arrowstyle": "->",
                              "color": "red",
                              "shrinkB": .01})

     plt.xlabel("Мембрана, нм")
     plt.ylabel("Концентрация иона \\(\\mathrm{Na}^+ \\),\
                 моль/\\(\\text{м}^3\\)")
     plt.title("Распределение ионов \\(\\mathrm{Na}^+ \\) внутри мембраны")
     plt.savefig("plots/" + name + ".pdf")
     plt.cla()

main("linear_field")
main("linear_conc")
