using PyPlot

rc("text", usetex=true)
rc("font", family="serif", size=10)
rc("text.latex", unicode=true)
rc("text.latex", preamble="\\usepackage[utf8x]{inputenc}
                           \\usepackage[T2A]{fontenc}
                           \\usepackage[russian]{babel}
                           \\usepackage{amsmath}
                           \\usepackage{pscyr}")
rc("figure", figsize = (4, 4))

c_out = 430
c_in = 50
thickness = 1e-8

nodes = 101
dx = thickness / (nodes - 1)

# устанавливаем начальное условие
c = linspace(0, 0, nodes)
c[1] = c_out
c[nodes] = c_in

x = linspace(0, 10, nodes)
# теперь очередь разностного оператора
#
# j+1 |     y'
#  j  | … x y z …
#
# он рассчитывает значение y'
r =  .25 # <.5 для сходимости схемы
D = 1.3e-9
q = 0.44
dt = r * dx * dx / D

function sub_operator(x, y, z)
    y * (1 - 2 * r) + z * (r - q * dt / 2 / dx) + x * (r + q * dt / 2 / dx)
end

# теперь считаем:
t = 0
step = 5e-9
new_c = c

for i in [1:4]
    while t < step
        for j in [2:nodes-1]
            new_c[j] = sub_operator(c[j-1], c[j], c[j+1])
        end
        c = new_c
        t += dt
    end
    plot(x, c, color="black", linewidth=.5, linestyle="-")
    if (i == 4)
        annotate(@sprintf("\\(t = %.0f~\\text{нс}\\)", step * i * 1e9),
            xy = (x[50], c[50]), xytext = (x[50], c[50] + 30),
            arrowprops = {"arrowstyle" => "->",
                          "color" => "black",
                          "shrinkB" => .01})
    end
    t = 0
end

function stat_conc(x)
    k = q * thickness / D
    num = c_out * e^k - c_in - (c_out - c_in) * e^(k * x / thickness)
    num / (e^k - 1)
end

s_conc = map(stat_conc, x / 10 * thickness)
plot(x, s_conc, color="red")
annotate("\\(t \\to \\infty \\)", xy = (x[70], s_conc[70]),
            xytext = (x[70], s_conc[70] + 30),
            arrowprops = {"arrowstyle" => "->",
                          "color" => "red",
                          "shrinkB" => .01})

xlabel("Мембрана, нм")
ylabel("Концентрация иона \\( \\mathrm{Na}^+ \\), моль/\\(\\text{м}^3\\)")
title("Распределение концентрации ионов \\( \\mathrm{Na}^+ \\)\\
       внутри мембраны")
savefig("plots/linear_field.pdf")

