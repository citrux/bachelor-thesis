using PyPlot

rc("text", usetex=true)
rc("font", family="serif", size=10)
rc("text.latex", unicode=true)
rc("text.latex", preamble="\\usepackage[russian]{babel}\\usepackage{amsmath}")
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
x = linspace(0, 10, nodes)
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
    annotate(string("t = ", step * i),
        xy = (x[50], c[50]))
    t = 0
end

xlabel("Мембрана, нм")
ylabel("Концентрация иона \\( \\mathrm{Na}^+ \\), моль/\\(\\text{м}^3\\)")
savefig("1.png")

