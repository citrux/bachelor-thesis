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

# внешние параметры и постоянные
q = 1.6e-19
k = 1.38e-23
T = 298

# параметры мембраны
c_out = 430
c_in = 50
thickness = 1e-8
phi = 0.1
# начальное распределение концентраций
start = x -> 0.0

# параметры ионов
u = 4.4e-8
z = 1
D = u * k * T / q / abs(z)

# параметры разностной схемы
nodes = 101
r = .25 # <.5 для сходимости схемы
xi = linspace(0, 1, nodes)
dxi = 1 / ( nodes - 1 )
tau = 0
dtau = r * dxi * dxi

t = tau -> tau * thickness^2 / D
x = thickness * xi
abscissa = 1e9 * x


# теперь очередь разностного оператора
#
# j+1 |           i
#  j  | …   i-1   i   i+1   …
#
# он рассчитывает значение
function linear_field_operator( c, i )
    w = q * abs(z) * phi / k / T
    s = .5 * w * dxi
    c[i] * (1 - 2 * r) + c[i + 1] * r * (1 - s) + c[i - 1] * r * (1 + s)
end

function linear_conc_operator( c, i )
    w = q * abs(z) * phi / k / T
    s = .5 * w * dxi
    b = w * dtau
    g(i) = 1 / (xi[i] - c_out / (c_out - c_in))
    (c[i] * (1 - 2 * r + b * g(i)^2) + c[i + 1] * r * (1 - s * g(i)) +
    c[i - 1] * r * (1 + s * g(i)))
end


function linear_field_stat_conc(x)
    a = z * q * phi / k / T
    num = c_out * e^a - c_in - (c_out - c_in) * e^(a * x / thickness)
    num / (e^a - 1)
end

function linear_conc_stat_conc(x)
    c_out - (c_out - c_in) * x / thickness
end

function main( name )
    if name == "linear_conc"
        stat_conc = linear_conc_stat_conc
        operator = linear_conc_operator
    end
    if name == "linear_field"
        stat_conc = linear_field_stat_conc
        operator = linear_field_operator
    end
    # устанавливаем начальное условие
    c = map( start, linspace( 0, thickness, nodes ) )
    c[1] = c_out
    c[nodes] = c_in
    tau = 0


    # теперь считаем:
    step = 5e-9 * D / thickness^2
    new_c = c

    for j in [1:4]
        while tau < step * j
            for i in [2:nodes-1]
                new_c[i] = operator( c, i )
            end
            c = new_c
            tau += dtau
        end
        plot(abscissa, c, color="black", linewidth=.5, linestyle="-")
        if j == 4
            annotate(@sprintf("\\(t = %.0f~\\text{нс}\\)", t(tau) * 1e9),
                xy = (abscissa[50], c[50]), xytext = (abscissa[50], c[50] + 30),
                arrowprops = {"arrowstyle" => "->",
                              "color" => "black",
                              "shrinkB" => .01})
        end
    end

    s_conc = map( stat_conc, x )
    plot(abscissa, s_conc, color="red")
    annotate("\\(t \\to \\infty \\)", xy = (abscissa[70], s_conc[70]),
                xytext = (abscissa[70], s_conc[70] + 30),
                arrowprops = {"arrowstyle" => "->",
                              "color" => "red",
                              "shrinkB" => .01})

    xlabel("Мембрана, нм")
    ylabel("Концентрация иона \\( \\mathrm{Na}^+ \\), моль/\\(\\text{м}^3\\)")
    title("Распределение ионов \\( \\mathrm{Na}^+ \\) внутри мембраны")
    savefig("plots/" * name * ".pdf")
    cla()
end

main("linear_field")
main("linear_conc")
