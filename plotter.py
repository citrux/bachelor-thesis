from data import *
from matplotlib import pyplot as plt


plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=14)
plt.rc("text.latex", unicode=True)
plt.rc("text.latex", preamble="\\usepackage[utf8x]{inputenc}\
                           \\usepackage[T2A]{fontenc}\
                           \\usepackage[russian]{babel}\
                           \\usepackage{amsmath}\
                           \\usepackage{pscyr}")

labels = ["$\mathrm{Na}^+$", "$\mathrm{K}^+$", "$\mathrm{Cl}^-$"]
for i, e in enumerate(E):
    plt.plot(abscissa, e, label=labels[i])

#E = [E0[i] + E1[i] + E2[i] for i in range(len(abscissa))]
#plt.plot(abscissa, E, "k-", label="ions' field")

plt.xlabel(r"$x,\ \text{нм}$")
plt.ylabel(r"$E,\ \frac{\text{В}}{\text{м}}$")
plt.grid(True)
plt.legend()
plt.savefig("plots/stat_field.pdf")
plt.cla()

for i, c in enumerate(C):
    plt.plot(abscissa, c, label=labels[i])
plt.xlabel(r"$x,\ \text{нм}$")
plt.ylabel(r"$C,\ \frac{\text{моль}}{\text{м}^3}$")
plt.grid(True)
plt.legend()
plt.savefig("plots/stat_conc.pdf")
plt.cla()
