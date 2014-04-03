PRODUCTION = True

from data import *
from matplotlib import pyplot as plt

if PRODUCTION:
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif", size=14)
    plt.rc("text.latex", unicode=True)
    plt.rc("text.latex", preamble="\\usepackage[utf8x]{inputenc}\
                               \\usepackage[T2A]{fontenc}\
                               \\usepackage[russian]{babel}\
                               \\usepackage{amsmath}\
                               \\usepackage{pscyr}")


for E, l in [(E0, "Na"), (E1, "K"), (E2, "Cl")]:
    plt.plot(abscissa, E, label=l)

#E = [E0[i] + E1[i] + E2[i] for i in range(len(abscissa))]
#plt.plot(abscissa, E, "k-", label="ions' field")

plt.xlabel(r"$x,\ \text{нм}$")
plt.ylabel(r"$E,\ 10^9\frac{\text{В}}{\text{м}}$")
plt.grid(True)
plt.legend()
plt.savefig("plots/stat_field.pdf")
plt.cla()

for C, l in [(C0, "Na"), (C1, "K"), (C2, "Cl")]:
    plt.plot(abscissa, C, label=l)
plt.xlabel(r"$x,\ \text{нм}$")
plt.ylabel(r"$C,\ \frac{\text{моль}}{\text{м}^3}$")
plt.grid(True)
plt.legend()
plt.savefig("plots/stat_conc.pdf")
plt.cla()
