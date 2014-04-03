from data import abscissa, E0, E1, E2
from matplotlib import pyplot as plt

for E, l in [(E0, "Na"), (E1, "K"), (E2, "Cl")]:
    plt.plot(abscissa, E, label=l)
plt.xlabel("x, nm")
plt.ylabel("E, GV/m")
plt.grid(True)
plt.legend()
plt.savefig("field.png")
plt.cla()

F = 96485 * 1e-9
eps = 1
eps0 = 8.85e-12

dx = (abscissa[1] - abscissa[0]) * 1e-9
#print(dx)
for E, l in [(E0, "Na"), (E1, "K"), (E2, "Cl")]:
    C = [0] * len(abscissa)
    C[0] = eps * eps0 / F * (E[1] - E[0]) / dx
    C[-1] = eps * eps0 / F * (E[-1] - E[-2]) / dx
    for i in range(1, len(abscissa) - 1):
        C[i] = eps * eps0 / F * (E[i+1] - E[i-1]) / 2 / dx
    plt.plot(abscissa, C, label=l)
plt.xlabel("x, nm")
plt.ylabel("c, mol/m^3")
plt.grid(True)
plt.legend()
plt.savefig("conc.png")

phi = -(sum(E) - .5 * (E[0] + E[-1])) * dx * 1e9
print(phi)
