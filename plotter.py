from data import abscissa, E
from matplotlib import pyplot as plt

plt.plot(abscissa, E)
plt.grid(True)
plt.savefig("field.png")
plt.cla()

F = 96485
eps = 1
eps0 = 8.85e-12

dx = (abscissa[1] - abscissa[0]) * 1e-9
print(dx)
C = [0] * len(abscissa)
C[0] = eps * eps0 / F * (E[1] - E[0]) / dx
C[-1] = eps * eps0 / F * (E[-1] - E[-2]) / dx
for i in range(1, len(abscissa) - 1):
    C[i] = eps * eps0 / F * (E[i+1] - E[i-1]) / 2 / dx
plt.plot(abscissa, C)
plt.grid(True)
plt.savefig("conc.png")

phi = (sum(E) - .5 * (E[0] + E[-1])) * dx
print(phi)
