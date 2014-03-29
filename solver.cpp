#include <cstdio>
#include <vector>

typedef std::vector<double> vd;

double field_operator(vd E, size_t i)
{
    double a = u * phi / 2 / D / dxi * dtau;
    double b = u * thickness / 2 / D / dxi * dtau;
    return E[i] * (1 - 2 * r - b * (E[i+1] - E[i-1])) +\
           E[i + 1] * (r - a) + E[i - 1] * (r + a);
}

int main(int argc, char const *argv[])
{
    // # устанавливаем начальное условие
    vd E(nodes);
    double dE = [](c) { return F / eps / eps0 * c; };
    E[0] = E[1] - dxi * thickness * dE(c_out);
    E[-1] = E[-2] + dxi * thickness * dE(c_in);
    double tau = 0;

    // # теперь считаем:
    // step = 5e-9 * D / thickness ** 2
    vd new_E = E;

    while (tau < 10)
    {
        for (int i = 1; i < nodes-1; ++i)
            new_E[i] = field_operator(E, i);
        new_E[0] = new_E[1] - dxi * thickness * dE(c_out);
        new_E[-1] = new_E[-2] + dxi * thickness * dE(c_in);
        E = new_E;
        tau += dtau;
    }
    return 0;
}
