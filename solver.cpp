#include <cstdio>
#include <vector>
#include <cmath>

#define forn(i, n)     for ( size_t i = 0; i < n; i++ )
#define formn(i, m, n) for ( size_t i = m; i < n; i++ )

using namespace std;

typedef vector<double> vd;
typedef vector<vd> vvd;

//# внешние параметры и постоянные
const double q = 1.6e-19;
const double F = 96485;
const double k = 1.38e-23;
const double T = 298;
const double eps0 = 8.85e-12;

//# параметры мембраны
//                Na   K   Cl
const vd c_out = {430,  20, 556};
const vd c_in =  { 50, 397,  40};
const double thickness = 1e-8;
const double phi = 0.0;
const double eps = 1;
//# начальное распределение концентраций
//start = x -> 0.0

//# параметры ионов
const vd u = {4.44e-8, 6.53e-8, 6.59e-8};
const vd z = {1, 1, -1};
vd D(3);
//# параметры разностной схемы
size_t nodes = 101;
double r = .25;  // # <.5 для сходимости схемы
double dxi = 1.0 / ( nodes - 1 );
double tau = 0;
double dtau = r * dxi * dxi;


double field_operator(vvd E, size_t j, size_t i)
{
    double a = D[j] / D[0] / dxi / dxi * dtau;
    double b = z[j] / abs(z[j]) * u[j] * thickness / 2 / D[0] / dxi * dtau;
    double sumE = 0;
    forn(k, 3)
        sumE += E[k][i];
    return E[j][i] * (1 - 2 * a) +
           E[j][i + 1] * (a - b * sumE) + E[j][i - 1] * (a + b * sumE);
}

int main(int argc, char const *argv[])
{
    forn(i, 3)
        D[i] = u[i] * k * T / q / abs(z[i]);

    auto t = []() {return tau * pow(thickness, 2) / D[0];};
    vd xi(nodes);
    for (size_t i = 0; i < nodes; i++)
        xi[i] = (double) i / (nodes - 1);
    vd x(nodes);
    for (size_t i = 0; i < xi.size(); i++)
        x[i] = thickness * xi[i];
    vd abscissa(nodes);
    for (size_t i = 0; i < xi.size(); i++)
        abscissa[i] = 10 * xi[i];

    // # устанавливаем начальное условие
    vvd E(3);
    auto dE = [](double c) { return F / eps / eps0 * c; };
    forn(i, 3)
    {
        E[i][0] = E[i][1] - dxi * thickness * dE(c_out[i]);
        E[i][E.size()-1] = E[i][E.size()-2] + dxi * thickness * dE(c_in[i]);
    }
    tau = 0;

    // # теперь считаем:
    // step = 5e-9 * D / thickness ** 2
    vvd new_E = E;

    while (tau < 1)
    {
        forn(j, 3)
        {
            formn (i, 1, nodes-1)
                new_E[j][i] = field_operator(E, j, i);
            new_E[j][0] = new_E[j][1] - dxi * thickness * dE(c_out[j]);
            new_E[j][new_E.size()-1] =
                new_E[j][new_E.size()-2] + dxi * thickness * dE(c_in[j]);
        }
        E = new_E;
        tau += dtau;
        //printf("%f\n", tau);
    }

    FILE *data;
    data = fopen("data.py", "w");
    fprintf(data, "abscissa = [ ");
    for (size_t i = 0; i < abscissa.size(); i++)
        fprintf(data, "%f, ", abscissa[i]);
    forn(j, 3)
    {
        fprintf(data, "]\nE%d = [ ", j);
        forn(i, E.size())
            fprintf(data, "%f, ", E[j][i]);
        fprintf(data, "]");
    }
    fclose(data);
    return 0;
}
