#include <cstdio>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> vd;

//# внешние параметры и постоянные
double q = 1.6e-19;
double F = 96485;
double k = 1.38e-23;
double T = 298;
double eps0 = 8.85e-12;

//# параметры мембраны
double c_out = 430;
double c_in = 50;
double thickness = 1e-8;
double phi = 0.0;
double eps = 1;
//# начальное распределение концентраций
//start = x -> 0.0

//# параметры ионов
double u = 4.4e-8;
double z = 1;
double D = u * k * T / q / abs(z);

//# параметры разностной схемы
size_t nodes = 101;
double r = .25;  // # <.5 для сходимости схемы
double dxi = 1.0 / ( nodes - 1 );
double tau = 0;
double dtau = r * dxi * dxi;


double field_operator(vd E, size_t i)
{
    double a = u * phi / 2 / D / dxi * dtau;
    double b = u * thickness / 2 / D / dxi * dtau;
    return E[i] * (1 - 2 * r - b * (E[i+1] - E[i-1])) +\
           E[i + 1] * (r - a) + E[i - 1] * (r + a);
}

int main(int argc, char const *argv[])
{
    auto t = []() {return tau * pow(thickness, 2) / D;};
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
    vd E(nodes);
    auto dE = [](double c) { return F / eps / eps0 * c; };
    E[0] = E[1] - dxi * thickness * dE(c_out);
    E[E.size()-1] = E[E.size()-2] + dxi * thickness * dE(c_in);
    tau = 0;

    // # теперь считаем:
    // step = 5e-9 * D / thickness ** 2
    vd new_E = E;

    while (tau < 1)
    {
        for (int i = 1; i < nodes-1; ++i)
            new_E[i] = field_operator(E, i);
        new_E[0] = new_E[1] - dxi * thickness * dE(c_out);
        new_E[new_E.size()-1] =
            new_E[new_E.size()-2] + dxi * thickness * dE(c_in);
        E = new_E;
        tau += dtau;
        //printf("%f\n", tau);
    }

    FILE *data;
    data = fopen("data.py", "w");
    fprintf(data, "abscissa = [ ");
    for (size_t i = 0; i < abscissa.size(); i++)
        fprintf(data, "%f, ", abscissa[i]);
    fprintf(data, "]\nE = [ ");
    for (size_t i = 0; i < E.size(); i++)
        fprintf(data, "%f, ", E[i]);
    fprintf(data, "]");
    fclose(data);
    return 0;
}
