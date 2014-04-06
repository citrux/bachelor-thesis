#include <cstdio>
#include <cstdlib>
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

const int types_used = 3;
//# параметры мембраны
//                 Na   K   Cl
const vd c_out = {430,  20, 556};
const vd c_in =  { 50, 397,  40};
const double thickness = 1e-8;
// const double phi = 0.0;
const double eps = 1;
//# начальное распределение концентраций
//start = x -> 0.0

//# параметры ионов
const vd u = {4.44e-8, 6.53e-8, 6.59e-8};
const vd z = {1, 1, -1};
vd D(types_used);
//# параметры разностной схемы
size_t nodes = 4001;
double r = 1e-2;  // # <.5 для сходимости схемы
double dxi = 1.0 / ( nodes - 1 );
double tau = 0;
double dtau = r * dxi * dxi;

vvd E(types_used);
vvd C(types_used);
vd xi(nodes);


double t()
{
	return tau * pow(thickness, 2) / D[0];
}

double field_operator(size_t type, size_t i)
{
    double a = D[type] / D[0] / dxi / dxi * dtau;
    double b = z[type] / abs(z[type]) * u[type] * thickness / 2 / D[0] / dxi * dtau;
    double sumE = 0;
    forn(k, types_used)
        sumE += E[k][i];
    return E[type][i] * (1 - 2 * a) +
           E[type][i + 1] * (a - b * sumE) + E[type][i - 1] * (a + b * sumE);
}

void border_condition()
{
    forn(type, types_used)
    {
        E[type][0] = E[type][1] - dxi * thickness * z[type] * F / eps / eps0 * c_out[type];
        E[type][nodes-1] = E[type][nodes-2] + dxi * thickness * z[type] * F / eps / eps0 * c_in[type];
        double correction = (E[type][0] + E[type][nodes-1]) / 2;
        forn(i, nodes)
            E[type][i] -= correction;
    }
}


void init()
{
    forn(i, nodes)
        xi[i] = (double) i / (nodes - 1);
	forn(type, types_used)
    {
        D[type] = u[type] * k * T / q / abs(z[type]);
        E[type] = vd(nodes);
        C[type] = vd(nodes);
    }
    border_condition();
    tau = 0;
}

void calculate_field(double time)
{
    vvd new_E = E;
    while (t() < time)
    {
        forn(type, types_used)
            formn (i, 1, nodes-1)
                new_E[type][i] = field_operator(type, i);
        E = new_E;
        border_condition();
        tau += dtau;
        printf("%.1f%% complete\n", t() / time * 100);
        printf("\033[F\033[J"); // update previous line
    }
}

void calculate_concentrations()
{
    forn(type, types_used)
    {
        C[type][0] = c_out[type];
        C[type][nodes-1] = c_in[type];
        formn(i, 1, nodes - 1)
            C[type][i] = eps * eps0 / z[type] / F * (E[type][i+1] - E[type][i-1]) / 2 / dxi / thickness;
    }
}

int main(int argc, char const *argv[])
{
    // устанавливаем начальное условие
    init();

    // теперь считаем:
    double t = 1e-9;
    if (argc > 1)
        t = atof(argv[1]);
    calculate_field(t);
    calculate_concentrations();

    printf("write results to data.py… ");
    FILE *data;
    data = fopen("data.py", "w");
    fprintf(data, "abscissa = [ ");
    forn(i, nodes)
        fprintf(data, "%f, ", xi[i] * 10);
    fprintf(data, "]\n");

    fprintf(data, "E = [\n");
    forn(type, types_used)
    {
        fprintf(data, "[");
        forn(i, nodes)
            fprintf(data, "%f, ", E[type][i]);
        fprintf(data, "],\n");
    }
    fprintf(data, "]\n");

    fprintf(data, "C = [\n");
    forn(type, types_used)
    {
        fprintf(data, "[");
        forn(i, nodes)
            fprintf(data, "%f, ", C[type][i]);
        fprintf(data, "],\n");
    }
    fprintf(data, "]\n");
    fclose(data);
    printf("done!\n");

    return 0;
}
