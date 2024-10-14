#include <stdio.h>
#include <math.h>


//Função para calcular circuito superarmortecido(sigma > omega)
void calc_superamortecido(double sigma, double omega, double resistencia, double indutancia, double capacitancia, double *A1, double *A2, double tensaoCapacitor, double correnteIndutor){
    //calcular s1 e s2
   double s1= -sigma + sqrt(pow(sigma, 2) - pow(omega, 2));
   double s2= -sigma - sqrt(pow(sigma, 2) - pow(omega, 2));
    *A1 = (correnteIndutor - capacitancia*tensaoCapacitor*s2);
}
void calc_criticamente(double sigma){
    double s1 = -sigma;
    double s2 = -sigma;
}

void calc_subamortecido(double sigma, double omega, double omegaD,double A1, double A2, double B1, double B2){
    B1 = A1 + A2;
    B2 = A1 - A2;

}

int main(){

    double resistencia, indutancia, capacitancia;
    double tensaoCapacitor, correnteIndutor;
    double sigma = 0;
    double omega = 0;
    double omegaD = 0;
    double tm, vtm;
    double A1, A2, B1, B2;

    sigma = resistencia /(2 * indutancia);
    omega = 1 / sqrt(indutancia * capacitancia);
    omegaD = sqrt(pow(omega, 2) - pow(sigma, 2)) ;
}
