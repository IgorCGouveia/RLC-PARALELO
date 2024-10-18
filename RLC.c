#include <stdio.h>
#include <math.h>

#define m  pow(10, -3);
#define Mi  pow(10, -6);
#define k  pow(10, 3);
#define Me  pow(10, 6);

double converter_unidade(double valor, char unidade) {
    switch (unidade) {
        case 'm': // mili (10^-3)
            return valor * pow(10, -3);
        case 'u': // micro (10^-6)
            return valor * pow(10, -6);
        case 'k': // kilo (10^3)
            return valor * pow(10, 3);
        case 'M': // Mega (10^6)
            return valor * pow(10, 6);
        default: // Caso o usuário não informe uma unidade válida, retornamos o valor original
            return valor;
    }}
//Função para calcular circuito superarmortecido(sigma > omega)
void calc_superamortecido(double R, double L, double C ,  double VC0 , double IL0, double sigma, double omega){

    //calcular s1 e s2
   double s1 = -sigma + sqrt(pow(sigma, 2) - pow(omega, 2));
   double s2 = -sigma - sqrt(pow(sigma, 2) - pow(omega, 2));
   
   double Corrente_R0 = VC0 / R;
   double corrente_C0 = IL0 - Corrente_R0;
   double var = corrente_C0  / C;

  
   
    double A1 = (var - VC0*s2)/(s1 - s2);
    double A2 = VC0 - A1;
    
    double tm = (log(fabs(s1*A1)) - log(fabs(s2*A2)))/(s2 - s1);
    double vtm = A1 * exp(s1 * tm) + A2 * exp(s2 * tm);
    printf(" Valor de A1 = %.3f \n Valor de A2 = %.3f\n sigma = %.3f(s^-1)\n omega = %.3f(rad/s)\n s1 = %.3f\n s2 = %.3f\n",A1,A2,sigma,omega,s1,s2);
    printf("tm= %.10lf(seg)\n VTM = %.10lf(V)\n", tm, vtm);
}



void calc_criticamente(double R, double L, double C ,  double VC0 , double IL0,double sigma, double omega){

    // Calcular s1 e s2 (s1 = s2 = -sigma)
    double s1 = -sigma;

    double Corrente_R0 = VC0 / R;
    double corrente_C0 = IL0 - Corrente_R0;
    double var = corrente_C0 / C;

    // Calcular A1 e A2
    double A1 = VC0; // A1 é o coeficiente constante
    double A2 = var + (A1 * s1); // A2 é o coeficiente que multiplica o termo t

    // Calcular o tempo onde a tensão é máxima
    double tm = -(A2 + A1 * s1) / (A2 * s1);
    double vtm = (A1 + A2 * tm) * exp(s1 * tm);

    printf("Valor de A1 = %.3f \nValor de A2 = %.3f\nSigma = %.3f(s^-1)\nOmega = %.3f(rad/s)\ns1 = %.3f\nTempo tm (tensão máxima) = %.10lf(seg)\nVtm = %.10lf(V)\n", A1, A2, sigma, omega, s1, tm , vtm);
}



// Função para calcular o circuito subamortecido (sigma < omega)
void calc_subamortecido(double R, double L, double C, double VC0, double IL0,double sigma, double omega) {
    
    // Calcular a frequência amortecida omegaD
    double omegaD = sqrt(pow(omega, 2) - pow(sigma, 2));
    
    // Calcular as constantes A1 e A2
    double Corrente_R0 = VC0 / R;
    double corrente_C0 = IL0 - Corrente_R0;
    double var = corrente_C0 / C;

    double B1 = VC0;
    double B2 = (var + sigma * VC0) / omegaD;

    // Calcular o tempo onde a tensão é máxima
    double tm = (atan((B2 * omegaD - B1 * sigma) / (B1 * omegaD + B2 * sigma)))/omegaD;
    double vtm = B1 * exp(-sigma * tm) * cos(omegaD * tm) + B2 * exp(-sigma * tm) * sin(omegaD * tm);
    
    printf("Valor de B1 = %.3f \n Valor de B2 = %.3f\n sigma = %.3f(s^-1)\n omega = %.3f(rad/s)\nomegaD = %.3f(rad/s)\n tm = %.10lf(seg)\nVtm = %.10lf(V)\n", B1, B2, sigma, omega, omegaD, tm, vtm);
}

int main(){

    double R, L, C;
    double VC0, IL0;
    char unidade_R = '\0', unidade_L = '\0', unidade_C = '\0', unidade_VC0 = '\0', unidade_IL0 = '\0';
    int cont = 1;
    while(cont == 1){
    printf("Digite o valor da Resistência: ");
    scanf("%lf", &R);
    printf("Digite a unidade da Resistência (m para milli, u para micro, k para kilo, M para Mega ou pressione Enter para nenhuma): ");
    getchar(); // Limpar o buffer do teclado
    unidade_R = getchar();

    printf("Digite o valor da Capacitância: ");
    scanf("%lf", &C);
    printf("Digite a unidade da Capacitância (m para milli, u para micro, k para kilo, M para Mega ou pressione Enter para nenhuma): ");
    getchar(); // Limpar o buffer do teclado
    unidade_C = getchar();

    printf("Digite o valor da Indutância: ");
    scanf("%lf", &L);
    printf("Digite a unidade da Indutância (m para milli, u para micro, k para kilo, M para Mega ou pressione Enter para nenhuma): ");
    getchar(); // Limpar o buffer do teclado
    unidade_L = getchar();

    printf("Digite o valor da Tensão Inicial (VC0): ");
    scanf("%lf", &VC0);
    printf("Digite a unidade da Tensão (m para milli, u para micro, k para kilo, M para Mega ou pressione Enter para nenhuma): ");
    getchar(); // Limpar o buffer do teclado
    unidade_VC0 = getchar();

    printf("Digite o valor da Corrente Inicial (IL0): ");
    scanf("%lf", &IL0);
    printf("Digite a unidade da Corrente (m para milli, u para micro, k para kilo, M para Mega ou pressione Enter para nenhuma): ");
    getchar(); // Limpar o buffer do teclado
    unidade_IL0 = getchar();

    // Converter todas as unidades para o padrão (base)
    R = converter_unidade(R, unidade_R);
    C = converter_unidade(C, unidade_C);
    L = converter_unidade(L, unidade_L);
    VC0 = converter_unidade(VC0, unidade_VC0);
    IL0 = converter_unidade(IL0, unidade_IL0);

    if (R <= 0 || L <= 0 || C <= 0) {
    printf("Erro: Resistência, Indutância e Capacitância devem ser maiores que zero.\n");
    return;
}
    double sigma = 0;
    double omega = 0;
    sigma = 1 / (2 * R * C);
    omega = 1 / sqrt(L * C);

    if(sigma > omega){
        calc_superamortecido(R,L,C,VC0,IL0,sigma,omega);
    } else if ( fabs(sigma - omega) < 1e-6){
         calc_criticamente(R,L,C,VC0,IL0,sigma,omega);
    } else if (sigma < omega){
        calc_subamortecido(R,L,C,VC0,IL0,sigma,omega);
    }
    printf("Voce quer colocar outros valores?(1- sim/2- nao): ");
    scanf("%d", &cont);
    }
}