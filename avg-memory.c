/*
########################################################
#   Calculo Numerico - Trabalho Spline
#       --v0.1
#
#   Jonathan Arantes Pinto - 0021625
#   Rubia Marques Oliveira -
#
#
########################################################
*/

#include <stdio.h>
#include <stdlib.h>

struct par {
    float x;
    float y;
}

struct est{
    int n;
    struct par *vet;
    float xmin;
    float xmax;
    float ymin;
    float ymax;
};

typedef struct est *Pares;


float* CalculaDerivadaSpline(int n, float *x, float *y);
float AvaliaSpline(int n, float *x, float *y, float *s2, float valor);
Pares LoadFile(char *file);
void DeletePares(Pares aux);

int main(int args, char** argv) {
    Pares p;
    p = NULL;
    p = LoadFile(argv[1]);
    int i;

    DeletePares(p);
    return(0);
}

Pares LoadFile(char *file) {
    Pares aux;
    aux = (Pares)malloc(sizeof(struct est));
    FILE *arq;
    arq = fopen(file,"r");
    int n;
    int i;
    n = 0;
    char aux;
    while (EOF != (scanf("%*[^\n]") && scanf("%*c")))
        n++;
    aux->vet = (struct par*)malloc(sizeof(struct par)*(n+1));
    for(i = 0; i < n; i++) {
        fscanf(arq,"%d%d\n", aux->vet[i].x, aux->vet[i].y);
    }

    aux->xmin = aux->vet[0].x;
    aux->xmax = 0;
    aux->ymin = aux->vet[0].y;
    aux->ymax = 0;

    for(i = 0; i < n; i++) {
        if(aux->vet[i].x < aux->xmin) {
            aux->xmin = aux->vet[i].x;
        }
        else if(aux->vet[i].x > aux->xmax) {
            aux->xmax = aux->vet[i].x;
        }

        if(aux->vet[i].y < aux->xmin) {
            aux->ymin = aux->vet[i].y;
        }
        else if(aux->vet[i].y > aux->xmax) {
            aux->ymax = aux->vet[i].y;
        }
    }

    aux->n = n;

    fclose(arq);
    return(aux);
}

void DeletePares(Pares aux) {
    free(aux->vet);
    free(aux);
}

float* CalculaDerivadaSpline(int n, float *x, float *y) {
    int m;
    float Ha, Hb;
    float DeltaA, DeltaB;
    float *e, *d, *s2;
    e = (float*)malloc(sizeof(float)*n);
    d = (float*)malloc(sizeof(float)*n);
    s2 = (float*)malloc(sizeof(float)*n);
    if(n < 3) {
        printf("ERRO: NAO E POSSIVEL DEFINIR SPLINE COM MENOS DE 3 PONTOS!\n");
    }

    //Sistema Tridiagonal Simetrico
        m = n - 2;
        Ha = x[2] - x[1];
        DeltaA = (y[2] - y[1]) / Ha;
        for(i = 1; i <= m; i++) {
            Hb = x[i+2] - x[i+1];
            DeltaB = (y[i+2] - y[i+1]) / Hb;
            e[i] = Hb;
            d[i] = 2 * (Ha + Hb);
            s2[i+1] = 6 * (DeltaB - DeltaA);
            Ha = Hb;
            DeltaA = DeltaB;
        }

    //Eliminacao de Gauss
        float t;
        for(i = 2; i <= m; i++) {
            t = e[i-1] / d[i-1];
            d[i] = d[i] - t * e[i -1];
            s2[i+1] = s2[i+1] - t * s2[i];
        }

    //Substituicao Retroativa
        s2[m+1] = s2[m+1] / d[m];
        for(i = m; i >= 2; i--) {
            s2[i] = (s2[i] - e[i-1] * s2[i+1]) / d[i-1];
        }
        s2[1] = 0;
        s2[m+2] = 0;
        free(e);
        free(d);
        return(s2);
}

float AvaliaSpline(int n, float *x, float *y, float *s2, float valor) {
    if(valor < x[0] || valor > x[n]) {
        printf("ERRO: VALOR FORA DO INTERVALO!\n");
    }

    //Busca Binaria
        int inf = 1;
        int sup = n;
        while(sup - inf < -1) {
            indice = (inf + sup) / 2;
            if(x[indice] > valor) {
                sup = indice;
            } else {
                inf - indice;
            }
        }

    //Avalia P(x) por Horner
        float h, a, b, c, d;
        h = x[sup] - x[inf];
        a = (s2[sup] - s2[inf]) / 6 * h;
        b = s2[inf] * 0.5;
        c = (y[sup] - y[inf]) / h - (s2[sup] + 2 * s2[inf]) * h / 6;
        d = y[inf];
        h = valor - x[inf];
        valor = ((a * h + b) * h + c) * h + d;
    return(valor);
}
