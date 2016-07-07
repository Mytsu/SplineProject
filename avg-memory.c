/*
#######################################################################
#   Calculo Numerico - Trabalho Spline
#       --v0.1
#
#   Jonathan Arantes Pinto - 0021625
#   Rubia Marques Oliveira - 0022430
#
#
#######################################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct par {
    float x;
    float y;
};

struct est{
    int n;
    struct par *vet;
    float xmin;
    float xmax;
    float ymin;
    float ymax;
};

typedef struct est *Pares;


float* CalculaDerivadaSpline(Pares p);
float AvaliaSpline(Pares p, float *s2, float valor);
Pares LoadFile(char *file);
void DeletePares(Pares aux);
void PrintExit(int n, double avg, char *saida);
float Random(float xmin, float xmax);
void RscriptGenerate(Pares p, double avg, char *saida, float *s2);
double IntegralMonteCarlo(Pares p, float *s2);

int main(int args, char** argv) {
    Pares p;
    float *s2;
    double avg;
    p = LoadFile(argv[1]);
    s2 = CalculaDerivadaSpline(p);
    avg = IntegralMonteCarlo(p, s2);
    RscriptGenerate(p, avg, argv[2], s2);
    PrintExit(p->n, avg, argv[2]);
    DeletePares(p);
    return(0);
}

/* Função que carrega o arquivo e preenche a estrutura com os valores informados */
Pares LoadFile(char *file) {
    Pares aux;
    aux = (Pares)malloc(sizeof(struct est));
    FILE *arq;
    arq = fopen(file,"r");
    int n;
    int i;
    n = 0;

    /* Contador de Linhas */
    while (fscanf(arq,"%*f %*f\n") != EOF) {
        n++;
    }
    fclose(arq);
    arq = fopen(file,"r");

    aux->vet = (struct par*)malloc(sizeof(struct par)*(n));
    for(i = 0; i < n; i++) {
        fscanf(arq,"%f %f\n", &aux->vet[i].x, &aux->vet[i].y);
    }

    /* Setando valores de X min e max, e Y min e max */
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

/* Funcao que mostra o resultado de saida do terminal para o usuario */
void PrintExit(int n, double avg, char *saida) {
    printf("Number of Samples\t: %d\n", n);
    printf("Average Memory Usage\t: %.3f\n", avg);
    printf("\n");
    printf("Run 'Rscript %s' to generate Average Memory Usage Chart\n", saida);
    return;
}

/* Funcao propria para limpar a estrutura */
void DeletePares(Pares aux) {
    free(aux->vet);
    free(aux);
    return;
}

/* Funcao U(xmin,xmax) */
float Random(float xmin, float xmax) {
    float random;
    srand((unsigned)time(NULL));
    random = ((rand())/(float)RAND_MAX)*(xmax-xmin)+xmin;
    return(random);
}

void RscriptGenerate(Pares p, double avg, char *saida, float *s2) {
    FILE *arq;
    int i;

    arq = fopen(saida,"w");
    fprintf(arq, "#\n");
    fprintf(arq, "# Generated automatically by 'avg-memory' application");
    fprintf(arq, "#\n");
    fprintf(arq, "\n");
    fprintf(arq, "# Original points (x coordinates)\n");
    fprintf(arq, "xorig <- c(\n");
    for(i = 0; i < p->n; i++) {
        fprintf(arq, "\t%.0f,\n", p->vet[i].x);
    }
    fprintf(arq, "\t%.0f\n", p->vet[p->n+1].x);
    fprintf(arq, ");\n");
    fprintf(arq, "\n");
    fprintf(arq, "# Original points (y coordinates)\n");
    fprintf(arq, "yorig <- c (\n");
    for(i = 0; i < p->n; i++) {
        fprintf(arq, "\t%.0f,\n", p->vet[i].y);
    }
    fprintf(arq, "\t%.0f\n", p->vet[p->n+1].y);
    fprintf(arq, ");\n");
    fprintf(arq, "\n");
    fprintf(arq, "# Spline points (x coordinates, sampling interval = 0.01)\n");
    fprintf(arq, "xspl <- c(\n");
    float intervalo;
    float xinc;
    intervalo = 0.01;
    xinc = p->xmin;
    while(xinc < p->xmax) {
        fprintf(arq, "\t%.2f,\n", xinc);
        xinc = xinc + intervalo;
    }
    fprintf(arq, "\t%.2f\n", xinc);
    fprintf(arq, ");\n");
    fprintf(arq, "\n");
    fprintf(arq, "# Spline points (y coordinates, sampling interval = 0.01)\n");
    fprintf(arq, "yspl <- c(\n");
    xinc = p->xmin;
    while(xinc < p->xmax) {
        fprintf(arq, "\t%.5f,\n", AvaliaSpline(p, s2, xinc));
        xinc = xinc + intervalo;
    }
    fprintf(arq, "\t%.5f\n", AvaliaSpline(p, s2, xinc));
    fprintf(arq, ");\n");
    fprintf(arq, "\n");
    fprintf(arq, "# Average Memory Usage\n");
    fprintf(arq, "AvgMemory <- %.3f\n", avg);
    fprintf(arq, "\n");
    fprintf(arq, "# Plot the values in .png file\n");
    fprintf(arq, "png(file='%s.png', width=1200);\n", saida);
    fprintf(arq, "title <- paste('AVG Memory Usage: %.3f Kb (%d Samples)');\n", avg, p->n);
    fprintf(arq, "plot(xspl, yspl, type='l', col='blue', main=title, xlab='Samples', ylab='Mem. Usage', lwd=3);\n");
    fprintf(arq, "points(xorig, yorig, pch=19, col='red');\n");
    fprintf(arq, "lines( c(min(xorig), max(xorig)), c(AvgMemory, AvgMemory), col='black', lty=2, lwd=3);\n");
    fprintf(arq, "dev.off();");

    fclose(arq);
}

/* Funcao Derivada por Spline Cubica informada no enunciado */
float* CalculaDerivadaSpline(Pares p) {
    int m;
    float Ha, Hb;
    float DeltaA, DeltaB;
    float *e, *d, *s2;
    int i;
    e = (float*)malloc(sizeof(float)*p->n);
    d = (float*)malloc(sizeof(float)*p->n);
    s2 = (float*)malloc(sizeof(float)*p->n);
    if(p->n < 3) {
        printf("ERRO: NAO E POSSIVEL DEFINIR SPLINE COM MENOS DE 3 PONTOS!\n");
        DeletePares(p);
        exit(0);
    }

    /* Sistema Tridiagonal Simetrico */
        m = p->n - 1;
        Ha = p->vet[1].x - p->vet[0].x;
        DeltaA = (p->vet[1].y - p->vet[0].y) / Ha;
        for(i = 0; i <= m; i++) {
            Hb = p->vet[i+1].x - p->vet[i].x;
            DeltaB = (p->vet[i+1].y - p->vet[i].y) / Hb;
            e[i] = Hb;
            d[i] = 2 * (Ha + Hb);
            s2[i+1] = 6 * (DeltaB - DeltaA);
            Ha = Hb;
            DeltaA = DeltaB;
        }

    /* Eliminacao de Gauss */
        float t;
        for(i = 1; i <= m; i++) {
            t = e[i-1] / d[i-1];
            d[i] = d[i] - t * e[i -1];
            s2[i+1] = s2[i+1] - t * s2[i];
        }

    /* Substituicao Retroativa */
        s2[m] = s2[m+1] / d[m];
        for(i = m; i >= 1; i--) {
            s2[i] = (s2[i] - e[i-1] * s2[i+1]) / d[i-1];
        }
        s2[1] = 0;
        s2[m+1] = 0;
        free(e);
        free(d);
        return(s2);
}

/* Funcao Avalia Spline Cubica ou f(x) informada no enunciado */
float AvaliaSpline(Pares p, float *s2, float valor) {
    int indice;
    if((valor < p->xmin ) || (valor > p->xmax)) {
        printf("Valor Aleatorio: %.2f\n", valor);
        printf("ERRO: VALOR FORA DO INTERVALO!\n");
    }

    /* Busca Binaria */
        int inf = 1;
        int sup = p->n;
        while(sup - inf > 1) {
            indice = (inf + sup) / 2;
            if(p->vet[indice].x > valor) {
                sup = indice;
            } else {
                inf = indice;
            }
        }

    /* Avalia P(x) por Horner */
        float h, a, b, c, d;
        double resultado;
        h = p->vet[sup].x - p->vet[inf].x;
        a = (s2[sup] - s2[inf]) / 6 * h;
        b = s2[inf] * 0.5;
        c = (p->vet[sup].y - p->vet[inf].y) / h - (s2[sup] + 2 * s2[inf]) * h / 6;
        d = p->vet[inf].y;
        h = valor - p->vet[inf].x;
        resultado = ((a * h + b) * h + c) * h + d;
    return(resultado);
}

/* Integração por Monte Carlo */

double IntegralMonteCarlo(Pares p, float *s2) {
    float num_abaixo , x, y;
    double AreaTotal, Area;
    int i;
    num_abaixo = 0.0;
     for (i = 1; i <= p->n; i++){ /* laço de 1 até o número menor ou igual a n */
        x = Random(p->xmin,p->xmax);  /* usando a função para gerar números pseudoaleatórios */
        y = Random(p->ymin,p->ymax);
         if (y <= AvaliaSpline(p, s2, x)){
            num_abaixo++; /* contador */
         }
     }
     AreaTotal = (p->xmax - p->xmin) * (p->ymax - p->ymin); /* cálculo da área total */
     Area = AreaTotal * (num_abaixo / p->n); /* cálculo da área */
    return (Area);
}
