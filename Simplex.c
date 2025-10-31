#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>


//sudo apt-get install libgsl-dev
/*gcc Simplex.c -o simplex -lgsl -lm */

// Observação: Foram usadas notações inspiradas no livro "Linear and Nonlinear Programming" de David G. Luenberger e Yinyu Ye.

typedef struct {
    int m;
    int n;
    double **A;
    double *c;
    double *b;
} ProgramaLinear;



typedef struct {
    int *ind_b;
    int *ind_nb;
} Base;

double* ResolveSistemaLU(double **A, double *b, int n) {
    gsl_matrix *gA = gsl_matrix_alloc(n, n);
    gsl_vector *gb = gsl_vector_alloc(n);
    gsl_vector *gx = gsl_vector_alloc(n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int sinal;

    for (int i = 0; i < n; i++) {
        gsl_vector_set(gb, i, b[i]);
        for (int j = 0; j < n; j++) {
            gsl_matrix_set(gA, i, j, A[i][j]);
        }
    }

    gsl_linalg_LU_decomp(gA, p, &sinal);

    gsl_linalg_LU_solve(gA,p, gb, gx);

    double *x = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        x[i] = gsl_vector_get(gx, i);
    }

    gsl_permutation_free(p);
    gsl_vector_free(gb);
    gsl_vector_free(gx);
    gsl_matrix_free(gA);

    return x;

}

void Imprime_vetor(double *v, int n, const char *msg) {
    printf("%s\n", msg);
    for(int i = 0; i < n; i++)
        printf("%lf ", v[i]);
    
    printf("\n");
}

void Imprime_matriz(double **A, int m, int n, const char *msg) {
    printf("%s\n", msg);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            printf("%8.3lf", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double Multiplica_vetores(double *x, double *y, int n) {
    double s = 0;
    for(int i = 0; i < n; i++) {
        s += x[i] * y[i];
    }
    return(s);
}

int indice_minimo(double *x, int n) {
    double p = DBL_MAX; int min;
    for(int i = 0; i < n; i++) {
       if(x[i] < p) {
        min = i;
        p = x[i];
       }
    }
    return(min);
}

int indice_minimo_proporcao(double *x, double *y, int m) {
    double l = DBL_MAX; double r; int min = -1;
    for(int i = 0; i < m; i++) {
        if(y[i] > 0) {
            r = x[i] / y[i];
            if(r < l) {
                l = r;
                min = i;
            }
        }
    }
    return(min);
}

double *Aloca_vetor(int n) {
    double *v = malloc(n * sizeof(double));
    if(!v) {
        printf("Erro ao alocar memoria para o vetor\n");
        exit(1);
    }
    return(v);
}

double **Aloca_Matriz(int m, int n) {
    double **A = (double**)malloc(m * sizeof(double*));
    for(int i = 0; i < m; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
    }
    if(!A) {
        printf("Erro ao alocar memoria para a matriz\n");
        exit(1);
    }
    return(A);
}

void Libera_vetor(double **v) {
    if((*v) != NULL) {
        free(*v);
        *v = NULL;
    }
}
void Libera_Matriz(double*** A, int m) {
    if((*A) != NULL) {
        for(int i = 0; i < m; i++) {
            free((*A)[i]);
        }
        free(*A);
        (*A) = NULL;
    }
}


ProgramaLinear *Aloca_PL(int m, int n) {
    ProgramaLinear *pl = malloc(sizeof(ProgramaLinear));
    pl->m = m; pl->n = n;
    double **A = Aloca_Matriz(m,n);
    pl->A = A;
    double *c = Aloca_vetor(n);
    double *b = Aloca_vetor(m);
    pl->c = c; pl->b = b;
    if(!pl) {
        printf("Erro ao alocar memoria para o PL\n");
        exit(1);
    }
    return(pl);
}

Base *Aloca_Base(int m, int n) {
    Base *B = malloc(sizeof(Base));
    int *ind_b = (int*)malloc(m * sizeof(int));
    int *ind_nb = (int*)malloc((n - m) * sizeof(int));
    B->ind_b = ind_b;
    B->ind_nb = ind_nb;
    if(!B) {
        printf("Erro ao alocar memoria para a Base\n");
        exit(1);
    }
    return(B); 
}

void Libera_Base(Base **Bs){

     if(!(*Bs)) 
     return;

     if ((*Bs)->ind_b) 
     free((*Bs)->ind_b);

     if ((*Bs)->ind_nb)
     free((*Bs)->ind_nb);

     free((*Bs));

}

void Libera_PL(ProgramaLinear **pl) {
      if (!pl) 
      return;
      Libera_Matriz(&((*pl)->A), (*pl)->m);
      if(!(*pl)->b) 
        free((*pl)->b);
      if(!(*pl)->c) 
        free((*pl)->c);
    
        free((*pl));
    
}

double **Transpoe_Matriz(double **A, int m, int n){
    double **A_t = Aloca_Matriz(n,m);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            A_t[i][j] = A[j][i];
        }
    }
    return(A_t);
} 

Base* Atualiza_Base(Base * Bs, int e, int o, int m, int n) {
    int k = -1; int l = -1; int i = 0;
    for(int i = 0; i < m; i++) {
        if(Bs->ind_b[i] == o) k = i;
    }
    for(int i = 0; i < (n - m); i++) {
        if(Bs->ind_nb[i] == e) l = i;
    }

    Bs->ind_b[k] = e;
    Bs->ind_nb[l] = o;
    return(Bs);
}

void Simplex(ProgramaLinear* PL, Base *Bs) { 

    int m = PL->m; int n = PL->n; int o = -1; int l = 0;
    double *a_0 = Aloca_vetor(m); double *y = Aloca_vetor(m); double*r_N = Aloca_vetor(n - m); double *A_e = Aloca_vetor(m); double *a_e = Aloca_vetor(m);
    double **B = Aloca_Matriz(m, m); 


    while(l < 20) { // Loop ate encontrar solucao otima ou ilimitada
        printf("Iteracao %d da fase 2 do Simplex:\n", l + 1);
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < m; j++) {
                B[i][j] = PL->A[i][Bs->ind_b[j]]; // Constroi a matriz base B
            }
        }

        Imprime_matriz(B, m, m, "Matriz B:");

        double **B_t = Transpoe_Matriz(B, m, m);
        double *c_b = Aloca_vetor(m);

        for(int i = 0; i < m; i++) {
            c_b[i] = PL->c[Bs->ind_b[i]]; // Determina vetor de custos basicos c_b
        }
        Imprime_vetor(c_b, m, "Vetor custos basicos:");

        a_0 = ResolveSistemaLU(B, PL->b, m); // Resolve B * a_0 = b (determina o valor das variaveis basicas iniciais)
        Imprime_vetor(a_0, m, "Variaveis basicas iniciais a0: ");

        y = ResolveSistemaLU(B_t, c_b, m); // Resolve B^T * y = c_b (determina os multiplicadores simplex)
        Imprime_vetor(y, m, "Multiplicadores Simplex: ");

        for(int k = 0; k < (n - m); k++) { 
            double s = 0;
            for(int l = 0; l < m; l++) {
                s += y[l] * PL->A[l][Bs->ind_nb[k]]; // Calcula produto y^T * N
            }
            r_N[k] = PL->c[Bs->ind_nb[k]] - s; // Calcula custos reduzidos r_N = c_N - N^T * y

        }
        Imprime_vetor(r_N, (n - m), "Custos reduzidos rn:");

        int e_idc = indice_minimo(r_N, (n - m));
        if(r_N[e_idc] >= 0) {
            Imprime_vetor(a_0, m, "Solucao otima encontrada: ");
            double f = Multiplica_vetores(a_0, c_b, m); 
            printf("Valor otimo: %lf", f); // Se todos os custos reduzidos forem nao-negativos, solucao otima foi encontrada
            break;
        }
        int e = Bs->ind_nb[e_idc]; 
        printf("Indice a entrar na base: %d \n", e);
        
        for(int i = 0; i < m; i++) {
            A_e[i] = PL->A[i][e]; // Define coluna A_e a entrar na base
        }
        Imprime_vetor(A_e, m, "Coluna a entrar na base: " );

        a_e = ResolveSistemaLU(B, A_e, m);
        Imprime_vetor(a_e, m, "Combinacao linear das colunas de B que formam coluna A_e");

        int o_idc = indice_minimo_proporcao(a_0, a_e, m); // Encontra indice da variavel a sair da base
        int o = Bs->ind_b[o_idc];
        printf("indice a sair base: %d \n", o);
         
        if(o == -1) {
            printf("Solucao Ilimitada \n");
            return(NULL);
        }
   
        Bs = Atualiza_Base(Bs, e, o, m, n);
        printf("Indices basicos encontrados:\n");
        for(int i = 0; i < m; i++){ 
            printf("%d \n", (Bs->ind_b)[i]);
        }
        printf("\n");
        Libera_vetor(&a_0);
        Libera_Matriz(&B_t, m);
        l++;

    }
    Libera_PL;
    return;
}

Base *Fase1(ProgramaLinear *pl) {

    int m = pl->m; int n = pl->n;

    Base *Bs_fase1 = Aloca_Base(m, (n + m));
    
    ProgramaLinear *pl_fase1 = Aloca_PL(m, (n + m));

    printf("Inicia fase 1:\n");

   for(int i = 0; i < m; i++) {
       Bs_fase1->ind_b[i] = i;
   }

   for(int j = 0; j < pl->n; j++) {
    Bs_fase1->ind_nb[j] = j + pl->m;
   }

   pl_fase1->m = pl->m; pl_fase1->n = (pl->m + pl->n); pl_fase1->b = pl->b;
   for(int j = 0; j < (pl->m + pl->n); j++) {
        for(int i = 0; i < pl->m; i++ ) {
            if(j < pl->m) {
                if(i == j) {
                    pl_fase1->A[i][j] = 1;
                } 
                else {
                    pl_fase1->A[i][j] = 0;
                }
            }
            else {
                pl_fase1->A[i][j] = pl->A[i][j - pl->m];
            }
        }
    }

    for(int i = 0; i < pl->n; i++) {
        pl_fase1->c[i] = 1; 
    }

    Base *Base_Simplex = Aloca_Base(m, n);
    Base_Simplex = Simplex(pl_fase1, Bs_fase1);
    
    Libera_Base(&Bs_fase1);
    Libera_PL(&pl_fase1);

    return(Base_Simplex);

}


int main() {

    printf("Digite m (numero de restricoes): ");

    int m; scanf("%d", &m);

    printf("Digite n (numero de variaveis): ");
    int n; scanf("%d", &n); 

    double *c = Aloca_vetor(n);
    printf("Digite o vetor c (coeficientes da funcao objetivo):\n");
    for(int i = 0; i < n; i++) {
        scanf("%lf", &c[i]);
    }
    
    double *b = Aloca_vetor(m);
    printf("Digite o vetor b (lado direito das restricoes):\n");
     for(int j = 0; j < m; j++) {
        scanf("%lf", &b[j]);
    }
    
    double **A = Aloca_Matriz(m,n);
    printf("Digite a matriz A (coeficientes das restricoes):\n");
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    ProgramaLinear *pl = Aloca_PL(m, n);

    Base *B = Aloca_Base(m, n);
    B->ind_b[0] = 0;
    B->ind_b[1] = 2;
    B->ind_nb[0] = 1;
    B->ind_nb[1] = 3;

    if(B == NULL) {
        printf("Nao foi possivel encontrar uma base viavel!\n");
        return(0);
    }

    printf("Base viavel encontrada pela Fase I:\n");
    printf("Indices basicos:\n");
    for(int i = 0; i < m; i++)
        printf("%d ", B->ind_b[i]);

    pl->m = m; pl->n = n; pl->c = c; pl->b = b; pl->c = c; pl->A = A;

    Simplex(pl, B);

    Libera_PL(&pl);

    printf("Indices basicos encontrados:\n");
    for(int i = 0; i < m; i++)
        printf("%d\n", B->ind_b[i]);

    x = Aloca_vetor(m);
    for(int i = 0; i < m; i++) {
        x[i] = 0;
    }
    

    return(1);
    

}




    









    
