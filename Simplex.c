#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

//sudo apt-get install libgsl-dev
/*gcc Simplex.c -o simplex -lgsl -lm */


double* ResolveSistemaLU(double **A, double *b) {

}

int indice_soma_coluna_max(double **A, int m, int n) {
    int s0, e; int s = INT_MIN;
    for(int j = 0; j < n; j++) {
        s0 = 0;
        for(int i = 0; i < m; i++) {
            s0 += A[i][j];
        }
        if(s < s0)  {
            s = s0;
            e = j;
        }
    }
    if (s <= 0) {
        printf("Solucao atual otima")
        return(0);
    }
    return(e);
}

int indice_minimo_proporcao(double *a_0, double *a_e, int m) {
    double p, p0 = FLT_MAX;
    int l;
    for(int i = 0; i < m; i++) {
        if(a_e[i] > 0) {
            p = (a_0[i]) / (a_e[i]);
        }
        if(p < p0) {
            p0 = p;
            l = i;
        }
    }
    return(l);
}


double* Fase1(double** A, double *b, int m, int n) {
    int e, int s;
    double *a_0 = (double*)malloc(n, sizeof(double));
    for(int i = 0; i < m; i++) {
        a_0[i] = b[i];
    }
    e = indice_soma_coluna_max(A, m, n);
    for(int i = 0; i < m; i++) {
        if(A[i][e] <= 0) {
            return("Solucao ilimitada");
        }
    }

    
}

double *Alocar_vetor(int n) {
    double *v = (double*)malloc(n, sizeof(double));
    return(v);
}

double **Alocar_Matriz(int m, int n) {
    double **A = (double**)malloc(m, sizeof(double*));
    for(int i = 0; i < n; i++) {
        A[i] = (double*)malloc(n, sizeof(double));
    }
    return(A);
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

double **Transpoe_Matriz(double **A, int m, int n){
    double **A_t = Alocar_Matriz(n,m);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            A_t[i][j] = A[j][i];
        }
    }
    return(A_t);
}


double* Simplex(double **A, double *c, double *b, int m, int n) {
    double **B; double* x_b, a_0, c_b, y;
    x_b = Alocar_vetor(m);
    a_0 = Alocar_vetor(m);
    c_b = Alocar_vetor(m);
    y = Alocar_vetor(m);
    B = Alocar_Matriz(m,m);
    int *indices_b = Fase1(A, b, m, n);

    a_0 = ResolveSistemaLU(B, b);
    x_b = a_0;
    for(int i = 0; i < m; i++) {
        c_b[i] = c[indices_b[i]];
    }
    
    y = ResolveSistemaLU(B_t, c_b);

    
}



int main() {
    int m; scanf("%d", &m);
    int n; scanf("%d", &n); 
    double *c = (double*)malloc(n, sizeof(double));
    double *b = (double*)malloc(m, sizeof(double));
    double **A = (double**)malloc(m, sizeof(double*));
    for(int i = 0; i < n; i++) {
        A[i] = (double*)malloc(n, sizeof(double));
    }

}




    









    
