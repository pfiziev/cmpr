#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
//#define ABS(a) (((a) >= 0) ? (a) : (-(a)))


#define TRUE 1
#define FALSE 0

#define E_UNDERFLOW 1
#define SUCCESS 0

typedef double ldouble;
//typedef long double ldouble;
#define INF HUGE_VAL

//typedef long double ldouble;
//#define INF HUGE_VALL

typedef int BOOL;
#define CACHE_NOT_FOUND INF

#define LEFT_BOUNDARY 0
#define RIGHT_BOUNDARY 1

#define FORWARD_DIRECTION 0
#define BACKWARD_DIRECTION 1

inline int ABS(int x) {
    return x >= 0 ? x : -x;
}

inline ldouble FABS(ldouble x) {
    return x >= 0 ? x : -x;
}

inline int sign(int x){
    return x / ABS(x);
}

ldouble ** matrix(int n_rows, int n_columns) {
    ldouble **M = (ldouble **) calloc(n_rows, sizeof (ldouble *));
    if (M == NULL) { fprintf(stderr, "Cannot allocate memory for matrix M!"); exit(1); }

    for (int i = 0; i < n_rows; i ++) {
        M[i] = (ldouble *) calloc(n_columns, sizeof(ldouble));
        if (M[i] == NULL) { fprintf(stderr, "Cannot allocate memory for matrix M!"); exit(1); }

    }
    return M;
}

ldouble * new_array(int n) {
    ldouble *a = (ldouble *) calloc(n, sizeof (ldouble));
    if (a == NULL) { fprintf(stderr, "Cannot allocate memory for array a!"); exit(1); }
    return a;
}

int ** integer_matrix(int n_rows, int n_columns) {
    int **M = (int **) calloc(n_rows, sizeof (int *));
    if (M == NULL) { fprintf(stderr, "Cannot allocate memory for matrix M!"); exit(1); }
    for (int i = 0; i < n_rows; i ++) {
        M[i] = (int *) calloc(n_columns, sizeof(int));
        if (M[i] == NULL) { fprintf(stderr, "Cannot allocate memory for matrix M!"); exit(1); }
    }
    return M;
}

void free_matrix(ldouble **M, int n) {
    for (int i = 0; i < n; i ++) {
        free(M[i]);
    }
    free(M);
}

void free_integer_matrix(int **M, int n) {
    for (int i = 0; i < n; i ++) {
        free(M[i]);
    }
    free(M);
}


void set_array(ldouble *array, int n, ldouble value) {
    for (int i = 0; i < n; i ++) {
        array[i] = value;
    }
}


void set_matrix(ldouble **M, int n_rows, int n_cols, ldouble value) {
    for (int i = 0; i < n_rows; i ++) {
        for (int j = 0; j < n_cols; j ++) {
            M[i][j] = value;
        }
    }
}


void set_integer_matrix(int **M, int n_rows, int n_cols, int value) {
    for (int i = 0; i < n_rows; i ++) {
        for (int j = 0; j < n_cols; j ++) {
            M[i][j] = value;
        }
    }
}


void print_array(ldouble *array, int n){
    printf("[ ");
    for (int i = 0; i < n; i ++) {
        printf("%.16le", array[i]);
        if (i < n - 1) {
            printf(", ");
        }
    }
    printf(" ]\n");
}

void print_array_int(int *array, int n){
    printf("[ ");
    for (int i = 0; i < n; i ++) {
        printf("%d", array[i]);
        if (i < n - 1) {
            printf(", ");
        }
    }
    printf(" ]\n");
}

void print_matrix(ldouble **M, int n_rows, int n_cols) {
    printf("[");
    for (int i = 0; i < n_rows; i ++) {
        printf("[");
        for (int j = 0; j < n_cols; j ++) {
            printf("%.16le", M[i][j]);
            if (j < n_cols - 1) { printf(", "); }
        }
        printf("]");
        if (i < n_rows - 1) { printf(",\n"); }

    }
    printf("]\n");
}

void print_matrix_int(int **M, int n_rows, int n_cols) {
    for (int i = 0; i < n_rows; i ++) {
        for (int j = 0; j < n_cols; j ++) {
            printf("%d\t", M[i][j]);
        }
        printf("\n");
    }
}


ldouble *** cube(int Z, int n_rows, int n_columns) {

    ldouble ***C = (ldouble ***) calloc(Z, sizeof (ldouble **));
    if (C == NULL) { fprintf(stderr, "Cannot allocate memory for cube C!"); exit(1); }

    for (int z = 0; z < Z; z ++) {
        C[z] = matrix(n_rows, n_columns);
    }
    return C;
}

int *** integer_cube(int Z, int n_rows, int n_columns) {
    int ***C = (int ***) calloc(Z, sizeof (int **));
    if (C == NULL) { fprintf(stderr, "Cannot allocate memory for cube C!"); exit(1); }

    for (int z = 0; z < Z; z ++) {
        C[z] = integer_matrix(n_rows, n_columns);
    }
    return C;
}

void free_cube(ldouble ***C, int Z, int n) {
    for (int z = 0; z < Z; z ++){
        for (int i = 0; i < n; i ++) {
            free(C[z][i]);
        }
        free(C[z]);
    }
    free(C);
}

void free_integer_cube(int ***C, int Z, int n) {
    for (int z = 0; z < Z; z ++){
        for (int i = 0; i < n; i ++) {
            free(C[z][i]);
        }
        free(C[z]);
    }
    free(C);
}


void set_cube(ldouble ***C, int Z, int n_rows, int n_cols, ldouble value) {
    for (int z = 0; z < Z; z ++) {
        for (int i = 0; i < n_rows; i ++) {
            for (int j = 0; j < n_cols; j ++) {
                C[z][i][j] = value;
            }
        }
    }
}

void set_integer_cube(int ***C, int Z, int n_rows, int n_cols, int value) {
    for (int z = 0; z < Z; z ++) {
        for (int i = 0; i < n_rows; i ++) {
            for (int j = 0; j < n_cols; j ++) {
                C[z][i][j] = value;
            }
        }
    }
}

void print_cube(ldouble ***C, int Z, int n_rows, int n_cols){

    for (int z = 0; z < Z; z ++) {
        printf("z = %d\n", z);
        printf("[");
        for (int i = 0; i < n_rows; i ++) {
            printf("[");
            for (int j = 0; j < n_cols; j ++) {
                printf("%.16le", C[z][i][j]);
                if (j < n_cols - 1) {printf(", ");}
            }
            printf("]");
            if (i < n_rows - 1) {printf(",\n");}

        }
        printf("]\n\n");
    }
}

void print_cube_int(int ***C, int Z, int n_rows, int n_cols){
    for (int z = 0; z < Z; z ++) {
        printf("z = %d\n", z);
        for (int i = 0; i < n_rows; i ++) {
            for (int j = 0; j < n_cols; j ++) {
                printf("%d\t", C[z][i][j]);
            }
            printf("\n");
        }
    }
}

ldouble * ldouble_array(int length) {
    ldouble * array = (ldouble *) calloc (length, sizeof (ldouble));
    if (array == NULL) { fprintf(stderr, "Cannot allocate memory for ldouble array with length %d", length); exit(1); }
    for (int i = 0; i < length; i++) {
        array[i] = 0;
    }
    return array;
}


void matcopy(ldouble **from_matrix, ldouble **to_matrix, int n_rows, int n_columns) {
    for (int i = 0; i < n_rows; i ++) {
        for (int j = 0; j < n_columns; j ++) {
            to_matrix[i][j] = from_matrix[i][j];
        }
    }
}

void array_copy(ldouble *from, ldouble *to, int n) {
    for (int i = 0; i < n; i ++) {
        to[i] = from[i];
    }
}


ldouble add_log_probs(ldouble log_X, ldouble log_Y) {
    if (log_X == -INF) {
        return log_Y;
    } else if (log_Y == -INF) {
        return log_X;
    }

    // swap them if log_Y is the bigger number
    if (log_X < log_Y) {
        ldouble _tmp = log_X;
        log_X = log_Y;
        log_Y = _tmp;
        }

    ldouble to_add = log(1 + exp(log_Y - log_X));
    if (to_add == -INF || to_add == INF) {
        return log_X;
    } else {
        return log_X + to_add;
    }
}


ldouble KL(ldouble * P, ldouble * Q, int n_probs) {
    ldouble kl = 0;

    for (int prob_idx = 0; prob_idx < n_probs; prob_idx ++ ) {
        kl += P[prob_idx] * log(P[prob_idx] / Q[prob_idx]) / log(2);
    }

    return (kl > 0) ? kl : 0; // protect from underflow errors that can yield a negative number
}


ldouble symmetric_KL_divergence(ldouble * probs_1, ldouble * probs_2, int n_probs) {
    return (KL(probs_1, probs_2, n_probs) + KL(probs_2, probs_1, n_probs)) / 2.;
}


//
//int main() {
//    for (long int x = 0; x < 2091259980; x++) {
//        //printf("%d %le\n", x, poisson_pmf(x, 5));
//        ldouble z = poisson_pmf(x % 1000, 50.);
//    }
//    return 0;
//}