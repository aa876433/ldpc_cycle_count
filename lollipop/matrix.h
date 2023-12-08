#ifndef MATRIX_H
#define MATRIX_H

typedef struct MATRIX_T
{
    int n;
    int m;
    int **data;
} MATRIX_T;

void matrix_show(MATRIX_T *mat);
void matrix_clear(MATRIX_T *mat);
MATRIX_T *matrix_alloc(int n, int m);
void matrix_free(MATRIX_T *mat);
void matrix_add(MATRIX_T *mat_a, MATRIX_T *mat_b, MATRIX_T *add);
void matrix_minus(MATRIX_T *mat_a, MATRIX_T *mat_b, MATRIX_T *minus);
void matrix_prod(MATRIX_T *mat_a, MATRIX_T *mat_b, MATRIX_T *prod);
void matrix_dot_coef(int coef_a, MATRIX_T *mat_a, int coef_b, MATRIX_T *mat_b, MATRIX_T *prod);
void matrix_z(MATRIX_T *mat);
void matrix_diag(MATRIX_T *mat);
int matrix_trace(MATRIX_T *mat);
void matrix_copy(MATRIX_T *dst, MATRIX_T *src);
void matrix_combination(int coef, MATRIX_T *src, int k, MATRIX_T *dst);

#endif // MATRIX_H
