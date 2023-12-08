#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include "matrix.h"

int find_matrix_max(const MATRIX_T *mat)
{
    int max = mat->data[0][0];
    for (int i = 0; i < mat->n; i++)
    {
        for (int j = 0; j < mat->m; j++)
        {
            if (mat->data[i][j] > max)
            {
                max = mat->data[i][j];
            }
        }
    }
    return max;
}

void matrix_show(MATRIX_T *mat)
{
    int max_width = snprintf(NULL, 0, "%d", find_matrix_max(mat));
    for (int i = 0; i < mat->n; i++)
    {
        for (int j = 0; j < mat->m; j++)
        {
            printf("%*d ", max_width, mat->data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

MATRIX_T *matrix_alloc(int n, int m)
{
    MATRIX_T *mat = malloc(sizeof(MATRIX_T));
    assert(mat != NULL);
    mat->n = n;
    mat->m = m;
    mat->data = malloc(n * sizeof(int *));
    assert(mat->data != NULL);
    for (int i = 0; i < n; i++)
    {
        mat->data[i] = malloc(m * sizeof(int));
        assert(mat->data[i] != NULL);
    }

    return mat;
}

void matrix_free(MATRIX_T *mat)
{
    for (int i = 0; i < mat->n; i++)
    {
        free(mat->data[i]);
    }
    free(mat->data);
    free(mat);
}

void matrix_clear(MATRIX_T *mat)
{
    for (int i = 0; i < mat->n; i++)
    {
        for (int j = 0; j < mat->m; j++)
        {
            mat->data[i][j] = 0;
        }
    }
}

void matrix_minus(MATRIX_T *mat_a, MATRIX_T *mat_b, MATRIX_T *minus)
{
    if (mat_a->n != mat_b->n || mat_a->m != mat_b->m)
    {
        assert(0);
        return;
    }

    for (int i = 0; i < mat_a->n; i++)
    {
        for (int j = 0; j < mat_a->m; j++)
        {
            minus->data[i][j] = mat_a->data[i][j] - mat_b->data[i][j];
        }
    }
}

void matrix_add(MATRIX_T *mat_a, MATRIX_T *mat_b, MATRIX_T *add)
{
    if (mat_a->n != mat_b->n || mat_a->m != mat_b->m)
    {
        assert(0);
        return;
    }

    for (int i = 0; i < mat_a->n; i++)
    {
        for (int j = 0; j < mat_a->m; j++)
        {
            add->data[i][j] = mat_a->data[i][j] + mat_b->data[i][j];
        }
    }
}

void matrix_prod(MATRIX_T *mat_a, MATRIX_T *mat_b, MATRIX_T *prod)
{
    if (mat_a->m != mat_b->n || (prod->n != mat_a->n || prod->m != mat_b->m))
    {
        assert(0);
        return;
    }

    int sum;
    for (int i = 0; i < mat_a->n; i++)
    {
        for (int j = 0; j < mat_b->m; j++)
        {
            sum = 0;
            for (int k = 0; k < mat_a->m; k++)
            {
                sum += mat_a->data[i][k] * mat_b->data[k][j];
            }
            prod->data[i][j] = sum;
        }
    }
}

void matrix_dot_coef(int coef_a, MATRIX_T *mat_a, int coef_b, MATRIX_T *mat_b, MATRIX_T *dot)
{
    if ((mat_a->n != mat_b->n || mat_a->m != mat_b->m) || (dot->n != mat_a->n || dot->m != mat_b->m))
    {
        assert(0);
        return;
    }

    for (int i = 0; i < mat_a->n; i++)
    {
        for (int j = 0; j < mat_b->m; j++)
        {
            dot->data[i][j] = (coef_a * mat_a->data[i][j]) * (coef_b * mat_b->data[i][j]);
        }
    }
}

int matrix_trace(MATRIX_T *mat)
{
    if (mat->n != mat->m)
    {
        assert(0);
        return -1;
    }

    int sum = 0;
    for (int i = 0; i < mat->n; i++)
    {
        sum += mat->data[i][i];
    }
    return sum;
}

void matrix_diag(MATRIX_T *mat)
{
    for (int i = 0; i < mat->n; i++)
    {
        for (int j = 0; j < mat->m; j++)
        {
            if (i != j)
            {
                mat->data[i][j] = 0;
            }
        }
    }
}

void matrix_z(MATRIX_T *mat)
{
    for (int i = 0; i < mat->n; i++)
    {
        mat->data[i][i] = 0;
    }
}

void matrix_copy(MATRIX_T *dst, MATRIX_T *src)
{
    if (dst->n != src->n || dst->m != src->m)
    {
        assert(0);
        return;
    }

    for (int i = 0; i < dst->n; i++)
    {
        for (int j = 0; j < dst->m; j++)
        {
            dst->data[i][j] = src->data[i][j];
        }
    }
}

int combination(int a, int b)
{
    if (a <= b)
    {
        return a == b;
    }

    if (b > a - b)
        b = a - b;

    int result = 1;
    for (int i = 0; i < b; ++i)
    {
        result *= (a - i);
        result /= (i + 1);
    }

    return result;
}

void matrix_combination(int coef, MATRIX_T *src, int k, MATRIX_T *dst)
{
    if (src->n != dst->n || src->m != dst->m)
    {
        assert(0);
        return;
    }

    for (int i = 0; i < src->n; i++)
    {
        for (int j = 0; j < src->m; j++)
        {
            dst->data[i][j] = coef * combination(src->data[i][j], k);
        }
    }
}