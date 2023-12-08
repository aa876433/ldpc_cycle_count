#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "matrix.h"
#include "calc_cycle.h"
#include "read_h.h"

#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define MIN(a, b) ((a) < (b)) ? (a) : (b)

typedef struct CYCLE_PRE_INFO_T
{
    CYCLE_INFO_T cycle;
    MATRIX_T *edge;
    MATRIX_T *edge_t;
    MATRIX_T *p_u_u;       // P^U_2
    MATRIX_T *p_u_w;       // P^U_3
    MATRIX_T *p_w_w;       // P^W_2
    MATRIX_T *p_w_u;       // P^W_3
    MATRIX_T *l_u_u;       // L^U_{2, 2}
    MATRIX_T *l_u_w;       // L^U_{1, 2}
    MATRIX_T *l_w_w;       // L^W_{2, 2}
    MATRIX_T *l_w_u;       // L^W_{1, 2}
    MATRIX_T *c_p_u_u;     // P^U_g-2
    MATRIX_T *c_p_u_w;     // P^U_g-1
    MATRIX_T *c_p_w_w;     // P^W_g-2
    MATRIX_T *c_p_w_u;     // P^W_g-1
    MATRIX_T *c_l_u_u;     // L^U_{g-2, 2}
    MATRIX_T *c_l_u_w;     // L^U_{g-3, 2}
    MATRIX_T *c_l_w_w;     // L^W_{g-2, 2}
    MATRIX_T *c_l_w_u;     // L^W_{g-3, 2}
    MATRIX_T *l_u_u_m;     // max(L^U_{0, 2} - 1, 0)
    MATRIX_T *l_w_w_m;     // max(L^W_{0, 2} - 1, 0)
    MATRIX_T *l_u_u_m2;    // max(L^U_{0, 2} - 2, 0)
    MATRIX_T *l_w_w_m2;    // max(L^W_{0, 2} - 2, 0)
    MATRIX_T *l_u_u_g;     // L^U_{0, g}
    MATRIX_T *l_w_w_g;     // L^W_{0, g}

} CYCLE_PRE_INFO_T;

void edge_init(PARITY_INFO_T *h, CYCLE_PRE_INFO_T *c_info)
{
    int n = h->blk * h->row;
    int m = h->blk * h->col;
    c_info->edge = matrix_alloc(n, m);
    c_info->edge_t = matrix_alloc(m, n);
    if (h->blk == 1)
    {
        for (int i = 0; i < c_info->edge->n; i++)
        {
            for (int j = 0; j < c_info->edge->m; j++)
            {
                c_info->edge->data[i][j] = h->build_matirx[i][j];
            }
        }
    }
    else
    {
        for (int i = 0; i < c_info->edge->n; i++)
        {
            for (int j = 0; j < c_info->edge->m; j++)
            {
                c_info->edge->data[i][j] = 0;
            }
        }

        int r, c;
        r = 0;
        for (int i = 0; i < h->row; i++)
        {
            c = 0;
            for (int j = 0; j < h->col; j++)
            {
                if (h->build_matirx[i][j] != -1)
                {
                    for (int k = 0; k < h->blk; k++)
                    {
                        int t = (k + h->build_matirx[i][j]) % h->blk;
                        c_info->edge->data[r + k][c + t] = 1;
//                        printf("%d, %d\n", r+k, c + t);
                    }
                }
                c += h->blk;
            }
            r += h->blk;
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            c_info->edge_t->data[j][i] = c_info->edge->data[i][j];
        }
    }
}

void update_lollipop_girth(CYCLE_PRE_INFO_T *c_info, MATRIX_T *p_u_w, MATRIX_T *p_w_u)
{
    // calc L^U_{0, g}
    matrix_prod(p_u_w, c_info->edge_t, c_info->l_u_u_g);
    matrix_diag(c_info->l_u_u_g);

    // calc L^W_{0, g}
    matrix_prod(p_w_u, c_info->edge, c_info->l_w_w_g);
    matrix_diag(c_info->l_w_w_g);
}

void l_p_init(CYCLE_PRE_INFO_T *c_info)
{
    int n = c_info->edge->n;
    int m = c_info->edge->m;
    MATRIX_T *buf_uu = matrix_alloc(n, n);
    MATRIX_T *buf_uw = matrix_alloc(n, m);
    MATRIX_T *buf_ww = matrix_alloc(m, m);
    MATRIX_T *buf_wu = matrix_alloc(m, n);

    // calc P^U_2, max(L^U_{0, 2} - 1, 0);
    matrix_prod(c_info->edge, c_info->edge_t, buf_uu);
    matrix_copy(c_info->l_u_u_m, buf_uu);
    matrix_diag(c_info->l_u_u_m);
    matrix_minus(buf_uu, c_info->l_u_u_m, c_info->p_u_u);
    for (int i = 0; i < c_info->l_u_u_m->n; i++)
    {
        c_info->l_u_u_m2->data[i][i] = MAX(c_info->l_u_u_m->data[i][i] - 2, 0);
        c_info->l_u_u_m->data[i][i] = MAX(c_info->l_u_u_m->data[i][i] - 1, 0);
    }

    // calc P^W_2, max(L^W_{0, 2} - 1, 0)
    matrix_prod(c_info->edge_t, c_info->edge, buf_ww);
    matrix_copy(c_info->l_w_w_m, buf_ww);
    matrix_diag(c_info->l_w_w_m);
    matrix_minus(buf_ww, c_info->l_w_w_m, c_info->p_w_w);
    for (int i = 0; i < c_info->l_w_w_m->n; i++)
    {
        c_info->l_w_w_m2->data[i][i] = MAX(c_info->l_w_w_m->data[i][i] - 2, 0);
        c_info->l_w_w_m->data[i][i] = MAX(c_info->l_w_w_m->data[i][i] - 1, 0);
    }

    // calc L^U_{1, 2}, P^U_3
    matrix_prod(c_info->edge, c_info->l_w_w_m, c_info->l_u_w);
    matrix_prod(c_info->p_u_u, c_info->edge, buf_uw);
    matrix_minus(buf_uw, c_info->l_u_w, c_info->p_u_w);

    // calc L^W_{1, 2}, P^W_3
    matrix_prod(c_info->edge_t, c_info->l_u_u_m, c_info->l_w_u);
    matrix_prod(c_info->p_w_w, c_info->edge_t, buf_wu);
    matrix_minus(buf_wu, c_info->l_w_u, c_info->p_w_u);

    // calc L^U_{2, 2}
    matrix_prod(c_info->edge, c_info->l_w_u, c_info->l_u_u);
    matrix_z(c_info->l_u_u);

    // calc L^W_{2, 2}
    matrix_prod(c_info->edge_t, c_info->l_u_w, c_info->l_w_w);
    matrix_z(c_info->l_w_w);

    // calc L^U_{0, g}, calc L^W_{0, g}
    update_lollipop_girth(c_info, c_info->p_u_w, c_info->p_w_u);

    matrix_copy(c_info->c_p_u_u, c_info->p_u_u);
    matrix_copy(c_info->c_p_u_w, c_info->p_u_w);
    matrix_copy(c_info->c_p_w_w, c_info->p_w_w);
    matrix_copy(c_info->c_p_w_u, c_info->p_w_u);

    matrix_copy(c_info->c_l_u_u, c_info->l_u_u);
    matrix_copy(c_info->c_l_u_w, c_info->l_u_w);
    matrix_copy(c_info->c_l_w_w, c_info->l_w_w);
    matrix_copy(c_info->c_l_w_u, c_info->l_w_u);

    matrix_free(buf_uu);
    matrix_free(buf_uw);
    matrix_free(buf_ww);
    matrix_free(buf_wu);
}

void *cycle_init(void *h_info)
{
    PARITY_INFO_T *h = (PARITY_INFO_T *) h_info;
    CYCLE_PRE_INFO_T *c_info = malloc(sizeof(CYCLE_PRE_INFO_T));
    edge_init(h, c_info);
    int n = c_info->edge->n;
    int m = c_info->edge->m;
    c_info->p_u_u = matrix_alloc(n, n);
    c_info->p_u_w = matrix_alloc(n, m);
    c_info->p_w_w = matrix_alloc(m, m);
    c_info->p_w_u = matrix_alloc(m, n);
    c_info->l_u_u = matrix_alloc(n, n);
    c_info->l_u_w = matrix_alloc(n, m);
    c_info->l_w_w = matrix_alloc(m, m);
    c_info->l_w_u = matrix_alloc(m, n);

    c_info->c_p_u_u = matrix_alloc(n, n);
    c_info->c_p_u_w = matrix_alloc(n, m);
    c_info->c_p_w_w = matrix_alloc(m, m);
    c_info->c_p_w_u = matrix_alloc(m, n);
    c_info->c_l_u_u = matrix_alloc(n, n);
    c_info->c_l_u_w = matrix_alloc(n, m);
    c_info->c_l_w_w = matrix_alloc(m, m);
    c_info->c_l_w_u = matrix_alloc(m, n);

    c_info->l_u_u_m = matrix_alloc(n, n);
    c_info->l_w_w_m = matrix_alloc(m, m);
    c_info->l_u_u_m2 = matrix_alloc(n, n);
    matrix_clear(c_info->l_u_u_m2);
    c_info->l_w_w_m2 = matrix_alloc(m, m);
    matrix_clear(c_info->l_w_w_m2);

    c_info->l_u_u_g = matrix_alloc(n, n);
    c_info->l_w_w_g = matrix_alloc(m, m);

    l_p_init(c_info);
    return c_info;
}

int calc_cycle_count(MATRIX_T *l_u_u_g, int k2)
{
    int count = matrix_trace(l_u_u_g) / (k2);
    return count;
}

void calc_girth(CYCLE_PRE_INFO_T *c_info)
{
    int n = c_info->edge->n;
    int m = c_info->edge->m;
    MATRIX_T *buf_uu = matrix_alloc(n, n);
    MATRIX_T *buf_uw = matrix_alloc(n, m);
    MATRIX_T *buf_ww = matrix_alloc(m, m);
    MATRIX_T *buf_wu = matrix_alloc(m, n);

    MATRIX_T *n_p_u_u = matrix_alloc(n, n);
    MATRIX_T *n_p_u_w = matrix_alloc(n, m);
    MATRIX_T *n_p_w_w = matrix_alloc(m, m);
    MATRIX_T *n_p_w_u = matrix_alloc(m, n);

    MATRIX_T *n_l_u_u = matrix_alloc(n, n);
    MATRIX_T *n_l_u_w = matrix_alloc(n, m);
    MATRIX_T *n_l_w_w = matrix_alloc(m, m);
    MATRIX_T *n_l_w_u = matrix_alloc(m, n);

    int girth = 4;
    int count = calc_cycle_count(c_info->l_u_u_g, girth);
    while (count == 0)
    {
        // calc L^U_{2k-1, 2}
        matrix_prod(c_info->edge, c_info->c_l_w_w, n_l_u_w);
        matrix_prod(c_info->l_u_u_m, c_info->c_l_u_w, buf_uw);
        matrix_minus(n_l_u_w, buf_uw, n_l_u_w);

        // calc L^W_{2k-1, 2}
        matrix_prod(c_info->edge_t, c_info->c_l_u_u, n_l_w_u);
        matrix_prod(c_info->l_w_w_m, c_info->c_l_w_u, buf_wu);
        matrix_minus(n_l_w_u, buf_wu, n_l_w_u);

        // calc L^U_{2k, 2}
        matrix_prod(c_info->edge, n_l_w_u, n_l_u_u);
        matrix_prod(c_info->l_u_u_m, c_info->c_l_u_u, buf_uu);
        matrix_minus(n_l_u_u, buf_uu, n_l_u_u);

        // calc L^W_{2k, 2}
        matrix_prod(c_info->edge_t, n_l_u_w, n_l_w_w);
        matrix_prod(c_info->l_w_w_m, c_info->c_l_w_w, buf_ww);
        matrix_minus(n_l_w_w, buf_ww, n_l_w_w);

        // calc P^U_2k
        matrix_prod(c_info->c_p_u_w, c_info->edge_t, buf_uu);
        matrix_minus(buf_uu, c_info->c_l_u_u, n_p_u_u);

        // calc P^U_2k+1
        matrix_prod(n_p_u_u, c_info->edge, buf_uw);
        matrix_minus(buf_uw, n_l_u_w, n_p_u_w);

        // calc P^W_2k
        matrix_prod(c_info->c_p_w_u, c_info->edge, buf_ww);
        matrix_minus(buf_ww, c_info->c_l_w_w, n_p_w_w);

        // calc P^W_2k+1
        matrix_prod(n_p_w_w, c_info->edge_t, buf_wu);
        matrix_minus(buf_wu, n_l_w_u, n_p_w_u);

        // calc L^U_{0, g}, L^W_{0, g}
        update_lollipop_girth(c_info, n_p_u_w, n_p_w_u);

        girth += 2;
        count = calc_cycle_count(c_info->l_u_u_g, girth);
        matrix_copy(c_info->c_p_u_u, n_p_u_u);
        matrix_copy(c_info->c_p_u_w, n_p_u_w);
        matrix_copy(c_info->c_p_w_w, n_p_w_w);
        matrix_copy(c_info->c_p_w_u, n_p_w_u);
        matrix_copy(c_info->c_l_u_u, n_l_u_u);
        matrix_copy(c_info->c_l_u_w, n_l_u_w);
        matrix_copy(c_info->c_l_w_w, n_l_w_w);
        matrix_copy(c_info->c_l_w_u, n_l_w_u);
    }

    c_info->cycle.girth = girth;
    c_info->cycle.no_g = count;

    matrix_free(buf_uu);
    matrix_free(buf_uw);
    matrix_free(buf_ww);
    matrix_free(buf_wu);
    matrix_free(n_p_u_u);
    matrix_free(n_p_u_w);
    matrix_free(n_p_w_w);
    matrix_free(n_p_w_u);
    matrix_free(n_l_u_u);
    matrix_free(n_l_u_w);
    matrix_free(n_l_w_w);
    matrix_free(n_l_w_u);
}

CYCLE_INFO_T calc_cycle(void *ctx, CALC_CYCLE_OPTION option)
{
    CYCLE_PRE_INFO_T *c_info = (CYCLE_PRE_INFO_T *) ctx;
    c_info->cycle = (CYCLE_INFO_T){-1, -1, -1, -1};

    if (option == CALC_NONE)
    {
        return c_info->cycle;
    }

    calc_girth(c_info);

    if (option == CALC_GIRTH)
    {
        goto EXIT_GIRTH;
    }


    int n = c_info->edge->n;
    int m = c_info->edge->m;
    int g = c_info->cycle.girth;

    MATRIX_T *buf_uu = matrix_alloc(n, n);
    MATRIX_T *buf_uw = matrix_alloc(n, m);
    MATRIX_T *buf_ww = matrix_alloc(m, m);
    MATRIX_T *buf_wu = matrix_alloc(m, n);

    MATRIX_T *p_u_g = matrix_alloc(n, n);
    MATRIX_T *p_u_g1 = matrix_alloc(n, m);
    MATRIX_T *p_w_g = matrix_alloc(m, m);
    MATRIX_T *p_w_g1 = matrix_alloc(m, n);
    MATRIX_T *l_u_1g_2 = matrix_alloc(n, m);
    MATRIX_T *l_w_1g_2 = matrix_alloc(m, n);
    MATRIX_T *l_u_1_g = matrix_alloc(n, m);
    MATRIX_T *l_w_1_g = matrix_alloc(m, n);
    MATRIX_T *l_u_0_g2 = matrix_alloc(n, n);
    MATRIX_T *l_w_0_g2 = matrix_alloc(m, m);

    // calc P^U_g
    matrix_prod(c_info->c_p_u_w, c_info->edge_t, p_u_g);
    matrix_minus(p_u_g, c_info->l_u_u_g, p_u_g);
    matrix_minus(p_u_g, c_info->c_l_u_u, p_u_g);

    // calc P_W_g
    matrix_prod(c_info->c_p_w_u, c_info->edge, p_w_g);
    matrix_minus(p_w_g, c_info->l_w_w_g, p_w_g);
    matrix_minus(p_w_g, c_info->c_l_w_w, p_w_g);

    // calc L^U_{1, g}
    matrix_prod(c_info->edge, c_info->l_w_w_g, l_u_1_g);
    matrix_dot_coef(2, c_info->c_p_u_w, 1, c_info->edge, buf_uw);
    matrix_minus(l_u_1_g, buf_uw, l_u_1_g);

    // calc L^W_{1, g}
    matrix_prod(c_info->edge_t, c_info->l_u_u_g, l_w_1_g);
    matrix_dot_coef(2, c_info->c_p_w_u, 1, c_info->edge_t, buf_wu);
    matrix_minus(l_w_1_g, buf_wu, l_w_1_g);

    // calc L^U_{g - 1, 2}
    matrix_prod(c_info->edge, c_info->c_l_w_w, l_u_1g_2);
    matrix_dot_coef(1, c_info->c_p_u_w, 1, c_info->edge, buf_uw);
    matrix_minus(l_u_1g_2, buf_uw, l_u_1g_2);
    matrix_prod(c_info->l_u_u_m, c_info->c_l_u_w, buf_uw);
    matrix_minus(l_u_1g_2, buf_uw, l_u_1g_2);

    // calc L^W_{g - 1, 2}
    matrix_prod(c_info->edge_t, c_info->c_l_u_u, l_w_1g_2);
    matrix_dot_coef(1, c_info->c_p_w_u, 1, c_info->edge_t, buf_wu);
    matrix_minus(l_w_1g_2, buf_wu, l_w_1g_2);
    matrix_prod(c_info->l_w_w_m, c_info->c_l_w_u, buf_wu);
    matrix_minus(l_w_1g_2, buf_wu, l_w_1g_2);

    // calc P^U_g+1
    matrix_prod(p_u_g, c_info->edge, p_u_g1);
    matrix_minus(p_u_g1, l_u_1_g, p_u_g1);
    matrix_minus(p_u_g1, l_u_1g_2, p_u_g1);

    // calc P^W_g+1
    matrix_prod(p_w_g, c_info->edge_t, p_w_g1);
    matrix_minus(p_w_g1, l_w_1_g, p_w_g1);
    matrix_minus(p_w_g1, l_w_1g_2, p_w_g1);

    // calc L^U_{0, g+2}
    matrix_prod(p_u_g1, c_info->edge_t, l_u_0_g2);
    matrix_diag(l_u_0_g2);

    // calc L^W_{0, g+2}
    matrix_prod(p_w_g1, c_info->edge, l_w_0_g2);
    matrix_diag(l_w_0_g2);

    c_info->cycle.no_g2 = calc_cycle_count(l_u_0_g2, g + 2);

    if (option == CALC_GIRTH_TO_2)
    {
        goto EXIT_GIRTH_2;
    }

    MATRIX_T *l_u_2_g = matrix_alloc(n, n);
    MATRIX_T *l_w_2_g = matrix_alloc(m, m);

    MATRIX_T *l_u_g_2 = matrix_alloc(n, n);
    MATRIX_T *l_w_g_2 = matrix_alloc(m, m);

    MATRIX_T *p_u_g2 = matrix_alloc(n, n);
    MATRIX_T *p_w_g2 = matrix_alloc(m, m);

    MATRIX_T *l_u_1_g2 = matrix_alloc(n, m);
    MATRIX_T *l_w_1_g2 = matrix_alloc(m, n);

    MATRIX_T *l_u_3_g = matrix_alloc(n, m);
    MATRIX_T *l_w_3_g = matrix_alloc(m, n);

    MATRIX_T *l_u_g1_2 = matrix_alloc(n, m);
    MATRIX_T *l_w_g1_2 = matrix_alloc(m, n);

    MATRIX_T *p_u_g3 = matrix_alloc(n, m);
    MATRIX_T *p_w_g3 = matrix_alloc(m, n);

    MATRIX_T *l_u_0_g4 = matrix_alloc(n, n);
    MATRIX_T *l_w_0_g4 = matrix_alloc(m, m);

    // calc L^U_{2, g}
    matrix_prod(c_info->edge, l_w_1_g, l_u_2_g);
    matrix_z(l_u_2_g);
    if (g == 4)
    {
        matrix_combination(6, c_info->p_u_u, 3, buf_uu);
        matrix_minus(l_u_2_g, buf_uu, l_u_2_g);
    }

    // calc L^W_{2, g}
    matrix_prod(c_info->edge_t, l_u_1_g, l_w_2_g);
    matrix_z(l_w_2_g);
    if (g == 4)
    {
        matrix_combination(6, c_info->p_w_w, 3, buf_ww);
        matrix_minus(l_w_2_g, buf_ww, l_w_2_g);
    }


    // calc L^U_{g, 2}
    matrix_prod(c_info->edge, l_w_1g_2, l_u_g_2);
    matrix_z(l_u_g_2);
    matrix_prod(c_info->l_u_u_m, c_info->c_l_u_u, buf_uu);
    matrix_minus(l_u_g_2, buf_uu, l_u_g_2);
    matrix_dot_coef(1, c_info->c_p_u_u, 1, c_info->p_u_u, buf_uu);
    matrix_add(l_u_g_2, buf_uu, l_u_g_2);
    if (g == 4)
    {
        matrix_minus(l_u_g_2, c_info->p_u_u, l_u_g_2);
    }

    // calc L_W_{g, 2}
    matrix_prod(c_info->edge_t, l_u_1g_2, l_w_g_2);
    matrix_z(l_w_g_2);
    matrix_prod(c_info->l_w_w_m, c_info->c_l_w_w, buf_ww);
    matrix_minus(l_w_g_2, buf_ww, l_w_g_2);
    matrix_dot_coef(1, c_info->c_p_w_w, 1, c_info->p_w_w, buf_ww);
    matrix_add(l_w_g_2, buf_ww, l_w_g_2);
    if (g == 4)
    {
        matrix_minus(l_w_g_2, c_info->p_w_w, l_w_g_2);
    }

    // calc P^U_g+2
    matrix_prod(p_u_g1, c_info->edge_t, p_u_g2);
    matrix_minus(p_u_g2, l_u_0_g2, p_u_g2);
    matrix_minus(p_u_g2, l_u_2_g, p_u_g2);
    matrix_minus(p_u_g2, l_u_g_2, p_u_g2);

    // calc P^W_g+2
    matrix_prod(p_w_g1, c_info->edge, p_w_g2);
    matrix_minus(p_w_g2, l_w_0_g2, p_w_g2);
    matrix_minus(p_w_g2, l_w_2_g, p_w_g2);
    matrix_minus(p_w_g2, l_w_g_2, p_w_g2);

    // calc L^U_{1, g+2}
    matrix_prod(c_info->edge, l_w_0_g2, l_u_1_g2);
    matrix_dot_coef(2, p_u_g1, 1, c_info->edge, buf_uw);
    matrix_minus(l_u_1_g2, buf_uw, l_u_1_g2);
    if (g == 4)
    {
        matrix_combination(1, c_info->p_u_w, 2, buf_uw);
        matrix_dot_coef(2, buf_uw, 1, c_info->edge, buf_uw);
        matrix_minus(l_u_1_g2, buf_uw, l_u_1_g2);

        matrix_combination(1, c_info->p_w_w, 2, buf_ww);
        matrix_prod(c_info->edge, buf_ww, buf_uw);
        matrix_minus(buf_uw, c_info->p_u_w, buf_uw);
        matrix_dot_coef(2, buf_uw, 1, c_info->edge, buf_uw);
        matrix_add(l_u_1_g2, buf_uw, l_u_1_g2);

        matrix_combination(1, c_info->p_u_u, 2, buf_uu);
        matrix_prod(buf_uu, c_info->edge, buf_uw);
        matrix_minus(buf_uw, c_info->p_u_w, buf_uw);
        matrix_dot_coef(2, buf_uw, 1, c_info->edge, buf_uw);
        matrix_add(l_u_1_g2, buf_uw, l_u_1_g2);
    }


    // calc L^W_{1, g+2}
    matrix_prod(c_info->edge_t, l_u_0_g2, l_w_1_g2);
    matrix_dot_coef(2, p_w_g1, 1, c_info->edge_t, buf_wu);
    matrix_minus(l_w_1_g2, buf_wu, l_w_1_g2);
    if (g == 4)
    {
        matrix_combination(1, c_info->p_w_u, 2, buf_wu);
        matrix_dot_coef(2, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_minus(l_w_1_g2, buf_wu, l_w_1_g2);

        matrix_combination(1, c_info->p_u_u, 2, buf_uu);
        matrix_prod(c_info->edge_t, buf_uu, buf_wu);
        matrix_minus(buf_wu, c_info->p_w_u, buf_wu);
        matrix_dot_coef(2, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_add(l_w_1_g2, buf_wu, l_w_1_g2);

        matrix_combination(1, c_info->p_w_w, 2, buf_ww);
        matrix_prod(buf_ww, c_info->edge_t, buf_wu);
        matrix_minus(buf_wu, c_info->p_w_u, buf_wu);
        matrix_dot_coef(2, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_add(l_w_1_g2, buf_wu, l_w_1_g2);
    }


    // calc L^U_{3, g}
    matrix_prod(c_info->edge, l_w_2_g, l_u_3_g);
    matrix_prod(c_info->l_u_u_m, l_u_1_g, buf_uw);
    matrix_minus(l_u_3_g, buf_uw, l_u_3_g);
    if (g == 6)
    {
        matrix_combination(6, c_info->p_u_w, 3, buf_uw);
        matrix_minus(l_u_3_g, buf_uw, l_u_3_g);
    }

    if (g == 4)
    {
        matrix_combination(4, c_info->p_u_w, 2, buf_uw);
        matrix_dot_coef(1, buf_uw, 1, c_info->edge, buf_uw);
        matrix_minus(l_u_3_g, buf_uw, l_u_3_g);

        matrix_combination(1, c_info->p_w_w, 2, buf_ww);
        matrix_prod(c_info->edge, buf_ww, buf_uw);
        matrix_minus(buf_uw, c_info->p_u_w, buf_uw);
        matrix_dot_coef(6, buf_uw, 1, c_info->edge, buf_uw);
        matrix_add(l_u_3_g, buf_uw, l_u_3_g);

        matrix_combination(1, c_info->p_u_u, 2, buf_uu);
        matrix_prod(buf_uu, c_info->edge, buf_uw);
        matrix_minus(buf_uw, c_info->p_u_w, buf_uw);
        matrix_dot_coef(4, buf_uw, 1, c_info->edge, buf_uw);
        matrix_add(l_u_3_g, buf_uw, l_u_3_g);
    }

    // calc L^W_{3, g}
    matrix_prod(c_info->edge_t, l_u_2_g, l_w_3_g);
    matrix_prod(c_info->l_w_w_m, l_w_1_g, buf_wu);
    matrix_minus(l_w_3_g, buf_wu, l_w_3_g);
    if (g == 6)
    {
        matrix_combination(6, c_info->p_w_u, 3, buf_wu);
        matrix_minus(l_w_3_g, buf_wu, l_w_3_g);
    }
    if (g == 4)
    {
        matrix_combination(4, c_info->p_w_u, 2, buf_wu);
        matrix_dot_coef(1, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_minus(l_w_3_g, buf_wu, l_w_3_g);

        matrix_combination(1, c_info->p_u_u, 2, buf_uu);
        matrix_prod(c_info->edge_t, buf_uu, buf_wu);
        matrix_minus(buf_wu, c_info->p_w_u, buf_wu);
        matrix_dot_coef(6, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_add(l_w_3_g, buf_wu, l_w_3_g);

        matrix_combination(1, c_info->p_w_w, 2, buf_ww);
        matrix_prod(buf_ww, c_info->edge_t, buf_wu);
        matrix_minus(buf_wu, c_info->p_w_u, buf_wu);
        matrix_dot_coef(4, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_add(l_w_3_g, buf_wu, l_w_3_g);
    }

    // calc L^U_{g+1, 2}
    matrix_prod(c_info->edge, l_w_g_2, l_u_g1_2);
    matrix_prod(c_info->l_u_u_g, c_info->l_u_w, buf_uw);
    matrix_minus(l_u_g1_2, buf_uw, l_u_g1_2);
    matrix_dot_coef(1, p_u_g1, 1, c_info->edge, buf_uw);
    matrix_minus(l_u_g1_2, buf_uw, l_u_g1_2);
    matrix_prod(c_info->l_u_u_m, l_u_1g_2, buf_uw);
    matrix_minus(l_u_g1_2, buf_uw, l_u_g1_2);
    matrix_prod(c_info->c_p_u_w, c_info->l_w_w_m2, buf_uw);
    matrix_dot_coef(2, buf_uw, 1, c_info->edge, buf_uw);
    matrix_add(l_u_g1_2, buf_uw, l_u_g1_2);
    matrix_dot_coef(1, l_u_1g_2, 1, c_info->edge, buf_uw);
    matrix_add(l_u_g1_2, buf_uw, l_u_g1_2);
    matrix_dot_coef(2, c_info->c_p_u_w, 1, c_info->edge, buf_uw);
    matrix_add(l_u_g1_2, buf_uw, l_u_g1_2);
    if (g == 4)
    {
        matrix_combination(1, c_info->p_u_u, 2, buf_uu);
        matrix_prod(buf_uu, c_info->edge, buf_uw);
        matrix_minus(buf_uw, c_info->p_u_w, buf_uw);
        matrix_dot_coef(2, buf_uw, 1, c_info->edge, buf_uw);
        matrix_add(l_u_g1_2, buf_uw, l_u_g1_2);
    }

    // calc L^W_{g+1, 2}
    matrix_prod(c_info->edge_t, l_u_g_2, l_w_g1_2);
    matrix_prod(c_info->l_w_w_g, c_info->l_w_u, buf_wu);
    matrix_minus(l_w_g1_2, buf_wu, l_w_g1_2);
    matrix_dot_coef(1, p_w_g1, 1, c_info->edge_t, buf_wu);
    matrix_minus(l_w_g1_2, buf_wu, l_w_g1_2);
    matrix_prod(c_info->l_w_w_m, l_w_1g_2, buf_wu);
    matrix_minus(l_w_g1_2, buf_wu, l_w_g1_2);
    matrix_prod(c_info->c_p_w_u, c_info->l_u_u_m2, buf_wu);
    matrix_dot_coef(2, buf_wu, 1, c_info->edge_t, buf_wu);
    matrix_add(l_w_g1_2, buf_wu, l_w_g1_2);
    matrix_dot_coef(1, l_w_1g_2, 1, c_info->edge_t, buf_wu);
    matrix_add(l_w_g1_2, buf_wu, l_w_g1_2);
    matrix_dot_coef(2, c_info->c_p_w_u, 1, c_info->edge_t, buf_wu);
    matrix_add(l_w_g1_2, buf_wu, l_w_g1_2);
    if (g == 4)
    {
        matrix_combination(1, c_info->p_w_w, 2, buf_ww);
        matrix_prod(buf_ww, c_info->edge_t, buf_wu);
        matrix_minus(buf_wu, c_info->p_w_u, buf_wu);
        matrix_dot_coef(2, buf_wu, 1, c_info->edge_t, buf_wu);
        matrix_add(l_w_g1_2, buf_wu, l_w_g1_2);
    }

    // calc P^U_g+3
    matrix_prod(p_u_g2, c_info->edge, p_u_g3);
    matrix_minus(p_u_g3, l_u_1_g2, p_u_g3);
    matrix_minus(p_u_g3, l_u_3_g, p_u_g3);
    matrix_minus(p_u_g3, l_u_g1_2, p_u_g3);

    // calc P^W_g+3
    matrix_prod(p_w_g2, c_info->edge_t, p_w_g3);
    matrix_minus(p_w_g3, l_w_1_g2, p_w_g3);
    matrix_minus(p_w_g3, l_w_3_g, p_w_g3);
    matrix_minus(p_w_g3, l_w_g1_2, p_w_g3);


    // calc L^U_{0, g+4}
    matrix_prod(p_u_g3, c_info->edge_t, l_u_0_g4);
    matrix_diag(l_u_0_g4);

    // calc L^W_{0, g+4}
    matrix_prod(p_w_g3, c_info->edge, l_w_0_g4);
    matrix_diag(l_w_0_g4);

    c_info->cycle.no_g4 = calc_cycle_count(l_u_0_g4, g + 4);

    matrix_free(l_u_2_g);
    matrix_free(l_w_2_g);
    matrix_free(l_u_g_2);
    matrix_free(l_w_g_2);
    matrix_free(p_u_g2);
    matrix_free(p_w_g2);
    matrix_free(l_u_1_g2);
    matrix_free(l_w_1_g2);
    matrix_free(l_u_3_g);
    matrix_free(l_w_3_g);
    matrix_free(l_u_g1_2);
    matrix_free(l_w_g1_2);
    matrix_free(p_u_g3);
    matrix_free(p_w_g3);
    matrix_free(l_u_0_g4);
    matrix_free(l_w_0_g4);

    EXIT_GIRTH_2:
    matrix_free(buf_uu);
    matrix_free(buf_uw);
    matrix_free(buf_ww);
    matrix_free(buf_wu);
    matrix_free(p_u_g);
    matrix_free(p_u_g1);
    matrix_free(p_w_g);
    matrix_free(p_w_g1);
    matrix_free(l_u_1g_2);
    matrix_free(l_w_1g_2);
    matrix_free(l_u_1_g);
    matrix_free(l_w_1_g);
    matrix_free(l_u_0_g2);
    matrix_free(l_w_0_g2);

    EXIT_GIRTH:

    return c_info->cycle;
}

void cycle_release(void *ctx)
{
    CYCLE_PRE_INFO_T *c_info = (CYCLE_PRE_INFO_T *)ctx;
    matrix_free(c_info->edge);
    matrix_free(c_info->edge_t);
    matrix_free(c_info->p_u_u);
    matrix_free(c_info->p_u_w);
    matrix_free(c_info->p_w_w);
    matrix_free(c_info->p_w_u);
    matrix_free(c_info->l_u_u);
    matrix_free(c_info->l_u_w);
    matrix_free(c_info->l_w_w);
    matrix_free(c_info->l_w_u);
    matrix_free(c_info->c_p_u_u);
    matrix_free(c_info->c_p_u_w);
    matrix_free(c_info->c_p_w_w);
    matrix_free(c_info->c_p_w_u);
    matrix_free(c_info->c_l_u_u);
    matrix_free(c_info->c_l_u_w);
    matrix_free(c_info->c_l_w_w);
    matrix_free(c_info->c_l_w_u);
    matrix_free(c_info->l_u_u_m);
    matrix_free(c_info->l_w_w_m);
    matrix_free(c_info->l_u_u_m2);
    matrix_free(c_info->l_w_w_m2);
    matrix_free(c_info->l_u_u_g);
    matrix_free(c_info->l_w_w_g);
}