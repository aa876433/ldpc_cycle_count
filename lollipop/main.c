#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "read_h.h"
#include "matrix.h"
#include "calc_cycle.h"

int main()
{
    PARITY_INFO_T *info = get_parity_matrix_info("h.txt");
    CALC_CYCLE_OPTION option = CALC_GIRTH_TO_4;
    void *ctx = cycle_init(info);
    CYCLE_INFO_T cycle = calc_cycle(ctx, option);
    printf("girth %d\n", cycle.girth);
    printf("cycle_count %d, %d, %d\n", cycle.no_g, cycle.no_g2, cycle.no_g4);
    cycle_release(ctx);
    return 0;
}
