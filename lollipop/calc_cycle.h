#ifndef CALC_CYCLE_H
#define CALC_CYCLE_H

typedef struct CYCLE_INFO_T
{
    int girth;
    int no_g;
    int no_g2;
    int no_g4;
} CYCLE_INFO_T;

typedef enum {
    CALC_NONE = 0,
    CALC_GIRTH,
    CALC_GIRTH_TO_2,
    CALC_GIRTH_TO_4
} CALC_CYCLE_OPTION;

void *cycle_init(void *p_info);
CYCLE_INFO_T calc_cycle(void *ctx, CALC_CYCLE_OPTION option);
void cycle_release(void *ctx);

#endif // CALC_CYCLE_H
