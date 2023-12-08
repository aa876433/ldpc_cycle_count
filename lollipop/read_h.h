#ifndef READ_H_H
#define READ_H_H

#include <stdint.h>

typedef struct PARITY_INFO_T_
{
    int blk;
    int row;
    int col;
    int **build_matirx;
} PARITY_INFO_T;

PARITY_INFO_T *get_parity_matrix_info(const char *file_name);
void release_parity_info(PARITY_INFO_T *parity_info);

#endif // READ_H_H
