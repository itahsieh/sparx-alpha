#ifndef __UNIT_H__
#define __UNIT_H__

#include "data_structs.h"


typedef enum {
    UNIT_K,
    UNIT_JYPX,
    UNIT_MKS,
    UNIT_CGS,
    UNIT_JYPC,
    UNIT_KPC,
    N_UNIT
}UNIT_TYPE;

static DatINode UNITS[] = {
    {"K", UNIT_K},
    {"JY/PIXEL", UNIT_JYPX},
    {"MKS", UNIT_MKS},
    {"CGS", UNIT_CGS},
    {"JY/PC", UNIT_JYPC},
    {"K/PC", UNIT_KPC},
    {0, 0}
};

#endif
