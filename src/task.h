#ifndef __TASK_H__
#define __TASK_H__

#include "data_structs.h"

typedef enum {
        TASK_AMC,
        TASK_LINE,
        TASK_ZEEMAN,
        TASK_CONT,
        TASK_COLDENS,
        N_TASK_TYPE
} TASK_TYPE;


static DatINode TASKS[] = {
        {"amc",TASK_AMC},
        {"line", TASK_LINE},
        {"zeeman", TASK_ZEEMAN},
        {"cont", TASK_CONT},
        {"coldens", TASK_COLDENS},
        {0, 0}
};

#endif