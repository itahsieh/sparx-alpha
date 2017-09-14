#ifndef __TASK_H__
#define __TASK_H__

#include "data_structs.h"

typedef enum {
        TASK_AMC,
        TASK_LINE,
        TASK_ZEEMAN,
        TASK_CONT,
        TASK_COLDENS,
        TASK_ZEEMANCTB,
        TASK_LINECTB,
        TASK_CONTCTB,
        TASK_MODEL2VTK,
        N_TASK_TYPE
} TASK_TYPE;


static DatINode TASKS[] = {
        {"amc",     TASK_AMC},
        {"line",    TASK_LINE},
        {"zeeman",  TASK_ZEEMAN},
        {"cont",    TASK_CONT},
        {"coldens", TASK_COLDENS},
        {"zeemanctb", TASK_ZEEMANCTB},
        {"linectb", TASK_LINECTB},
        {"contctb", TASK_CONTCTB},
        {"model2vtk",     TASK_MODEL2VTK},
        {0, 0}
};

#endif
