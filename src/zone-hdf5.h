#ifndef __ZONE_HDF5_H__
#define __ZONE_HDF5_H__

#include <hdf5.h>

#define ZoneH5_KAPPLEN ((size_t)64)

/* This structure is used only by this function,
 * so make it private */


typedef struct ZoneH5_Record_Zone {
        int level;
        double max[3], min[3], cen[3];
        long nchildren;
        long naxes[3];
} ZoneH5_Record_Zone;

typedef struct ZoneH5_Record_Grid {
        int level;
        long pos;
        double max[3], min[3], cen[3];
        double n_H2, T_k, X_pH2, X_oH2, X_e, X_H, X_He, V_t;
        double v_cen[3];
        long nchildren;
        long naxes[3];
} ZoneH5_Record_Grid;

typedef struct ZoneH5_Record_Molec {
        double X_mol;
} ZoneH5_Record_Molec;

typedef struct ZoneH5_Record_Dust {
        double T_d, dust_to_gas;
        char kapp_d[ZoneH5_KAPPLEN];
} ZoneH5_Record_Dust;

typedef struct ZoneH5_Record_Polariz {
        double B_cen[3];
        double alpha;
} ZoneH5_Record_Polariz;

int ZoneH5_FwriteTable_Zone(hid_t h5f_id, const ZoneH5_Record_Zone *records);
int ZoneH5_FwriteTable_Grid(hid_t h5f_id, const ZoneH5_Record_Grid *records, size_t nrecord);
int ZoneH5_FwriteTable_Molec(hid_t h5f_id, const ZoneH5_Record_Molec *records, size_t nrecord);
int ZoneH5_FwriteTable_Dust(hid_t h5f_id, const ZoneH5_Record_Dust *records, size_t nrecord);
int ZoneH5_FwriteTable_Polariz(hid_t h5f_id, const ZoneH5_Record_Polariz *records, size_t nrecord);

int ZoneH5_FreadTable_Zone(     hid_t h5f_id, ZoneH5_Record_Zone        *records);
int ZoneH5_FreadTable_Grid(     hid_t h5f_id, ZoneH5_Record_Grid        *records);
int ZoneH5_FreadTable_Molec(    hid_t h5f_id, ZoneH5_Record_Molec       *records);
int ZoneH5_FreadTable_Dust(     hid_t h5f_id, ZoneH5_Record_Dust        *records);
int ZoneH5_FreadTable_Polariz(  hid_t h5f_id, ZoneH5_Record_Polariz     *records);
#endif
