#ifndef __ZONE_HDF5_H__
#define __ZONE_HDF5_H__

#include <hdf5.h>

#define ZoneH5_GEOMLEN ((size_t)6)
#define ZoneH5_KAPPLEN ((size_t)64)

/* This structure is used only by this function,
 * so make it private */
typedef struct ZoneH5_Record {
	int level;
	long pos;
	char geom[ZoneH5_GEOMLEN];
	double max[3], min[3], cen[3];
	double n_H2, T_k, X_mol, X_pH2, X_oH2, X_e, X_H, X_He, V_t;
	double vedge[6][3];
	double v_cen[3],b_cen[3];
	double ds; /* Average path length */
	long nchildren;
	long naxes[3];
	double T_d;
	char kapp_d[ZoneH5_KAPPLEN];
	double T_ff;
	char kapp_ff[ZoneH5_KAPPLEN];
	double T_bb;
} ZoneH5_Record;

int ZoneH5_FwriteTable(hid_t h5f_id, const char *name, const ZoneH5_Record *records, size_t nrecord);
int ZoneH5_FreadTable(hid_t h5f_id, const char *name, ZoneH5_Record *records);

#endif