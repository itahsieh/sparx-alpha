#include <assert.h>
#include "memory.h"
#include "error.h"
#include "zone.h"
#include "zone-hdf5.h"
#include "hdf5_hl.h"

/* Size of record */
static size_t record_size = sizeof(ZoneH5_Record);

/* CHK when format change */
#define NFIELDS ((hsize_t)26)

/* CHK when format change */
/* Offsets of each field */
static size_t field_offsets[NFIELDS] = {
	HOFFSET(ZoneH5_Record, level),
	HOFFSET(ZoneH5_Record, pos),
	HOFFSET(ZoneH5_Record, geom),
	HOFFSET(ZoneH5_Record, max),
	HOFFSET(ZoneH5_Record, min),
	HOFFSET(ZoneH5_Record, cen),
	HOFFSET(ZoneH5_Record, n_H2),
	HOFFSET(ZoneH5_Record, T_k),
	HOFFSET(ZoneH5_Record, X_mol),
	HOFFSET(ZoneH5_Record, X_pH2),
	HOFFSET(ZoneH5_Record, X_oH2),
	HOFFSET(ZoneH5_Record, X_e),
	HOFFSET(ZoneH5_Record, X_H),
	HOFFSET(ZoneH5_Record, X_He),
	HOFFSET(ZoneH5_Record, V_t),
	HOFFSET(ZoneH5_Record, vedge),
	HOFFSET(ZoneH5_Record, v_cen),
	HOFFSET(ZoneH5_Record, b_cen),
	HOFFSET(ZoneH5_Record, ds),
	HOFFSET(ZoneH5_Record, nchildren),
	HOFFSET(ZoneH5_Record, naxes),
	HOFFSET(ZoneH5_Record, T_d),
	HOFFSET(ZoneH5_Record, kapp_d),
	HOFFSET(ZoneH5_Record, T_ff),
	HOFFSET(ZoneH5_Record, kapp_ff),
	HOFFSET(ZoneH5_Record, T_bb)
};

/* Dummy record for calculating sizes */
static ZoneH5_Record dummy_record;

/* Sizes of each field */
static size_t field_sizes[NFIELDS] = {
	sizeof(dummy_record.level),
	sizeof(dummy_record.pos),
	sizeof(dummy_record.geom),
	sizeof(dummy_record.max),
	sizeof(dummy_record.min),
	sizeof(dummy_record.cen),
	sizeof(dummy_record.n_H2),
	sizeof(dummy_record.T_k),
	sizeof(dummy_record.X_mol),
	sizeof(dummy_record.X_pH2),
	sizeof(dummy_record.X_oH2),
	sizeof(dummy_record.X_e),
	sizeof(dummy_record.X_H),
	sizeof(dummy_record.X_He),
	sizeof(dummy_record.V_t),
	sizeof(dummy_record.vedge),
	sizeof(dummy_record.v_cen),
	sizeof(dummy_record.b_cen),
	sizeof(dummy_record.ds),
	sizeof(dummy_record.nchildren),
	sizeof(dummy_record.naxes),
	sizeof(dummy_record.T_d),
	sizeof(dummy_record.kapp_d),
	sizeof(dummy_record.T_ff),
	sizeof(dummy_record.kapp_ff),
	sizeof(dummy_record.T_bb)
};

/* Field names */
static const char *field_names[NFIELDS] = {
	"LEVEL",
	"POS",
	"geom",
	"X_max",
	"X_min",
	"X_cen",
	"n_H2",
	"T_k",
	"X_mol",
	"X_pH2",
	"X_oH2",
	"X_e",
	"X_H",
	"X_He",
	"V_t",
	"V_edge",
	"V_cen",
	"B_cen",
	"ds",
	"NCHILDREN",
	"NAXES",
	"T_d",
	"kapp_d",
	"T_ff",
	"kapp_ff",
	"T_bb"
};

/*----------------------------------------------------------------------------*/

int ZoneH5_FwriteTable(hid_t h5f_id, const char *name, const ZoneH5_Record *records, size_t nrecord)
{
	herr_t hstatus;

	hid_t field_types[NFIELDS], vec3d_type, vec3s_type, vel_type, geom_type, kapp_type;
	hsize_t chunk_size = 10, vec3_size = 3, vel_size[2] = {6, 3};
	int *fill_data = NULL, compress = 0;

	/* Create array data type */
	vec3d_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &vec3_size);
	vel_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, vel_size);
	vec3s_type = H5Tarray_create(H5T_NATIVE_LONG, 1, &vec3_size);

	/* geom column */
	geom_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(geom_type, (size_t)ZoneH5_GEOMLEN);

	/* dust column */
	kapp_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(kapp_type, (size_t)ZoneH5_KAPPLEN);

	field_types[0] = H5T_NATIVE_INT;
	field_types[1] = H5T_NATIVE_LONG;
	field_types[2] = geom_type;
	field_types[3] = vec3d_type;
	field_types[4] = vec3d_type;
	field_types[5] = vec3d_type;
	field_types[6] = H5T_NATIVE_DOUBLE;
	field_types[7] = H5T_NATIVE_DOUBLE;
	field_types[8] = H5T_NATIVE_DOUBLE;
	field_types[9] = H5T_NATIVE_DOUBLE;
	field_types[10] = H5T_NATIVE_DOUBLE;
	field_types[11] = H5T_NATIVE_DOUBLE;
	field_types[12] = H5T_NATIVE_DOUBLE;
	field_types[13] = H5T_NATIVE_DOUBLE;
	field_types[14] = H5T_NATIVE_DOUBLE;
	field_types[15] = vel_type;
	field_types[16] = vec3d_type;
	field_types[17] = vec3d_type;
	field_types[18] = H5T_NATIVE_DOUBLE;
	field_types[19] = H5T_NATIVE_LONG;
	field_types[20] = vec3s_type;
	field_types[21] = H5T_NATIVE_DOUBLE;
	field_types[22] = kapp_type;
	field_types[23] = H5T_NATIVE_DOUBLE;
	field_types[24] = kapp_type;
	field_types[25] = H5T_NATIVE_DOUBLE;

// 	field_types[17] = H5T_NATIVE_DOUBLE;
// 	field_types[18] = H5T_NATIVE_LONG;
// 	field_types[19] = vec3s_type;
// 	field_types[20] = H5T_NATIVE_DOUBLE;
// 	field_types[21] = kapp_type;
// 	field_types[22] = H5T_NATIVE_DOUBLE;
// 	field_types[23] = kapp_type;
// 	field_types[24] = H5T_NATIVE_DOUBLE;

	/* Make the table */
	hstatus = H5TBmake_table(
		"Grid table",
		h5f_id,
		name,
		NFIELDS,
		(hsize_t)nrecord,
		record_size,
		field_names,
		field_offsets,
		field_types,
		chunk_size,
		fill_data,
		compress,
		records
	);

	H5Tclose(vec3d_type);
	H5Tclose(vec3s_type);
	H5Tclose(vel_type);
	H5Tclose(geom_type);
	H5Tclose(kapp_type);

	if(hstatus < 0)
		return Err_SETSTRING("Error reading HDF5 `%s' table", name);
	else
		return 0;
}

/*----------------------------------------------------------------------------*/

int ZoneH5_FreadTable(hid_t h5f_id, const char *name, ZoneH5_Record *records)
{
	herr_t hstatus;

	/* read the table */
	hstatus = H5TBread_table(h5f_id, name, record_size, field_offsets, field_sizes, records);

	if(hstatus < 0)
		return Err_SETSTRING("Error reading HDF5 `%s' table", name);
	else
		return 0;
}

#undef NFIELDS

/*----------------------------------------------------------------------------*/
