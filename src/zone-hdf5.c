#include <assert.h>
#include "memory.h"
#include "error.h"
#include "zone.h"
#include "zone-hdf5.h"
#include "hdf5_hl.h"

/* Size of record */
static size_t record_size_zone = sizeof(ZoneH5_Record_Zone);
static size_t record_size_grid = sizeof(ZoneH5_Record_Grid);
static size_t record_size_molec = sizeof(ZoneH5_Record_Molec);
static size_t record_size_dust = sizeof(ZoneH5_Record_Dust);
static size_t record_size_polariz = sizeof(ZoneH5_Record_Polariz);

/* CHK when format change */

#define NFIELDS_ZONE ((hsize_t)6)
#define NFIELDS_GRID ((hsize_t)16)
#define NFIELDS_MOLEC ((hsize_t)1)
#define NFIELDS_DUST ((hsize_t)3)
#define NFIELDS_POLARIZ ((hsize_t)2)


/* CHK when format change */
/* Offsets of each field */
static size_t field_offsets_zone[NFIELDS_ZONE] = {
        HOFFSET(ZoneH5_Record_Zone, level),
        HOFFSET(ZoneH5_Record_Zone, max),
        HOFFSET(ZoneH5_Record_Zone, min),
        HOFFSET(ZoneH5_Record_Zone, cen),
        HOFFSET(ZoneH5_Record_Zone, nchildren),
        HOFFSET(ZoneH5_Record_Zone, naxes)
};
static size_t field_offsets_grid[NFIELDS_GRID] = {
	HOFFSET(ZoneH5_Record_Grid, level),
	HOFFSET(ZoneH5_Record_Grid, pos),
	HOFFSET(ZoneH5_Record_Grid, max),
	HOFFSET(ZoneH5_Record_Grid, min),
	HOFFSET(ZoneH5_Record_Grid, cen),
	HOFFSET(ZoneH5_Record_Grid, n_H2),
	HOFFSET(ZoneH5_Record_Grid, T_k),
	HOFFSET(ZoneH5_Record_Grid, X_pH2),
	HOFFSET(ZoneH5_Record_Grid, X_oH2),
	HOFFSET(ZoneH5_Record_Grid, X_e),
	HOFFSET(ZoneH5_Record_Grid, X_H),
	HOFFSET(ZoneH5_Record_Grid, X_He),
	HOFFSET(ZoneH5_Record_Grid, V_t),
	HOFFSET(ZoneH5_Record_Grid, v_cen),
	HOFFSET(ZoneH5_Record_Grid, nchildren),
	HOFFSET(ZoneH5_Record_Grid, naxes)
};
static size_t field_offsets_molec[NFIELDS_MOLEC] = {
        HOFFSET(ZoneH5_Record_Molec, X_mol)
};
static size_t field_offsets_dust[NFIELDS_DUST] = {
        HOFFSET(ZoneH5_Record_Dust, T_d),
        HOFFSET(ZoneH5_Record_Dust, kapp_d),
        HOFFSET(ZoneH5_Record_Dust, dust_to_gas)
};
static size_t field_offsets_polariz[NFIELDS_POLARIZ] = {
        HOFFSET(ZoneH5_Record_Polariz, B_cen),
        HOFFSET(ZoneH5_Record_Polariz, alpha)
};

/* Dummy record for calculating sizes */
static ZoneH5_Record_Zone       dummy_record_zone;
static ZoneH5_Record_Grid       dummy_record_grid;
static ZoneH5_Record_Molec      dummy_record_molec;
static ZoneH5_Record_Dust       dummy_record_dust;
static ZoneH5_Record_Polariz    dummy_record_polariz;

/* Sizes of each field */
static size_t field_sizes_zone[NFIELDS_ZONE] = {
        sizeof(dummy_record_zone.level),
        sizeof(dummy_record_zone.max),
        sizeof(dummy_record_zone.min),
        sizeof(dummy_record_zone.cen),
        sizeof(dummy_record_zone.nchildren),
        sizeof(dummy_record_zone.naxes)
};
/* Sizes of each field */
static size_t field_sizes_grid[NFIELDS_GRID] = {
	sizeof(dummy_record_grid.level),
	sizeof(dummy_record_grid.pos),
	sizeof(dummy_record_grid.max),
	sizeof(dummy_record_grid.min),
	sizeof(dummy_record_grid.cen),
	sizeof(dummy_record_grid.n_H2),
	sizeof(dummy_record_grid.T_k),
	sizeof(dummy_record_grid.X_pH2),
	sizeof(dummy_record_grid.X_oH2),
	sizeof(dummy_record_grid.X_e),
	sizeof(dummy_record_grid.X_H),
	sizeof(dummy_record_grid.X_He),
	sizeof(dummy_record_grid.V_t),
	sizeof(dummy_record_grid.v_cen),
	sizeof(dummy_record_grid.nchildren),
	sizeof(dummy_record_grid.naxes)
};
static size_t field_sizes_molec[NFIELDS_MOLEC] = {
        sizeof(dummy_record_molec.X_mol)
};
static size_t field_sizes_dust[NFIELDS_DUST] = {
        sizeof(dummy_record_dust.T_d),
        sizeof(dummy_record_dust.kapp_d),
        sizeof(dummy_record_dust.dust_to_gas)
};
static size_t field_sizes_polariz[NFIELDS_POLARIZ] = {
        sizeof(dummy_record_polariz.B_cen),
        sizeof(dummy_record_polariz.alpha)
};


/* Field names */
static const char *field_names_zone[NFIELDS_ZONE] = {
        "LEVEL",
        "X_max",
        "X_min",
        "X_cen",
        "NCHILDREN",
        "NAXES"
};
static const char *field_names_grid[NFIELDS_GRID] = {
        "LEVEL",
        "POS",
        "X_max",
        "X_min",
        "X_cen",
        "n_H2",
        "T_k",
        "X_pH2",
        "X_oH2",
        "X_e",
        "X_H",
        "X_He",
        "V_t",
        "V_cen",
        "NCHILDREN",
        "NAXES"
};
static const char *field_names_molec[NFIELDS_MOLEC] = {
        "X_mol"
};
static const char *field_names_dust[NFIELDS_DUST] = {
        "T_d",
        "kapp_d",
        "dust_to_gas"
};
static const char *field_names_polariz[NFIELDS_POLARIZ] = {
        "B_cen",
        "alpha"
};

/*----------------------------------------------------------------------------*/

int ZoneH5_FwriteTable_Zone(hid_t h5f_id, const ZoneH5_Record_Zone *records)
{
        herr_t hstatus;

        hid_t field_types[NFIELDS_ZONE];
        hid_t vec3d_type, vec3s_type, vel_type;
        hsize_t chunk_size = 10, vec3_size = 3, vel_size[2] = {6, 3};
        int *fill_data = NULL, compress = 0;

        /* Create array data type */
        vec3d_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &vec3_size);
        vel_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, vel_size);
        vec3s_type = H5Tarray_create(H5T_NATIVE_LONG, 1, &vec3_size);

        field_types[0] = H5T_NATIVE_INT;
        field_types[1] = vec3d_type;
        field_types[2] = vec3d_type;
        field_types[3] = vec3d_type;
        field_types[4] = H5T_NATIVE_LONG;
        field_types[5] = vec3s_type;

        /* Make the table */
        hstatus = H5TBmake_table(
                "Zone table",
                h5f_id,
                "ZONE",
                NFIELDS_ZONE,
                (hsize_t)1,
                record_size_zone,
                field_names_zone,
                field_offsets_zone,
                field_types,
                chunk_size,
                fill_data,
                compress,
                records
        );

        H5Tclose(vec3d_type);
        H5Tclose(vec3s_type);

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 GRID table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/

int ZoneH5_FwriteTable_Grid(hid_t h5f_id, const ZoneH5_Record_Grid *records, size_t nrecord)
{
        herr_t hstatus;

        hid_t field_types[NFIELDS_GRID];
        hid_t vec3d_type, vec3s_type, vel_type;
        hsize_t chunk_size = 10, vec3_size = 3, vel_size[2] = {6, 3};
        int *fill_data = NULL, compress = 0;

        /* Create array data type */
        vec3d_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &vec3_size);
        vel_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, vel_size);
        vec3s_type = H5Tarray_create(H5T_NATIVE_LONG, 1, &vec3_size);

        field_types[0] = H5T_NATIVE_INT;
        field_types[1] = H5T_NATIVE_LONG;
        field_types[2] = vec3d_type;
        field_types[3] = vec3d_type;
        field_types[4] = vec3d_type;
        field_types[5] = H5T_NATIVE_DOUBLE;
        field_types[6] = H5T_NATIVE_DOUBLE;
        field_types[7] = H5T_NATIVE_DOUBLE;
        field_types[8] = H5T_NATIVE_DOUBLE;
        field_types[9] = H5T_NATIVE_DOUBLE;
        field_types[10] = H5T_NATIVE_DOUBLE;
        field_types[11] = H5T_NATIVE_DOUBLE;
        field_types[12] = vel_type;
        field_types[13] = vec3d_type;     
        field_types[14] = H5T_NATIVE_LONG;
        field_types[15] = vec3s_type;


        /* Make the table */
        hstatus = H5TBmake_table(
                "Grid table",
                h5f_id,
                "GRID",
                NFIELDS_GRID,
                (hsize_t)nrecord,
                record_size_grid,
                field_names_grid,
                field_offsets_grid,
                field_types,
                chunk_size,
                fill_data,
                compress,
                records
        );

        H5Tclose(vec3d_type);
        H5Tclose(vec3s_type);
        H5Tclose(vel_type);

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 GRID table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/

int ZoneH5_FwriteTable_Molec(hid_t h5f_id, const ZoneH5_Record_Molec *records, size_t nrecord)
{
        herr_t hstatus;

        hid_t field_types[NFIELDS_MOLEC];
        hsize_t chunk_size = 10;
        int *fill_data = NULL, compress = 0;

        field_types[0] = H5T_NATIVE_DOUBLE;
        
        /* Make the table */
        hstatus = H5TBmake_table(
                "Molec table",
                h5f_id,
                "MOLEC",
                NFIELDS_MOLEC,
                (hsize_t)nrecord,
                record_size_molec,
                field_names_molec,
                field_offsets_molec,
                field_types,
                chunk_size,
                fill_data,
                compress,
                records
        );

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 GRID table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/

int ZoneH5_FwriteTable_Dust(hid_t h5f_id, const ZoneH5_Record_Dust *records, size_t nrecord)
{
        herr_t hstatus;

        hid_t field_types[NFIELDS_DUST];
        hid_t kapp_type;
        hsize_t chunk_size = 10;
        int *fill_data = NULL, compress = 0;

        /* dust column */
        kapp_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(kapp_type, (size_t)ZoneH5_KAPPLEN);

        field_types[0] = H5T_NATIVE_DOUBLE;
        field_types[1] = kapp_type;
        field_types[2] = H5T_NATIVE_DOUBLE;

        /* Make the table */
        hstatus = H5TBmake_table(
                "Dust table",
                h5f_id,
                "DUST",
                NFIELDS_DUST,
                (hsize_t)nrecord,
                record_size_dust,
                field_names_dust,
                field_offsets_dust,
                field_types,
                chunk_size,
                fill_data,
                compress,
                records
        );
        
        H5Tclose(kapp_type);

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 GRID table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/

int ZoneH5_FwriteTable_Polariz(hid_t h5f_id, const ZoneH5_Record_Polariz *records, size_t nrecord)
{
        herr_t hstatus;

        hid_t field_types[NFIELDS_POLARIZ];
        hid_t vec3d_type;
        hsize_t chunk_size = 10, vec3_size = 3;
        int *fill_data = NULL, compress = 0;

        /* Create array data type */
        vec3d_type = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &vec3_size);

        field_types[0] = vec3d_type;
        field_types[1] = H5T_NATIVE_DOUBLE;

        /* Make the table */
        hstatus = H5TBmake_table(
                "Polariz table",
                h5f_id,
                "POLARIZ",
                NFIELDS_POLARIZ,
                (hsize_t)nrecord,
                record_size_polariz,
                field_names_polariz,
                field_offsets_polariz,
                field_types,
                chunk_size,
                fill_data,
                compress,
                records
        );

        H5Tclose(vec3d_type);

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 GRID table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/
int ZoneH5_FreadTable_Zone(hid_t h5f_id, ZoneH5_Record_Zone *records)
{
        herr_t hstatus;

        /* read the table */
        hstatus = H5TBread_table(
                h5f_id, 
                "ZONE", 
                record_size_zone, 
                field_offsets_zone, 
                field_sizes_zone, 
                records
        );

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 ZONE table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/
int ZoneH5_FreadTable_Grid(hid_t h5f_id, ZoneH5_Record_Grid *records)
{
        herr_t hstatus;

        /* read the table */
        hstatus = H5TBread_table(
                h5f_id, 
                "GRID", 
                record_size_grid, 
                field_offsets_grid, 
                field_sizes_grid, 
                records
        );

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 GRID table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/
int ZoneH5_FreadTable_Molec(hid_t h5f_id, ZoneH5_Record_Molec *records)
{
        herr_t hstatus;

        /* read the table */
        hstatus = H5TBread_table(
                h5f_id, 
                "MOLEC", 
                record_size_molec, 
                field_offsets_molec, 
                field_sizes_molec, 
                records
        );

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 MOLEC table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/
int ZoneH5_FreadTable_Dust(hid_t h5f_id, ZoneH5_Record_Dust *records)
{
        herr_t hstatus;

        /* read the table */
        hstatus = H5TBread_table(
                h5f_id, 
                "DUST", 
                record_size_dust, 
                field_offsets_dust, 
                field_sizes_dust, 
                records
        );

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 DUST table");
        else
                return 0;
}

/*----------------------------------------------------------------------------*/
int ZoneH5_FreadTable_Polariz(hid_t h5f_id, ZoneH5_Record_Polariz *records)
{
        herr_t hstatus;

        /* read the table */
        hstatus = H5TBread_table(
                h5f_id, 
                "POLARIZ", 
                record_size_polariz, 
                field_offsets_polariz, 
                field_sizes_polariz, 
                records
        );

        if(hstatus < 0)
                return Err_SETSTRING("Error reading HDF5 POLARIZ table");
        else
                return 0;
}
