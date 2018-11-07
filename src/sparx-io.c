#include "sparx.h"

/*----------------------------------------------------------------------------*/

void SpIO_SetTaskName(const char *str)
/* Point Sp_parm.task to a static buffer within this function */
{
	static char buffer[BUFSIZ];

	strncpy(buffer, str, (size_t)BUFSIZ);

	Sp_parm.task = buffer;

	return;
}

/*----------------------------------------------------------------------------*/

void SpIO_Print(const char *file, int line, const char *func, const char *format, ...)
{
	va_list ap;

	#ifdef HAVE_MPI
	if(Sp_MPIRANK != 0)
		return;
	#endif

	if(Sp_parm.debug)
		fprintf(Sp_parm.out_fp, "%s:%d In function `%s':\n", file, line, func);

	if(Sp_parm.task)
		fprintf(Sp_parm.out_fp, "%s: ", Sp_parm.task);
	else
		fprintf(Sp_parm.out_fp, "%s: ", Sp_parm.prog);

	va_start(ap, format);
	vfprintf(Sp_parm.out_fp, format, ap);
	fflush(Sp_parm.out_fp);
	va_end(ap);

	return;
}

/*----------------------------------------------------------------------------*/

void SpIO_Printf(const char *file, int line, const char *func, const char *format, ...)
{
	va_list ap;

	#ifdef HAVE_MPI
	if(Sp_MPIRANK != 0)
		return;
	#endif

	if(Sp_parm.debug)
		fprintf(Sp_parm.out_fp, "%s:%d In function `%s':\n", file, line, func);

	va_start(ap, format);
	vfprintf(Sp_parm.out_fp, format, ap);
	fflush(Sp_parm.out_fp);
	va_end(ap);

	return;
}

/*----------------------------------------------------------------------------*/

void SpIO_Pwarn(const char *file, int line, const char *func, const char *format, ...)
{
	va_list ap;

	if(Sp_parm.debug)
		fprintf(Sp_parm.err_fp, "%s:%d In function `%s':\n", file, line, func);

	if(Sp_parm.task)
		fprintf(Sp_parm.err_fp, "%s warning: ", Sp_parm.task);
	else
		fprintf(Sp_parm.err_fp, "%s warning: ", Sp_parm.prog);

	va_start(ap, format);
	vfprintf(Sp_parm.err_fp, format, ap);
	fflush(Sp_parm.err_fp);
	va_end(ap);

	return;
}

/*----------------------------------------------------------------------------*/

void SpIO_Perror(const char *file, int line, const char *func, const char *format, ...)
{
	va_list ap;

	if(Sp_parm.debug)
		fprintf(Sp_parm.err_fp, "%s:%d In function `%s':\n", file, line, func);

	if(Sp_parm.task)
		fprintf(Sp_parm.err_fp, "%s error: ", Sp_parm.task);
	else
		fprintf(Sp_parm.err_fp, "%s error: ", Sp_parm.prog);

	va_start(ap, format);
	vfprintf(Sp_parm.err_fp, format, ap);
	fflush(Sp_parm.err_fp);
	va_end(ap);

	return;
}

/*----------------------------------------------------------------------------*/

SpFile *SpIO_OpenFile(const char *fname, int mode)
{
	SpFile *sfp = 0;
	hid_t file_id = 0;
	
	switch(mode) {
		case Sp_NEW:
			file_id = H5Fcreate(fname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
			break;

		case Sp_OLD:
			file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
			break;

		case Sp_TRUNC:
			file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			break;

		default: /* No such mode */
			Deb_ASSERT(0);
	}
	
	if(file_id >= 0) {
		sfp = Mem_CALLOC(1, sfp);
		sfp->name = Mem_STRDUP(fname);
		sfp->h5f_id = file_id;
	}
	else {
		printf("Error opening file `%s'\n", fname);
	}

	return sfp;
}

/*----------------------------------------------------------------------------*/

int SpIO_OpenFile2(const char *fname, int mode, SpFile **fp)
{
	int sts = 0;

	*fp = SpIO_OpenFile(fname, mode);

	if(!(*fp)) {
		printf( "Error opening SPARX file '%s'\n", fname);
		sts = 1;
	}

	return sts;
}

/*----------------------------------------------------------------------------*/

void SpIO_CloseFile(SpFile *sfp)
{
	if(sfp->name)
		free(sfp->name);

	if(sfp->h5f_id >= 0)
		H5Fclose(sfp->h5f_id);

	free(sfp);

	return;
}

/*----------------------------------------------------------------------------*/

int SpIO_FwriteModel(SpFile *sfp, SpModel model)
{
	int status = 0;
	
	herr_t hstatus;

        /* Write the version of format */
        if(!status) {
                hstatus = H5LTset_attribute_string(sfp->h5f_id, "/", "format", "SPARX format v3");
                if(hstatus < 0)
                        status = 1;
        }

        /* Write the type of coordinat system */
        if(!status) {
                DatINode *geom = Dat_IList_IdxLookup(GEOM_TYPES, model.parms.geom);
                Deb_ASSERT(geom != NULL);
                hstatus = H5LTset_attribute_string(sfp->h5f_id, "/", "geom", geom->name);
                if(hstatus < 0)
                        status = 1;
        }

	/* Write T_cmb */
	if(!status) {
		hstatus = H5LTset_attribute_double(sfp->h5f_id, "/", "T_cmb", &model.parms.T_cmb, (size_t)1);
		if(hstatus < 0)
			status = 1;
	}

	/* Write T_in */
        if(!status) {
                hstatus = H5LTset_attribute_double(sfp->h5f_id, "/", "T_in", &model.parms.T_in, (size_t)1);
                if(hstatus < 0)
                        status = 1;
        }
        

        /* Write molecule name */
        if(!status) {
                Molec *mol = model.parms.mol;
                hstatus = H5LTset_attribute_string(sfp->h5f_id, "/", "molec", mol ? mol->name : "");
                if(hstatus < 0)
                        status = 1;
        }
	
	/* Write pops-switch */
        if(!status) {
                hstatus = H5LTset_attribute_int(sfp->h5f_id, "/", "pops", &model.parms.pops, (size_t)1);
                if(hstatus < 0)
                        status = 1;
        }

        /* Write dust-switch */
        if(!status) {
                hstatus = H5LTset_attribute_int(sfp->h5f_id, "/", "dust", &model.parms.dust, (size_t)1);
                if(hstatus < 0)
                        status = 1;
        }
        /* Write polariz-switch */
        if(!status) {
                hstatus = H5LTset_attribute_int(sfp->h5f_id, "/", "polariz", &model.parms.polariz, (size_t)1);
                
                
                if(hstatus < 0)
                        status = 1;
        }
        /* Write polariz-switch */
        if(!status) {
                hstatus = H5LTset_attribute_double(sfp->h5f_id, "/", "z", &model.parms.z, (size_t)1);
                if(hstatus < 0)
                        status = 1;
        }

        
        /* Write Outer_Source */
        if(!status) {
            int nSource = model.parms.Outer_Source;
            hstatus = H5LTset_attribute_int(sfp->h5f_id, "/", "Outer_Source", &nSource, (size_t)1);
            if(hstatus < 0)
                status = 1;
        
            if (nSource != 0){
                ZoneH5_Record_Source *data; 
                /* Allocate table data */ 
                data = Mem_CALLOC( (size_t)nSource, data); 
                /* Load data values */ 
                for(size_t i = 0; i < (size_t)nSource; i++) 
                    data[i] = SpIO_SourceToH5Record( &model.parms.source[i] ); 
                /* Create and write grid data */ 
                if(!status) 
                    status = ZoneH5_FwriteTable_Source( sfp->h5f_id, data, (size_t)nSource); 
                /* Free table data */ 
                free(data); 
            }
        }

	/* Write grid */
	if(!status)
		status = SpIO_H5WriteGrid(sfp->h5f_id, model.grid, &model.parms);        
        
	if(status)
		status = printf("Error writing model to file `%s'\n", sfp->name);

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_FreadModel(const SpFile *sfp, const SpFile *popsfp, SpModel *model, int *read_pops)
{
	int status = 0;

	Deb_ASSERT(model->parms.mol == NULL);
	Deb_ASSERT(model->grid == NULL);
#if 1
        /* pre-check format version */
        if(!status){
                const char *default_format = "SPARX format v3";
                char *format = NULL;
                
                status = SpIO_H5GetAttribute_string(sfp->h5f_id, "/", "format", &format);

                if(strncmp(format, default_format, strlen(default_format))){
                        printf("The format is not '%s' (%s)\n", default_format, format);
                        status = 1;
                }
                free(format);
        }
#endif 

        /* Read molecule name */
        if(!status){
            char *mol_name = NULL;
            status = SpIO_H5GetAttribute_string(sfp->h5f_id, "/", "molec", &mol_name);

            /* Load molecule if present */
            if(strlen(mol_name) > 0) {
                if(!(model->parms.mol = SpIO_FreadMolec(mol_name)))
                    PyWrErr_SetString(PyExc_Exception, "Molecular file name not found '%s'", mol_name);
                    status = 1;
            }
            free(mol_name);
        }

        /* Read coordinate name */
        if(!status){
                char *coordinate = NULL; 
                status = SpIO_H5GetAttribute_string(sfp->h5f_id, "/", "geom", &coordinate);

                /* Load molecule if present */
                if(strlen(coordinate) > 0) {

                    DatINode *geom = Dat_IList_NameLookup(GEOM_TYPES, coordinate);
     
                    model->parms.geom = geom->idx;
                }
                free(coordinate);
        }
        
        herr_t hstatus;

        /* Read T_cmb */
        if(!status) {
                hstatus = H5LTget_attribute_double(sfp->h5f_id, "/", "T_cmb", &model->parms.T_cmb);
                if(hstatus < 0)
                        status = 1;
        }
        /* Read T_in */
        if(!status) {
                hstatus = H5LTget_attribute_double(sfp->h5f_id, "/", "T_in", &model->parms.T_in);
                if(hstatus < 0)
                        status = 1;
        }
	
     

        /* Read number of outer source */
        if(!status) {
            hstatus = H5LTget_attribute_int(sfp->h5f_id, "/", "Outer_Source", &model->parms.Outer_Source);
            
            if(hstatus < 0)
                status = 1;
            size_t nSource = (size_t) model->parms.Outer_Source;
            if (nSource != 0){
                ZoneH5_Record_Source *data; 
                data = Mem_CALLOC( nSource, data);

                if(!status) 
                    status = ZoneH5_FreadTable_Source(sfp->h5f_id, data); 

                model->parms.source = Mem_CALLOC( nSource, model->parms.source); 
                if(!status) { 
                    for(size_t i = 0; i < nSource; i++)
//                         printf("OK %e \n", model->parms.source[i].temperature);
                        SpIO_SourceFromH5Record( &model->parms.source[i], data[i]); 
                } 
                free(data); 
            }
            
        }
        
        /* Read pops-switch */
        if(!status) {
                hstatus = H5LTget_attribute_int(sfp->h5f_id, "/", "pops", &model->parms.pops);
                if(hstatus < 0)
                        status = 1;
        }
        /* Read dust-switch */
        if(!status) {
                hstatus = H5LTget_attribute_int(sfp->h5f_id, "/", "dust", &model->parms.dust);
                if(hstatus < 0)
                        status = 1;
        }
        /* Read polariz-switch */
        if(!status) {
                hstatus = H5LTget_attribute_int(sfp->h5f_id, "/", "polariz", &model->parms.polariz);
                if(hstatus < 0)
                        status = 1;
        }
        /* Read z (zeeman splitting) factor */
        if(!status) {
                hstatus = H5LTget_attribute_double(sfp->h5f_id, "/", "z", &model->parms.z);
                if(hstatus < 0)
                        status = 1;
        }


#if 0 
	/* Set velocity field: the "velfield" attribute is expected
	 * to contain a Python function named "Vgas((x1,x2,x3),(min1,min2,min3),(max1,max2,max3)" where
	 * (x1,x2,x3) are the three coordinates of the ray position, (min1,min2,min3) are the lower
	 * bounds of the bounding box, and (max1,max2,max3) are the upper baounds of the bounding
	 * box. An example Vgas function would be:
	 * def Vgas((x1, x2, x3), (min1, min2, min3), (max1, max2, max3)):"
	 *	return [0,0,0]"
	 */
	if(!status) {
                char *velfield = NULL;
		//if(!strncmp(velfield, "grid", strlen(velfield))) {
		if(!strncmp(velfield, "grid", strlen("grid"))) {
			model->parms.velfield = NULL;
		}
		else {
			status = PyRun_SimpleString(velfield);
			PyObject *__main__ = NULL, *ret = NULL;
			__main__ = PyImport_AddModule("__main__");
			model->parms.velfield = PyObject_GetAttrString(__main__, "Vgas");
			Py_INCREF(model->parms.velfield);
			/* Try calling the function and make sure it returns a sequence of length 3 */
                        char format[] = "(d,d,d),(d,d,d),(d,d,d)";
			ret = PyObject_CallFunction(model->parms.velfield, format, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			if(!PySequence_Check(ret)) {
				printf( "Vgas() does not return a vector!\n");
				status = 1;
			}
		}
		
		free(velfield);
	}
#endif

	/* Load grid */	
	if(!status)
		status = 
		SpIO_H5ReadGrid( 
                        sfp->h5f_id, 
                        (*read_pops) ? popsfp->h5f_id : NULL, 
                        &model->grid, 
                        &model->parms, 
                        read_pops
                );
	/* Cleanup */	
        

	if(status)
		printf( "Error reading model from file `%s'\n", sfp->name);
        
	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_OpenModel(const char *sourcefname, const char *popsfname, SpModel *model, int *read_pops)
{
	int status = 0;

	SpFile *sfp = SpIO_OpenFile(sourcefname, Sp_OLD);
        if( !sfp )
                status = 1;

        SpFile *popsfp;
        if (*read_pops){
                popsfp = SpIO_OpenFile(popsfname, Sp_OLD);
                if( !popsfp )
                        status = 1;
        }

        if(!status)
		status = SpIO_FreadModel(sfp, popsfp, model, read_pops);
        

        /* clean-up file */
        if(sfp)
                SpIO_CloseFile(sfp);

        if (*read_pops)
                SpIO_CloseFile(popsfp);
        
	return status;
}

/*----------------------------------------------------------------------------*/

Molec *SpIO_FreadMolec(const char *molname)
{
	FILE *fp;
	char *path;
	Molec *mol = 0;

	path = Mem_Sprintf("%s/%s.dat", Sp_parm.molec_path, molname);
	fp = fopen(path, "r");
        
        
	if(!fp) {
		printf("File `%s' could not be opened\n", path);
        }
	else {
                
		mol = Mol_ReadLamda(fp, path, molname);

		if(!mol)
			printf("Error loading molecule from `%s'\n", path);
		else
			SpPhys_ProcLamda(mol);	
	}
    
	if(fp)
		fclose(fp);
	free(path);

	return mol;
}


/*----------------------------------------------------------------------------*/
Kappa *SpIO_FreadKappa(const char *name)
{
	char *path;
	Kappa *kap = 0;
	FILE *fp;

	path = Mem_Sprintf("%s/%s.tab", Sp_parm.kappa_path, name);

	fp = fopen(path, "r");

	if(!fp) {
		printf("File `%s' could not be opened\n", path);
		goto cleanup;
	}
	else {
		kap = Kap_New_Table(name, path, fp);
		if(!kap)
			printf("Error loading opacity from `%s'\n", path);
	}

	cleanup:
	if(fp)
		fclose(fp);
	free(path);

	return kap;
}

/*----------------------------------------------------------------------------*/

Kappa *SpIO_LoadKappa(const char *string)
{
	char arg1[BUFSIZ], arg2[BUFSIZ], arg3[BUFSIZ], arg4[BUFSIZ];

	sscanf(string, "%[^,],%[^,],%[^,],%[^,]", arg1, arg2, arg3, arg4);

	if(!strcmp(arg1, "powerlaw"))
		return Kap_New_Powerlaw(atof(arg2), atof(arg3), atof(arg4));
	else if(!strcmp(arg1, "table"))
		return SpIO_FreadKappa(arg2);
	else
		Deb_ASSERT(0);

	return NULL;
}

/*----------------------------------------------------------------------------*/

ZoneH5_Record_Zone SpIO_ZoneToH5Record(const Zone *zone)
{
	ZoneH5_Record_Zone zone_data;
        Mem_BZERO(&zone_data);
	zone_data.level = zone->level;

	Mem_MEMCPY(zone_data.max, zone->voxel.max.x, 3);
	Mem_MEMCPY(zone_data.min, zone->voxel.min.x, 3);
	Mem_MEMCPY(zone_data.cen, zone->voxel.cen.x, 3);

	zone_data.nchildren = (long int)zone->nchildren;
	for(size_t i = 0; i < 3; i++) {
                zone_data.naxes[i] = (long int)GeVec3_X(zone->naxes, i);
	}

	return zone_data;
}

/*----------------------------------------------------------------------------*/

ZoneH5_Record_Grid SpIO_GridToH5Record(const Zone *zone)
{
        ZoneH5_Record_Grid grid_data;
        Mem_BZERO(&grid_data);
        
        SpPhys *pp = zone->data;
        grid_data.level = zone->level;
        grid_data.pos = (long int)zone->pos;
        grid_data.nchildren = (long int)zone->nchildren;

        Mem_MEMCPY(grid_data.max, zone->voxel.max.x, 3);
        Mem_MEMCPY(grid_data.min, zone->voxel.min.x, 3);
        Mem_MEMCPY(grid_data.cen, zone->voxel.cen.x, 3);
        #define CPYPP(attr) (grid_data.attr = pp->attr)
        CPYPP(n_H2);
        CPYPP(T_k);
        CPYPP(X_pH2);
        CPYPP(X_oH2);
        CPYPP(X_e);
        CPYPP(X_H);
        CPYPP(X_He);
        CPYPP(V_t);
        #undef CPYPP

        for(size_t i = 0; i < 3; i++) {
                grid_data.v_cen[i] = GeVec3_X(pp->v_cen, i);
        }

        return grid_data;
}

/*----------------------------------------------------------------------------*/

ZoneH5_Record_Molec SpIO_MolecToH5Record(const Zone *zone)
{
        ZoneH5_Record_Molec molec_data;
        Mem_BZERO(&molec_data);
        
        SpPhys *pp = zone->data;

        #define CPYPP(attr) (molec_data.attr = pp->attr)
        CPYPP(X_mol);
        #undef CPYPP

        return molec_data;
}

/*----------------------------------------------------------------------------*/

ZoneH5_Record_Dust SpIO_DustToH5Record(const Zone *zone)
{
        ZoneH5_Record_Dust dust_data;
        Mem_BZERO(&dust_data);
        
        SpPhys *pp = zone->data;

        #define CPYPP(attr) (dust_data.attr = pp->attr)
        CPYPP(T_d);
        CPYPP(dust_to_gas);
        #undef CPYPP

        strncpy(dust_data.kapp_d, pp->kapp_d, ZoneH5_KAPPLEN);

        return dust_data;
}

/*----------------------------------------------------------------------------*/

ZoneH5_Record_Polariz SpIO_PolarizToH5Record(const Zone *zone)
{
        ZoneH5_Record_Polariz polariz_data;
        Mem_BZERO(&polariz_data);
        
        SpPhys *pp = zone->data;

        #define CPYPP(attr) (polariz_data.attr = pp->attr)
        CPYPP(alpha);
        #undef CPYPP

        for(size_t i = 0; i < 3; i++) {
                polariz_data.b_cen[i] = GeVec3_X(pp->b_cen, i);
        }

        return polariz_data;
}

/*----------------------------------------------------------------------------*/

ZoneH5_Record_Source SpIO_SourceToH5Record(const SourceData *source)
{
    ZoneH5_Record_Source source_data;
    Mem_BZERO(&source_data);
    
    #define CPYPP(attr) (source_data.attr = source->attr)
    CPYPP(temperature);
    CPYPP(theta);
    CPYPP(phi);
    CPYPP(radius);
    CPYPP(distance);
    #undef CPYPP
    
    return source_data;
}

/*----------------------------------------------------------------------------*/

void SpIO_ZoneFromH5Record(Zone *zone, ZoneH5_Record_Zone record)
{
	zone->level = record.level;
        
        GeVec3_d min, max;
	Mem_MEMCPY(max.x, record.max, 3);
	Mem_MEMCPY(min.x, record.min, 3);
        zone->voxel = GeVox_Init2(zone->voxel.geom, min, max);
        
	zone->nchildren = (size_t)record.nchildren;
	for( size_t i = 0; i < 3; i++) {
		GeVec3_X(zone->naxes, i) = (size_t)record.naxes[i];
	}

	return;
}
/*----------------------------------------------------------------------------*/

void SpIO_GridFromH5Record(Zone *zone, ZoneH5_Record_Grid record)
{
        SpPhys *pp = zone->data;
        zone->level = record.level;
        zone->pos = (size_t)record.pos;
      
        GeVec3_d min, max;
        Mem_MEMCPY(max.x, record.max, 3);
        Mem_MEMCPY(min.x, record.min, 3);
        zone->voxel = GeVox_Init2(zone->parent->voxel.geom, min, max);
        
        #define CPYPP(attr) (pp->attr = record.attr)
        CPYPP(n_H2);
        CPYPP(T_k);
        CPYPP(X_pH2);
        CPYPP(X_oH2);
        CPYPP(X_e);
        CPYPP(X_H);
        CPYPP(X_He);
        CPYPP(V_t);
        #undef CPYPP

        Mem_MEMCPY(pp->v_cen.x, record.v_cen, 3);

        zone->nchildren = (size_t)record.nchildren;
        for(size_t i = 0; i < 3; i++) {
                GeVec3_X(zone->naxes, i) = (size_t)record.naxes[i];
        }

        return;
}

/*----------------------------------------------------------------------------*/
void SpIO_MolecFromH5Record(Zone *zone, ZoneH5_Record_Molec record)
{
        SpPhys *pp = zone->data;

        #define CPYPP(attr) (pp->attr = record.attr)
        CPYPP(X_mol);
        #undef CPYPP

        return;
}

/*----------------------------------------------------------------------------*/
void SpIO_DustFromH5Record(Zone *zone, ZoneH5_Record_Dust record)
{
        SpPhys *pp = zone->data;

        #define CPYPP(attr) (pp->attr = record.attr)
        CPYPP(T_d);
        CPYPP(dust_to_gas);
        #undef CPYPP
        
        strncpy(pp->kapp_d, record.kapp_d, ZoneH5_KAPPLEN);

        return;
}

/*----------------------------------------------------------------------------*/
void SpIO_PolarizFromH5Record(Zone *zone, ZoneH5_Record_Polariz record)
{
        SpPhys *pp = zone->data;

        #define CPYPP(attr) (pp->attr = record.attr)
        CPYPP(alpha);
        #undef CPYPP
        
        Mem_MEMCPY(pp->b_cen.x, record.b_cen, 3);

        return;
}

/*----------------------------------------------------------------------------*/
void SpIO_SourceFromH5Record(SourceData *source, ZoneH5_Record_Source record)
{
    #define CPYPP(attr) (source->attr = record.attr)
    CPYPP(temperature);
    CPYPP(theta);
    CPYPP(phi);
    CPYPP(radius);
    CPYPP(distance);
    #undef CPYPP    
    return;
}
/*----------------------------------------------------------------------------*/

int SpIO_H5WriteGrid(hid_t h5f_id, const Zone *zone, const SpPhysParm *parms)
/* Write data of all children to an HDF5 table */
{
        ZoneH5_Record_Zone this_zone;
	/* Write properties of this zone */
	this_zone = SpIO_ZoneToH5Record(zone);
	int status = ZoneH5_FwriteTable_Zone(h5f_id, &this_zone);

        #define COPY_TYPE_TO_RECORD(DataType, ZoneH5_FwriteTable_Type, SpIO_Type ) \
            {   DataType *data; \
                /* Allocate table data */ \
                data = Mem_CALLOC(zone->nchildren, data); \
                /* Load data values */ \
                for(size_t i = 0; i < zone->nchildren; i++) \
                        data[i] = SpIO_Type(zone->children[i]); \
                /* Create and write grid data */ \
                if(!status) \
                        status = ZoneH5_FwriteTable_Type(h5f_id, data, zone->nchildren); \
                /* Free table data */ \
                free(data); \
            }        
        
        
        COPY_TYPE_TO_RECORD(ZoneH5_Record_Grid, ZoneH5_FwriteTable_Grid, SpIO_GridToH5Record);
        /* if any children of the zone is leaf (zp->nchildren = 0) then write all radiative-relatied attributes */
        int write_attr = 0;
        for(size_t i = 0; i < zone->nchildren; i++) {
		if(zone->children[i]->nchildren == 0) {
                    write_attr = 1;
                    break;
		}
	}
        
        if (write_attr){
            /* Create and write pops if present */
            if(parms->mol){
                    COPY_TYPE_TO_RECORD(ZoneH5_Record_Molec, ZoneH5_FwriteTable_Molec, SpIO_MolecToH5Record);
                    SpIO_H5WritePops(h5f_id, zone);
            }
            if(parms->dust){
                    COPY_TYPE_TO_RECORD(ZoneH5_Record_Dust, ZoneH5_FwriteTable_Dust, SpIO_DustToH5Record);
            }
            if(parms->polariz){
                    COPY_TYPE_TO_RECORD(ZoneH5_Record_Polariz, ZoneH5_FwriteTable_Polariz, SpIO_PolarizToH5Record);
            }
        }
        #undef COPY_TYPE_TO_RECORD
	/* Create and write tau if present */
        //SpPhys *pp = zone->data;
	//if(pp->mol)
		//SpIO_H5WriteTau(h5f_id, zone);

	/* Loop through children and write sub-grids if present */
	for(size_t i = 0; i < zone->nchildren; i++) {
		/* Pointer to child */
		Zone *zp = zone->children[i];

		if(zp->nchildren > 0) {
			/* Sub-grid present, create group */
			Deb_ASSERT(zp->children != NULL);

			/* Build group name string: POS */
			//strp = Mem_Sprintf("grid%lu", (unsigned long)zone->pos);
			char *strp = Mem_Sprintf("grid%lu", i);

			/* Create a new group for this grid */
			hid_t group_id = H5Gcreate(h5f_id, strp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			/* Recursively call SpIO_H5WriteGrid() to write
			 * sub-grids */
			SpIO_H5WriteGrid(group_id, zp, parms);

			/* Close group */
			H5Gclose(group_id);
		}
	}

	return 0;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5ReadGrid(hid_t h5f_id, hid_t popsh5f_id, Zone **zone, SpPhysParm *parms, int *read_pops)
/* Write data of all children to an HDF5 table */
{
	int status = 0;

	if(!(*zone)){
                // root zone initialization
		*zone = SpZone_ALLOC(parms);
		//*zone = Zone_Alloc(pos, parent, SpPhys_Alloc, parms);
                
                // root zone's coordinate is the definition in parms
                (*zone)->voxel.geom = parms->geom;
        }
	/* Read data for this zone */
        ZoneH5_Record_Zone this_zone;
	status = ZoneH5_FreadTable_Zone(h5f_id, &this_zone);

	if(!status) {
		SpIO_ZoneFromH5Record((*zone), this_zone);

		/* Grow grid */
		SpZone_GROW(*zone, (*zone)->naxes, parms);
	}


        

	#define COPY_TYPE_FROM_RECORD(DataType, ZoneH5_FreadTable_Type, SpIO_Type ) \
	{   DataType *data; \
            data = Mem_CALLOC((*zone)->nchildren, data); \
            if(!status) \
                status = ZoneH5_FreadTable_Type(h5f_id, data); \
                if(!status) { \
                    for(size_t i = 0; i < (*zone)->nchildren; i++) \
                        SpIO_Type((*zone)->children[i], data[i]); \
                } \
                free(data); \
        }
        COPY_TYPE_FROM_RECORD(ZoneH5_Record_Grid, ZoneH5_FreadTable_Grid, SpIO_GridFromH5Record);

        int read_leaf = 0;
        for(size_t i = 0; i < (*zone)->nchildren; i++){
            if ( (*zone)->children[i]->nchildren == 0 ){
                read_leaf = 1;
                break;
            }
        }
        
        if (read_leaf){
            if(parms->mol)
                COPY_TYPE_FROM_RECORD(ZoneH5_Record_Molec, ZoneH5_FreadTable_Molec, SpIO_MolecFromH5Record);
            if(parms->dust)
                COPY_TYPE_FROM_RECORD(ZoneH5_Record_Dust, ZoneH5_FreadTable_Dust, SpIO_DustFromH5Record);
            if(parms->polariz)
                COPY_TYPE_FROM_RECORD(ZoneH5_Record_Polariz, ZoneH5_FreadTable_Polariz, SpIO_PolarizFromH5Record);
            /* Read pops if present */
            if(!status && parms)
		if( parms->mol && (*read_pops) )
                    status = SpIO_H5ReadPops(popsh5f_id, *zone);
        }
        #undef COPY_TYPE_FROM_RECORD

        /* Recursively read subgrids */
        for(size_t i = 0; i < (*zone)->nchildren; i++) {
            if(!status) {
                if((*zone)->children[i]->nchildren > 0){
                    char *grid_idx = Mem_Sprintf("grid%lu", (unsigned long)i);

                    hid_t group_id = H5Gopen(h5f_id, grid_idx, H5P_DEFAULT);
                    hid_t popsgroup_id = (*read_pops) ? H5Gopen(popsh5f_id, grid_idx, H5P_DEFAULT) : NULL;
                    status = SpIO_H5ReadGrid(
                        group_id, 
                        popsgroup_id, 
                        &(*zone)->children[i], 
                        parms, 
                        read_pops);

                }
            }
        }
        
	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5WritePops(hid_t h5f_id, const Zone *zone)
{
	SpPhys *pp = zone->data;

	/* Just in case the programmer did something stupid */
	Deb_ASSERT(pp->mol != NULL);
        Deb_ASSERT(pp->pops_preserve != NULL);

	herr_t hstatus = 0;
	int status = 0;
	size_t
		i, j,
		nlev = pp->mol->nlev,
		record_size =  sizeof(double) * nlev,
		field_offset[nlev];
	const char **field_names = Mem_CALLOC(nlev, field_names);
	char **level_names = Mem_CALLOC(nlev, level_names);
	hid_t field_type[nlev];
	hsize_t chunk_size = 10;
	int *fill_data = NULL, compress  = 0;
	double *pops = Mem_CALLOC(zone->nchildren * nlev, pops);

	/* Init fields */
	for(i = 0; i < nlev; i++) {
		field_offset[i] = i * sizeof(double);
		field_type[i] = H5T_NATIVE_DOUBLE;
		level_names[i] = Mem_Sprintf("lev%lu", (unsigned long)i);
		field_names[i] = level_names[i];
	}

	/* Load data */
	#define POPS(i, j)\
		pops[(j) + nlev * (i)]

	for(i = 0; i < zone->nchildren; i++) {
		pp = zone->children[i]->data;
		for(j = 0; j < nlev; j++) {
                    POPS(i, j) = pp->pops_preserve[j];
		}
	}

	#undef POPS

	/* Write table */
	hstatus = H5TBmake_table(
		"Level populations",
		h5f_id,
		"POPS",
		(hsize_t)nlev,
		(hsize_t)zone->nchildren,
		record_size,
		field_names,
		field_offset,
		field_type,
		chunk_size,
		fill_data,
		compress,
		pops
	);

	/* Cleanup */
	for(i = 0; i < nlev; i++)
		free(level_names[i]);
	free(level_names);
	free(field_names);
	free(pops);

	if(hstatus < 0)
		status = printf("Error writing `POPS' table\n");

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5ReadPops(hid_t h5f_id, Zone *zone)
/* Write data of all children to an HDF5 table */
{
	SpPhys *pp = zone->data;

	/* Just in case the programmer did something stupid */
	Deb_ASSERT(pp->mol != NULL);
        Deb_ASSERT(pp->pops_preserve != NULL);

	int status = 0;
	size_t
		i, j,
		nlev = pp->mol->nlev,
		record_size =  sizeof(double) * nlev,
		field_sizes[nlev],
		field_offsets[nlev];
	double *pops = Mem_CALLOC(zone->nchildren * nlev, pops);
	herr_t hstatus;

	/* Init fields */
	for(i = 0; i < nlev; i++) {
		field_offsets[i] = i * sizeof(double);
		field_sizes[i] = sizeof(double);
	}
	
	/* Read pops table */
	hstatus = H5TBread_table(h5f_id, "POPS", record_size, field_offsets, field_sizes, pops);
	if(hstatus < 0)
		status = printf("Error reading HDF5 `%s' table\n", "POPS");
	#define POPS(i, j)\
		pops[(j) + nlev * (i)]

	for(i = 0; i < zone->nchildren; i++) {
		pp = zone->children[i]->data;

		for(j = 0; j < nlev; j++) {
                    pp->pops_preserve[j] = POPS(i, j);
		}
	}

	#undef POPS

	free(pops);

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5WriteTau(hid_t h5f_id, const Zone *zone)
{
	SpPhys *pp = zone->data;

	/* Just in case the programmer did something stupid */
	Deb_ASSERT(pp->mol != NULL);
	Deb_ASSERT(pp->tau != NULL);

	herr_t hstatus = 0;
	int status = 0;
	size_t
		i, j,
		nrad = pp->mol->nrad,
		record_size =  sizeof(double) * nrad,
		field_offset[nrad];
	const char **field_names = Mem_CALLOC(nrad, field_names);
	char **level_names = Mem_CALLOC(nrad, level_names);
	hid_t field_type[nrad];
	hsize_t chunk_size = 10;
	int *fill_data = NULL, compress  = 0;
	double *tau = Mem_CALLOC(zone->nchildren * nrad,  tau);

	/* Init fields */
	for(i = 0; i < nrad; i++) {
		field_offset[i] = i * sizeof(double);
		field_type[i] = H5T_NATIVE_DOUBLE;
		level_names[i] = Mem_Sprintf("line%lu", (unsigned long)i);
		field_names[i] = level_names[i];
	}

	/* Load data */
	#define TAU(i, j)\
		tau[(j) + nrad * (i)]

	for(i = 0; i < zone->nchildren; i++) {
		pp = zone->children[i]->data;
		for(j = 0; j < nrad; j++) {
			TAU(i, j) = pp->tau[j];
		}
	}

	#undef TAU

	/* Write table */
	hstatus = H5TBmake_table(
		"Level populations",
		h5f_id,
		"TAU",
		(hsize_t)nrad,
		(hsize_t)zone->nchildren,
		record_size,
		field_names,
		field_offset,
		field_type,
		chunk_size,
		fill_data,
		compress,
		tau
	);

	/* Cleanup */
	for(i = 0; i < nrad; i++)
		free(level_names[i]);
	free(level_names);
	free(field_names);
	free(tau);

	if(hstatus < 0)
		status = printf("Error writing `TAU' table\n");

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5ReadTau(hid_t h5f_id, Zone *zone)
/* Write data of all children to an HDF5 table */
{
	SpPhys *pp = zone->data;

	/* Just in case the programmer did something stupid */
	Deb_ASSERT(pp->mol != NULL);
	Deb_ASSERT(pp->tau != NULL);

	int status = 0;
	size_t
		i, j,
		nrad = pp->mol->nrad,
		record_size =  sizeof(double) * nrad,
		field_sizes[nrad],
		field_offsets[nrad];
	double *tau = Mem_CALLOC(zone->nchildren * nrad, tau);
	herr_t hstatus;

	/* Init fields */
	for(i = 0; i < nrad; i++) {
		field_offsets[i] = i * sizeof(double);
		field_sizes[i] = sizeof(double);
	}

	/* Read tau table */
	hstatus = H5TBread_table(h5f_id, "TAU", record_size, field_offsets, field_sizes, tau);
	if(hstatus < 0)
		status = printf("Error reading HDF5 `%s' table\n", "TAU");

	#define TAU(i, j)\
		tau[(j) + nrad * (i)]

	for(i = 0; i < zone->nchildren; i++) {
		pp = zone->children[i]->data;

		for(j = 0; j < nrad; j++) {
			pp->tau[j] = TAU(i, j);
		}
	}

	#undef TAU

	free(tau);

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5GetAttribute_string(hid_t h5f_id, const char *obj_name, const char *attr_name, char **attribute)
{
	int status = 0;
	herr_t hstatus = 0;
	hsize_t dims;
	H5T_class_t class;
	size_t size;

	if(hstatus >= 0)
		hstatus = H5LTget_attribute_info(h5f_id, obj_name, attr_name, &dims, &class, &size);
        
	if(hstatus >= 0) {
		*attribute = Mem_CALLOC(size+1, *attribute);
		hstatus = H5LTget_attribute_string(h5f_id, obj_name, attr_name, *attribute);
	}
	
        
	if(hstatus < 0) {
		printf( "Error getting attribute '%s' from '%s'\n", attr_name, obj_name);
		status = 1;
	}

	return status;
}





