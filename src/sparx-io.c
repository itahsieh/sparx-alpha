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
		Err_SETSTRING("Error opening file `%s'", fname);
	}

	return sfp;
}

/*----------------------------------------------------------------------------*/

int SpIO_OpenFile2(const char *fname, int mode, SpFile **fp)
{
	int sts = 0;

	*fp = SpIO_OpenFile(fname, mode);

	if(!(*fp)) {
		PyWrErr_SetString(PyExc_Exception, "Error opening SPARX file '%s'", fname);
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
	Molec *mol = model.parms.mol;
	const char *mol_name = mol ? mol->name : "";
	herr_t hstatus;

	/* Write molecule name */
	if(!status) {
		hstatus = H5LTset_attribute_string(sfp->h5f_id, "/", "molec", mol_name);
		if(hstatus < 0)
			status = 1;
	}

	/* Write T_cmb */
	if(!status) {
		hstatus = H5LTset_attribute_double(sfp->h5f_id, "/", "T_cmb", &model.parms.T_cmb, (size_t)1);
		if(hstatus < 0)
			status = 1;
	}

	/* Write gas_to_dust */
	if(!status) {
		hstatus = H5LTset_attribute_double(sfp->h5f_id, "/", "gas_to_dust", &model.parms.gas_to_dust, (size_t)1);
		if(hstatus < 0)
			status = 1;
	}

	/* Write velocity field info */
	//debug
	if(!status) {
		#if 0
		hstatus = H5LTset_attribute_string(sfp->h5f_id, "/", "velfield",
			"def Vgas(x, min, max):\n"
			"	return [0,0,0]\n"
			);
		#endif
		hstatus = H5LTset_attribute_string(sfp->h5f_id, "/", "velfield", "grid");
		if(hstatus < 0)
			status = 1;
	}

	/* Write grid */
	if(!status)
		status = SpIO_H5WriteGrid(sfp->h5f_id, model.grid);

	if(status)
		status = Err_SETSTRING("Error writing model to file `%s'", sfp->name);

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_FreadModel(const SpFile *sfp, const SpFile *popsfp, SpModel *model)
{
	int status = 0;
	herr_t hstatus;
	//char mol_name[BUFSIZ], velfield[BUFSIZ];
	char *mol_name = NULL, *velfield = NULL, format[] = "(d,d,d),(d,d,d),(d,d,d)";
	PyObject *__main__ = NULL, *ret = NULL;

	Deb_ASSERT(model->parms.mol == NULL);
	Deb_ASSERT(model->grid == NULL);

	/* Read molecule name */
	if(!status)
		status = SpIO_H5GetAttribute_string(popsfp->h5f_id, "/", "molec", &mol_name);
	/* Read velocity field info */
	if(!status)
		status = SpIO_H5GetAttribute_string(sfp->h5f_id, "/", "velfield", &velfield);
	#if 0
	hstatus = H5LTget_attribute_string(sfp->h5f_id, "/", "molec", mol_name);
	if(hstatus < 0)
		status = 1;

	hstatus = H5LTget_attribute_string(sfp->h5f_id, "/", "velfield", velfield);
	if(hstatus < 0)
		status = 1;
	#endif

	/* Read T_cmb */
	if(!status) {
		hstatus = H5LTget_attribute_double(sfp->h5f_id, "/", "T_cmb", &model->parms.T_cmb);
		if(hstatus < 0)
			status = 1;
	}
	/* Read gas_to_dust */
	if(!status) {
		hstatus = H5LTget_attribute_double(sfp->h5f_id, "/", "gas_to_dust", &model->parms.gas_to_dust);
		if(hstatus < 0)
			status = 1;
	}
	/* Load molecule if present */
	if(!status && strlen(mol_name) > 0) {
		if(!(model->parms.mol = SpIO_FreadMolec(mol_name)))
		status = 1;
	}
	/* Set velocity field: the "velfield" attribute is expected
	 * to contain a Python function named "Vgas((x1,x2,x3),(min1,min2,min3),(max1,max2,max3)" where
	 * (x1,x2,x3) are the three coordinates of the ray position, (min1,min2,min3) are the lower
	 * bounds of the bounding box, and (max1,max2,max3) are the upper baounds of the bounding
	 * box. An example Vgas function would be:
	 * def Vgas((x1, x2, x3), (min1, min2, min3), (max1, max2, max3)):"
	 *	return [0,0,0]"
	 */
	if(!status) {
		//if(!strncmp(velfield, "grid", strlen(velfield))) {
		if(!strncmp(velfield, "grid", strlen("grid"))) {
			model->parms.velfield = NULL;
		}
		else {
			status = PyRun_SimpleString(velfield);
			__main__ = PyImport_AddModule("__main__");
			model->parms.velfield = PyObject_GetAttrString(__main__, "Vgas");
			Py_INCREF(model->parms.velfield);
			/* Try calling the function and make sure it returns a sequence of length 3 */
			ret = PyObject_CallFunction(model->parms.velfield, format, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			if(!PySequence_Check(ret)) {
				PyWrErr_SetString(PyExc_Exception, "Vgas() does not return a vector!");
				status = 1;
			}
		}
	}
	/* Load grid */	
	if(!status)
		status = SpIO_H5ReadGrid(sfp->h5f_id, popsfp->h5f_id, &model->grid, &model->parms);
	/* Cleanup */	
	free(mol_name);
	free(velfield);

	if(status)
		PyWrErr_SetString(PyExc_Exception, "Error reading model from file `%s'", sfp->name);

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_OpenModel(const char *fname, const char *popsfname, SpModel *model)
{
	int status = 0;
	SpFile *sfp, *popsfp;

	sfp = SpIO_OpenFile(fname, Sp_OLD);
	popsfp = SpIO_OpenFile(popsfname, Sp_OLD);

	if(!sfp)
		status = 1;

	if(!status)
		status = SpIO_FreadModel(sfp, popsfp, model);


	if(sfp)
		SpIO_CloseFile(sfp);
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
		Err_SETSTRING("File `%s' could not be opened", path);
		goto cleanup;
		
	}
	else {
		mol = Mol_ReadLamda(fp, path, molname);
		if(!mol)
			Err_SETSTRING("Error loading molecule from `%s'", path);
		else
			SpPhys_ProcLamda(mol);		
	}

	cleanup:
	if(fp)
		fclose(fp);
	free(path);

	return mol;
}

/*----------------------------------------------------------------------------*/
Molec *SpIO_FreadMolec_hyper(const char *molname)
{
	FILE *fp;
	char *path;
	Molec *mol = 0;

	path = Mem_Sprintf("%s/%s.dat", Sp_parm.molec_path, molname);
	fp = fopen(path, "r");

	if(!fp) {
		Err_SETSTRING("File `%s' could not be opened", path);
		goto cleanup;
		
	}
	else {
		mol = Mol_ReadLamda_hyper(fp, path, molname);
		if(!mol)
			Err_SETSTRING("Error loading molecule from `%s'", path);
		else
			SpPhys_ProcLamda(mol);		
	}

	cleanup:
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
		Err_SETSTRING("File `%s' could not be opened", path);
		goto cleanup;
	}
	else {
		kap = Kap_New_Table(name, path, fp);
		if(!kap)
			Err_SETSTRING("Error loading opacity from `%s'", path);
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

ZoneH5_Record SpIO_ZoneToH5Record(const Zone *zone)
{
	size_t i;
	ZoneH5_Record zone_data;
	SpPhys *pp;
	DatINode *geom;

	Mem_BZERO(&zone_data);

	pp = zone->data;
	zone_data.level = zone->level;
	zone_data.pos = (long int)zone->pos;
	geom = Dat_IList_IdxLookup(GEOM_TYPES, zone->voxel.geom);
	Deb_ASSERT(geom != NULL);
	strncpy(zone_data.geom, geom->name, ZoneH5_GEOMLEN);
	Mem_MEMCPY(zone_data.max, zone->voxel.max.x, 3);
	Mem_MEMCPY(zone_data.min, zone->voxel.min.x, 3);
	Mem_MEMCPY(zone_data.cen, zone->voxel.cen.x, 3);
	#define CPYPP(attr) (zone_data.attr = pp->attr)
	CPYPP(n_H2);
	CPYPP(T_k);
	CPYPP(T_d);
	CPYPP(T_ff);
	CPYPP(T_bb);
	CPYPP(X_mol);
	CPYPP(X_pH2);
	CPYPP(X_oH2);
        
	CPYPP(X_e);
	CPYPP(X_H);
	CPYPP(X_He);
	CPYPP(V_t);
	CPYPP(ds);
	#undef CPYPP

	strncpy(zone_data.kapp_d, pp->kapp_d, ZoneH5_KAPPLEN);
	strncpy(zone_data.kapp_ff, pp->kapp_ff, ZoneH5_KAPPLEN);

	for(i = 0; i < 6; i++) {
		Mem_MEMCPY(&zone_data.vedge[i][0], pp->v_edge[i].x, 3);
	}
	zone_data.nchildren = (long int)zone->nchildren;
	for(i = 0; i < 3; i++) {
                zone_data.v_cen[i] = GeVec3_X(pp->v_cen, i);
                zone_data.b_cen[i] = GeVec3_X(pp->b_cen, i);
                zone_data.naxes[i] = (long int)GeVec3_X(zone->naxes, i);
	}

	return zone_data;
}

/*----------------------------------------------------------------------------*/

void SpIO_ZoneFromH5Record(Zone *zone, ZoneH5_Record record)
{
	size_t i;
	SpPhys *pp;
	GeVec3_d min, max;
	DatINode *geom;

	pp = zone->data;
	zone->level = record.level;
	zone->pos = (size_t)record.pos;
	Mem_MEMCPY(max.x, record.max, 3);
	Mem_MEMCPY(min.x, record.min, 3);
	geom = Dat_IList_NameLookup(GEOM_TYPES, record.geom);
	Deb_ASSERT(geom != NULL);
	zone->voxel = GeVox_Init2(geom->idx, min, max);
	#define CPYPP(attr) (pp->attr = record.attr)
	CPYPP(n_H2);
	CPYPP(T_k);
	CPYPP(T_d);
	CPYPP(T_ff);
	CPYPP(T_bb);
	CPYPP(X_mol);
	CPYPP(X_pH2);
	CPYPP(X_oH2);
	CPYPP(X_e);
	CPYPP(X_H);
	CPYPP(X_He);
	CPYPP(V_t);
	CPYPP(ds);
	#undef CPYPP

	strncpy(pp->kapp_d, record.kapp_d, ZoneH5_KAPPLEN);
	strncpy(pp->kapp_ff, record.kapp_ff, ZoneH5_KAPPLEN);

	for(i = 0; i < 6; i++) {
		Mem_MEMCPY(pp->v_edge[i].x, &record.vedge[i][0], 3);
	}
	Mem_MEMCPY(pp->v_cen.x, record.v_cen, 3);
	Mem_MEMCPY(pp->b_cen.x, record.b_cen, 3);
	zone->nchildren = (size_t)record.nchildren;
	for(i = 0; i < 3; i++) {
		GeVec3_X(zone->naxes, i) = (size_t)record.naxes[i];
	}

	return;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5WriteGrid(hid_t h5f_id, const Zone *zone)
/* Write data of all children to an HDF5 table */
{
	int status = 0;
	SpPhys *pp = zone->data;
	char *strp;
	size_t i;
	Zone *zp;
	ZoneH5_Record *zone_data, this_zone;
	hid_t group_id;

	/* Write properties of this zone */
	this_zone = SpIO_ZoneToH5Record(zone);
	status = ZoneH5_FwriteTable(h5f_id, "ZONE", &this_zone, (size_t)1);

	/* Allocate table data */
	zone_data = Mem_CALLOC(zone->nchildren, zone_data);

	/* Load data values */
	for(i = 0; i < zone->nchildren; i++)
		zone_data[i] = SpIO_ZoneToH5Record(zone->children[i]);

	/* Create and write grid data */
	if(!status)
		status = ZoneH5_FwriteTable(h5f_id, "GRID", zone_data, zone->nchildren);

	/* Create and write pops if present */
	if(pp->mol)
		SpIO_H5WritePops(h5f_id, zone);

	/* Create and write tau if present */
	if(pp->mol)
		SpIO_H5WriteTau(h5f_id, zone);

	/* Free table data */
	free(zone_data);

	/* Loop through children and write sub-grids if present */
	for(i = 0; i < zone->nchildren; i++) {
		/* Pointer to child */
		zp = zone->children[i];

		if(zp->nchildren > 0) {
			/* Sub-grid present, create group */
			Deb_ASSERT(zp->children != NULL);

			/* Build group name string: POS */
			//strp = Mem_Sprintf("grid%lu", (unsigned long)zone->pos);
			strp = Mem_Sprintf("grid%lu", i);

			/* Create a new group for this grid */
			group_id = H5Gcreate(h5f_id, strp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			/* Recursively call SpIO_H5WriteGrid() to write
			 * sub-grids */
			SpIO_H5WriteGrid(group_id, zp);

			/* Close group */
			H5Gclose(group_id);
		}
	}

	return 0;
}

/*----------------------------------------------------------------------------*/

//int SpIO_H5ReadGrid(hid_t h5f_id, Zone **zone, SpPhysParm *parms, size_t pos, Zone *parent)
int SpIO_H5ReadGrid(hid_t h5f_id, hid_t popsh5f_id, Zone **zone, SpPhysParm *parms)
/* Write data of all children to an HDF5 table */
{
	int status = 0;
	ZoneH5_Record *zone_data, this_zone;
	size_t i;
	hid_t group_id, popsgroup_id;
	char *strp;

	if(!(*zone))
		*zone = SpZone_ALLOC(parms);
		//*zone = Zone_Alloc(pos, parent, SpPhys_Alloc, parms);

	/* Read data for this zone */
	status = ZoneH5_FreadTable(h5f_id, "ZONE", &this_zone);

	if(!status) {
		SpIO_ZoneFromH5Record((*zone), this_zone);

		/* Grow grid */
		SpZone_GROW(*zone, (*zone)->naxes, parms);
	}

	/* Read pops if present */
	if(!status && parms) {
		//printf("%d\n",*parms->mol);
		if(parms->mol){
			// modified by I-Ta 2012.10.26
			//status = SpIO_H5ReadPops(h5f_id, *zone);
			status = SpIO_H5ReadPops(popsh5f_id, *zone);
		}
	}

	/* Read tau if present */
	if(!status && parms) {
		
		if(parms->mol)
			// modified by I-Ta 2012.10.26
			//status = SpIO_H5ReadTau(h5f_id, *zone);
			status = SpIO_H5ReadTau(popsh5f_id, *zone);
			
	}
	zone_data = Mem_CALLOC((*zone)->nchildren, zone_data);

	/* Read grid */
	if(!status)
		status = ZoneH5_FreadTable(h5f_id, "GRID", zone_data);
	if(!status) {
		for(i = 0; i < (*zone)->nchildren; i++)
			SpIO_ZoneFromH5Record((*zone)->children[i], zone_data[i]);
			//printf("level=%d pos= nchild=%lu \n", (*zone)->level+1/*, (*zone)->children[i]->pos*/, (*zone)->children[i]->nchildren);
	}

	/* Recursively read subgrids */
	for(i = 0; i < (*zone)->nchildren; i++) {
		if(!status) {
			if((*zone)->children[i]->nchildren > 0){
				strp = Mem_Sprintf("grid%lu", (unsigned long)i);
				group_id = H5Gopen(h5f_id, strp,H5P_DEFAULT);
				popsgroup_id = H5Gopen(popsh5f_id, strp,H5P_DEFAULT);
				status = SpIO_H5ReadGrid(group_id, popsgroup_id, &(*zone)->children[i], parms);
			}
		}
	}

	free(zone_data);
	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5WritePops(hid_t h5f_id, const Zone *zone)
{
	SpPhys *pp = zone->data;

	/* Just in case the programmer did something stupid */
	Deb_ASSERT(pp->mol != NULL);
	Deb_ASSERT(pp->pops[0] != NULL);

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
			POPS(i, j) = pp->pops[0][j];
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
		status = Err_SETSTRING("Error writing `POPS' table");

	return status;
}

/*----------------------------------------------------------------------------*/

int SpIO_H5ReadPops(hid_t h5f_id, Zone *zone)
/* Write data of all children to an HDF5 table */
{
	SpPhys *pp = zone->data;

	/* Just in case the programmer did something stupid */
	Deb_ASSERT(pp->mol != NULL);
	Deb_ASSERT(pp->pops[0] != NULL);

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
		status = Err_SETSTRING("Error reading HDF5 `%s' table", "POPS");
	#define POPS(i, j)\
		pops[(j) + nlev * (i)]

	for(i = 0; i < zone->nchildren; i++) {
		pp = zone->children[i]->data;

		for(j = 0; j < nlev; j++) {
			pp->pops[0][j] = POPS(i, j);
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
		status = Err_SETSTRING("Error writing `TAU' table");

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
		status = Err_SETSTRING("Error reading HDF5 `%s' table", "TAU");

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
		*attribute = Mem_CALLOC(size, *attribute);
		hstatus = H5LTget_attribute_string(h5f_id, obj_name, attr_name, *attribute);
	}

	if(hstatus < 0) {
		PyWrErr_SetString(PyExc_Exception, "Error getting attribute '%s' from '%s'", attr_name, obj_name);
		status = 1;
	}

	return status;
}





