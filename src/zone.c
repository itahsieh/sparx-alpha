#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_interp.h>

#include "memory.h"
#include "zone.h"
#include "debug.h"

/*----------------------------------------------------------------------------*/

Zone *_Zone_Alloc(size_t pos, Zone *parent, void *(*DataAlloc)(const Zone *, const void *), const void *data_parms)
{
	Zone *zp = Mem_CALLOC(1, zp);

	/* level is parent's + 1, and -1 for root */
	if(!parent) {
		/* This is the root zone, level = -1 */
		zp->level = -1;

		/* Root zone is of course this zone */
		zp->root = zp;
	}
	else {
		/* Level of this zone is parent level + 1 */
		zp->level = parent->level + 1;

		/* Inherit root zone pointer from parent */
		zp->root = parent->root;
	}

	/* pos is relative to the parent->children array */
	zp->pos = pos;

	/* Parent node */
	zp->parent = parent;

	/* Allocate data if constructor is specified */
	if(DataAlloc)
		zp->data = DataAlloc(zp, data_parms);


	return zp;
}

/*----------------------------------------------------------------------------*/

void Zone_Free(void *ptr, void (*DataFree)(void *ptr))
{
	size_t i;
	Zone *zone = ptr;

	if(DataFree) {
		Deb_ASSERT(zone->data != NULL);

		/* Free data */
		DataFree(zone->data);
	}

	if(zone->children) {
		/* Free children */
		for(i = 0; i < zone->nchildren; i++) {
			Zone_Free(zone->children[i], DataFree);
		}

		/* Free pointer to cildren */
		free(zone->children);
	}

	free(zone);
}

/*----------------------------------------------------------------------------*/

void Zone_Fprintf(FILE *fp, const Zone *zone, void (*DataFprintf)(void *data, FILE *fp))
{
	size_t nindent = (size_t)(4 * (zone->level + 1));
	size_t i;
	DatINode *geom;

	#define PINDENT\
		{for(i = 0; i < nindent; i++) fprintf(fp, " ");}

	if(!zone->parent) {
		fprintf(fp, "== ROOT Zone ==\n");
	}
	else {
		PINDENT; fprintf(fp, "Zone[%lu] ", (unsigned long)zone->pos);
		GeVec3_PRINT(fp, zone->index);
		fprintf(fp, "\n");
	}
	geom = Dat_IList_IdxLookup(GEOM_TYPES, zone->voxel.geom);
	Deb_ASSERT(geom != NULL);
	PINDENT; fprintf(fp, "Geometry: "); fprintf(fp, "%s\n", geom->name);
	PINDENT; fprintf(fp, "Min: "); GeVec3_PRINT(fp, zone->voxel.min); fprintf(fp, "\n");
	PINDENT; fprintf(fp, "Cen: "); GeVec3_PRINT(fp, zone->voxel.cen); fprintf(fp, "\n");
	PINDENT; fprintf(fp, "Max: "); GeVec3_PRINT(fp, zone->voxel.max); fprintf(fp, "\n");

	if(DataFprintf) {
		PINDENT; fprintf(fp, "Data: ");
		DataFprintf(zone->data, fp);
	}

	fprintf(fp, "\n");

	if(zone->children) {
		PINDENT;
		fprintf(fp, "=== L%d Grid ", zone->level + 1);
		GeVec3_PRINT(fp, zone->naxes);
		fprintf(fp, "===\n");

		for(i = 0; i < zone->nchildren; i++)
			Zone_Fprintf(fp, zone->children[i], DataFprintf);
	}

	#undef PINDENT

	return;
}

/*----------------------------------------------------------------------------*/

void Zone_GrowChildren(Zone *zone, GeVec3_s naxes, void *(*DataAlloc)(const Zone *, const void *), const void *data_parms)
/* Grow 'children' from the current zone according
   to the number of divisions specified */
{
	size_t i;

	/* Calculate total number of children */
	zone->nchildren = naxes.x[0] * naxes.x[1] * naxes.x[2];
	Deb_ASSERT(zone->nchildren >= 1);
	zone->naxes = naxes;

	/* Allocate array of pointers for children */
	zone->children = Mem_CALLOC(zone->nchildren, zone->children);

	/* Loop through pointers and allocate each child */
	for(i = 0; i < zone->nchildren; i++) {
		/* Allocate sub-zone */
		zone->children[i] = Zone_Alloc(i, zone, DataAlloc, data_parms);

		/* Set sub-zone index */
		zone->children[i]->index = Ge_IelemToIndex(i, &naxes);

		/* Set sub-zone voxel */
		Zone_Set_child_voxel(zone, i);
	}

	return;
}

/*----------------------------------------------------------------------------*/

void Zone_Set_child_voxel(Zone *zone, size_t pos)
/* Calculated voxel according to parent zone */
{
	Zone *child = zone->children[pos];

	child->voxel = GeVox_GetSubVox(&zone->voxel, &child->index, &zone->naxes);

	return;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetMinLeaf(Zone *zone)
{
	if(zone->children)
		return Zone_GetMinLeaf(zone->children[0]);
	else
		return zone;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetMaxLeaf(Zone *zone)
{
	if(zone->children)
		return Zone_GetMaxLeaf(zone->children[zone->nchildren - 1]);
	else
		return zone;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetInner(Zone *zone)
{
	Zone *parent = zone->parent;

	if(parent) {
		if(zone->pos > 0) {
			return Zone_GetMaxLeaf(parent->children[zone->pos - 1]);
		}
		else if(parent->parent) {
			return Zone_GetInner(parent);
		}
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetOuter(Zone *zone)
{
	Zone *parent = zone->parent;

	if(parent) {
		if(zone->pos < parent->nchildren - 1) {
			return Zone_GetMinLeaf(parent->children[zone->pos + 1]);
		}
		else if(parent->parent) {
			return Zone_GetOuter(parent);
		}
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_AscendTree(Zone *zone)
{
	Zone *parent = zone->parent;

	if(parent) {
		if(zone->pos < parent->nchildren - 1)
			return Zone_GetMinLeaf(parent->children[zone->pos + 1]);
		else
			return parent;
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetLeaf(Zone *zone, size_t side, const GeVec3_d *pt, const GeRay *ray)
{
	Zone *leaf = NULL;

	switch(zone->voxel.geom) {
	  case GEOM_SPH1D:
	  	leaf = Zone_GetLeaf_sph1d(zone, side);
	  	break;
          
	  case GEOM_SPH3D:
	  	leaf = Zone_GetLeaf_sph3d(zone, side, pt, ray);
	  	break;
	  
	  case GEOM_REC3D:
	  	leaf = Zone_GetLeaf_rec3d(zone, side, pt);
	  	break;
	  	
	  case GEOM_CYL3D:
	  	leaf = Zone_GetLeaf_cyl3d(zone, side, pt, ray);
	  	break;
          
	  default: /* Shouldn't happen */
	  	Deb_ASSERT(0);
	}

	return leaf;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetLeaf_sph1d(Zone *zone, size_t side)
/* Locate leaf zone nearest to side */
{
	Zone *leaf = 0;

	if(side == 0) { /* Entering from inner sphere */
		leaf = Zone_GetMinLeaf(zone);
	}
	else if(side == 1) { /* Entering from outer sphere */
		leaf = Zone_GetMaxLeaf(zone);
	}
	else { /* Shouldn't happen */
		Deb_ASSERT(0);
	}

	return leaf;
}

/*----------------------------------------------------------------------------*/
#define RTHRESHOLD 1e-6
Zone *Zone_GetLeaf_sph3d(Zone *zone, size_t side, const GeVec3_d *pt, const GeRay *ray)
/* Locate leaf zone on side containing pt */
{
	/* If zone does not have children, zone is the leaf */
	if(!zone->children)
		return zone;

        GeVec3_d * SphPos = GeVec3_Cart2Sph(pt);
        double R     = SphPos->x[0];
        double theta = SphPos->x[1];
        double phi   = SphPos->x[2];
        GeVec3_s pos = GeVec3_INIT(0, 0, 0);
        
	/* Otherwise search for child containing pt */
	for(size_t i = 0; i < 3; i++) {
		size_t n = GeVec3_X(zone->naxes, i);
                size_t axis = side / 2;

		if(i == axis) {
			/* We already know the position on the driving axis -- this saves
			 * some time */
			GeVec3_X(pos, i) = (side % 2 == 0) ? 0 : n - 1;
		}
		else {
			/* Allocate array for binary search and load coordinate boundaries */
			double * array = Mem_CALLOC(n + 1, array);

			/* The Oth element should contain the lower bound */
                        GeVec3_s idx = GeVec3_INIT(0, 0, 0);
			Zone * child = Zone_CHILD(zone, idx);
			array[0] = GeVec3_X(child->voxel.min, i);

			/* All other elements are upper bounds */
			for(size_t j = 0; j < n; j++) {
				GeVec3_X(idx, i) = j;
				child = Zone_CHILD(zone, idx);
				array[j + 1] = GeVec3_X(child->voxel.max, i);
			}

			/* Do a binary search -- there MUST always be a match */
			if(i==0){ // Radius search
// 				Deb_PRINT("R=%g\n",R);
				GeVec3_X(pos, i) = 
					gsl_interp_bsearch(array, R, (size_t)0, n);
			}
			else if(i==1){ // Theta search
				if(R < RTHRESHOLD)
					theta = acos(GeRay_D(*ray,2));
// 				Deb_PRINT("R=%g theta=%g\n",R,theta);
				GeVec3_X(pos, i) = 
					gsl_interp_bsearch(array,theta, (size_t)0, n);
				//Deb_PRINT("solve: low=%g up=%g\n", array[GeVec3_X(pos, i)],array[GeVec3_X(pos, i)+1] );
			}
			else if(i==2){ // Phi search
				double Rc = R * sin(theta);
				if(Rc < RTHRESHOLD){
					GeVec3_d * SphDir = GeVec3_Cart2Sph(&ray->d);
                                        phi = SphDir->x[2];
				}
// 				Deb_PRINT("Rc=%g phi=%g\n",Rc,phi);
				GeVec3_X(pos, i) = 
					gsl_interp_bsearch(array,phi,(size_t)0, n);
			}
	
			/* Cleanup */
			free(array);
		}
	}

	/* pos should be at the appropriate child zone by now, recursively descend to
	 * leaf zone */
	return Zone_GetLeaf_sph3d(Zone_CHILD(zone, pos), side, pt, ray);
}

/*----------------------------------------------------------------------------*/
Zone *Zone_GetLeaf_rec3d(Zone *zone, size_t side, const GeVec3_d *pt)
/* Locate leaf zone on side containing pt */
{
	size_t i, j;
	GeVec3_s idx = GeVec3_INIT(0, 0, 0), pos = GeVec3_INIT(0, 0, 0);
	size_t axis = side / 2, n;
	Zone *child = 0;
	double *array;

	/* If zone does not have children, zone is the leaf */
	if(!zone->children)
		return zone;

	/* Otherwise search for child containing pt */
	for(i = 0; i < 3; i++) {
		n = GeVec3_X(zone->naxes, i);

		if(i == axis) {
			/* We already know the position on the driving axis -- this saves
			 * some time */
			GeVec3_X(pos, i) = (side % 2 == 0) ? 0 : n - 1;
		}
		else {
			/* Allocate array for binary search and load coordinate boundaries */
			array = Mem_CALLOC(n + 1, array);

			/* The Oth element should contain the lower bound */
			child = Zone_CHILD(zone, idx);
			array[0] = GeVec3_X(child->voxel.min, i);

			/* All other elements are upper bounds */
			for(j = 0; j < n; j++) {
				GeVec3_X(idx, i) = j;
				child = Zone_CHILD(zone, idx);
				array[j + 1] = GeVec3_X(child->voxel.max, i);
			}

			/* Do a binary search -- there MUST always be a match */
			GeVec3_X(pos, i) = gsl_interp_bsearch(array, GeVec3_X(*pt, i), (size_t)0, n);

			/* Cleanup */
			free(array);
		}
	}

	/* pos should be at the appropriate child zone by now, recursively descend to
	 * leaf zone */
	return Zone_GetLeaf_rec3d(Zone_CHILD(zone, pos), side, pt);
}

/*----------------------------------------------------------------------------*/
Zone *Zone_GetLeaf_cyl3d(Zone *zone, size_t side, const GeVec3_d *pt, const GeRay *ray)
/* Locate leaf zone on side containing pt */
{
	/* If zone does not have children, zone is the leaf */
	if(!zone->children)
		return zone;

        GeVec3_s pos = GeVec3_INIT(0, 0, 0);
        GeVec3_d * CylPos = GeVec3_Cart2Cyl(pt);
        double Rc  = CylPos->x[0];
        double phi = CylPos->x[1];
        double Hz  = CylPos->x[2];
        /* Otherwise search for child containing pt */
        for(size_t i = 0; i < 3; i++) {
                size_t n = GeVec3_X(zone->naxes, i);

                size_t axis = side / 2;
                if(i == axis) {
                        /* We already know the position on the driving axis -- this saves
                        * some time */
                        GeVec3_X(pos, i) = (side % 2 == 0) ? 0 : n - 1;
                }
                else {
                        /* Allocate array for binary search and load coordinate boundaries */
                        double *array;
                        array = Mem_CALLOC(n + 1, array);
                        
                        /* The Oth element should contain the lower bound */
                        GeVec3_s idx = GeVec3_INIT(0, 0, 0);
                        Zone * child = Zone_CHILD(zone, idx);
                        array[0] = GeVec3_X(child->voxel.min, i);

                        /* All other elements are upper bounds */
                        for(size_t j = 0; j < n; j++) {
                                GeVec3_X(idx, i) = j;
                                child = Zone_CHILD(zone, idx);
                                array[j + 1] = GeVec3_X(child->voxel.max, i);
                        }
                        
 
                        
                        /* Do a binary search -- there MUST always be a match */
                        switch (i){
                          case 0: // Rc search
                                  //Deb_PRINT("Rc=%g\n",Rc);
                                  GeVec3_X(pos, i) = 
                                          gsl_interp_bsearch(array, Rc, (size_t)0, n);
                                  break;
                        
                          case 1: // phi search
                                  //Deb_PRINT("Rc=%g phi=%g\n", Rc, phi);
                                  if ( Rc < RTHRESHOLD){
                                          GeVec3_d  * CylDirec = GeVec3_Cart2Cyl(&ray->d);
                                          phi = CylDirec->x[1];
                                  }
                                  GeVec3_X(pos, i) = 
                                          gsl_interp_bsearch(array, phi, (size_t)0, n);
                                  //Deb_PRINT("solve: low=%g up=%g\n", array[GeVec3_X(pos, i)],array[GeVec3_X(pos, i)+1] );
                                  break;
                        
                          case 2: // vertical height search
                                  //Deb_PRINT("Rc=%g phi=%g\n", Rc, phi);
                                  GeVec3_X(pos, i) = 
                                          gsl_interp_bsearch(array, Hz, (size_t)0, n);
                                  //Deb_PRINT("solve: low=%g up=%g\n", array[GeVec3_X(pos, i)],array[GeVec3_X(pos, i)+1] );
                                  break;
                          default:
                                  /* Not a valid axis */
                                  Deb_ASSERT(0);
                        }
                        
                        /* Cleanup */
                        free(array);
                }
        }

        /* pos should be at the appropriate child zone by now, 
           recursively descend to leaf zone                     */
        return Zone_GetLeaf_cyl3d(Zone_CHILD(zone, pos), side, pt, ray);
}
#undef RTHRESHOLD
/*----------------------------------------------------------------------------*/

Zone *Zone_GetNext_sph1d(Zone *zone, size_t *side)
{
	Zone *next = NULL;

	if( *side == 0) { /* Going inward */
		next = Zone_GetInner(zone);
                *side = 1;
	}
	else if( *side == 1) { /* Going outward */
		next = Zone_GetOuter(zone);
                *side = 0;
	}
	else { /* Shouldn't happen */
		Deb_ASSERT(0);
	}

	return next;
}

/*----------------------------------------------------------------------------*/
Zone *Zone_GetNext_sph3d(Zone *zone, size_t *side, const GeVec3_d *pt, const GeRay *ray)
/* Given face and point of intersection, locate next zone to enter */
{
	size_t axis = *side / 2; /* driving axis */
	GeVec3_s idx = zone->index;
	Zone *parent = zone->parent;
	int step, pos;	

// 	Deb_PRINT("checkpoint: Zone_GetNext_sph3d\n");
	/* Check for next-zones if not at root level */
	if(parent) {
		/* step is the increment/decrement on the driving axis */
		step = (*side % 2 == 0) ? -1 : 1;

		/* pos is the anticipated position on the driving axis,
		* which could be off the grid */
		pos = (int)GeVec3_X(idx, axis) + step;
// 		Deb_PRINT("pos=%d\n",pos);
		if((pos >= 0) && (pos < (int)GeVec3_X(parent->naxes, axis))) {
			/* If anticipated position is within current grid, descend
			 * to leaf of child zone at pos */
			GeVec3_X(idx, axis) = (size_t)pos;
			*side = (*side % 2 == 0) ? *side + 1 : *side - 1;
			return Zone_GetLeaf_sph3d( Zone_CHILD(parent, idx), *side, pt, ray );
		}
		else {
			/* Edge of current level reached, go to parent level */
// 			Deb_PRINT("getting upper level!\n");
			return Zone_GetNext_sph3d(parent, side, pt, ray);
		}
	}
	// no parent, root zone
	else if(axis==0){ 
		//if(side==0)
		//	return Zone_GetLeaf_sph3d(zone, 0, pt, ray );
		//else if(side==1)
			/* Outermost zone reached, no next-zone */
// 			Deb_PRINT("getting out!\n");
			return NULL;
	}
	else if(axis==1){
		return Zone_GetLeaf_sph3d(zone, *side, pt, ray );
	}
	else if(axis==2){
		*side = (*side % 2 == 0) ? *side + 1 : *side - 1;
		return Zone_GetLeaf_sph3d(zone, *side , pt, ray );
	}
	else
		/* Shouldn't happen */
		Deb_ASSERT(0);	

	#if 0
	Deb_PRINT("solve: newside=%d\n", *side);
	Deb_PAUSE();
	#endif
	
	
}

/*----------------------------------------------------------------------------*/
Zone *Zone_GetNext_rec3d(Zone *zone, size_t side, const GeVec3_d *pt)
/* Given face and point of intersection, locate next zone to enter */
{
	size_t axis = side / 2; /* driving axis */
	GeVec3_s idx = zone->index;
	Zone *parent = zone->parent;
	int step, pos;
	
	/* step is the increment/decrement on the driving axis */
	step = (side % 2 == 0) ? -1 : 1;

	/* pos is the anticipated position on the driving axis,
	 * which could be off the grid */
	pos = (int)GeVec3_X(idx, axis) + step;

	/* Check for next-zones if not at root level */
	if(parent) {
		if((pos >= 0) && (pos < (int)GeVec3_X(parent->naxes, axis))) {
			/* If anticipated position is within current grid, descend
			 * to leaf of child zone at pos */
			GeVec3_X(idx, axis) = (size_t)pos;
			return Zone_GetLeaf_rec3d(Zone_CHILD(parent, idx), side % 2 == 0 ? side + 1 : side - 1, pt);
		}
		else {
			/* Edge of current level reached, go to parent level */
			return Zone_GetNext_rec3d(parent, side, pt);
		}
	}

	/* Outermost zone reached, no next-zone */
	return NULL;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_GetNext_cyl3d(Zone *zone, size_t *side, const GeVec3_d *pt, const GeRay *ray)
/* Given face and point of intersection, locate next zone to enter */
{
	size_t axis = *side / 2; /* driving axis */
	Zone *parent = zone->parent;	

// 	Deb_PRINT("checkpoint: Zone_GetNext_sph3d\n");
	/* Check for next-zones if not at root level */
	if(parent) {
		/* step is the increment/decrement on the driving axis */
		int step = (*side % 2 == 0) ? -1 : 1;

		/* pos is the anticipated position on the driving axis,
		* which could be off the grid */
                GeVec3_s idx = zone->index;
		int pos = (int)GeVec3_X(idx, axis) + step;
// 		Deb_PRINT("pos=%d\n",pos);
		if( 0 <= pos < (int)GeVec3_X(parent->naxes, axis)) {
			/* If anticipated position is within current grid, descend
			 * to leaf of child zone at pos */
			GeVec3_X(idx, axis) = (size_t)pos;
			*side = (size_t)((int) *side - step);
			return Zone_GetLeaf_cyl3d( Zone_CHILD(parent, idx), *side, pt, ray );
		}
		else {
			/* Edge of current level reached, go to parent level */
// 			Deb_PRINT("getting upper level!\n");
			return Zone_GetNext_cyl3d(parent, side, pt, ray);
		}
	}
	// no parent, root zone
	else if( axis == 0 ){ 
		switch (*side){
                  case 0 :
                          return Zone_GetLeaf_cyl3d(zone, *side, pt, ray );
                  case 1 :
                          /* Outermost zone reached, no next-zone */
                          //Deb_PRINT("getting out!\n");
                          return NULL;
                }
        }
	else if( axis == 1 ){
		*side = (*side % 2 == 0) ? *side + 1 : *side - 1;
		return Zone_GetLeaf_cyl3d(zone, *side, pt, ray );
	}
	else if( axis == 2 ){
		return NULL;
	}

	else
		/* Shouldn't happen */
		Deb_ASSERT(0);		
}
/*----------------------------------------------------------------------------*/
#if Sp_MIRSUPPORT
void Zone_Cpgplot(Zone *zone, GeCam *cam)
{
	size_t i;

	/* plot current zone */
	GeVox_Cpgplot(&zone->voxel, cam);

	/* plot children */
	for(i = 0; i < zone->nchildren; i++)
		Zone_Cpgplot(zone->children[i], cam);
}
#endif
/*----------------------------------------------------------------------------*/

Zone *Zone_GetNext(Zone *zone, size_t *side, const GeRay *ray)
{
	Zone *next = NULL;

	switch(zone->voxel.geom) {
	  case GEOM_SPH1D:
	  	next = Zone_GetNext_sph1d(zone, side);
	  	break;
          
	  case GEOM_SPH3D:
	  	next = Zone_GetNext_sph3d(zone, side, &(ray->e), ray);
	  	break;
	  
	  case GEOM_REC3D:
	  	next = Zone_GetNext_rec3d(zone, *side, &(ray->e));
	  	break;
	  	
	  case GEOM_CYL3D:
	  	next = Zone_GetNext_cyl3d(zone, side, &(ray->e), ray);
	  	break;
          
	  default: /* Shouldn't happen */
	  	Deb_ASSERT(0);
	}

	return next;
}

/* esc 09May06: Deprecate the following code */
/*----------------------------------------------------------------------------*/

size_t Zone_Fwrite(const Zone *zone, size_t (*DataFwrite)(void *data, FILE *fp), FILE *fp)
{
	size_t i, nbytes = 0;

	nbytes += Mem_FWRITE(zone, 1, fp);

	if(DataFwrite) {
		nbytes += DataFwrite(zone->data, fp);
	}

	if(zone->nchildren > 0)
		Deb_ASSERT(zone->children != NULL);

	for(i = 0; i < zone->nchildren; i++) {
		nbytes += Zone_Fwrite(zone->children[i], DataFwrite, fp);
	}

	return nbytes;
}

/*----------------------------------------------------------------------------*/

Zone *Zone_Fread(void *(*DataAlloc)(const void *data_parms), const void *data_parms, 
	size_t (*DataFread)(void *data, FILE *fp), FILE *fp)
{
	size_t i;
	Zone *zone;

	zone = Mem_CALLOC(1, zone);
	Mem_FREAD(zone, 1, fp);
	
	if(zone->data) {
		Deb_ASSERT(DataFread != NULL);

		zone->data = DataAlloc(data_parms);
		DataFread(zone->data, fp);
	}

	if(zone->nchildren > 0) {
		zone->children = Mem_CALLOC(zone->nchildren, zone->children);

		for(i = 0; i < zone->nchildren; i++) {
			zone->children[i] = Zone_Fread(DataAlloc, data_parms, DataFread, fp);
			zone->children[i]->parent = zone;
		}
	}

	return zone;
}




