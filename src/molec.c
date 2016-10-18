#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "error.h"
#include "memory.h"
#include "molec.h"
#include "debug.h"

static char *_get_nline_rel(FILE *fp, size_t nline);
static char *_get_nline(FILE *fp, size_t nline);

#define get_nline_rel(fp, nline)\
	_get_nline_rel((fp), (size_t)nline)

#define get_nline(fp, nline)\
	_get_nline((fp), (size_t)nline)

const char *Mol_species[6] = {
	"H2", "p-H2", "o-H2", "e", "H", "He"
};

const char **Mol_species_p = Mol_species;

/*----------------------------------------------------------------------------*/

static char *_get_nline_rel(FILE *fp, size_t nline)
{
  static char buffer[BUFSIZ] = "";
  char *sptr = buffer;
  size_t i;

  for(i = 0; sptr && i < nline; i++)
    sptr = fgets(buffer, BUFSIZ, fp);

  return sptr;
}

/*----------------------------------------------------------------------------*/

static char *_get_nline(FILE *fp, size_t nline)
{
  rewind(fp);

  return get_nline_rel(fp, nline);
}

/*----------------------------------------------------------------------------*/

Molec *Mol_Alloc(size_t nlev)
/* Allocate a bare-bones molecule that contains only level information */
{
	size_t i;
	Molec *mp = Mem_CALLOC(1, mp);

	mp->nlev = nlev;
	mp->lev = Mem_CALLOC(nlev, mp->lev);

	for(i = 0; i < nlev; i++) {
		mp->lev[i] = Mem_CALLOC(1, mp->lev[i]);
	}

	return mp;
}

/*----------------------------------------------------------------------------*/

void Mol_AllocRad(Molec *molec, size_t nrad)
/* Allocate radiative transitions for molec */
{
	size_t i;
	Deb_ASSERT(molec->rad == NULL);

	molec->nrad = nrad;
	molec->rad = Mem_CALLOC(nrad, molec->rad);

	for(i = 0; i < nrad; i++) {
		molec->rad[i] = Mem_CALLOC(1, molec->rad[i]);
	}

	return;
}

/*----------------------------------------------------------------------------*/

void Mol_AllocCol(Molec *molec, size_t ntmp, size_t ntr)
/* 'Add' a collisional partner to molec */
{
	size_t i, icol = molec->ncol;
	MolColPart *cp;

	/* Resize collisional partner array */
	molec->ncol += 1;
	molec->col = Mem_REALLOC(molec->ncol, molec->col);

	/* Allocate collisional partner */
	cp = molec->col[icol] = Mem_CALLOC(1, molec->col[icol]);

	/* Allocate temperature array */
	cp->ntmp = ntmp;
	cp->tmp = Mem_CALLOC(ntmp, cp->tmp);

	/* Allocate transitions array */
	cp->ntr = ntr;
	cp->tr = Mem_CALLOC(ntr, cp->tr);

	/* Allocate downward rates array for each transition */
	for(i = 0; i < ntr; i++) {
		cp->tr[i] = Mem_CALLOC(1, cp->tr[i]);
		cp->tr[i]->K_ul = Mem_CALLOC(ntmp, cp->tr[i]->K_ul);
	}

	return;
}

/*----------------------------------------------------------------------------*/

void Mol_Free(void *ptr)
{
	Molec *mp = ptr;
	size_t i, j;

	#define FREE(ptr)\
		{if(ptr) free(ptr);}

	if(mp->name)
		free(mp->name);

	if(mp->chemname)
		free(mp->chemname);

	if(mp->qnum)
		free(mp->qnum);

	/* Free levels */
	if(mp->lev) {
		for(i = 0; i < mp->nlev; i++) {
			if(mp->lev[i]) {
				FREE(mp->lev[i]->qstate);
				free(mp->lev[i]);
			}
		}
		free(mp->lev);
	}

	/* Free radiative transitions */
	if(mp->rad) {
		for(i = 0; i < mp->nrad; i++) {
			free(mp->rad[i]);
		}
		free(mp->rad);
	}

	/* Free collisional partners */
	if(mp->col) {
		for(i = 0; i < mp->ncol; i++) {
			FREE(mp->col[i]->ref);
			FREE(mp->col[i]->tmp);
			for(j = 0; j < mp->col[i]->ntr; j++) {
				free(mp->col[i]->tr[j]->K_ul);
				free(mp->col[i]->tr[j]);
			}
			free(mp->col[i]->tr);
			free(mp->col[i]);
		}
		free(mp->col);
	}

	free(mp);

	#undef FREE

	return;
}

/*----------------------------------------------------------------------------*/

void Mol_Fprintf(FILE *fp, const Molec *mol)
{
	size_t i, j, k;
	#define PRINT(...)\
		fprintf(fp, __VA_ARGS__);

	PRINT("Name: `%s'\n", mol->name);
	PRINT("Chemical name: `%s'\n", mol->chemname);
	PRINT("Quantum numbers: `%s'\n", mol->qnum);
	PRINT("Molecular weight: %g\n", mol->weight);

	/* Print energy levels */
	PRINT("Number of energy levels: %lu\n", (unsigned long)mol->nlev);
	PRINT("%10s|%15s|%15s|%15s\n", "LEVEL", "ENERGIES", "WEIGHT", mol->qnum);
	PRINT("----------|---------------|---------------|----------------\n");
	for(i = 0; i < mol->nlev; i++) {
		PRINT("%10lu|%15g|%15g|%15s\n", (unsigned long)mol->lev[i]->id, mol->lev[i]->E, mol->lev[i]->g, mol->lev[i]->qstate);
	}

	/* Print radiative transitions */
	if(mol->rad) {
		PRINT("Number of radiative transitions: %lu\n", (unsigned long)mol->nrad);
		PRINT("%10s|%10s|%10s|%15s|%15s|%15s|%15s\n", "TRANS", "UP", "LOW", "A_ul", "B_lu", "B_ul", "FREQ");
		PRINT("----------|----------|----------|---------------|---------------|---------------|----------------\n");
		for(i = 0; i < mol->nrad; i++) {
			PRINT("%10lu|%10lu|%10lu|%15g|%15g|%15g|%15g\n",
				(unsigned long)mol->rad[i]->id,
				(unsigned long)mol->rad[i]->up,
				(unsigned long)mol->rad[i]->lo,
				mol->rad[i]->A_ul,
				mol->rad[i]->B_lu,
				mol->rad[i]->B_ul,
				mol->rad[i]->freq);
		}
	}

	/* Print collissional partners */
	if(mol->col) {
		PRINT("Number of collissional partners: %lu\n", (unsigned long)mol->ncol);

		for(i = 0; i < mol->ncol; i++) {
			PRINT("Partner %lu:\n", (unsigned long)(i + 1));
			PRINT("  Species: %s\n", Mol_species[mol->col[i]->species - 1]);
			PRINT("  Reference: %s\n", mol->col[i]->ref);
			PRINT("  Number of temperatures: %lu\n", (unsigned long)mol->col[i]->ntmp);
			PRINT("  Number of collisional transitions: %lu\n", (unsigned long)mol->col[i]->ntr);

			PRINT("%5s|%5s|%5s", "TRANS", "UP", "LOW");
			for(j = 0; j < mol->col[i]->ntmp; j++)
				PRINT("|%9.2gK", mol->col[i]->tmp[j]);
			PRINT("\n");

			PRINT("-----|-----|-----");
			for(j = 0; j < mol->col[i]->ntmp; j++)
				PRINT("|----------");
			PRINT("\n");

			for(j = 0; j < mol->col[i]->ntr; j++) {
				PRINT("%5lu|%5lu|%5lu",
					(unsigned long)j,
					(unsigned long)mol->col[i]->tr[j]->up,
					(unsigned long)mol->col[i]->tr[j]->lo);
				for(k = 0; k < mol->col[i]->ntmp; k++) {
					PRINT("|%10.2e", mol->col[i]->tr[j]->K_ul[k])
				}
				PRINT("\n");
			}
		}
	}

	#undef PRINT
	
}

/*----------------------------------------------------------------------------*/

Molec *Mol_ReadLamda(FILE *fp, const char *fname, const char *name)
/* Read molecular parameters from a LAMDA format data file */
{
	int status = 0, species = 0;
	Molec *mp = 0;
	char *sptr = 0, *sptr2 = 0, *ref = 0, **list = 0;
	size_t n, ncol = 0, ntmp = 0, ntr = 0, i, j, k;

	/* Get number of levels from line 6 and allocate molecule */
	if(!status && !(sptr = get_nline(fp, 6)))
		status = Err_SETSTRING("In file `%s': Number of levels not found", fname);

	if(!status) {
		n = (size_t)atoi(sptr);
		if(n == 0)
			status = Err_SETSTRING("In file `%s': Number of levels must be > 0", fname);
		else
			mp = Mol_Alloc(n);
	}

	/* Assign name */
	if(!status)
		mp->name = Mem_STRDUP(name);

	/* Get chemical name from line 2 */
	if(!status && !(sptr = get_nline(fp, 2)))
		status = Err_SETSTRING("In file `%s': Chemical name not found", fname);

	if(!status)
		mp->chemname = Mem_STRDUP(MemStr_Stripwhite(sptr));

	/* Get molecular weight from line 4 */
	if(!status && !(sptr = get_nline(fp, 4)))
		status = Err_SETSTRING("In file `%s': Molecular weight not found", fname);

	if(!status) {
		mp->weight = atof(sptr);

		if(mp->weight <= 0.0)
			status = Err_SETSTRING("In file `%s': Molecular weight must be > 0", fname);
	}

	/* Get quantum numbers from last column of line 7 */
	if(!status && !(sptr = get_nline(fp, 7)))
		status = Err_SETSTRING("In file `%s': Quantum numbers not found", fname);

	for(i = 0; !status && i < 4; i++) {
		if(i == 0)
			sptr2 = strtok(sptr, "+");
		else if(i == 3)
			sptr2 = strtok(NULL, "\n");
		else
			sptr2 = strtok(NULL, "+");

		if(!sptr2)
			status = Err_SETSTRING("In file `%s': Quantum numbers not found", fname);
	}

	if(!status)
		mp->qnum = Mem_STRDUP(MemStr_Stripwhite(sptr2));

	/* get level info from nlev lines after line 7 */
	if(!status && !(sptr = get_nline(fp, 7)))
		status = Err_SETSTRING("In file `%s': Energy levels not found", fname);

	for(i = 0; i < mp->nlev; i++) {
		if(!status && !(sptr = get_nline_rel(fp, 1)))
			status = Err_SETSTRING("In file `%s': Anticipated %d rows in levels table, but got %d", fname, mp->nlev, i);

		if(!status) { /* Get first column (index) */
			if((sptr2 = strtok(sptr, "\t ")))
				mp->lev[i]->id = i;
			else
				status = Err_SETSTRING("In file `%s': Error reading level index at row %d", fname, i + 1);
		}

		if(!status) { /* Get second column (energy) */
			if((sptr2 = strtok(NULL, "\t ")))
				mp->lev[i]->E = atof(sptr2);
			else
				status = Err_SETSTRING("In file `%s': Error reading level energy at row %d", fname, i + 1);
		}

		if(!status) { /* Get third column (statistical weight) */
			if((sptr2 = strtok(NULL, "\t ")))
				mp->lev[i]->g = atof(sptr2);
			else
				status = Err_SETSTRING("In file `%s': Error reading level statistical weight at row %d", fname, i + 1);
		}

		if(!status) { /* Get fourth column (quantum state) */
			if((sptr2 = strtok(NULL, "\0")))
				mp->lev[i]->qstate = Mem_STRDUP(MemStr_Stripwhite(sptr2));
			else
				status = Err_SETSTRING("In file `%s': Error reading level quantum state at row %d", fname, i + 1);
		}
	}

	/* Get number of radiative transitions from line 9 + nlev */
	if(!status) {
		if(!(sptr = get_nline(fp, 9 + mp->nlev)))
			status = Err_SETSTRING("In file `%s': Number of radiative transitions not found", fname);

		if(!status) {
			n = (size_t)atoi(sptr);
			if(n == 0)
				status = Err_SETSTRING("In file `%s': Number of radiative transitions must be > 0", fname);
			else
				Mol_AllocRad(mp, n);
		}
	}

	/* load line info from nline lines after line 10 + nlev */
	if(!status) {
		if(!(sptr = get_nline(fp, 10 + mp->nlev)))
			status = Err_SETSTRING("In file `%s': Radiative transitions not found", fname);

		if(!status) {
			for(i = 0; i < mp->nrad; i++) {
				if(!status && !(sptr = get_nline_rel(fp, 1)))
					status = Err_SETSTRING("In file `%s': Anticipated %d rows in radiative transitions table, but got %d", fname, mp->nrad, i);

				if(!status) {
					list = MemStr_Split(MemStr_Stripwhite(sptr), "\t ", &n);
					if(n != 6)
						status = Err_SETSTRING("In file `%s': Anticipated %d columns in radiative transitions table, but got %d", fname, 6, n);
				}

				if(!status) {
					mp->rad[i]->id = i;
					mp->rad[i]->up = (size_t)(atoi(list[1]) - 1);
					mp->rad[i]->lo = (size_t)(atoi(list[2]) - 1);
					mp->rad[i]->A_ul = atof(list[3]);
					mp->rad[i]->freq = atof(list[4]);
                                        
				}
				MemStr_FreeList(list, n);
			}
		}
	}

	/* Get number of collisonal partners from 12 + nlev + nrad */
	if(!status) {
		if(!(sptr = get_nline(fp, 12 + mp->nlev + mp->nrad)))
			status = Err_SETSTRING("In file `%s': Number of collisonal partners not found", fname);

		if(!status) {
			ncol = (size_t)atoi(sptr);
			if(ncol == 0)
				status = Err_SETSTRING("In file `%s': Number of collisional partners must be > 0", fname);
		}
	}

	/* Loop through collisional partners */
	for(i = 0; i < ncol; i++) {
		/* Get species & reference from 2 lines below */
		if(!status && !(sptr = get_nline_rel(fp, 2))) {
			status = Err_SETSTRING("In file `%s': Collisional partner reference not found", fname);
		}

		if(!status) {
			ref = Mem_STRDUP(MemStr_Stripwhite(sptr));

			/* Extract species code */
			sptr2 = strtok(sptr, " ");
			if(sptr2)
				species = atoi(sptr2);
			if(species < 1 || species > 6)
				status = Err_SETSTRING("In file `%s': Collisional partner species must be >= 1 and <= 6 (string=%s)", fname, sptr2);
		}

		if(!status) {
			/* Extrance reference */
			sptr = strtok(NULL, "\0");
			free(ref);
			ref = Mem_STRDUP(MemStr_Stripwhite(sptr));
		}

		/* Get number of transitions from 2 lines below */
		if(!status && !(sptr = get_nline_rel(fp, 2))) {
			status = Err_SETSTRING("In file `%s': Number of collisional transitions not found", fname);
		}

		if(!status) {
			ntr = (size_t)atoi(sptr);
			if(ntr == 0)
				status = Err_SETSTRING("In file `%s': Number of collisional transitions must be > 0", fname);
		}

		/* Get number of temperatures from 2 lines below */
		if(!status && !(sptr = get_nline_rel(fp, 2))) {
			status = Err_SETSTRING("In file `%s': Collisional partner reference not found", fname);
		}

		if(!status) {
			ntmp = (size_t)atoi(sptr);
			if(ntmp == 0)
				status = Err_SETSTRING("In file `%s': Number of collisional transitions must be > 0", fname);
		}

		/* Allocate collisional partner */
		if(!status) {
			Mol_AllocCol(mp, ntmp, ntr);

			/* Load reference */
			mp->col[i]->ref = ref;
			ref = 0;
			mp->col[i]->species = species;
		}

		/* Get list of temperatures 2 lines below */
		if(!status && !(sptr = get_nline_rel(fp, 2))) {
			status = Err_SETSTRING("In file `%s': Collisional partner reference not found", fname);
		}

		if(!status) {
			list = MemStr_Split(MemStr_Stripwhite(sptr), "\t ", &n);
			if(n != ntmp)
				status = Err_SETSTRING("In file `%s': Wrong number of temperatures (should be %d, only %d found)", fname, ntmp, n);

			if(!status) {
				for(j = 0; !status && j < ntmp; j++) {
					mp->col[i]->tmp[j] = atof(list[j]);
					if(mp->col[i]->tmp[j] <= 0)
						status = Err_SETSTRING("In file `%s': Temperature must be > 0", fname);
				}
			}
			MemStr_FreeList(list, n);
		}

		/* Skip one line and start getting K_ul ntr lines below */
		if(!status && !(sptr = get_nline_rel(fp, 1)))
			status = Err_SETSTRING("In file `%s': Can't find collisional rate coefficients table", fname);

		for(j = 0; j < ntr; j++) {
			if(!status && !(sptr = get_nline_rel(fp, 1)))
				status = Err_SETSTRING("In file `%s': Anticipated %d rows of collisional rate coeff table, but got %d", fname, ntr, j);

			if(!status) {
				list = MemStr_Split(MemStr_Stripwhite(sptr), "\t ", &n);
				if(n != ntmp + 3)
					status = Err_SETSTRING("In file `%s': Anticipated %d columns in collisional rate coeff table, but got %d", fname, ntmp + 3, n);
			}

			if(!status) {
				mp->col[i]->tr[j]->up = (size_t)atoi(list[1]) - 1;
				mp->col[i]->tr[j]->lo = (size_t)atoi(list[2]) - 1;
				for(k = 0; k < ntmp; k++)
					mp->col[i]->tr[j]->K_ul[k] = atof(list[k + 3]);
			}

			if(list) {
				MemStr_FreeList(list, n);
				list = 0;
			}
		}
	}

	if(ref)
		free(ref);

	/* Cleanup if things went wrong */
	if(status) {
		Mol_Free(mp);
		mp = NULL;
	}

	return mp;
}


/*----------------------------------------------------------------------------*/
void Mol_FwriteBinary(const Molec *mol, FILE *fp)
/* Write the Molec structure to binary file */
{
	size_t i;

	/* Write base structure */
	Mem_FWRITE(mol, 1, fp);
	Mem_FputStr(mol->name, fp);
	Mem_FputStr(mol->chemname, fp);
	Mem_FputStr(mol->qnum, fp);

	/* Write levels */
	for(i = 0; i < mol->nlev; i++) {
		Mem_FWRITE(mol->lev[i], 1, fp);
		Mem_FputStr(mol->lev[i]->qstate, fp);
	}

	/* Write number of radiative transitions */
	Mem_FwriteSize_t(mol->nrad, fp);

	/* Write radiative transitions */
	for(i = 0; i < mol->nrad; i++) {
		Mem_FWRITE(mol->rad[i], 1, fp);
	}

	return;
}

/*----------------------------------------------------------------------------*/

Molec *Mol_FreadBinary(FILE *fp)
/* Read the Molec structure from binary file */
{
	Molec *tmp = Mem_CALLOC(1, tmp), *mol;
	size_t i;

	/* Allocate and read base structure */
	Mem_FREAD(tmp, 1, fp);
	mol = Mol_Alloc(tmp->nlev);
	mol->name = Mem_FgetStr(fp);
	mol->chemname = Mem_FgetStr(fp);
	mol->qnum = Mem_FgetStr(fp);
	mol->weight = tmp->weight;
	free(tmp);

	/* Read levels */
	for(i = 0; i < mol->nlev; i++) {
		Mem_FREAD(mol->lev[i], 1, fp);
		mol->lev[i]->qstate = Mem_FgetStr(fp);
	}

	/* Read number of radiative transitions */
	Mem_FreadSize_t(&mol->nrad, fp);

	/* Allocate radiative transitions */
	Mol_AllocRad(mol, mol->nrad);

	/* Read radiative transitions */
	for(i = 0; i < mol->nrad; i++) {
		Mem_FREAD(mol->rad[i], 1, fp);
	}

	return mol;
}

/*----------------------------------------------------------------------------*/
















