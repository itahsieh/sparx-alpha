#ifndef __MOLEC_H__
#define __MOLEC_H__

/* This is part of the readout API for molecular data files
 * following the LAMDA format.
 */

#define COLL_H2 1
#define COLL_PH2 2
#define COLL_OH2 3
#define COLL_E 4
#define COLL_H 5
#define COLL_HE 6

typedef struct {
  size_t id;
  double E, g;
  char *qstate,*F1,*F; /* Quantum state */
} MolLev;

typedef struct {
  size_t id;
  size_t up, lo;
  double freq, A_ul, B_ul, B_lu;

} MolTrRad;

typedef struct {
  size_t id;
  size_t up, lo;
  int overlap;
  double RelativeVel;
} MolOverlap;

typedef struct {
  size_t id;
  size_t up, lo;
  double *K_ul;
} MolTrCol;

typedef struct {
  char *ref;
  int species;
  size_t ntr, ntmp;
  MolTrCol **tr;
  double *tmp;
} MolColPart;

typedef struct {
  char *name, /* File name */
       *chemname, /* Chemical name */
       *qnum; /* Quantum numbers */
  double weight; /* Molecular weight */
  size_t nlev, nrad, ncol;
  MolLev **lev; /* Array of molecular levels */
  MolTrRad **rad; /* Array of radiative transitions */
  MolColPart **col; /* Array of collissional transitions */
  MolOverlap **OL; /* Array of overlapping */
} Molec;

Molec *Mol_Alloc(size_t nlev);
void Mol_AllocRad(Molec *molec, size_t nrad);
void Mol_AllocCol(Molec *molec, size_t ntmp, size_t ntr);
void Mol_Free(void *ptr);
void Mol_Fprintf(FILE *fp, const Molec *mol);

Molec *Mol_ReadLamda(FILE *fp, const char *fname, const char *name);
void Mol_FwriteBinary(const Molec *mol, FILE *fp);
Molec *Mol_FreadBinary(FILE *fp);


#endif










