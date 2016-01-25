#ifndef __DATA_STRUCTS_H__
#define __DATA_STRUCTS_H__

/* Node structure for singly & doubly-linked lists */
struct LNode {
	char *name;
	void *value;
	void (*freef)(void *);
	struct LNode *next, *prev;
};

typedef struct DatINode {
	const char *name;
	int idx;
} DatINode;

typedef struct LNode LNode;

DatINode *Dat_IList_NameLookup(DatINode *list, const char *name);
DatINode *Dat_IList_IdxLookup(DatINode *list, int idx);
LNode *Dat_Llst_Lookup(LNode *list, const char *name);
LNode *Dat_Llst_Push(LNode *list, void *value, void (*freef)(void *));
LNode *Dat_Llst_Insert(LNode *list, const char *name, void *value, void (*freef)(void *));
void Dat_LNode_Free(void *node);
void Dat_Llst_Free(void *list);

#endif

