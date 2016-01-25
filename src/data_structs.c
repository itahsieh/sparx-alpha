#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_structs.h"
#include "memory.h"
#include "debug.h"

/*----------------------------------------------------------------------------*/

DatINode *Dat_IList_NameLookup(DatINode *list, const char *name)
{
	size_t i;

	for(i = 0; list[i].name; i++) {
		if(!strcmp(list[i].name, name))
			return &list[i];
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

DatINode *Dat_IList_IdxLookup(DatINode *list, int idx)
{
	size_t i;

	for(i = 0; list[i].name; i++) {
		if(list[i].idx == idx)
			return &list[i];
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

LNode *Dat_Llst_Lookup(LNode *list, const char *name)
/* Lookup LNode associated with name in a linked list of LNodes */
{
	LNode *np;

	for(np = list; np; np = np->next) {
		if(!strcmp(np->name, name))
			return np;
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

LNode *Dat_Llst_Push(LNode *list, void *value, void (*freef)(void *))
/* Insert value associated with name in a linked list of LNodes */ 
{
	LNode *np;

	/* Name not found, create new node */
	np = Mem_CALLOC(1, np);

	/* Assign new value and freef */
	np->value = value;
	np->freef = freef;

	/* If list is not empty, prepend np to list */
	if(list) {
		np->next = list;
		np->next->prev = np;
	}

	return np;
}

/*----------------------------------------------------------------------------*/

LNode *Dat_Llst_Insert(LNode *list, const char *name, void *value, void (*freef)(void *))
/* Insert value associated with name in a linked list of LNodes */ 
{
	LNode *np, *newlst;

	/* Check if name already exists in list */
	np = Dat_Llst_Lookup(list, name);

	if(!np) {
		/* Name not found, create new node */
		np = Mem_CALLOC(1, np);

		/* Assign name */
		np->name = Mem_STRDUP(name);

		/* If list is not empty, prepend np to list */
		if(list) {
			np->next = list;
			np->next->prev = np;
		}

		/* list now begins with np */
		newlst = np;
	}
	else {
		/* Node with name already exists, free associated value if np->freef is specified */
		if(np->freef)
			np->freef(np->value);

		/* newlst still begins with list */
		newlst = list;
	}

	/* Assign new value and freef */
	np->value = value;
	np->freef = freef;

	return newlst;
}

/*----------------------------------------------------------------------------*/

void Dat_LNode_Free(void *node)
/* Free a single LNode */
{
	LNode *np = node;

	free(np->name);
	if(np->freef)
		np->freef(np->value);
	free(np);

	return;
}

/*----------------------------------------------------------------------------*/

void Dat_Llst_Free(void *list)
/* Free a linked list of LNodes */
{
	LNode *lp = list, *np;

	if(lp->next) {
		/* Free current node and move on to next node */
		np = lp;
		lp = lp->next;
		Dat_LNode_Free(np);
		Dat_Llst_Free(lp);
	}
	else {
		/* Last node of list reached, free lp and that's it */
		Dat_LNode_Free(lp);
	}

	return;
}











