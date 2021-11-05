#ifndef GLIST_H
#define GLIST_H

typedef struct node
{
   void *data;
   struct node *prev;
   struct node *next;
} node;
typedef struct node node;

typedef struct
{
   int size;
   node *first;
   node *last;
}glist;

glist *glist_create();
void glist_free(glist *list);
void glist_push(glist *list, void *data);
void glist_print(glist *list, void (*print)(void *data));


#endif
