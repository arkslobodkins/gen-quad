#include <stdlib.h>
#include <stdio.h>
#include "glist.h"

glist *glist_create()
{
   glist *new_list = (glist *)malloc(sizeof (glist));
   new_list->size = 0;
   new_list->first = NULL;
   new_list->last = NULL;
   return new_list;
}


void glist_free(glist *list)
{
   struct node *curr = list->first;
   struct node *next;

   while (curr != NULL)
   {
      next = curr->next;
      free(curr->data);
      free(curr);
      curr = next;
   }

   free(list);
}


void glist_push(glist *list, void *data)
{
   if (list == NULL)
   {
      fprintf(stderr, "glist_push: list is null\n");
      return;
   }

   struct node *new_node = (node*)malloc(sizeof(node));
   new_node->data = data;
   new_node->next = NULL;
   if(list->size == 0)
   {
      list->first = new_node;
      list->first->prev = NULL;
      list->first->next = NULL;
      list->last = list->first;
   }
   else
   {
      list->last->next = new_node;
      list->last->prev = list->last;
      list->last = new_node;
   }
   ++list->size;

}


void glist_print(glist *list, void (*print)(void *))
{
   if(list == NULL || list->first == NULL)
   {

      fprintf(stderr, "glist_print: list is null\n");
      return;
   }

   if(list->size == 0)
   {
      fprintf(stderr, "glist_print: list is empty\n");
      return;
   }

   struct node *curr = list->first;
   for(int i = 0; i < list->size; ++i)
   {
      print(curr->data);
      curr = curr->next;
   }
}

