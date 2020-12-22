#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "allocations.h"

#define SIZE 262139 //it must be a prime number!

typedef struct alloc alloc;

struct alloc{
   char* name;   
   void* address;
   alloc *next, *prev;
};

typedef struct {
    int n;
    alloc *head, *tail;
} alloc_list;


alloc_list** hashArray = NULL;

unsigned long int hashCode(unsigned long int key) {

  //int div = 251;

  /*size_t shift = (size_t)log2(1 + sizeof(key));
  return (size_t)(key) >> shift;*/

  //if ((menor = 0) || (key > menor))
  //  menor = key;

  //printf("%lu\n",key);

  //key = (key/* - 25886736*/)*div;

  /*key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);*/

  return key % SIZE;
}

void insert_alloc_list(alloc_list** list, alloc*item) {
	if ((*list) == NULL) {
		(*list) = (alloc_list*)malloc(sizeof(alloc_list));
		(*list)->n = 1;
		(*list)->head = (*list)->tail = item;
	}
	else if ((*list)->n == 0) {
		(*list)->head = item;
		(*list)->tail = item;
		(*list)->n=1;
	}
	else {
		item->prev = (*list)->tail;
		(*list)->tail->next = item;
		(*list)->tail = item;
		(*list)->n++;
	}
}

void remove_alloc_list(alloc_list* list, void* address) {
	//printf("[removing %lu]\n",(unsigned long int)address);
	alloc*aux;
	if(list==NULL){
		//printf("[Free in non-allocated memory!]\n");
	}
	else if (list->n == 0) {
		//printf("[Free in non-allocated memory!]\n");
	}
	else if (list->n == 1) {
		if (list->head->address == address) {
			//free(list->head->address);
			free(list->head->name);
			free(list->head);
			list->head = list->tail = NULL;
			list->n = 0;
		}
		//else
			//printf("[Free in non-allocated memory!]\n");
	}
	else if (list->head->address == address) {
		aux = list->head;
		list->head = aux->next;
		list->head->prev = NULL;
		//free(aux->address);
		free(aux->name);
		free(aux);
		list->n--;
	}
	else if (list->tail->address == address) {
		aux = list->tail;
		list->tail = aux->prev;
		list->tail->next = NULL;
		//free(aux->address);
		free(aux->name);
		free(aux);
		list->n--;
	}
	else {
		for(aux = list->head; aux; aux=aux->next) {
			if(aux->address == address) {
				aux->prev->next = aux->next;
				aux->next->prev = aux->prev;
				//free(aux->address);
				free(aux->name);
				free(aux);
				list->n--;
				return;
			}
		}
		//printf("[Free in non-allocated memory!]\n");
	}

}

void insert_hash_table(void*address,char*name){
	if(hashArray == NULL)
		hashArray = (alloc_list**)calloc(SIZE,sizeof(alloc_list*));
	alloc* item = (alloc*)malloc(sizeof(alloc));
	item->address = address;
	item->name = (char*)malloc((strlen(name)+1)*sizeof(char));
	strcpy(item->name,name);
	item->next = item->prev = NULL;
	unsigned long int ind = hashCode((unsigned long int)address);
	insert_alloc_list(&hashArray[ind],item);
}

void remove_hash_table(void*address) {
	unsigned long int ind = hashCode((unsigned long int)address);
	remove_alloc_list(hashArray[ind],address);
}

void print_alloc_list(alloc_list* list) {
	for(alloc* aux = list->head;aux;aux = aux->next)
		printf("[Memory leak at %s]\n",aux->name);
}

void list_leaks() {
	for(int i = 0; i < SIZE; i++) {
		if(hashArray[i]) {
			print_alloc_list(hashArray[i]);
		}
	}
}

void *mycalloc(char *var_name, int n, int struct_size)
{
	void *ptr;

	ptr = (void *) calloc(n,struct_size);
	
	if (ptr == NULL){
		printf("Memory allocation error in %s!\n", var_name);
		exit(1);
	}

	#ifdef check_memory_leak
		insert_hash_table(ptr,var_name);
	#endif
	
	return ptr;	
}

void myfree(void *ptr) {
	#ifdef check_memory_leak
		remove_hash_table(ptr);
	#endif
	free(ptr);
}
