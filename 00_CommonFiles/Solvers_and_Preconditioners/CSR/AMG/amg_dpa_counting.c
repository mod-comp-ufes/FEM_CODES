#include <stdlib.h>
#include "amg_dpa_counting.h"
# include "../../../Allocation_Operations/allocations.h"

/*----------------------------------------------------------------*/

counting * create_c(int n);

void destroy_c(counting * c);

Nodec * insert_c(counting *c, int i, int m) ;

void remove_c(counting *c, int i);

Nodec* move_c(counting *c, int i, int new_m);

/*----------------------------------------------------------------*/

// list of pointers for counting array elements
listc * create_lc() {
    listc *l = (listc *) mycalloc("l of create_lc",1,sizeof (listc));
    l->n = 0;
    l->head = l->tail = NULL;
    return l;
};

void destroy_lc(listc *l) {
    Nodec *h;
    while (l->head) {
        h = l->head->next;
        myfree(l->head);
        l->head = h;
    }
    myfree(l);
};

/*----------------------------------------------------------------*/

counting * create_c(int n) {
    counting * c = (counting*) mycalloc("c of create_c",1,sizeof (counting));
    c->count = (listc**) mycalloc("c->count of create_c",2 * n , sizeof (listc*));
    for (int i = 0; i < 2 * n; i++) c->count[i] = create_lc();
    c->nodes = (Nodec**) mycalloc("c->nodes of create_c",n , sizeof (Nodec*));
    for (int i = 0; i < n; i++) c->nodes[i] = NULL;
    c->m = (int*) mycalloc("m of create_c",n, sizeof (int));
    c->n = n;
    c->min_m = 2 * n;
    return c;
}

void destroy_c(counting * c) {
    for (int i = 0; i < c->n*2; i++) destroy_lc(c->count[i]);
    myfree(c->count);
    myfree(c->nodes);
    myfree(c->m);
    myfree(c);
}

Nodec * insert_c(counting *c, int i, int m) {
    listc *l = c->count[m];
    Nodec *t = (Nodec *) mycalloc("t of insert_c",1,sizeof (Nodec));
    c->nodes[i] = t;
    t->i = i;
    t->m = m;
    t->next = NULL;
    t->prev = l->tail;
    if (l->tail) {
        l->tail->next = t;
        l->tail = t;
    } else {
        l->head = l->tail = t;
    }
    l->n++;
    c->m[i] = m;
    if (c->m[i] < c->min_m)
        c->min_m = c->m[i];

    return t;
};

void remove_c(counting *c, int i) {
    Nodec *n = c->nodes[i];
    int m = n->m;

    if (n->prev) {
        n->prev->next = n->next;
    } else {
        c->count[m]->head = n->next;
    }

    if (n->next) {
        n->next->prev = n->prev;
    } else {
        c->count[m]->tail = n->prev;
    }

    myfree(n);
    c->nodes[i] = NULL;
    c->count[m]->n--;

    if (c->count[c->min_m]->n == 0) {
        for (c->min_m = m + 1; ((c->min_m < 2*(c->n)) && (c->count[c->min_m]->n == 0)); c->min_m++) {
        }
    }
}

Nodec* move_c(counting *c, int i, int new_m) {
    /* declarations */
    listc *l = c->count[new_m]; // new list
    Nodec *t = c->nodes[i]; // moving node
    int m = t->m; // old t

    // update m[] array
    c->m[i] = new_m;

    // remove t from old array
    if (t->prev) {
        t->prev->next = t->next;
    } else {
        c->count[m]->head = t->next;
    }
    if (t->next) {
        t->next->prev = t->prev;
    } else {
        c->count[m]->tail = t->prev;
    }

    // add t to new array
    t->m = new_m;
    t->next = NULL;
    t->prev = l->tail;
    if (l->tail) {
        l->tail->next = t;
        l->tail = t;
    } else {
        l->head = l->tail = t;
    }
    l->n++;

    // update nodes[] array
    c->nodes[i] = t;
    c->count[m]->n--;

    // update min_m
    if (new_m < c->min_m)
        c->min_m = new_m;
    else {
        if (c->count[c->min_m]->n == 0) {
            for (c->min_m = m + 1; ((c->min_m < 2*(c->n)) && (c->count[c->min_m]->n == 0)); c->min_m++) {
            }
        }
    }

    return t;
}
