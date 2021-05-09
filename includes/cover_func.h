#ifndef COVER_FUNC_H
#define COVER_FUNC_H

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <err.h>
#include <getopt.h>
#include <sys/time.h>


double start;

char * in_filename;			// nom du fichier contenant la matrice
bool print_solutions;		// affiche chaque solution
long long report_delta;		// affiche un rapport tous les ... noeuds
long long next_report;		// prochain rapport affiché au noeud...
long long max_solutions;	// stop après ... solutions

static const char DIGITS[62];


struct instance_t {
	int n_items;
	int n_primary;
	int n_options;
	char ** item_name;	// potentiellement NULL, sinon de taille n_items
	int * options;		// l'option i contient les objets options[ptr[i]:ptr[i+1]]
	int * ptr;			// taille n_options + 1
};

struct sparse_array_t {
	int len;			// nombre d'éléments stockés
	int capacity;		// taille maximale
	int * p;			// contenu de l'ensemble = p[0:len] 
	int * q;			// taille capacity (tout comme p)
};

struct context_t {
	struct sparse_array_t * active_items;		// objets actifs
	struct sparse_array_t ** active_options;	// options actives contenant l'objet i
	int * chosen_options;						// options choisies à ce stade
	int * child_num;							// numéro du fils exploré
	int * num_children;							// nombre de fils à explorer
	int level;									// nombre d'options choisies
	long long nodes;							// nombre de noeuds explorés
	long long solutions;						// nombre de solutions trouvées 
};


double wtime();
void usage(char **argv);
bool item_is_primary(const struct instance_t * instance, int item);
void print_option(const struct instance_t * instance, int option);

bool sparse_array_membership(const struct sparse_array_t * S, int x);
bool sparse_array_empty(const struct sparse_array_t * S);
void sparse_array_add(struct sparse_array_t * S, int x);
void sparse_array_remove(struct sparse_array_t * S, int x);
void sparse_array_unremove(struct sparse_array_t * S);
void sparse_array_unadd(struct sparse_array_t * S);

bool item_is_active(const struct context_t * ctx, int item);

void solution_found(const struct instance_t * instance, struct context_t * ctx);

void deactivate(const struct instance_t * instance, struct context_t * ctx, int option, int covered_item);
void reactivate(const struct instance_t * instance, struct context_t * ctx, int option, int uncovered_item);

void cover(const struct instance_t * instance, struct context_t * ctx, int item);
void uncover(const struct instance_t * instance, struct context_t * ctx, int item);

void progress_report(const struct context_t *ctx);

struct instance_t * load_matrix(const char * filename);

#endif
