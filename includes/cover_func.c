#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <err.h>
#include <getopt.h>
#include <sys/time.h>

#include "cover_func.h"


double start = 0.0;

char * in_filename = NULL;						// nom du fichier contenant la matrice
bool print_solutions = false;					// affiche chaque solution
long long report_delta = 1e6;					// affiche un rapport tous les ... noeuds
long long next_report;							// prochain rapport affiché au noeud...
long long max_solutions = 0x7fffffffffffffff;	// stop après ... solutions


static const char DIGITS[62] = {
	'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
	'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 
	'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
	'u', 'v', 'w', 'x', 'y', 'z',
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
	'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
	'U', 'V', 'W', 'X', 'Y', 'Z'
};


double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

void usage(char **argv)
{
	printf("%s --in FILENAME [OPTIONS]\n\n", argv[0]);
	printf("Options:\n");
	printf(" --progress-report N   display a message every N nodes (0 to disable)\n");
	printf(" --print-solutions     display solutions when they are found\n");
	printf(" --stop-after N        stop the search once N solutions are found\n");
	exit(0);
}

bool item_is_primary(const struct instance_t * instance, int item)
{
	return item < instance->n_primary;
}

void print_option(const struct instance_t * instance, int option)
{
	if (instance->item_name == NULL)
		errx(1, "tentative d'affichage sans noms d'objet");
	for (int p = instance->ptr[option]; p < instance->ptr[option + 1]; p++) {
		int item = instance->options[p];
		printf("%s ", instance->item_name[item]);
	}
	printf("\n");
}


bool sparse_array_membership(const struct sparse_array_t * S, int x)
{
	return (S->q[x] < S->len);
}

bool sparse_array_empty(const struct sparse_array_t * S)
{
	return (S->len == 0);
}

void sparse_array_add(struct sparse_array_t * S, int x)
{
	int i = S->len;
	S->p[i] = x;
	S->q[x] = i;
	S->len = i + 1;
}

void sparse_array_remove(struct sparse_array_t * S, int x)
{
	int j = S->q[x];
	int n = S->len - 1;
	// échange p[j] et p[n] 
	int y = S->p[n];
	S->p[n] = x;
	S->p[j] = y;
	// met q à jour
	S->q[x] = n;
	S->q[y] = j;
	S->len = n;
}

void sparse_array_unremove(struct sparse_array_t * S)
{
	S->len++;
}

void sparse_array_unadd(struct sparse_array_t * S)
{
	S->len--;
}


bool item_is_active(const struct context_t * ctx, int item)
{
	return sparse_array_membership(ctx->active_items, item);
}

void solution_found(const struct instance_t * instance, struct context_t * ctx)
{
	ctx->solutions++;
	if (!print_solutions)
		return;
	printf("Trouvé une nouvelle solution au niveau %d après %lld noeuds\n", ctx->level, ctx->nodes);
	printf("Options : \n");
	for (int i = 0; i < ctx->level; i++) {
		int option = ctx->chosen_options[i];
		printf("+ %d : ", option);
		print_option(instance, option);
	}
	printf("\n");
	printf("----------------------------------------------------\n");
}


void deactivate(const struct instance_t * instance, struct context_t * ctx, int option, int covered_item)
{
	for (int k = instance->ptr[option]; k < instance->ptr[option+1]; ++k) {
		int item = instance->options[k];
		if (item == covered_item)
			continue;
		sparse_array_remove(ctx->active_options[item], option);
	}
}

void reactivate(const struct instance_t * instance, struct context_t * ctx, int option, int uncovered_item)
{
	for (int k = instance->ptr[option + 1] - 1; k >= instance->ptr[option]; --k) {
		int item = instance->options[k];
		if (item == uncovered_item)
			continue;
		sparse_array_unremove(ctx->active_options[item]);
	}
}


void cover(const struct instance_t * instance, struct context_t * ctx, int item)
{
	if (item_is_primary(instance, item))
		sparse_array_remove(ctx->active_items, item);
	struct sparse_array_t * active_options = ctx->active_options[item];
	for (int i = 0; i < active_options->len; ++i) {
		int option = active_options->p[i];
		deactivate(instance, ctx, option, item);
	}
}

void uncover(const struct instance_t * instance, struct context_t * ctx, int item)
{
	struct sparse_array_t * active_options = ctx->active_options[item];
	for (int i = active_options->len - 1; i >= 0; --i) {
		int option = active_options->p[i];
		reactivate(instance, ctx, option, item);
	}
	if (item_is_primary(instance, item))
		sparse_array_unremove(ctx->active_items);
}


void progress_report(const struct context_t *ctx)
{
	double now = wtime();
	printf("Exploré %lld noeuds, trouvé %lld solutions, temps écoulé %.1fs. ", 
			ctx->nodes, ctx->solutions, now - start);
	int i = 0;
	for (int k = 0; k < ctx->level; ++k) {
		if (i > 44)
			break;
		int n = ctx->child_num[k];
		int m = ctx->num_children[k];
		if (m == 1)
			continue;
		printf("%c%c ", (n < 62) ? DIGITS[n] : '*', (m < 62) ? DIGITS[m] : '*');
		++i;
	}
	printf("\n"),
	next_report += report_delta;
}


struct instance_t * load_matrix(const char * filename)
{
	struct instance_t * instance = malloc(sizeof(*instance));
	if (instance == NULL)
		err(1, "Impossible d'allouer l'instance");
	FILE * in = fopen(filename, "r");
	if (in == NULL)
		err(1, "Impossible d'ouvrir %s en lecture", filename);
	int n_it, n_op;
	if (fscanf(in, "%d %d\n", &n_it, &n_op) != 2)
		errx(1, "Erreur de lecture de la taille du problème\n");
	if (n_it == 0 || n_op == 0)
		errx(1, "Impossible d'avoir 0 objets ou 0 options");
	instance->n_items = n_it;
	instance->n_options = n_op;
	instance->item_name = malloc(n_it * sizeof(char *));
	instance->ptr = malloc((n_op + 1) * sizeof(int));
	instance->options = malloc(n_it * n_op * sizeof(int));         // surallocation massive
	if (instance->item_name == NULL || instance->ptr == NULL || instance->options == NULL)
		err(1, "Impossible d'allouer la mémoire pour stocker la matrice");


	enum state_t {START, ID, WHITESPACE, BAR, ENDLINE, ENDFILE};
	enum state_t state = START;

	char buffer[256];
	int i = 0;     // prochain octet disponible du buffer
	int n = 0;     // dernier octet disponible du buffer

	char id[65];
	id[64] = 0;    // sentinelle à la fin, quoi qu'il arrive
	int j = 0;     // longueur de l'identifiant en cours de lecture

	int current_item = 0;
	while (state != ENDLINE) {
		enum state_t prev_state = state;
		if (i >= n) {
			n = fread(buffer, 1, 256, in);
			if (n == 0) {
				if (feof(in))
					state = ENDFILE;
				if (ferror(in))
					err(1, "erreur lors de la lecture de %s", in_filename);
			}
			i = 0;
		}

		if (state == ENDFILE)
			NULL; 	// don't examine buffer[i]
		else if (buffer[i] == '\n')
			state = ENDLINE;
		else if (buffer[i] == '|')
			state = BAR;
		else if (isspace(buffer[i]))
			state = WHITESPACE;
		else
			state = ID;

		// traite le caractère lu
		if (state == ID) {
			if (j == 64)
				errx(1, "nom d'objet trop long : %s", id);
			id[j++] = buffer[i];
		}
		if (prev_state == ID && state != ID) {
			id[j] = '\0';
			if (current_item == instance->n_items)
				errx(1, "Objet excedentaire : %s", id);
			for (int k = 0; k < current_item; ++k)
				if (strcmp(id, instance->item_name[k]) == 0)
					errx(1, "Nom d'objets dupliqué : %s", id);
			instance->item_name[current_item] = malloc(j+1);
			strcpy(instance->item_name[current_item++], id);
			j = 0;
		}

		if (state == BAR)
			instance->n_primary = current_item;
		if (state == ENDFILE)
			errx(1, "Fin de fichier prématurée");
		++i;	// passe au prochain caractère
	}

	if (current_item != instance->n_items)
		errx(1, "Incohérence : %d objets attendus mais seulement %d fournis\n", instance->n_items, current_item);
	if (instance->n_primary == 0)
		instance->n_primary = instance->n_items;

	int current_option = 0;
	int p = 0;       // pointeur courant dans instance->options
	instance->ptr[0] = p;
	bool has_primary = false;
	while (state != ENDFILE) {
		enum state_t prev_state = state;
		if (i >= n) {
			n = fread(buffer, 1, 256, in);
			if (n == 0) {
				if (feof(in))
					state = ENDFILE;
				if (ferror(in))
					err(1, "erreur lors de la lecture de %s", in_filename);
			}
			i = 0;
		}

		if (state == ENDFILE)
			NULL; 	// don't examine buffer[i]
		else if (buffer[i] == '\n')
			state = ENDLINE;
		else if (buffer[i] == '|')
			state = BAR;
		else if (isspace(buffer[i]))
			state = WHITESPACE;
		else
			state = ID;

		// traite le caractère lu
		if (state == ID) {
			if (j == 64)
				errx(1, "nom d'objet trop long : %s", id);
			id[j++] = buffer[i];
		}
		if (prev_state == ID && state != ID) {
			id[j] = '\0';
			// identifie le numéro de l'objet en question
			int item_number = -1;
			for (int k = 0; k < instance->n_items; ++k)
				if (strcmp(id, instance->item_name[k]) == 0) {
					item_number = k;
					break;
				}
			if (item_number == -1)
				errx(1, "Objet %s inconnu dans l'option #%d", id, current_option);
			// détecte les objets répétés
			for (int k = instance->ptr[current_option]; k < p; ++k)
				if (item_number == instance->options[k])
					errx(1, "Objet %s répété dans l'option %d\n", 
							instance->item_name[item_number], current_option);
			instance->options[p++] = item_number;
			has_primary |= item_is_primary(instance, item_number);
			j = 0;
		}

		if (state == BAR) {
			errx(1, "Trouvé | dans une option.");
		}
		if ((state == ENDLINE || state == ENDFILE)) {
			// esquive les lignes vides
			if (p > instance->ptr[current_option]) {
				if (current_option == instance->n_options)
					errx(1, "Option excédentaire");
				if (!has_primary)
					errx(1, "Option %d sans objet primaire\n", current_option);
				instance->ptr[++current_option] = p;
				has_primary = false;
			}
		}
		// passe au prochain caractère
		++i;
	}
	if (current_option != instance->n_options)
		errx(1, "Incohérence : %d options attendues mais seulement %d fournies\n", instance->n_options, current_option);
	fclose(in);
	fprintf(stderr, "Lu %d objets (%d principaux) et %d options\n", instance->n_items, instance->n_primary, instance->n_options);
	return instance;
}
