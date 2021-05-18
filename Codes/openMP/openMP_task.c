#ifdef _OPENMP
#include <omp.h>
#endif

#include <cover_functions.h>


unsigned short int nb_thread;
unsigned long long nb_sol = 0;


void solve_OMP(const struct instance_t * instance, struct context_t * ctx)
{
	ctx->nodes++;
	if (ctx->nodes == next_report)
		progress_report(ctx);
	if (sparse_array_empty(ctx->active_items)) {
		solution_found(instance, ctx);
		return;		/* succès : plus d'objet actif */
	}
	int chosen_item = choose_next_item(ctx);
	struct sparse_array_t *active_options = ctx->active_options[chosen_item];
	if (sparse_array_empty(active_options))
		return;		/* échec : impossible de couvrir chosen_item */
	cover(instance, ctx, chosen_item);
	ctx->num_children[ctx->level] = active_options->len;

	if (!(ctx->level % 8))
		for (int k = 0; k < active_options->len; k++) {
			int option = active_options->p[k];
			ctx->child_num[ctx->level] = k;
			choose_option(instance, ctx, option, chosen_item);

			#pragma omp task
			solve_OMP(instance, ctx);

			printf("Created task for k = %d\n", k);

			if (ctx->solutions >= max_solutions)
				exit(0);
			unchoose_option(instance, ctx, option, chosen_item);
		}
	else
		for (int k = 0; k < active_options->len; k++) {
			int option = active_options->p[k];
			ctx->child_num[ctx->level] = k;
			choose_option(instance, ctx, option, chosen_item);
			solve_OMP(instance, ctx);

			if (ctx->solutions >= max_solutions)
				exit(0);
			unchoose_option(instance, ctx, option, chosen_item);
		}
	uncover(instance, ctx, chosen_item);	/* backtrack */
}


void init_solve_OMP(const struct instance_t * instance, struct context_t * ctx)
{
	ctx->nodes++;
	if (ctx->nodes == next_report)
		progress_report(ctx);
	if (sparse_array_empty(ctx->active_items)) {
		solution_found(instance, ctx);
		return;		/* succès : plus d'objet actif */
	}
	int chosen_item = choose_next_item(ctx);
	struct sparse_array_t *active_options = ctx->active_options[chosen_item];
	if (sparse_array_empty(active_options))
		return;		/* échec : impossible de couvrir chosen_item */
	cover(instance, ctx, chosen_item);
	ctx->num_children[ctx->level] = active_options->len;

	#pragma omp parallel
	#pragma omp single
	for (int k = 0; k < active_options->len; k++) {
		int option = active_options->p[k];
		ctx->child_num[ctx->level] = k;
		choose_option(instance, ctx, option, chosen_item);

		#pragma omp task
		solve_OMP(instance, ctx);

		printf("Created task for k = %d\n", k);

		if (ctx->solutions >= max_solutions)
			exit(0);
		unchoose_option(instance, ctx, option, chosen_item);
	}

	uncover(instance, ctx, chosen_item);	/* backtrack */
}


int main(int argc, char **argv)
{
	option_setup(argc, argv);

	struct instance_t * instance = load_matrix(in_filename);
	struct context_t * ctx = backtracking_setup(instance);
	start = wtime();
	init_solve_OMP(instance, ctx);
	printf("FINI. Trouvé %lld solutions en %.2fs\n", ctx->solutions, wtime() - start);
	exit(EXIT_SUCCESS);
}
