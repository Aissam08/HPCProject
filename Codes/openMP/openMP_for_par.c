#ifdef _OPENMP
#include <omp.h>
#endif

#include <cover_functions.h>


unsigned short int nb_thread;
unsigned long long nb_sol = 0;


void solve_OMP(const struct instance_t * instance, struct context_t ** ctxs)
{
	int i, chosen_item;
	struct sparse_array_t * active_options;

	chosen_item = choose_next_item(ctxs[0]);
	active_options = ctxs[0]->active_options[chosen_item];

	#pragma omp parallel for schedule(static, 1)
	for (i = 0; i < nb_thread; ++i)
	{
		cover(instance, ctxs[i], chosen_item);
		ctxs[i]->num_children[0] = active_options->len;
	}

	#pragma omp parallel for schedule(dynamic)
	for (i = 0; i < active_options->len; ++i)
	{
		unsigned int th_num = omp_get_thread_num();
		int option = active_options->p[i];
		ctxs[th_num]->child_num[ctxs[th_num]->level] = i;
		choose_option(instance, ctxs[th_num], option, chosen_item);
		solve(instance, ctxs[th_num]);

		if (ctxs[th_num]->solutions >= max_solutions)
			exit(0);
		unchoose_option(instance, ctxs[th_num], option, chosen_item);
	}

	#pragma omp parallel for schedule(static, 1) reduction(+:nb_sol)
	for (i = 0; i < nb_thread; ++i) {
		uncover(instance, ctxs[i], chosen_item);		/* backtrack */
		nb_sol += ctxs[i]->solutions;
	}
}


int main(int argc, char **argv)
{
	option_setup(argc, argv);

	nb_thread = omp_get_max_threads();

	double stop;
	struct instance_t * instance = load_matrix(in_filename);
	struct context_t ** ctxs = (struct context_t **) malloc(nb_thread * sizeof(* ctxs));

	#pragma omp parallel for schedule(static, 1)
	for (unsigned short int i = 0; i < nb_thread; ++i)
	{
		ctxs[i] = backtracking_setup(instance);
		ctxs[i]->nodes = 1;
	}

	start = wtime();
	solve_OMP(instance, ctxs);
	stop = wtime() - start;

	printf("FINI. Trouvé %lld solutions en %.2fs\n", nb_sol, stop);
	exit(EXIT_SUCCESS);
}
