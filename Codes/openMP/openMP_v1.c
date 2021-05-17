#ifdef _OPENMP
#include <omp.h>
#endif

#include <cover_functions.h>


unsigned short int nb_thread;


void solve_OMP(const struct instance_t * instance, struct context_t ** ctxs)
{
	int chosen_item;
	struct sparse_array_t * active_options;

	chosen_item = choose_next_item(ctxs[0]);
	active_options = ctxs[0]->active_options[chosen_item];

	#pragma omp parallel for schedule(static, 1)
	for (unsigned short int i = 0; i < nb_thread; ++i)
	{
		cover(instance, ctxs[i], chosen_item);
		ctxs[i]->num_children[0] = active_options->len;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < active_options->len; k++)
	{
		unsigned int th_num = omp_get_thread_num();
		int option = active_options->p[k];
		ctxs[th_num]->child_num[ctxs[th_num]->level] = k;
		choose_option(instance, ctxs[th_num], option, chosen_item);
		solve(instance, ctxs[th_num]);

		if (ctxs[th_num]->solutions >= max_solutions)
			exit(0);
		unchoose_option(instance, ctxs[th_num], option, chosen_item);
	}

	#pragma omp parallel for schedule(static, 1)
	for (unsigned short int i = 0; i < nb_thread; ++i)
		uncover(instance, ctxs[i], chosen_item);		/* backtrack */
}


int main(int argc, char **argv)
{
	struct option longopts[5] = {
		{"in", required_argument, NULL, 'i'},
		{"progress-report", required_argument, NULL, 'v'},
		{"print-solutions", no_argument, NULL, 'p'},
		{"stop-after", required_argument, NULL, 's'},
		{NULL, 0, NULL, 0}
	};

	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'i':
			in_filename = optarg;
			break;
		case 'p':
			print_solutions = true;
			break;
		case 's':
			max_solutions = atoll(optarg);
			break;
		case 'v':
			report_delta = atoll(optarg);
			break;          
		default:
			errx(1, "Unknown option\n");
		}
	}
	if (in_filename == NULL)
		usage(argv);
	next_report = report_delta;

	nb_thread = omp_get_max_threads();

	double stop;
	unsigned long long nb_sol = 0;
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

	#pragma omp parallel for reduction(+:nb_sol)
	for (unsigned short int i = 0; i < nb_thread; ++i)
		nb_sol += ctxs[i]->solutions;

	printf("FINI. Trouvé %lld solutions en %.1fs\n", nb_sol, stop);
	exit(EXIT_SUCCESS);
}
