#ifdef _OPENMP
#include <omp.h>
#endif

#include <cover_functions.h>


void solve_OMP(const struct instance_t * instance)
{
	unsigned short int nb_thread = omp_get_max_threads();
	int chosen_item;
	struct sparse_array_t * active_options;
	struct context_t ** ctxs = (struct context_t **) malloc(nb_thread * sizeof(* ctxs));
	
	// #pragma omp parallel for schedule(static, 1)
	for (unsigned short int i = 0; i < nb_thread; ++i)
	{
		ctxs[i] = backtracking_setup(instance);
		ctxs[i]->nodes++;
		if (!i) {
			chosen_item = choose_next_item(ctxs[0]);
			active_options = ctxs[0]->active_options[chosen_item];
		}
		cover(instance, ctxs[i], chosen_item);
		ctxs[i]->num_children[0] = active_options->len;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < active_options->len; k++) {

		// printf("Thread %d got k = %d\n", omp_get_thread_num(), k);

		int option = active_options->p[k];

		// printf("Thread %d DEBUG 1\n", omp_get_thread_num());

		ctxs[omp_get_thread_num()]->child_num[ctxs[omp_get_thread_num()]->level] = k;

		// printf("Thread %d DEBUG 2\n", omp_get_thread_num());

		choose_option(instance, ctxs[omp_get_thread_num()], option, chosen_item);

		// printf("Thread %d DEBUG 3\n", omp_get_thread_num());

		solve(instance, ctxs[omp_get_thread_num()]);

		// printf("Thread %d DEBUG 4\n", omp_get_thread_num());

		if (ctxs[omp_get_thread_num()]->solutions >= max_solutions)
			exit(0);

		// printf("Thread %d DEBUG 5\n", omp_get_thread_num());

		unchoose_option(instance, ctxs[omp_get_thread_num()], option, chosen_item);
	}
	uncover(instance, ctxs[omp_get_thread_num()], chosen_item);		/* backtrack */
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

	struct instance_t * instance = load_matrix(in_filename);
	// struct context_t * ctx = backtracking_setup(instance);
	start = wtime();
	solve_OMP(instance);
	printf("FINI. Trouvé %lld solutions en %.1fs\n", 0LL,//ctxs[0]->solutions, 
			wtime() - start);
	exit(EXIT_SUCCESS);
}
