#ifdef _OPENMP
#include <omp.h>
#endif

#include <cover_functions.h>


// void solve_OMP(const struct instance_t * instance)
// {
// 	struct context_t * ctx = backtracking_setup(instance);
// 	// ctx->nodes++;
// 	// if (ctx->nodes == next_report)
// 	// 	progress_report(ctx);
// 	// if (sparse_array_empty(ctx->active_items)) {
// 	// 	solution_found(instance, ctx);
// 	// 	return;			/* succès : plus d'objet actif */
// 	// }
// 	int chosen_item = choose_next_item(ctx);
// 	struct sparse_array_t *active_options = ctx->active_options[chosen_item];
// 	// if (sparse_array_empty(active_options))
// 	// 	return;			/* échec : impossible de couvrir chosen_item */
// 	cover(instance, ctx, chosen_item);
// 	ctx->num_children[ctx->level] = active_options->len;

// 	#pragma omp parallel for schedule(dynamic, 100)
// 	for (int k = 0; k < active_options->len; k++) {

// 		printf("Thread %d got k = %d\n", omp_get_thread_num(), k);

// 		int option = active_options->p[k];

// 		printf("Thread %d DEBUG 1\n", omp_get_thread_num());

// 		ctx->child_num[ctx->level] = k;

// 		printf("Thread %d DEBUG 2\n", omp_get_thread_num());

// 		choose_option(instance, ctx, option, chosen_item);

// 		printf("Thread %d DEBUG 3\n", omp_get_thread_num());

// 		solve(instance, ctx);

// 		printf("Thread %d DEBUG 4\n", omp_get_thread_num());

// 		if (ctx->solutions >= max_solutions)
// 			exit(0);

// 		printf("Thread %d DEBUG 5\n", omp_get_thread_num());

// 		unchoose_option(instance, ctx, option, chosen_item);
// 	}
// 	uncover(instance, ctx, chosen_item);		/* backtrack */
// }


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
	printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions, 
			wtime() - start);
	exit(EXIT_SUCCESS);
}
