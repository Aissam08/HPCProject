#include <mpi.h>

#include <cover_functions.h>


void solve_MPI(const struct instance_t *instance, struct context_t *ctx, MPI_Status status, int nb_proc, int my_rank)
{       
	ctx->nodes++;
	if (ctx->nodes == next_report)
		progress_report(ctx);
	if (sparse_array_empty(ctx->active_items)) {
		solution_found(instance, ctx);
		return;			/* succès : plus d'objet actif */
	}
	int chosen_item = choose_next_item(ctx);
	struct sparse_array_t *active_options = ctx->active_options[chosen_item];
	if (sparse_array_empty(active_options))
		return;			/* échec : impossible de couvrir chosen_item */
	cover(instance, ctx, chosen_item);
	ctx->num_children[ctx->level] = active_options->len;
	for (int k = 0; k < active_options->len; k++) {
		int option = active_options->p[k];
		ctx->child_num[ctx->level] = k;
		choose_option(instance, ctx, option, chosen_item);
		solve_MPI(instance, ctx, status, nb_proc, my_rank);
		if (ctx->solutions >= max_solutions)
			return;
		unchoose_option(instance, ctx, option, chosen_item);
	}
	uncover(instance, ctx, chosen_item);	/* backtrack */
}

int main(int argc, char **argv)
{
	MPI_Status status;
	int nb_proc, my_rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	printf("First step done\n");

	struct option longopts[5] = {
		{"in", required_argument, NULL, 'i'},
		{"progress-report", required_argument, NULL, 'v'},
		{"print-solutions", no_argument, NULL, 'p'},
		{"stop-after", required_argument, NULL, 's'},
		{NULL, 0, NULL, 0}
	};

	char ch;
	printf("Second step done\n");
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

	printf("NB itemq%d\n",instance->n_items);
	struct context_t * ctx = backtracking_setup(instance);

	printf("Third step done\n");
	start = wtime();
	solve_MPI(instance, ctx, status, nb_proc, my_rank);
	printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions, 
			wtime() - start);
	printf("Finish step done\n");
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
