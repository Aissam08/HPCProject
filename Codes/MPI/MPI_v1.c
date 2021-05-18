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

/*	MPI_Type_struct(3,B,D,type1,newtype)
	--> { (2,1,2), (0,16,26), (MPI_FLOAT, {(double,0), (char,8) } , MPI_CHAR), name_of_type }
	 - Il y aura 3 composantes dans la structure
	   - longueur des données (2 / 1 / 3)
	   - La structure contient 2 FLOAT, 1 map et 3 char
	   - le map est composé d'un double et d'un char 
	   - longueur d'écarts (0,16,26)
	   - le MPI_FLOAT est à 0 / 
	   - le type1 (map) commence à 16 (après 2 float de 8 octet)
	   - le MPI_CHAR commence à 26 octet (16 + 9)
	*/


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
	if (my_rank == 0)
	
	{
		struct instance_t * instance = load_matrix(in_filename);

		int to_send[] = {instance->n_items,instance->n_primary,instance->n_options};

		MPI_Send(to_send,3,MPI_INT, 1, 
			10,MPI_COMM_WORLD);

		MPI_Send(instance->options,instance->n_options ,MPI_INT, 1,
			10,MPI_COMM_WORLD);

		MPI_Send(instance->ptr,instance->n_items ,MPI_INT, 1,
			9,MPI_COMM_WORLD);
		/*
		printf("Nombre d'items : %d\n",instance->n_items);
		printf("Nombre primary : %d\n",instance->n_primary);
		printf("Nombre d'options : %d\n",instance->n_options);
		printf("Taille tableau options : %d\n",instance->n_items * instance->n_options);
		struct context_t * ctx = backtracking_setup(instance);
		start = wtime();
		solve_MPI(instance, ctx, status, nb_proc, my_rank);
		printf("FINI. Trouvé %lld solutions en %.1fs\n\n", ctx->solutions, 
				wtime() - start);*/

	}
	if (my_rank == 1)
	{
		int *to_recv = malloc(2*sizeof(int));

		MPI_Recv(to_recv,3,MPI_INT,0,
			10,MPI_COMM_WORLD,&status);

		printf("Reçu : %d / %d / %d\n",to_recv[0], to_recv[1], to_recv[2]);
		int *recv_opt = malloc(to_recv[2]*sizeof(int));
		int *recv_items = malloc(to_recv[1]*sizeof(int));

		MPI_Recv(recv_opt,to_recv[2],MPI_INT, 0,
			10,MPI_COMM_WORLD, &status);

		MPI_Recv(recv_items,to_recv[1],MPI_INT, 0,
			9,MPI_COMM_WORLD, &status);

	}
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
