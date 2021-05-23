#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>
#include <stddef.h>

#include <cover_functions.h>


int nb_MPI_proc, my_MPI_rank, nb_OMP_thread, my_OMP_thread;


void receive_pb_data(struct instance_t ** instance, struct context_t *** ctx, int * common_item)
{
	*instance = (struct instance_t *) malloc(sizeof(struct instance_t));
	*ctxs = (struct context_t **) malloc(nb_thread * sizeof(struct context_t *));
	
	int buffer[4];
	MPI_Bcast(buffer, 4, MPI_INT, 0, MPI_COMM_WORLD);
	(*instance)->n_items = buffer[0];
	(*instance)->n_primary = buffer[1];
	(*instance)->n_options = buffer[2];
	*common_item = buffer[3];

	(*instance)->item_name = NULL;
	(*instance)->options = (int *) malloc((*instance)->n_options * (*instance)->n_items * sizeof(int));
	MPI_Bcast((*instance)->options, (*instance)->n_options * (*instance)->n_items, MPI_INT, 0, MPI_COMM_WORLD);

	(*instance)->ptr = (int *) malloc(((*instance)->n_options + 1) * sizeof(int));
	MPI_Bcast((*instance)->ptr, (*instance)->n_options + 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("On process %d : received %d items & %d options.\n", my_MPI_rank, (*instance)->n_items, (*instance)->n_options);
	#pragma omp parallel for schedule(static, 1)
	for (unsigned short int i = 0; i < nb_OMP_thread; ++i)
	{
		*ctxs[i] = backtracking_setup(*instance);
		*ctxs[i]->nodes = 1;
	}
	printf("On process %d : %d contexts created.\n", my_MPI_rank, nb_OMP_thread);
}

void send_pb_data(struct instance_t ** instance, struct context_t ** ctx, int * common_item)
{
	*instance = load_matrix(in_filename);
	*ctx = backtracking_setup(*instance);
	*common_item = choose_next_item(*ctx);

	printf("From Master : got %d items, %d primary & %d options to broadcast.\n",
		(*instance)->n_items, (*instance)->n_primary, (*instance)->n_options);

	int buffer[4] = {(*instance)->n_items, (*instance)->n_primary, (*instance)->n_options, *common_item};
	MPI_Bcast(buffer, 4, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast((*instance)->options, (*instance)->n_options * (*instance)->n_items, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast((*instance)->ptr, (*instance)->n_options + 1, MPI_INT, 0, MPI_COMM_WORLD);
}


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
		my_thread = omp_get_thread_num();
		int option = active_options->p[i];
		ctxs[my_thread]->child_num[ctxs[my_thread]->level] = i;
		choose_option(instance, ctxs[my_thread], option, chosen_item);
		solve(instance, ctxs[my_thread]);

		if (ctxs[my_thread]->solutions >= max_solutions)
			exit(0);
		unchoose_option(instance, ctxs[my_thread], option, chosen_item);
	}

	#pragma omp parallel for schedule(static, 1) reduction(+:nb_sol)
	for (i = 0; i < nb_thread; ++i) {
		uncover(instance, ctxs[i], chosen_item);		/* backtrack */
		nb_sol += ctxs[i]->solutions;
	}
}

void solve_MPI(const struct instance_t * instance, struct context_t * ctx, int chosen_item, int option)
{
	if (ctx->nodes == next_report)
		progress_report(ctx);
	if (sparse_array_empty(ctx->active_items)) {
		solution_found(instance, ctx);
		return;		/* succès : plus d'objet actif */
	}
	choose_option(instance, ctx, option, chosen_item);
	solve(instance, ctx);
	if (ctx->solutions >= max_solutions)
		return;
	unchoose_option(instance, ctx, option, chosen_item);
}

void comm_solve_MPI(const struct instance_t * instance, struct context_t * ctx, int chosen_item, bool debug)
{
	MPI_Status status;
	bool end = false;
	int i;
	long long int com_buffer[3] = {0LL, 0LL, 0LL};		// node_id / nb_nodes / nb_sol
	long long int nb_nodes = 0LL, nb_sol = 0LL;

	if (!my_MPI_rank) {
		for (i = 0; i < nb_MPI_proc - 1; ++i) {
			com_buffer[0] = i;
			MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, i + 1, 1, MPI_COMM_WORLD);
		}
		while (i < ctx->active_options[chosen_item]->len && !end) {
			MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			switch (status.MPI_TAG)
			{
				case 2 :
					if (debug)
						printf("Received branch %llu with %llu new solutions from process %d.\n",
							com_buffer[0], com_buffer[2], status.MPI_SOURCE);
					nb_nodes += com_buffer[1];
					nb_sol += com_buffer[2];
					com_buffer[0] = i++;
					com_buffer[1] = nb_nodes;
					com_buffer[2] = nb_sol;
					MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
					break;

				case 99 :
					nb_nodes += com_buffer[1];
					nb_sol += com_buffer[2];
					end = true;
					com_buffer[1] = nb_nodes;
					com_buffer[2] = nb_sol;
					MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 99, MPI_COMM_WORLD);
					break;

				default:
					com_buffer[1] = nb_nodes;
					com_buffer[2] = nb_sol;
					printf("Received unknown type communication. Resetting process %d.\n", status.MPI_SOURCE);
					MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
					break;
			}
		}
		for (i = 0; i < (end ? nb_MPI_proc - 2 : nb_MPI_proc - 1); ++i) {
			MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			nb_nodes += com_buffer[1];
			nb_sol += com_buffer[2];
			MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 99, MPI_COMM_WORLD);
		}
		ctx->solutions = nb_sol;
		ctx->nodes = nb_nodes;
	} else {
		struct sparse_array_t * active_options = ctx->active_options[chosen_item];
		if (sparse_array_empty(active_options))
			return;		/* échec : impossible de couvrir chosen_item */
		cover(instance, ctx, chosen_item);
		ctx->num_children[ctx->level] = active_options->len;
		
		while (!end) {
			MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			switch (status.MPI_TAG)
			{
				case 1 :
					ctx->nodes = nb_nodes = com_buffer[1] + 1;
					ctx->solutions = nb_sol = com_buffer[2];
					ctx->child_num[ctx->level] = com_buffer[0];
					// if (debug)
					// 	printf("Node %d received branch %llu to solve.\n", my_MPI_rank, com_buffer[0]);
					solve_MPI(instance, ctx, chosen_item, active_options->p[com_buffer[0]]);

					com_buffer[1] = ctx->nodes - nb_nodes;
					com_buffer[2] = ctx->solutions - nb_sol;
					if (ctx->solutions >= max_solutions)
						MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, 0, 99, MPI_COMM_WORLD);
					else
						MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, 0, 2, MPI_COMM_WORLD);
					break;

				case 99 :
					if (debug)
						printf("Stop received on node %d.\n", my_MPI_rank);
					end = true;
					break;
			}
		}
		uncover(instance, ctx, chosen_item);	/* backtrack */
	}
}


int main(int argc, char **argv)
{
	option_setup(argc, argv);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_MPI_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_MPI_rank);
	printf("Hello from process %d! I have %d threads available.\n", my_MPI_rank, nb_OMP_thread);

	nb_OMP_thread = omp_get_max_threads();

	struct instance_t * instance = NULL;
	struct context_t ** ctxs = NULL;
	int common_item;

	if (my_MPI_rank)
		receive_pb_data(&instance, &ctxs, &common_item);
	else
		send_pb_data(&instance, &ctxs, &common_item);

	MPI_Barrier(MPI_COMM_WORLD);
	start = wtime();
	comm_solve_MPI(instance, ctxs, common_item, false);
	MPI_Barrier(MPI_COMM_WORLD);

	if (!my_MPI_rank)
		printf("FINI. Trouvé %lld solutions en %.2fs\n", ctx->solutions, wtime() - start);

	free(instance);
	free(ctx);
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
