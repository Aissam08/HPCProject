#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>
#include <stddef.h>

#include <cover_functions.h>

#define FACTOR 10


int nb_MPI_proc, my_MPI_rank, nb_OMP_thread;
unsigned char * threads_per_proc;


void receive_pb_data(struct instance_t ** instance, struct context_t *** ctx, int * common_item)
{
	*instance = (struct instance_t *) malloc(sizeof(struct instance_t));
	*ctxs = (struct context_t **) malloc(nb_OMP_thread * sizeof(struct context_t *));

	MPI_Gather(nb_OMP_thread, 1, MPI_UNSIGNED_CHAR, NULL, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	
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

	threads_per_proc = (char *) malloc(nb_MPI_proc * sizeof(char));
	MPI_Gather(NULL, 1, MPI_UNSIGNED_CHAR, threads_per_proc, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

	int buffer[4] = {(*instance)->n_items, (*instance)->n_primary, (*instance)->n_options, *common_item};
	MPI_Bcast(buffer, 4, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast((*instance)->options, (*instance)->n_options * (*instance)->n_items, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast((*instance)->ptr, (*instance)->n_options + 1, MPI_INT, 0, MPI_COMM_WORLD);
}


void solve_OMP(const struct instance_t * instance, struct context_t ** ctxs, int chosen_item, int * buffer)
{
	unsigned int i;
	struct sparse_array_t * active_options = ctxs[0]->active_options[chosen_item];
	nb_nodes = com_buffer[1] + 1;
	nb_sol = com_buffer[2];


	#pragma omp parallel for schedule(static, 1) reduction(+:nb_sol)
	for (i = 0; i < nb_OMP_thread; ++i)
	{
		ctxs[i]->nodes = nb_nodes;
		ctxs[i]->solutions = nb_sol;
	}

	unsigned int min_bound = com_buffer[0] + nb_OMP_thread * FACTOR < active_options->len ?
		com_buffer[0] + nb_OMP_thread * FACTOR : active_options->len;

	#pragma omp parallel for schedule(dynamic)
	for (i = com_buffer[0]; i < min_bound; ++i)
	{
		unsigned short int my_thread = omp_get_thread_num();
		int option = active_options->p[i];
		ctxs[my_thread]->child_num[0] = i;
		choose_option(instance, ctxs[my_thread], option, chosen_item);
		solve(instance, ctxs[my_thread]);

		if (ctxs[my_thread]->solutions >= max_solutions)
			exit(EXIT_SUCCESS);
		unchoose_option(instance, ctxs[my_thread], option, chosen_item);
	}

	com_buffer[1] = 0; com_buffer[2] = 0;
	#pragma omp parallel for schedule(static, 1) reduction(+:com_buffer[1]) reduction(+:com_buffer[2])
	for (i = 0; i < nb_OMP_thread; ++i) {
		com_buffer[1] = ctxs[i]->nodes - nb_nodes;
		com_buffer[2] = ctxs[i]->solutions - nb_sol;
	}
}


void solve_OMPI_master(const struct instance_t * instance, struct context_t * ctx, int chosen_item, bool debug)
{
	MPI_Status status;
	bool end = false;
	long long int i = 0LL, com_buffer[3] = {0LL, 0LL, 0LL};		// node_id / nb_nodes / nb_sol
	long long int nb_nodes = 0LL, nb_sol = 0LL;

	#pragma omp for parallel reduction(+:i)
	for (unsigned short int k = 1; k < nb_MPI_proc; ++k)
		i += (long long int) threads_per_proc[k] * FACTOR;

	while (i < ctx->active_options[chosen_item]->len && !end) {
		MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		switch (status.MPI_TAG)
		{
			case 2 :
				if (debug)
					printf("Received branches %llu:%llu with %llu new solutions from process %d.\n", com_buffer[0],
						com_buffer[0] + threads_per_proc[status.MPI_SOURCE] * FACTOR - 1, com_buffer[2], status.MPI_SOURCE);
				nb_nodes += com_buffer[1];
				nb_sol += com_buffer[2];
				com_buffer[0] = i;
				com_buffer[1] = nb_nodes;
				com_buffer[2] = nb_sol;
				i += (long long int) threads_per_proc[status.MPI_SOURCE] * FACTOR;
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
}


void solve_OMPI_slave(const struct instance_t * instance, struct context_t ** ctxs, int chosen_item, bool debug)
{
	MPI_Status status;
	bool end = false;
	long long int i = 0LL, com_buffer[3] = {0LL, 0LL, 0LL};		// node_id / nb_nodes / nb_sol

	struct sparse_array_t * active_options = ctxs[0]->active_options[chosen_item];
	if (sparse_array_empty(active_options))
		return;		/* échec : impossible de couvrir chosen_item */

	#pragma omp parallel for schedule(static, 1)
	for (i = 0; i < nb_OMP_thread; ++i)
	{
		cover(instance, ctxs[i], chosen_item);
		ctxs[i]->num_children[0] = active_options->len;
	}
	
	while (!end) {
		MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		switch (status.MPI_TAG)
		{
			case 1 :
				if (debug)
					printf("Node %d received branch %llu to solve.\n", my_MPI_rank, com_buffer[0]);
				solve_OMP(instance, ctx, chosen_item, com_buffer);

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

	#pragma omp parallel for schedule(static, 1)
	for (i = 0; i < nb_OMP_thread; ++i)
		uncover(instance, ctxs[i], chosen_item);		/* backtrack */
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
		send_pb_data(&instance, ctxs, &common_item);

	MPI_Barrier(MPI_COMM_WORLD);
	start = wtime();
	if (my_MPI_rank)
		solve_OMPI_slave(instance, ctxs, common_item, false);
	else
		solve_OMPI_master(instance, *ctxs, common_item, true);
	MPI_Barrier(MPI_COMM_WORLD);

	if (!my_MPI_rank)
		printf("FINI. Trouvé %lld solutions en %.2fs\n", ctxs[0]->solutions, wtime() - start);


	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
