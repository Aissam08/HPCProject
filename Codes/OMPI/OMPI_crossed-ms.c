#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>

#include <cover_functions.h>

#define FACTOR 8


unsigned int nb_MPI_proc, my_MPI_rank, nb_OMP_thread;
unsigned int * thread_per_proc, * jobs;


long long int min_job()
{	
	unsigned int min_val = jobs[0];
	#pragma omp parallel for reduction(min:min_val)
	for (unsigned short int k = 1; k < (unsigned short int) nb_MPI_proc; ++k)
		min_val = jobs[k];
	return (long long int) min_val;
}

void receive_pb_data(struct instance_t ** instance, struct context_t *** ctxs, int * common_item)
{
	*instance = (struct instance_t *) malloc(sizeof(struct instance_t));
	*ctxs = (struct context_t **) malloc(nb_OMP_thread * sizeof(struct context_t *));

	MPI_Gather(&nb_OMP_thread, 1, MPI_UNSIGNED, NULL, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

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

	printf("On process %u : received %d items & %d options\n", my_MPI_rank, (*instance)->n_items, (*instance)->n_options);
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < nb_OMP_thread; ++i)
	{
		(*ctxs)[i] = backtracking_setup(*instance);
		(*ctxs)[i]->nodes = 1;
	}
	printf("On process %u : %u contexts created\n", my_MPI_rank, nb_OMP_thread);
}

void send_pb_data(struct instance_t ** instance, struct context_t ** ctx, int * common_item)
{
	*instance = load_matrix(in_filename);
	*ctx = backtracking_setup(*instance);
	*common_item = choose_next_item(*ctx);

	printf("From Master : got %d items, %d primary & %d options to broadcast\n",
		(*instance)->n_items, (*instance)->n_primary, (*instance)->n_options);

	thread_per_proc = (unsigned int *) malloc(nb_MPI_proc * sizeof(unsigned int));
	MPI_Gather(&nb_OMP_thread, 1, MPI_UNSIGNED, thread_per_proc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	jobs = (unsigned int *) malloc(nb_MPI_proc * sizeof(unsigned int));
	#pragma omp parallel for
	for (unsigned int k = 0; k < nb_MPI_proc - 1; ++k)
		jobs[k] = k + thread_per_proc[k] * nb_MPI_proc * FACTOR;
	jobs[nb_MPI_proc - 1] = nb_MPI_proc - 1;

	int buffer[4] = {(*instance)->n_items, (*instance)->n_primary, (*instance)->n_options, *common_item};
	MPI_Bcast(buffer, 4, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast((*instance)->options, (*instance)->n_options * (*instance)->n_items, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast((*instance)->ptr, (*instance)->n_options + 1, MPI_INT, 0, MPI_COMM_WORLD);
}


void solve_OMP(const struct instance_t * instance, struct context_t ** ctxs, int chosen_item, long long int * buffer)
{
	unsigned int i;
	struct sparse_array_t * active_options = ctxs[0]->active_options[chosen_item];
	long long int nb_nodes = buffer[1] + 1;
	long long int nb_sol = buffer[2];

	#pragma omp parallel for schedule(static,1)
	for (i = 0; i < nb_OMP_thread; ++i)
	{
		ctxs[i]->nodes = nb_nodes;
		ctxs[i]->solutions = nb_sol;
	}

	unsigned int bound = buffer[0] + nb_OMP_thread * nb_MPI_proc * FACTOR;
	if ((int) bound > active_options->len)
		bound = (unsigned int) active_options->len;

	#pragma omp parallel for schedule(dynamic)
	for (i = buffer[0]; i < bound; i += nb_MPI_proc)
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

	buffer[1] = 0LL; buffer[2] = 0LL;
	#pragma omp parallel for schedule(static,1) reduction(+:buffer[1])
	for (i = 0; i < nb_OMP_thread; ++i)
		buffer[1] += ctxs[i]->nodes;
	buffer[1] -= nb_nodes * nb_OMP_thread;

	#pragma omp parallel for schedule(static,1) reduction(+:buffer[2])
	for (i = 0; i < nb_OMP_thread; ++i)
		buffer[2] += ctxs[i]->solutions;
	buffer[2] -= nb_sol * nb_OMP_thread;
}


void solve_OMPI_master(struct context_t * ctx, int chosen_item, bool debug)
{
	MPI_Status status;
	bool end = false;
	long long int com_buffer[3] = {0LL, 0LL, 0LL};		// node_id / nb_nodes / nb_sol
	long long int nb_nodes = 0LL, nb_sol = 0LL, leap;

	while (!end) {
		MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		switch (status.MPI_TAG)
		{
			case 2 :
				leap = thread_per_proc[status.MPI_SOURCE] * FACTOR * nb_MPI_proc;
				nb_nodes += com_buffer[1];
				nb_sol += com_buffer[2];

				if (debug)
					printf("Received branches %lld:%lld step %u with %lld new solutions from process %d\n", com_buffer[0],
						com_buffer[0] + leap - nb_MPI_proc, nb_MPI_proc, com_buffer[2], status.MPI_SOURCE);

				com_buffer[0] = min_job();
				jobs[com_buffer[0] % nb_MPI_proc] = com_buffer[0] + leap;

				com_buffer[1] = nb_nodes;
				com_buffer[2] = nb_sol;
				if (nb_sol < max_solutions && com_buffer[0] < (unsigned int) ctx->active_options[chosen_item]->len)
					MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
				else {
					end = true;
					MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 99, MPI_COMM_WORLD);
				}
				break;

			default:
				com_buffer[1] = nb_nodes;
				com_buffer[2] = nb_sol;
				printf("Received unknown type communication. Resetting process %d\n", status.MPI_SOURCE);
				MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
				break;
		}
	}
	for (unsigned int i = 1; i < nb_MPI_proc - 1; ++i) {
		MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (debug)
			printf("Received branches %lld:%lld step %u with %lld new solutions from process %d\n", com_buffer[0],
				com_buffer[0] + thread_per_proc[status.MPI_SOURCE] * FACTOR * nb_MPI_proc - nb_MPI_proc, nb_MPI_proc, com_buffer[2], status.MPI_SOURCE);
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
	unsigned int i = 0;
	long long int com_buffer[3] = {(long long int) my_MPI_rank - 1, 0LL, 0LL};		// node_id / nb_nodes / nb_sol

	struct sparse_array_t * active_options = ctxs[0]->active_options[chosen_item];
	if (sparse_array_empty(active_options))
		return;		/* échec : impossible de couvrir chosen_item */

	#pragma omp parallel for schedule(static, 1)
	for (i = 0; i < nb_OMP_thread; ++i)
	{
		cover(instance, ctxs[i], chosen_item);
		ctxs[i]->num_children[0] = active_options->len;
	}

	solve_OMP(instance, ctxs, chosen_item, com_buffer);
	MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, 0, 2, MPI_COMM_WORLD);

	while (!end) {
		MPI_Recv(com_buffer, 3, MPI_LONG_LONG_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		switch (status.MPI_TAG)
		{
			case 1 :
				if (debug)
					printf("Node %u received branch %llu to solve\n", my_MPI_rank, com_buffer[0]);
				solve_OMP(instance, ctxs, chosen_item, com_buffer);
				MPI_Send(com_buffer, 3, MPI_LONG_LONG_INT, 0, 2, MPI_COMM_WORLD);
				break;

			case 99 :
				if (debug)
					printf("Stop received on node %u\n", my_MPI_rank);
				end = true;
				break;
		}
	}

	#pragma omp parallel for schedule(static,1)
	for (i = 0; i < nb_OMP_thread; ++i)
		uncover(instance, ctxs[i], chosen_item);		/* backtrack */
}


int main(int argc, char **argv)
{
	option_setup(argc, argv);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, (int *) &nb_MPI_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, (int *) &my_MPI_rank);
	#ifdef _OPENMP
	nb_OMP_thread = (unsigned int) omp_get_max_threads();
	#endif

	printf("Hello from process %u! I have %u threads available.\n", my_MPI_rank, nb_OMP_thread);

	struct instance_t * instance = NULL;
	int common_item;

	if (my_MPI_rank) {
		struct context_t ** ctxs = NULL;
		receive_pb_data(&instance, &ctxs, &common_item);

		MPI_Barrier(MPI_COMM_WORLD);

		start = wtime();
		solve_OMPI_slave(instance, ctxs, common_item, false);

		MPI_Barrier(MPI_COMM_WORLD);

	} else {
		struct context_t * ctx = NULL;
		send_pb_data(&instance, &ctx, &common_item);

		MPI_Barrier(MPI_COMM_WORLD);

		start = wtime();
		solve_OMPI_master(ctx, common_item, false);

		MPI_Barrier(MPI_COMM_WORLD);

		printf("FINI. Trouvé %lld solutions en %.2fs, %lld noeud parcouru\n", ctx->solutions, wtime() - start, ctx->nodes);
	}

	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
