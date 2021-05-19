#include <mpi.h>
#include <stddef.h>

#include <cover_functions.h>


int nb_MPI_proc, my_MPI_rank;


void receive_pb_data(MPI_Datatype mpi_type, struct instance_t * instance, struct context_t * ctx, int * common_item)
{
	int *to_recv = malloc(2 * sizeof(int));

	MPI_Bcast(to_recv, 3, MPI_INT, 0, MPI_COMM_WORLD);

	printf("Reçu : %d / %d / %d\n", to_recv[0], to_recv[1], to_recv[2]);
	int *recv_opt = malloc(to_recv[2] * sizeof(int));
	int *recv_items = malloc(to_recv[1] * sizeof(int));

	MPI_Bcast(recv_opt, to_recv[2], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(recv_items, to_recv[1], MPI_INT, 0, MPI_COMM_WORLD);
}

void send_pb_data(MPI_Datatype mpi_type, struct instance_t * instance, struct context_t * ctx, int * common_item)
{
	instance = load_matrix(in_filename);
	ctx = backtracking_setup(instance);
	*common_item = choose_next_item(ctx);

	struct instance_t *to_send = malloc(sizeof(struct instance_t*));
	to_send->n_items = instance->n_items;
	to_send->n_primary = instance->n_primary;
	to_send->n_options = instance->n_options;
	to_send->item_name = instance->item_name;
	to_send->options = instance->options;
	to_send->ptr = instance->ptr;

	printf("Nombre d'items : %d\n", instance->n_items);
	printf("Nombre primary : %d\n", instance->n_primary);
	printf("Nombre d'options : %d\n", instance->n_options);
	printf("Taille tableau options : %d\n", instance->n_items * instance->n_options);

	MPI_Send(&to_send,1,mpi_type,1,0, MPI_COMM_WORLD);
	//MPI_Bcast(to_send, 3, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(instance->options, instance->n_options, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(instance->ptr, instance->n_items, MPI_INT, 0, MPI_COMM_WORLD);
}


void solve_MPI(const struct instance_t * instance, struct context_t * ctx, int chosen_item, int option)
{
	choose_option(instance, ctx, option, chosen_item);
	solve(instance, ctx);
	if (ctx->solutions >= max_solutions)
		return;
	unchoose_option(instance, ctx, option, chosen_item);
}

void comm_solve_MPI(const struct instance_t * instance, struct context_t * ctx, int chosen_item)
{
	MPI_Status status;
	bool end = false;
	int i, com_buffer[3];		// node_id / nb_nodes / nb_sol
	unsigned long long nb_nodes = 0, nb_sol = 0;

	if (!my_MPI_rank) {
		for (i = 0; i < nb_MPI_proc - 1; ++i) {
			com_buffer[0] = i;
			MPI_Send(com_buffer, 3, MPI_INT, i + 1, 1, MPI_COMM_WORLD);
		}
		while (i < ctx->active_options[chosen_item]->len && !end) {
			MPI_Recv(com_buffer, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			switch (status.MPI_TAG)
			{
				case 2 :
					nb_nodes += com_buffer[1];
					nb_sol += com_buffer[2];
					com_buffer[0] = i;
					com_buffer[1] = nb_nodes;
					com_buffer[2] = nb_sol;
					printf("Received branch %d with %d solutions from process %d.\n", com_buffer[0], com_buffer[2], status.MPI_SOURCE);
					MPI_Send(com_buffer, 3, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
					break;

				case 99 :
					end = true;
					break;

				default:
					--i;
					com_buffer[1] = nb_nodes;
					com_buffer[2] = nb_sol;
					printf("Received unknown type communication. Resetting process %d.\n", status.MPI_SOURCE);
					MPI_Send(com_buffer, 3, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
					break;
			}

			for (i = 1; i < nb_MPI_proc; ++i) {
				MPI_Recv(com_buffer, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				nb_nodes += com_buffer[1];
				nb_sol += com_buffer[2];
				MPI_Send(com_buffer, 3, MPI_INT, status.MPI_SOURCE, 99, MPI_COMM_WORLD);
			}
		}
	} else {
		int k;
		if (ctx->nodes == next_report)
			progress_report(ctx);
		if (sparse_array_empty(ctx->active_items)) {
			solution_found(instance, ctx);
			return;		/* succès : plus d'objet actif */
		}

		struct sparse_array_t * active_options = ctx->active_options[chosen_item];
		if (sparse_array_empty(active_options))
			return;		/* échec : impossible de couvrir chosen_item */
		cover(instance, ctx, chosen_item);
		ctx->num_children[ctx->level] = active_options->len;
		
		while (!end) {
			MPI_Recv(com_buffer, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			switch (status.MPI_TAG)
			{
				case 1 :
					k = com_buffer[0];
					ctx->nodes = com_buffer[1] + 1;
					ctx->solutions = com_buffer[2];
					ctx->child_num[ctx->level] = k;
					solve_MPI(instance, ctx, chosen_item, active_options->p[k]);
					
					com_buffer[0] = k;
					com_buffer[1] = ctx->nodes;
					com_buffer[2] = ctx->solutions;
					MPI_Ssend(com_buffer, 3, MPI_INT, 0, 2, MPI_COMM_WORLD);

				case 99 :
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

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_MPI_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_MPI_rank);

	MPI_Datatype mpi_instance; 	/*Création d'un nouveau type (pointeur vers structure) */
	MPI_Type_commit(&mpi_instance);	/*Le déclarer pour MPI*/
	MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_CHAR, MPI_INT, MPI_INT};
	/* Tableau des types de variables (un tableau est un MPI_INT de grande taille) */

	MPI_Aint displacements[6] = {
		offsetof(struct instance_t, n_items),
		offsetof(struct instance_t, n_primary),
		offsetof(struct instance_t, n_options),
		offsetof(struct instance_t, item_name),
		offsetof(struct instance_t, options),
		offsetof(struct instance_t, ptr)
	}; /*displacement prend en compte la taille de chaque élément de la structure de données avec offserof */

	/* La longueur de chaque block (c'est là qu'on met la taille des tableaux)*/
	const int block_length[] = {1, 1, 1, instance->n_items, instance->n_options, instance->n_options + 1};

	/*Creation de la structure MPI */
	MPI_Type_create_struct(6, block_length, displacements, types, &mpi_instance);

	struct instance_t * instance = NULL;
	struct context_t * ctx = NULL;
	int * common_item = (int *) malloc(sizeof(int));

	if (my_MPI_rank)
		receive_pb_data(instance, ctx, common_item);
	else
		send_pb_data(instance, ctx, common_item);

	// comm_solve_MPI(instance, ctx, *common_item);

	MPI_Type_free(&mpi_instance);
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
