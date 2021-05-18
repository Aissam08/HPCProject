#include <mpi.h>

#include <cover_functions.h>


int nb_MPI_proc, my_MPI_rank;


void receive_pb_data(struct instance_t * instance, struct context_t * ctx) 
{
	int *to_recv = malloc(2 * sizeof(int));

	MPI_Bcast(to_recv, 3, MPI_INT, 0, MPI_COMM_WORLD);

	printf("Reçu : %d / %d / %d\n", to_recv[0], to_recv[1], to_recv[2]);
	int *recv_opt = malloc(to_recv[2] * sizeof(int));
	int *recv_items = malloc(to_recv[1] * sizeof(int));

	MPI_Bcast(recv_opt, to_recv[2], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(recv_items, to_recv[1], MPI_INT, 0, MPI_COMM_WORLD);
}

void send_pb_data(struct instance_t * instance)
{
	instance = load_matrix(in_filename);

	printf("Nombre d'items : %d\n", instance->n_items);
	printf("Nombre primary : %d\n", instance->n_primary);
	printf("Nombre d'options : %d\n", instance->n_options);
	printf("Taille tableau options : %d\n", instance->n_items * instance->n_options);

	int to_send[] = {instance->n_items, instance->n_primary, instance->n_options};

	MPI_Bcast(to_send, 3, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(instance->options, instance->n_options, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(instance->ptr, instance->n_items, MPI_INT, 0, MPI_COMM_WORLD);
}


void solve_MPI(const struct instance_t * instance, struct context_t * ctx, MPI_Status status)
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
		solve(instance, ctx);
		if (ctx->solutions >= max_solutions)
			return;
		unchoose_option(instance, ctx, option, chosen_item);
	}
	uncover(instance, ctx, chosen_item);	/* backtrack */
}


int main(int argc, char **argv)
{
	option_setup(argc, argv);

	// MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_MPI_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_MPI_rank);

	// MPI_Type_struct(3, B, D, type1, newtype)
	// --> { (2,1,2), (0,16,26), (MPI_FLOAT, {(double,0), (char,8) } , MPI_CHAR), name_of_type }
		// - Il y aura 3 composantes dans la structure
		// 	- longueur des données (2 / 1 / 3)
		// 	- La structure contient 2 FLOAT, 1 map et 3 char
		// 	- le map est composé d'un double et d'un char 
		// 	- longueur d'écarts (0,16,26)
		// 	- le MPI_FLOAT est à 0 / 
		// 	- le type1 (map) commence à 16 (après 2 float de 8 octet)
		// 	- le MPI_CHAR commence à 26 octet (16 + 9)

	struct instance_t * instance = NULL;
	struct context_t * ctx = NULL;

	if (my_MPI_rank == 0)
		send_pb_data(instance);
	else
		receive_pb_data(instance, ctx);




	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
