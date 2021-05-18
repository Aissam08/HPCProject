#include <cover_functions.h>



int main(int argc, char **argv)
{
	option_setup(argc, argv);

	struct instance_t * instance = load_matrix(in_filename);
	struct context_t * ctx = backtracking_setup(instance);
	start = wtime();
	solve(instance, ctx);
	printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions, wtime() - start);
	exit(EXIT_SUCCESS);
}
