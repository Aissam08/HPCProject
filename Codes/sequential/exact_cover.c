#include <cover_functions.h>


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
	struct context_t * ctx = backtracking_setup(instance);
	start = wtime();
	solve(instance, ctx);
	printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions, 
			wtime() - start);
	exit(EXIT_SUCCESS);
}
