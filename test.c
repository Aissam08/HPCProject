#include <stdio.h>
#include <stdlib.h>





void save_struct(struct context_t * ctx,const  struct instance_t * instance)
{
    FILE* f = fopen("context.txt","w");
    fprintf(f, "%d\n",ctx->level);
    fprintf(f,"%lld\n" ,ctx->nodes);
    fprintf(f, "%lld\n",ctx->solutions);
    fprintf(f,"%ls\n", ctx->child_num);
    fprintf(f,"%d\n", instance->n_items);
    fprintf(f,"%d\n", instance->n_options);

    fclose(f);
}



int get_val(FILE * f) /* retourn la dernière ligne du fichier */
{
    int ret;
    while(!feof(f))
        fscanf(f,"%d", &ret);

    return ret;
}

struct context_t * load_struct(const char * filename)
{
    struct context_t * ctx = malloc(sizeof(*ctx));
    if (ctx == NULL)
        err(1, "impossible d'allouer un contexte");
    FILE* f = fopen(filename,"r");
    /*fscanf(f, "%d", &ctx->level);
    fscanf(f,"%lld", &ctx->nodes);
    fscanf(f,"%lld", &ctx->solutions);
    fscanf(f,"%d", ctx->child_num);*/
    fclose (f);   
    return ctx;
}

void read_ints (const char* file_name)
{
    FILE* file = fopen (file_name, "r");
    int *tab = malloc(sizeof(int*));

    int i = 0;

    for(int n=0; n< 4 ; n++)
    {
        fscanf (file, "%d", &i);
        printf("read %d\n",i);
    }
    fclose (file);        
}


int main()
{

    FILE * f = fopen("file/options.txt","r");

    int k = get_val(f);

    printf("%d\n",k);

    return 0;
}