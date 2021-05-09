CC := mpicc
CCFLAGS := -Wall -O3
LDFLAGS := -I../../includes -L../../lib -lcoverfunc
LFLAGS := 

SRC := $(wildcard *.c)
EXEC := $(SRC:%.c=%)

SAVES := 


all:
	@make $(EXEC)

%: %.c
	$(CC) $(CCFLAGS) $^ -o $@ $(LDFLAGS) $(LFLAGS)


reserve_%:
	oarsub -l nodes=$(@:reserve%=%),walltime=1 -I

run_%:
	@make $(@:run_%=%)
	mpiexec --mca btl_base_warn_component_unused 0 --hostfile '$$OAR_NODEFILE' --map-by ppr:1:node $(@:run_%=%)

run: $(EXEC)
	@$(foreach exe, $(EXEC), make run_$(exe);)


clean:
	@rm -f $(EXEC)
	@echo "Removed executable files"

distclean: clean
	@rm -f $(SAVES)
	@echo "Removed save files"

.PHONY: clean distclean