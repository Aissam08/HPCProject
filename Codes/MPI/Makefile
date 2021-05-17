CC := mpicc
CCFLAGS := -Wall -Wextra -g -O3
INFLAGS := -I../../includes
LFLAGS := 

SRC := $(wildcard *.c)
EXEC := $(SRC:%.c=%)

SAVES := 


all:
	@make $(EXEC)

%: %.c
	$(CC) $(CCFLAGS) $^ -o $@ $(INFLAGS) $(LFLAGS)


reserve_%:
	oarsub -l nodes=$(@:reserve%=%),walltime=1 -I

run_%:
	@make $(@:run_%=%)
	mpiexec --mca btl_base_warn_component_unused 0 --hostfile "$$OAR_NODEFILE" --map-by ppr:1:node $(@:run_%=%)

run: $(EXEC)
	@$(foreach exe, $(EXEC), make run_$(exe);)


clean:
	@rm -f $(EXEC)
	@echo "Removed executable files"

distclean: clean
	@rm -f $(SAVES)
	@echo "Removed save files"

.PHONY: clean distclean