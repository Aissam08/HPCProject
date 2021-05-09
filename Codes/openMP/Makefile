CC := gcc-10
CCFLAGS := -Wall -O3
LDFLAGS := -I../../includes -L../../lib -lcoverfunc
LFLAGS := -fopenmp

SRC := $(wildcard *.c)
EXEC := $(SRC:%.c=%)

SAVES := 


all: 
	@make $(EXEC)

%: %.c
	$(CC) $(CCFLAGS) $^ -o $@ $(LDFLAGS) $(LFLAGS)


clean:
	@rm -f $(EXEC)
	@echo "Removed executable files"

distclean: clean
	@rm -f $(SAVES)
	@echo "Removed save files"

.PHONY: clean distclean