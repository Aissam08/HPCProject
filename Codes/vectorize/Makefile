CC := gcc-10
CCFLAGS := -Wall -Wextra -g -O3
INFLAGS := -I../../includes
LFLAGS := -fopenmp

SRC := $(wildcard *.c)
EXEC := $(SRC:%.c=%)

SAVES := 


all: 
	@make $(EXEC)

%: %.c
	$(CC) $(CCFLAGS) $^ -o $@ $(INFLAGS) $(LFLAGS)


clean:
	@rm -f $(EXEC)
	@echo "Removed executable files"

distclean: clean
	@rm -f $(SAVES)
	@echo "Removed save files"

.PHONY: clean distclean