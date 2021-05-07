CC = gcc-10
CCFLAGS := -Wall -g -O3
LDFLAGS := -I../includes
LFLAGS := 

SRC := $(wildcard *.c)
EXEC := $(SRC:%.c=%)

SAVES := 


all: 
	@make $(EXEC)

%: %.c
	$(CC) $(CCFLAGS) $^ -o $@ $(LDFLAGS)


clean:
	@rm -f $(EXEC)
	@echo "Removed executable files"

distclean: clean
	@rm -f $(SAVES)
	@echo "Removed save files"

.PHONY: clean distclean