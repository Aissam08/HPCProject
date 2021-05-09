CC := gcc-10
CCFLAGS := -Wall -O3

SRC := $(wildcard *.c)
OBJ := $(SRC:%.c=%.o)

LIB := libcoverfunc.a


all: $(LIB)

lib%.a: $(OBJ)
	ar rcs $@ $^
	@make move
	@make clean

%.o: %.c
	$(CC) $(CCFLAGS) -o $@ -c $<


clean:
	@rm -f $(OBJ)
	@echo "Removed object files"

move:
	@mv lib*.a ../lib/
	@echo "Moved library files into it's folder"

.PHONY: clean move