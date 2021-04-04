CC = gcc-10
CCFLAGS = -Wall -g -O3
LDFLAGS = -fopenmp

DIR = Codes

SRC = $(wildcard $(DIR)/*.c)
OBJ = $(SRC:%.c=%.o)

CMD = $(CC) $(CCFLAGS) $(filter-out $(EXEC:%=$(DIR)/%.o), $(OBJ)) $(exe:%=%.o) -o $(exe) $(LDFLAGS)


EXEC = exact_cover


all: $(EXEC:%=$(DIR)/%)

$(EXEC:%=$(DIR)/%): $(OBJ)
	$(foreach exe, $@, $(CMD))

$(DIR)/%.o: $(DIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<


clean:
	@rm -f $(OBJ) $(EXEC:%=$(DIR)/%)
	@echo "Removed exe & object files"


.PHONY: clean
