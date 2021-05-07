APP = openMP MPI seq vectorize

all: $(APP)

openMP:
	@cd $@; make

MPI:
	@cd $@; make

seq:
	@cd $@; make


vectorize:
	@cd $@; make

clean:
	@$(foreach app, $(APP), cd $(app); make clean; cd ..;)

distclean: 
	@$(foreach app, $(APP), cd $(app); make distclean; cd ..;)

.PHONY: openMP MPI seq vectorize clean distclean
