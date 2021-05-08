DIRS := $(wildcard */)
APP := $(DIRS:%/=%)


all: $(APP)

*:
	@cd $@; make


clean:
	@$(foreach app, $(APP), cd $(app); make clean; cd ..;)

distclean: 
	@$(foreach app, $(APP), cd $(app); make distclean; cd ..;)


.PHONY: $(APP) clean distclean
