DIRS := $(wildcard */)
APP := $(DIRS:%/=%)


all: $(APP)
	@$(foreach app, $(APP), cd $(app); make;)

*:
	@cd $@; make


clean:
	@$(foreach app, $(APP), cd $(app); make clean; cd ..;)

distclean: 
	@$(foreach app, $(APP), cd $(app); make distclean; cd ..;)


.PHONY: $(APP) clean distclean
