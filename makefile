all: 
	+@$(MAKE) -C src 
	+@$(MAKE) -C grid_gen

clean: 
	+@$(MAKE) clean -C src
	+@$(MAKE) clean -C grid_gen
