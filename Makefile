PRECOMPDIR = src
BINDIR = bin

default: all

help:
	@echo " "
	@echo "Use one of the targets below"
	@echo "	all:		build all programs"
	@echo "	pre:		build precompute code"
	@echo "	clean:		clean all builds"

all: pre

pre:
	@cd $(PRECOMPDIR); make; mv correlate ../$(BINDIR)/;

clean: clean_pre clean_bin

clean_pre:
	@cd $(PRECOMPDIR); make clean;

clean_bin:
	@cd $(BINDIR);rm *;
