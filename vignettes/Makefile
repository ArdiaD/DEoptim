ifndef R_HOME
	R_HOME=$(shell R RHOME)
endif
ifndef GSCMD
        GSCMD=$(shell which gs || which gswin32c)
endif

all: DEoptim.pdf DEoptimPortfolioOptimization.pdf 
	@echo ${R_HOME} "\n"
	@echo ${R_GSCMD} "\n"

DEoptim.pdf: 
	rm -f DEoptim.tex

DEoptimPortfolioOptimization.pdf: 
	rm -f DEoptimPortfolioOptimization.tex
