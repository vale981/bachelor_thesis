SHERPA:=/home/hiro/src/sherpa_rel_2_2_9/build/install/bin/Sherpa

sherpa_xs: Run.dat
	$(SHERPA) -e 0 | sed 's/\x1b\[[0-9;]*m//g' | perl -ne 'while(/.*\s:\s([0-9]\.[0-9]+(?:e-?[0-9]+)?)\spb\s\+-\s\(\s([0-9]\.[0-9]+(?:e-?[0-9]+)?).*$$/g){print "$$1\n$$2";}'  > sherpa_xs

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
analysis: Run.dat
	RIVET_ANALYSIS_PATH="$(ROOT_DIR)/../../analysis/qqgg_proton" mpirun --use-hwthread-cpus --use-hwthread-cpus $(SHERPA) ANALYSIS=Rivet
	yodamerge *.yoda -o Analysis_all.yoda
	rivet-mkhtml Analysis_all.yoda -o analysis
	mv Analysis_all.yoda analysis/Analysis.yoda

.PHONY: clean
clean:
	find * ! -name 'flat_pdf.C' ! -name 'qqgg_simple' ! -name 'analysis' ! -name 'Run.dat' ! -name 'sherpa_xs' ! -name 'makefile' -exec rm -rf {} +
