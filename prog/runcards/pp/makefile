SHERPA:=Sherpa

sherpa_xs: Sherpa.yaml
	$(SHERPA) -e 0 -a "" 'GENERATE_RESULT_DIRECTORY: false' | sed 's/\x1b\[[0-9;]*m//g' | perl -ne 'while(/.*\s:\s([0-9]+\.[0-9]+(?:e-?[0-9]+)?)\spb\s\+-\s\(\s([0-9]+\.[0-9]+(?:e-?[0-9]+)?).*$$/g){print "$$1\n$$2\n";}' | tee sherpa_xs

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
analysis: Sherpa.yaml
	RIVET_ANALYSIS_PATH="$(ROOT_DIR)/../../analysis/qqgg_proton" mpirun --use-hwthread-cpus --use-hwthread-cpus Sherpa "OUTPUT: 0"
	yodamerge *.yoda -o Analysis_all.yoda
	rivet-mkhtml Analysis_all.yoda -o analysis
	mv Analysis_all.yoda analysis/Analysis.yoda

.PHONY: clean
clean:
	find * ! -name 'flat_pdf.C' ! -name 'qqgg_simple' ! -name 'analysis' ! -name 'Sherpa.yaml' ! -name 'sherpa_xs' ! -name 'makefile' -exec rm -rf {} +
