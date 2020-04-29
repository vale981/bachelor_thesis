sherpa_xs: Sherpa.yaml
	Sherpa -e 0 -a "" 'GENERATE_RESULT_DIRECTORY: false' | sed 's/\x1b\[[0-9;]*m//g' | sed -n  's/2_2__d__db__P__P\s:\s\([0-9]*\.[0-9]*\)\spb\s+-\s(\s\([0-9]*\.[0-9]*e-*[0-9]*\).*/\\\sigma = \\\SI{\1\\\pm \2}{\\\pico\\\barn}/p' > sherpa_xs

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
analysis: Sherpa.yaml
	RIVET_ANALYSIS_PATH="$(ROOT_DIR)/../../analysis/qqgg_proton" mpirun --use-hwthread-cpus --use-hwthread-cpus Sherpa "OUTPUT: 0"
	yodamerge *.yoda -o Analysis_all.yoda
	rivet-mkhtml Analysis_all.yoda -o analysis
	mv Analysis_all.yoda analysis/Analysis.yoda

.PHONY: clean
clean:
	find * ! -name 'qqgg_simple' ! -name 'analysis' ! -name 'Sherpa.yaml' ! -name 'sherpa_xs' ! -name 'makefile' -exec rm -rf {} +
