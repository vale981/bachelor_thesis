sherpa_xs: Sherpa.yaml
	Sherpa  'EVENTS: 0' "ANALYSIS: None" 'GENERATE_RESULT_DIRECTORY: false' | sed 's/\x1b\[[0-9;]*m//g' | sed -n  's/2_2__d__db__P__P\s:\s\([0-9]*\.[0-9]*\)\spb\s+-\s(\s\([0-9]*\.[0-9]*e-*[0-9]*\).*/\1, \2/p' | tee sherpa_xs

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
analysis: Sherpa.yaml
	RIVET_ANALYSIS_PATH="$(ROOT_DIR)/../../analysis/qqgg_simple/" Sherpa "OUTPUT: 2"
	rivet-mkhtml Analysis.yoda -o analysis
	mv Analysis.yoda analysis/Analysis.yoda

.PHONY: clean
clean:
	find * ! -name 'qqgg_simple' ! -name 'analysis' ! -name 'Sherpa.yaml' ! -name 'sherpa_xs' ! -name 'makefile' -exec rm -rf {} +
