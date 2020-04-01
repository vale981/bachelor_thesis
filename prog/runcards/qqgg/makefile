sherpa: Sherpa.yaml
	Sherpa

.PHONY: clean
clean:
	find * ! -name 'Sherpa.yaml' ! -name 'makefile' -exec rm -rf {} +

sherpa_xs: Sherpa.yaml
	Sherpa | sed 's/\x1b\[[0-9;]*m//g' | sed -n  's/2_2__d__db__P__P\s:\s\([0-9]*\.[0-9]*\)\spb\s+-\s(\s\([0-9]*\.[0-9]*e-*[0-9]*\).*/\\\sigma = \\\SI{\1\\\pm \2}{\\\pico\\\barn}/p' > sherpa_xs
