SHERPA:=Sherpa
ANALYSIS_DIR:=qqgg_proton
ANALYSIS_NAME:=MC_DIPHOTON_PROTON
ANALYSIS_SO_NAME:=RivetMcDiphotonProton
OUT_DIR:=out

THIS_DIR:=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SET_RIVET:=RIVET_ANALYSIS_PATH="$(THIS_DIR)$(ANALYSIS_DIR)"
SHERPA_CMD:=$(SHERPA) "RESULT_DIRECTORY: Results" "ANALYSIS_OUTPUT: analysis"

###############################################################################
#                               Integrat the XS                               #
###############################################################################

sherpa_xs: $(OUT_DIR)/sherpa_xs
$(OUT_DIR)/sherpa_xs: Sherpa.yaml
	cd $(OUT_DIR) && $(SHERPA_CMD) -e 0 -a "" \
		| sed 's/\x1b\[[0-9;]*m//g' \
		| perl -ne 'while(/.*\s:\s([0-9]+\.[0-9]+(?:e-?[0-9]+)?)\spb\s\+-\s\(\s([0-9]+\.[0-9]+(?:e-?[0-9]+)?).*$$/g){print "$$1\n$$2\n";}' \
		| tee sherpa_xs

###############################################################################
#                             Build Rivet Analysis                            #
###############################################################################

analysis_cpp: $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so

$(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so: $(ANALYSIS_DIR)/$(ANALYSIS_NAME).cc
	cd $(ANALYSIS_DIR) \
	&& rivet-build $(ANALYSIS_SO_NAME).so $(ANALYSIS_NAME).cc

###############################################################################
#                                 Run Analysis                                #
###############################################################################

analysis: $(OUT_DIR)/analysis/Analysis.yoda
$(OUT_DIR)/analysis/Analysis.yoda: sherpa_xs $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so
	cd $(OUT_DIR) && \
		$(SET_RIVET) mpirun --use-hwthread-cpus \
			--use-hwthread-cpus $(SHERPA_CMD) "OUTPUT: 0"
	cd $(OUT_DIR) && yodamerge *.yoda -o Analysis_all.yoda
	cd $(OUT_DIR) && $(SET_RIVET) rivet-mkhtml Analysis_all.yoda -o analysis
	cd $(OUT_DIR) && mv Analysis_all.yoda analysis/Analysis.yoda && rm *.yoda

.PHONY: clean
clean:
	rm -rf $(OUT_DIR)
	rm -f $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so


###############################################################################
#                   Create Directories and copy Sherpa.yaml                   #
###############################################################################

$(info $(shell mkdir -p $(OUT_DIR) && echo "Created output directory: $(OUT_DIR)"))
$(info $(shell cp Sherpa.yaml $(OUT_DIR) && echo "Copied Sherpa.yaml"))
