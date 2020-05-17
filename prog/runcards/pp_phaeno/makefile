# Sherpa Executable
SHERPA = Sherpa

# Runcard Directory
RUNCARD_DIR = runcards

# Analysis Directory
ANALYSIS_DIR = qqgg_proton

# Name of the Analysis
ANALYSIS_NAME = MC_DIPHOTON_PROTON
ANALYSIS_SO_NAME = RivetMcDiphotonProton

# Output Directory
OUT_DIR = out

# Event Count
EVENT_COUNT = 1000

# List of Runcards, see RUNCARD_FOLDER
RUNCARDS = basic test

###############################################################################
#                                    Macros                                   #
###############################################################################

THIS_DIR := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SET_RIVET := RIVET_ANALYSIS_PATH="$(THIS_DIR)$(ANALYSIS_DIR)"
SHERPA_CMD := $(SHERPA) -e $(EVENT_COUNT) "RESULT_DIRECTORY: Results" "ANALYSIS_OUTPUT: analysis"
OUT_YODA_FILES := $(foreach card,$(RUNCARDS),$(OUT_DIR)/$(card).yoda)
all: analysis_cpp analyses histograms

###############################################################################
#                             Build Rivet Analysis                            #
###############################################################################

analysis_cpp: $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so

$(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so: $(ANALYSIS_DIR)/$(ANALYSIS_NAME).cc
	cd $(ANALYSIS_DIR) \
	&& rivet-build $(ANALYSIS_SO_NAME).so $(ANALYSIS_NAME).cc

###############################################################################
#                                 Run Analyses                                #
###############################################################################

.PHONY: analyses
analyses: $(OUT_YODA_FILES)
$(OUT_YODA_FILES): $(OUT_DIR)/%.yoda: $(RUNCARD_DIR)/%/Sherpa.yaml $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so
	$(eval R_DIR:=$(dir $<))
	$(eval R_NAME:=$(shell basename $(R_DIR)))
	$(eval YODA_FILE:=$(R_NAME).yoda)
	cd $(R_DIR) && $(SHERPA_CMD) -e 0
	cd $(R_DIR) && \
		$(SET_RIVET) mpirun --use-hwthread-cpus \
			--use-hwthread-cpus $(SHERPA_CMD) "EVT_OUTPUT: 0"
	cd $(R_DIR) && yodamerge *.yoda -o $(YODA_FILE) \
		&& mv $(YODA_FILE) $(THIS_DIR)$(OUT_DIR)/

###############################################################################
#                               Make Histograms                               #
###############################################################################

.PHONY: histograms
histograms: $(OUT_DIR)/analysis/index.html
$(OUT_DIR)/analysis/index.html: $(OUT_YODA_FILES)
	cd $(OUT_DIR) && $(SET_RIVET) rivet-mkhtml *.yoda -o analysis

.PHONY: clean
clean:
	rm -rf $(OUT_DIR)
	rm -f $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so
	cd $(RUNCARD_DIR) && git clean -d -f -X


$(info $(shell mkdir -p $(OUT_DIR) && echo "Created output directory: $(OUT_DIR)"))
