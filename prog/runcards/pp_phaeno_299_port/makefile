# Sherpa Executable
SHERPA = /home/hiro/src/sherpa_rel_2_2_9/build/install/bin/Sherpa

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
EVENT_COUNT = 100000

# List of Runcards, see RUNCARD_FOLDER
RUNCARDS = with_pT_and_fragmentation_and_mi basic with_jets with_jets_and_pT with_pT_and_fragmentation

# Runcard Names
BASE_RUNCARD = Sherpa_Base.dat
SHERPA_EXT_RUNCARD = Sherpa_Ext.dat

# Event Lutput Level
EVT_OUTPUT = 2

###############################################################################
#                                    Macros                                   #
###############################################################################

NUM_CPUS = $(shell nproc)
EVENT_COUNT_CPU = $(shell echo "$(EVENT_COUNT) / $(NUM_CPUS)" | bc)
THIS_DIR := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SET_RIVET := RIVET_ANALYSIS_PATH="$(THIS_DIR)$(ANALYSIS_DIR)"
SHERPA_ARGS := -e $(EVENT_COUNT_CPU) "RESULT_DIRECTORY:=Results" "ANALYSIS_OUTPUT:=analysis" -a Rivet
OUT_YODA_FILE_NAMES := $(foreach card,$(RUNCARDS),$(card).yoda)
OUT_YODA_FILES := $(foreach name,$(OUT_YODA_FILE_NAMES),$(OUT_DIR)/$(name))
OUT_YODA_FILE_NAMES_WITH_FLAGS := $(foreach card,$(RUNCARDS),$(card).yoda:'$(shell cat $(RUNCARD_DIR)/$(card)/plot.flags)')
# check if outdir is list
ifneq ($(words $(OUT_DIR)),1)
$(error "OUT_DIR must be a single directory!")
endif

# check if outdir is /
ifeq ($(patsubst /%,,$(OUT_DIR)),)
$(error "OUT_DIR must be relative to the project root!")
endif

###############################################################################
#                             Build Rivet Analysis                            #
###############################################################################

all: analysis_cpp analyses histograms
analysis_cpp: $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so

$(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so: $(ANALYSIS_DIR)/$(ANALYSIS_NAME).cc
	cd $(ANALYSIS_DIR) \
	&& rivet-build $(ANALYSIS_SO_NAME).so $(ANALYSIS_NAME).cc -std=c++17

###############################################################################
#                                 Run Analyses                                #
###############################################################################

.PHONY: analyses
analyses: $(OUT_YODA_FILES)
$(OUT_YODA_FILES): $(OUT_DIR)/%.yoda: $(RUNCARD_DIR)/%/$(SHERPA_EXT_RUNCARD) $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so
	$(eval R_DIR:=$(dir $<))
	$(eval R_NAME:=$(shell basename $(R_DIR)))
	$(eval YODA_FILE:=$(R_NAME).yoda)

	@printf "\033[1;31mRunning Sherpa Config: $(R_NAME)\033[0m\n"
	@printf "\033[1;32mGenerating Run.dat\033[0m\n"
	sed -e '/% TEMPLATE/r $(R_DIR)$(SHERPA_EXT_RUNCARD)' $(BASE_RUNCARD) > $(R_DIR)Run.dat

	@printf "\033[1;32mRunning Sherpa\033[0m\n"
	cd $(R_DIR) && $(SHERPA) $(SHERPA_ARGS) -e 0
	cd $(R_DIR) && \
		$(SET_RIVET) $(THIS_DIR)/run_n_sherpas.fish $(NUM_CPUS) $(SHERPA_ARGS) "EVT_OUTPUT:=$(EVT_OUTPUT)"
	cd $(R_DIR) && yodamerge *.yoda -o $(YODA_FILE) \
		&& mv $(YODA_FILE) $(THIS_DIR)$(OUT_DIR)/

###############################################################################
#                               Make Histograms                               #
###############################################################################

.PHONY: histograms
histograms: $(OUT_DIR)/analysis/index.html
$(OUT_DIR)/analysis/index.html: $(OUT_YODA_FILES) $(ANALYSIS_DIR)/$(ANALYSIS_NAME).plot $(ANALYSIS_DIR)/$(ANALYSIS_NAME).info
	cd $(OUT_DIR) && $(SET_RIVET) rivet-mkhtml $(OUT_YODA_FILE_NAMES_WITH_FLAGS) -o analysis

.PHONY: clean
clean:
	rm -rf $(THIS_DIR)$(OUT_DIR)
	rm -f $(ANALYSIS_DIR)/$(ANALYSIS_SO_NAME).so
	cd $(RUNCARD_DIR) && git clean -d -f -X


# Create OUT_DIR if it is not present
ifeq ("$(wildcard $(OUT_DIR))", "")
$(info $(shell mkdir -p $(OUT_DIR) && echo "Created output directory: $(OUT_DIR)"))
endif
