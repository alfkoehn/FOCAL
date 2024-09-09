# this is a makefile
# a set of rules are defined below, where each rule consists of 3 parts:
# 	target: pre-requisites
# 		command1				<-- needs to be indented by 1 tab
# 		command2
# having an "@" before the command, will hide it from the console
# running make without arguments, starts the target "all"
#
# note that a pipe is used to separate normal and order-only pre-requisites
#   target: normal pre-requisites | order-only pre-requisites
# --> target will not be updated, if order-only pre-requisites is changed  
# --> typical use-case for this is when target are to be placed in separate
#     directory, which might not exist before make is run
#
# automatic variables:
# 	$<	1st pre-requisite filename
# 	$@	target filename
# 	$^	all pre-requisites' filenames, separated by spaces
#
# start with defining some variables to make life a bit easier

# filebase of all files that are going to be compiled
PROG_FILE_BASE_NAME = main focal antenna grid_io background_profiles power_calc save_data cJSON initialize_grid boundaries antenna_detector

# name of output (executable) file
EXE_NAME = exe

# define all relevant directories
BUILD_DIR  = build
HEADER_DIR = include
SOURCE_DIR = src
BIN_DIR    = bin

# list with source and object files, using previously defines filebase-names
SOURCES := $(addprefix $(SOURCE_DIR)/, $(addsuffix .c, $(PROG_FILE_BASE_NAME) ) )
OBJECTS := $(addprefix $(BUILD_DIR)/,  $(addsuffix .o, $(PROG_FILE_BASE_NAME) ) )

# set compiler and flags used during compilation
# NOTE: might be better to only use h5cc when required
# 		==> set 2 sets of compiler and flags to always use minimal required options
CC = h5cc
CFLAGS = -Wall -O2 -fopenmp -I$(HEADER_DIR)

# set command for creating folder
MKDIR_P	= mkdir -p

all: $(EXE_NAME)

# build executable
$(EXE_NAME): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$@ $^

# build object files
$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# create build directory if it does not exist
# (basically checks if target (build-dir) exists, if not, then run command)
$(BUILD_DIR):
	$(MKDIR_P) $@

# create bin directory if it does not exist
$(BIN_DIR):
	$(MKDIR_P) $@
	
# for debugging purposes
show: $(SOURCES) $(OBJECTS)
	echo $^

# clean build and bin directories
clean:
	rm -f $(BUILD_DIR)/*.o $(BIN_DIR)/$(EXE_NAME)

