# this is a makefile
# a set of rules are defined below, where each rule consists of 3 parts:
# 	target: pre-requisites
# 		command					<-- needs to be indented by 1 tab
# running make without arguments, starts the target "all"
#
# automatic variables:
# 	$<	1st pre-requisite filename
# 	$@	target filename
# 	$^	all pre-requisites' filenames, separated by spaces
#
# start with defining some variables to make life a bit easier

# filebase of all files that are going to be compiled
PROG_FILE_BASE_NAME = main focal antenna

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

all: $(EXE_NAME)

# build executable
$(EXE_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$@ $^

# build object files
$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# for debugging purposes
show: $(SOURCES) $(OBJECTS)
	echo $^

clean:
	rm -f $(BUILD_DIR)/*.o $(BIN_DIR)/$(EXE_NAME)


