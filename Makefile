#-------------------------------------------------------------------------------
# executable name
TARGET := libpfate

# files
SRCFILES  :=  $(filter-out src/RcppExports.cpp src/r_interface.cpp $(wildcard src/pybind*.cpp), $(wildcard src/*.cpp))
PYBINDFILES := $(wildcard src/pybind*.cpp)
HEADERS := $(wildcard src/*.tpp) $(wildcard include/*.h) $(wildcard tests/*.h)
# ------------------------------------------------------------------------------

# paths
#CUDA_INSTALL_PATH ?= /usr/local/cuda#-5.0

#ROOT_DIR := /home/jjoshi/codes

# ROOT_DIR := ${shell dirname ${shell pwd}}
EXTERNAL_DIR := external
# ^ Do NOT put trailing whitespaces or comments after the above line

# include and lib dirs (esp for cuda)
INC_PATH :=  -I./inst/include #-I./CppNumericalSolvers-1.0.0
INC_PATH +=  -I./src # This is to allow inclusion of .tpp files in headers
INC_PATH += -I$(EXTERNAL_DIR)/phydro/inst/include \
            -I$(EXTERNAL_DIR)/libpspm/include \
			-I$(EXTERNAL_DIR)/flare/include 
LIB_PATH := -L$(EXTERNAL_DIR)/libpspm/lib -L./lib
PTH_PATH := $(shell python3 -m pybind11 --includes)

# flags
PROFILING_FLAGS = -g -pg
CPPFLAGS = -O3 -std=c++17 -Wall -Wextra -DPHYDRO_ANALYTICAL_ONLY $(PROFILING_FLAGS)
LDFLAGS =  $(PROFILING_FLAGS)

CPPFLAGS +=    \
-pedantic-errors  -Wcast-align \
-Wcast-qual \
-Wdisabled-optimization \
-Wformat=2 \
-Wformat-nonliteral -Wformat-security  \
-Wformat-y2k \
-Wimport  -Winit-self   \
-Winvalid-pch   \
-Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn \
-Wpacked   -Wpointer-arith \
-Wredundant-decls \
-Wstack-protector \
-Wstrict-aliasing=2 \
-Wswitch-enum \
-Wunreachable-code \
-Wvariadic-macros \
-Wwrite-strings \
# -Wswitch-default \
# -Wunused \
# -Wunused-parameter \
# -Waggregate-return -Wpadded -Wfloat-equal \
# -Wlong-long \
# -Wshadow \
# -Winline \
# -Wconversion \
# -Weffc++ \


CPPFLAGS += -Wno-sign-compare -Wno-unused-variable \
-Wno-unused-but-set-variable -Wno-float-conversion \
-Wno-unused-parameter

# libs
AR = ar
LIBS = 	 -lpspm	# additional libs
#LIBS = -lcudart 			# cuda libs

# files
OBJECTS = $(patsubst src/%.cpp, build/%.o, $(SRCFILES))

all: dir external_libs $(TARGET) apps

# $(PYBINDFILES): build/%.o : src/%.cpp
# 	g++ -c $(CPPFLAGS) $(PTH_PATH) $(INC_PATH) $< -o $@

external_libs:
	(cd $(EXTERNAL_DIR)/libpspm && $(MAKE))


python: dir $(TARGET) # $(PYBINDFILES)
	pip3 install .

dir:
	mkdir -p lib build tests/build bin

hi:
	echo $(SRCFILES)

$(TARGET): $(OBJECTS)
	$(AR) rcs lib/$(TARGET).a $(OBJECTS) $(LIBS)	
# g++ $(LDFLAGS) -o $(TARGET) $(LIB_PATH) $(OBJECTS) $(LIBS)

$(OBJECTS): build/%.o : src/%.cpp $(HEADERS)
	g++ -c $(CPPFLAGS) $(INC_PATH) $< -o $@

libclean:
	rm -f $(TARGET) build/*.o lib/*.a src/*.o bin/* log.txt gmon.out

extclean:
	(cd $(EXTERNAL_DIR)/libpspm && $(MAKE) clean)

re: clean all

clean: libclean testclean


## EXECUTABLES (APPS) ##

APP_FILES = $(wildcard apps/*.cpp)
APP_OBJECTS = $(patsubst apps/%.cpp, build/%.o, $(APP_FILES))
APP_TARGETS = $(patsubst apps/%.cpp, bin/%, $(APP_FILES))

apps: dir $(TARGET) compile_apps
	@echo $(APP_TARGETS)

compile_apps: $(APP_TARGETS)

$(APP_TARGETS): bin/% : apps/%.cpp $(HEADERS)
	g++ $(LDFLAGS) $(CPPFLAGS) $(INC_PATH) $(LIB_PATH) -o $@ $(OBJECTS) $< $(LIBS) -lpfate 


## TESTING SUITE ##

TEST_FILES = tests/lho.cpp tests/pf_test.cpp #$(wildcard tests/*.cpp)
TEST_OBJECTS = $(patsubst tests/%.cpp, tests/%.o, $(TEST_FILES))
TEST_TARGETS = $(patsubst tests/%.cpp, tests/%.test, $(TEST_FILES))
TEST_RUNS = $(patsubst tests/%.cpp, tests/%.run, $(TEST_FILES))
ADD_OBJECTS =

check: dir external_libs $(TARGET) compile_tests clean_log run_tests

compile_tests: $(TEST_TARGETS)
	
clean_log:
	@rm -f log.txt

run_tests: $(TEST_RUNS)
	
$(TEST_RUNS): tests/%.run : tests/%.test
	@echo "~~~~~~~~~~~~~~~ $< ~~~~~~~~~~~~~~~~" >> log.txt
	@time ./$< && \
		printf "%b" "\033[0;32m[PASS]\033[m" ": $* \n"  || \
		printf "%b" "\033[1;31m[FAIL]\033[m" ": $* \n"

$(TEST_OBJECTS): tests/%.o : tests/%.cpp $(HEADERS)
	g++ -c $(CPPFLAGS) $(INC_PATH) $< -o $@

$(TEST_TARGETS): tests/%.test : tests/%.o $(HEADERS)
	g++ $(LDFLAGS) -o $@ $(LIB_PATH) $(OBJECTS) $(ADD_OBJECTS) $< $(LIBS) -lpfate 

testclean:
	rm -f tests/*.o tests/*.test

recheck: testclean check

.PHONY: $(TEST_RUNS) run_tests clean testclean
# ------------------------------------------------------------------------------


website:
	R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/resources/app/bin/quarto/bin/tools'); print(Sys.getenv('RSTUDIO_PANDOC')); pkgdown::clean_site(); pkgdown::init_site(); pkgdown::build_home(); pkgdown::build_articles(); pkgdown::build_tutorials(); pkgdown::build_news()"

api:
	doxygen doxygen/Doxyfile

clean_api:
	rm -rf docs/html

#-gencode=arch=compute_10,code=\"sm_10,compute_10\"  -gencode=arch=compute_20,code=\"sm_20,compute_20\"  -gencode=arch=compute_30,code=\"sm_30,compute_30\"

#-W -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function
#-m64 -fno-strict-aliasing
#-I. -I/usr/local/cuda/include -I../../common/inc -I../../../shared//inc
#-DUNIX -O2


#g++ -fPIC   -m64 -o ../../bin/linux/release/swarming_chasing_predator obj/x86_64/release/genmtrand.cpp.o  obj/x86_64/release/simpleGL.cu.o  -L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart
#-lGL -lGLU -lX11 -lXi -lXmu -lGLEW_x86_64 -L/usr/X11R6/lib64 -lGLEW_x86_64 -L/usr/X11R6/lib64 -lglut
#-L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart
#-L/usr/lib -lgsl -lgslcblas
#-lcutil_x86_64  -lshrutil_x86_64




#CXXWARN_FLAGS := \
#	-W -Wall \
#	-Wimplicit \
#	-Wswitch \
#	-Wformat \
#	-Wchar-subscripts \
#	-Wparentheses \
#	-Wmultichar \
#	-Wtrigraphs \
#	-Wpointer-arith \
#	-Wcast-align \
#	-Wreturn-type \
#	-Wno-unused-function \
#	$(SPACE)

#CWARN_FLAGS := $(CXXWARN_FLAGS) \
#	-Wstrict-prototypes \
#	-Wmissing-prototypes \
#	-Wmissing-declarations \
#	-Wnested-externs \
#	-Wmain \
#
	
#HEADERS  := $(wildcard *.h)
	

