include ../../conf/make_ext_config

DIRS  = $(patsubst %/, %, $(filter-out obj/,$(filter %/,$(shell ls -F))))
SRCS  = $(wildcard *.cc)
SRCS += $(foreach dir,$(DIRS),$(patsubst $(dir)/%.cc,%.cc,$(wildcard $(dir)/*.cc)))

OBJS = $(patsubst %.cc,./obj/%.o,$(SRCS))

# Include files could be everywhere, we use only two levels
UPPER_DIRS = $(filter-out test%,$(patsubst %/, %,$(filter %/,$(shell ls -F $(ORBIT_ROOT)/src))))
LOWER_DIRS = $(foreach dir,$(UPPER_DIRS),$(patsubst %/, $(ORBIT_ROOT)/src/$(dir)/%,$(filter %/,$(shell ls -F $(ORBIT_ROOT)/src/$(dir)))))

INCLUDES_LOCAL = $(patsubst %, -I$(ORBIT_ROOT)/src/%, $(UPPER_DIRS))
INCLUDES_LOCAL += $(filter-out %obj,$(patsubst %, -I%, $(LOWER_DIRS)))
INCLUDES_LOCAL += $(patsubst %, -I./%, $(filter %/,$(shell ls -F ./)))
INCLUDES_LOCAL += -I./

INC  = $(wildcard *.hh)
INC += $(wildcard *.h)
INC += $(foreach dir,$(DIRS),$(wildcard ./$(dir)/*.hh))
INC += $(foreach dir,$(DIRS),$(wildcard ./$(dir)/*.h))

# Add the GNU Scientific library
#LIBS += -lgsl -lgslcblas
#INCLUDES += -I$($(CONDA_PREFIX)/include)

# Wrappers CC FLAGS
WRAPPER_FLAGS = -fno-strict-aliasing

# CXXFLAGS - for shared library implementation
CXXFLAGS += -fPIC -std=c++11

# Monitor shared library
monitor_lib = monitor.so

#========rules=========================
compile: $(OBJS_WRAP) $(OBJS) $(INC)
	$(CXX)  $(LINKFLAGS) -o $(ORBIT_ROOT)/lib/$(monitor_lib) $(OBJS) $(LIBS)

./obj/wrap_%.o : wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/wrap_%.o : ./*/wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/%.o : %.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;
	
./obj/%.o : ./*/%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

clean:
	rm -rf ./obj/*.o
	rm -rf $(ORBIT_ROOT)/lib/$(monitor_lib)
