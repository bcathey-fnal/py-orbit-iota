include ../../conf/make_root_config

DIRS  = $(patsubst %/, %, $(filter-out obj/,$(filter %/,$(shell ls -F))))
SRCS  = $(wildcard *.cc)
SRCS += $(foreach dir,$(DIRS),$(patsubst $(dir)/%.cc,%.cc,$(wildcard $(dir)/*.cc)))

OBJS = $(patsubst %.cc,./obj/%.o,$(SRCS))

#include files could be everywhere, we use only two levels
UPPER_DIRS = $(filter-out test%,$(patsubst %/, %,$(filter %/,$(shell ls -F ..))))
LOWER_DIRS = $(foreach dir,$(UPPER_DIRS),$(patsubst %/, ../$(dir)/%,$(filter %/,$(shell ls -F ../$(dir)))))

INCLUDES_LOCAL = $(patsubst %, -I../%, $(UPPER_DIRS))
INCLUDES_LOCAL += $(filter-out %obj,$(patsubst %, -I%, $(LOWER_DIRS)))

# Include files from extensions if present
EXT_OUTER_DIRS = $(patsubst %/, ../../ext/%, $(filter %/,$(shell ls -F ../../ext/)))
EXT_INNER_DIRS = $(foreach dir,$(EXT_OUTER_DIRS),$(patsubst %/, $(dir)/%, $(filter-out obj/, $(filter %/,$(shell ls -F $(dir))))))

INCLUDES_EXT = -I../../ext
INCLUDES_EXT += $(patsubst %, -I%, $(EXT_OUTER_DIRS))
INCLUDES_EXT += $(patsubst %, -I%, $(EXT_INNER_DIRS))

INC  = $(wildcard *.hh)
INC += $(wildcard *.h)
INC += $(foreach dir,$(DIRS),$(wildcard ./$(dir)/*.hh))
INC += $(foreach dir,$(DIRS),$(wildcard ./$(dir)/*.h))

#wrappers CC FLAGS
WRAPPER_FLAGS = -fno-strict-aliasing

#========rules=========================
compile: $(OBJS_WRAP) $(OBJS) $(INC)

./obj/wrap_%.o : wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/wrap_%.o : ./*/wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/%.o : %.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES_EXT) $(INCLUDES) -c $< -o $@;
	
./obj/%.o : ./*/%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

clean:
	rm -rf ./obj/*.o
