MODULES  := src test
CXXFLAGS := $(patsubst %, -I%, $(MODULES))
FFLAGS   := $(patsubst %, -I%, $(MODULES))
LIBS     := -lboost_unit_test_framework
SRC      :=
CXX      := g++
FC       := gfortran

include $(patsubst %, %/module.mk, $(MODULES))

OBJ      := \
	$(patsubst %.cpp, %.o, $(filter %.cpp, $(SRC))) \
	$(patsubst %.f, %.o, $(filter %.f, $(SRC)))

DEP      := $(OBJ:.o=.d)

test_two_scale_solver: $(OBJ)
	$(CXX) -o $@ $^ $(LIBS)

ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEP)
endif

%.d: %.cpp
# -MT '$*.o' ensures that the target contains the full path
	$(CXX) $(CXXFLAGS) -MM -MP -MG -o $@ -MT '$*.o' $^

%.d: %.f
# the sed script ensures that the target contains the full path
	$(FC) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

clean:
	rm -f $(OBJ)
	rm -f $(DEP)
