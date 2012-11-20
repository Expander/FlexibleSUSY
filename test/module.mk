DIR                       := test
MODNAME                   := test

TEST_TWO_SCALE_SOLVER_SRC := \
		$(DIR)/test_two_scale_solver.cpp

TEST_TWO_SCALE_SOLVER_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_TWO_SCALE_SOLVER_SRC)))

TEST_TWO_SCALE_SOLVER_DEP := \
		$(TEST_TWO_SCALE_SOLVER_OBJ:.o=.d)

TEST_TWO_SCALE_SOLVER     := \
		$(DIR)/test_two_scale_solver

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(TEST_TWO_SCALE_SOLVER)

clean-$(MODNAME):
		rm -rf $(TEST_TWO_SCALE_SOLVER_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(TEST_TWO_SCALE_SOLVER_DEP)
		rm -rf $(TEST_TWO_SCALE_SOLVER)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(TEST_TWO_SCALE_SOLVER): $(TEST_TWO_SCALE_SOLVER_OBJ) $(ALLLIB)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

ALLDEP += $(TEST_TWO_SCALE_SOLVER_DEP)
ALLEXE += $(TEST_TWO_SCALE_SOLVER)
