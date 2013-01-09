DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_mssm_solver.cpp \
		$(DIR)/test_sm_smcw_two_scale_integration.cpp \
		$(DIR)/test_sm_two_scale.cpp \
		$(DIR)/test_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_SCRIPT := \
		$(DIR)/execute_test.sh

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(TEST_EXE)

clean-$(MODNAME):
		rm -rf $(TEST_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(TEST_DEP)
		rm -rf $(TEST_EXE)

execute-tests:  all-$(MODNAME)
		@echo "executing all tests ..."
		@./$(TEST_SCRIPT) $(TEST_EXE)
		@echo "all tests finished"

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(BOOSTLIBS)

$(DIR)/test_mssm_solver.x: $(DIR)/test_mssm_solver.o $(LIBMSSM) $(LIBFLEXI)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(FLIBS) $(BOOSTLIBS)

$(DIR)/test_running_precision.x: $(DIR)/test_running_precision.o $(LIBFLEXI)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(BOOSTLIBS)

$(DIR)/test_sm_smcw_two_scale_integration.x: $(DIR)/test_sm_smcw_two_scale_integration.o $(LIBSMCW) $(LIBSM) $(LIBFLEXI)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(BOOSTLIBS)

$(DIR)/test_sm_two_scale.x: $(DIR)/test_sm_two_scale.o $(LIBSM) $(LIBFLEXI)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(BOOSTLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(BOOSTLIBS)

%.x: %.o $(ALLLIB)
		$(CXX) $(BOOSTFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ \
		$^ $(BOOSTLIBS)

ALLDEP += $(TEST_DEP)
ALLTST += $(TEST_EXE)
