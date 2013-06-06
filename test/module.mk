DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_logger.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
TEST_SRC += \
		$(DIR)/test_two_scale_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp
ifeq ($(shell $(FSCONFIG) --with-smssm),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_mssm_solver.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-smssm --with-MSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_MSSM.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-sm),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-sm --with-smcw),yes yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm_smcw_integration.cpp
endif
endif

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
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_mssm_solver.x: $(DIR)/test_two_scale_mssm_solver.o $(LIBSMSSM) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_running_precision.x: $(DIR)/test_two_scale_running_precision.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_sm_smcw_integration.x: $(DIR)/test_two_scale_sm_smcw_integration.o $(LIBSMCW) $(LIBSM) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_sm.x: $(DIR)/test_two_scale_sm.o $(LIBSM) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_MSSM.x: $(DIR)/test_MSSM.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)

%.x: %.o $(ALLLIB)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

# add boost flags for the test object files only
$(TEST_OBJ): CPPFLAGS += $(BOOSTFLAGS)

ALLDEP += $(TEST_DEP)
ALLTST += $(TEST_EXE)
