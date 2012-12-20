DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_mssm_solver.cpp \
		$(DIR)/test_sm_smcw_two_scale_integration.cpp \
		$(DIR)/test_sm_two_scale.cpp \
		$(DIR)/test_two_scale_solver.cpp

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(TEST_EXE)

clean-$(MODNAME):
		rm -rf $(TEST_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(TEST_DEP)
		rm -rf $(TEST_EXE)

execute-tests:  all-$(MODNAME)
		@echo "executing all tests ..."
		@rm -f test.log
		@for x in $(TEST_EXE); do      \
			echo -n "executing test: $$x ... ";              \
			echo "**************************" >> test.log;   \
			echo "* executing test: $$x ... " >> test.log;   \
			echo "**************************" >> test.log;   \
			./$$x --log_level=test_suite >> test.log 2>> test.log; \
			if [ $$? = 0 ] ; then  \
				echo "OK";     \
			else                   \
				echo "FAILED"; \
			fi                     \
		done
		@echo "all tests finished"

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_mssm_solver.x: $(DIR)/test_mssm_solver.o $(LIBMSSM) $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_sm_smcw_two_scale_integration.x: $(DIR)/test_sm_smcw_two_scale_integration.o $(LIBSMCW) $(LIBSM) $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_sm_two_scale.x: $(DIR)/test_sm_two_scale.o $(LIBSM) $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.x: %.o $(ALLLIB)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

ALLDEP += $(TEST_DEP)
ALLTST += $(TEST_EXE)
