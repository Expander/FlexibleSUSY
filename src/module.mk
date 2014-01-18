DIR          := src
MODNAME      := libflexisusy

LIBFLEXI_SRC := \
		$(DIR)/betafunction.cpp \
		$(DIR)/command_line_options.cpp \
		$(DIR)/def.cpp \
		$(DIR)/dilog.f \
		$(DIR)/error.cpp \
		$(DIR)/gsl_utils.cpp \
		$(DIR)/linalg.cpp \
		$(DIR)/lowe.cpp \
		$(DIR)/numerics.cpp \
		$(DIR)/program_options.cpp \
		$(DIR)/rge.cpp \
		$(DIR)/rk.cpp \
		$(DIR)/scan.cpp \
		$(DIR)/slha_io.cpp \
		$(DIR)/stopwatch.cpp \
		$(DIR)/utils.cpp \
		$(DIR)/wrappers.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBFLEXI_SRC += \
		$(DIR)/two_scale_composite_convergence_tester.cpp \
		$(DIR)/two_scale_convergence_tester.cpp \
		$(DIR)/two_scale_running_precision.cpp \
		$(DIR)/two_scale_solver.cpp
endif

ifneq ($(findstring lattice,$(ALGORITHMS)),)
LIBFLEXI_SRC += \
		$(DIR)/lattice_model.cpp \
		$(DIR)/lattice_constraint.cpp \
		$(DIR)/lattice_numerical_constraint.cpp \
		$(DIR)/lattice_solver.cpp \
		$(DIR)/SM.cpp \
		$(DIR)/fortran_utils.f
endif

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFLEXI)

clean-$(MODNAME):
		-rm -f $(LIBFLEXI_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(LIBFLEXI_DEP)
		-rm -f $(LIBFLEXI)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(EIGENFLAGS)

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(GSLFLAGS) $(BOOSTFLAGS)
endif

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MAKELIB) $@ $^
else
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MAKELIB) $@ $^ $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(LAPACKLIBS)
endif

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
