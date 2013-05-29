DIR          := models/fmssm
MODNAME      := libfmssm

LIBFMSSM_SRC  :=

ifneq ($(findstring lattice,$(ALGORITHMS)),)
LIBFMSSM_SRC  += \
		$(DIR)/fmssm_lattice.cpp \
		$(DIR)/fmssm_lattice_mz_constraint.cpp \
		$(DIR)/fmssm_lattice_numerical_mz_constraint.cpp \
		$(DIR)/fmssm_lattice_msusy_constraint.cpp \
		$(DIR)/fmssm_lattice_numerical_msusy_constraint.cpp \
		$(DIR)/fmssm_lattice_mx_constraint.cpp \
		$(DIR)/fmssm_oneloop.cpp \
		$(DIR)/fmssm_lattice_rge.f \
		$(DIR)/fmssm_lattice_constraints.f \
		$(DIR)/fmssm_lattice_numerical_constraints.cpp \
		$(DIR)/fmssm_lattice_numerical_constraints_dependence.cpp \
		$(DIR)/fmssm_lattice_numerical_constraints_functions.f
endif

LIBFMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFMSSM_SRC)))

LIBFMSSM_DEP  := \
		$(LIBFMSSM_OBJ:.o=.d)

LIBFMSSM      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFMSSM)

clean-$(MODNAME):
		rm -rf $(LIBFMSSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBFMSSM_DEP)
		rm -rf $(LIBFMSSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/%.cpp : $(DIR)/%.cpp.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.f : $(DIR)/%.f.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.inc : $(DIR)/%.inc.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LIBFMSSM_DEP) $(LIBFMSSM_OBJ): CPPFLAGS += $(TVMETFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)
endif

$(LIBFMSSM): $(LIBFMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBFMSSM_DEP)
ALLLIB += $(LIBFMSSM)
