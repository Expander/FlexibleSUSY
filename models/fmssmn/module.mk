DIR          := models/fmssmn
MODNAME      := libfmssmn

LIBFMSSMN_SRC  :=

ifneq ($(findstring lattice,$(ALGORITHMS)),)
LIBFMSSMN_SRC  += \
		$(DIR)/fmssmn_lattice.cpp \
		$(DIR)/fmssmn_lattice_mz_constraint.cpp \
		$(DIR)/fmssmn_lattice_msusy_constraint.cpp \
		$(DIR)/fmssmn_lattice_mx_constraint.cpp \
		$(DIR)/fmssmn_lattice_rge.f \
		$(DIR)/fmssmn_lattice_constraints.f \
		$(DIR)/fmssmn_lattice_numerical_constraints.cpp \
		$(DIR)/fmssmn_lattice_numerical_constraints_dependence.cpp \
		$(DIR)/fmssmn_lattice_numerical_constraints_functions.f \
		$(DIR)/fmssm_fmssmn_lattice_matchings.f \
		$(DIR)/fmssm_fmssmn_lattice_numerical_matchings.cpp \
		$(DIR)/fmssm_fmssmn_lattice_numerical_matchings_dependence.cpp \
		$(DIR)/fmssm_fmssmn_lattice_numerical_matchings_functions.f
endif

LIBFMSSMN_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFMSSMN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFMSSMN_SRC)))

LIBFMSSMN_DEP  := \
		$(LIBFMSSMN_OBJ:.o=.d)

LIBFMSSMN      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFMSSMN)

clean-$(MODNAME):
		rm -rf $(LIBFMSSMN_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBFMSSMN_DEP)
		rm -rf $(LIBFMSSMN)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/%.cpp : $(DIR)/%.cpp.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.f : $(DIR)/%.f.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.inc : $(DIR)/%.inc.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LIBFMSSMN_OBJ): CPPFLAGS += $(TVMETFLAGS) $(GSLFLAGS)
endif

$(LIBFMSSMN): $(LIBFMSSMN_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBFMSSMN_DEP)
ALLLIB += $(LIBFMSSMN)
