DIR          := templates
MODNAME      := templates

TEMPLATES    := \
		$(DIR)/convergence_tester.hpp.in \
		$(DIR)/high_scale_constraint.hpp.in \
		$(DIR)/initial_guesser.hpp.in \
		$(DIR)/info.hpp.in \
		$(DIR)/info.cpp.in \
		$(DIR)/input_parameters.hpp.in \
		$(DIR)/low_scale_constraint.hpp.in \
		$(DIR)/model.hpp.in \
		$(DIR)/physical.hpp.in \
		$(DIR)/physical.cpp.in \
		$(DIR)/run.cpp.in \
		$(DIR)/scan.cpp.in \
		$(DIR)/slha_io.hpp.in \
		$(DIR)/slha_io.cpp.in \
		$(DIR)/susy_scale_constraint.hpp.in \
		$(DIR)/lattice_convergence_tester.hpp.in \
		$(DIR)/lattice_convergence_tester.cpp.in \
		$(DIR)/lattice_high_scale_constraint.hpp.in \
		$(DIR)/lattice_high_scale_constraint.cpp.in \
		$(DIR)/lattice_initial_guesser.hpp.in \
		$(DIR)/lattice_initial_guesser.cpp.in \
		$(DIR)/lattice_initial_guesser_low_scale_model.hpp.in \
		$(DIR)/lattice_initial_guesser_low_scale_model.cpp.in \
		$(DIR)/lattice_low_scale_constraint.hpp.in \
		$(DIR)/lattice_low_scale_constraint.cpp.in \
		$(DIR)/lattice_model.hpp.in \
		$(DIR)/lattice_model.cpp.in \
		$(DIR)/lattice_susy_scale_constraint.hpp.in \
		$(DIR)/lattice_susy_scale_constraint.cpp.in \
		$(DIR)/spectrum_generator.hpp.in \
		$(DIR)/low_scale_spectrum_generator.hpp.in \
		$(DIR)/two_scale_convergence_tester.hpp.in \
		$(DIR)/two_scale_convergence_tester.cpp.in \
		$(DIR)/two_scale_high_scale_constraint.hpp.in \
		$(DIR)/two_scale_high_scale_constraint.cpp.in \
		$(DIR)/two_scale_initial_guesser.hpp.in \
		$(DIR)/two_scale_initial_guesser.cpp.in \
		$(DIR)/two_scale_initial_guesser_low_scale_model.hpp.in \
		$(DIR)/two_scale_initial_guesser_low_scale_model.cpp.in \
		$(DIR)/two_scale_low_scale_constraint.hpp.in \
		$(DIR)/two_scale_low_scale_constraint.cpp.in \
		$(DIR)/two_scale_model.hpp.in \
		$(DIR)/two_scale_model.cpp.in \
		$(DIR)/two_scale_soft_beta_.cpp.in \
		$(DIR)/two_scale_soft_parameters.hpp.in \
		$(DIR)/two_scale_soft_parameters.cpp.in \
		$(DIR)/two_scale_susy_beta_.cpp.in \
		$(DIR)/two_scale_susy_parameters.hpp.in \
		$(DIR)/two_scale_susy_parameters.cpp.in \
		$(DIR)/two_scale_susy_scale_constraint.hpp.in \
		$(DIR)/two_scale_susy_scale_constraint.cpp.in \
		$(DIR)/utilities.hpp.in \
		$(DIR)/utilities.cpp.in

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
