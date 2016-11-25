DIR          := templates
MODNAME      := templates

BASE_TEMPLATES := \
		$(DIR)/a_muon.hpp.in \
		$(DIR)/a_muon.cpp.in \
		$(DIR)/effective_couplings.hpp.in \
		$(DIR)/effective_couplings.cpp.in \
		$(DIR)/info.hpp.in \
		$(DIR)/info.cpp.in \
		$(DIR)/input_parameters.hpp.in \
		$(DIR)/input_parameters.cpp.in \
		$(DIR)/mass_eigenstates.hpp.in \
		$(DIR)/mass_eigenstates.cpp.in \
		$(DIR)/model.hpp.in \
		$(DIR)/model_slha.hpp.in \
		$(DIR)/observables.hpp.in \
		$(DIR)/observables.cpp.in \
		$(DIR)/physical.hpp.in \
		$(DIR)/physical.cpp.in \
		$(DIR)/plot_rgflow.gnuplot.in \
		$(DIR)/plot_spectrum.gnuplot.in \
		$(DIR)/soft_beta_.cpp.in \
		$(DIR)/soft_parameters.hpp.in \
		$(DIR)/soft_parameters.cpp.in \
		$(DIR)/spectrum_generator.hpp.in \
		$(DIR)/spectrum_generator_interface.hpp.in \
		$(DIR)/susy_beta_.cpp.in \
		$(DIR)/susy_parameters.hpp.in \
		$(DIR)/susy_parameters.cpp.in \
		$(DIR)/utilities.hpp.in \
		$(DIR)/utilities.cpp.in

TWO_SCALE_TEMPLATES := \
		$(DIR)/convergence_tester.hpp.in \
		$(DIR)/high_scale_constraint.hpp.in \
		$(DIR)/high_scale_spectrum_generator.hpp.in \
		$(DIR)/high_scale_spectrum_generator.cpp.in \
		$(DIR)/initial_guesser.hpp.in \
		$(DIR)/librarylink.cpp.in \
		$(DIR)/librarylink.m.in \
		$(DIR)/low_scale_constraint.hpp.in \
		$(DIR)/low_scale_spectrum_generator.hpp.in \
		$(DIR)/low_scale_spectrum_generator.cpp.in \
		$(DIR)/run.cpp.in \
		$(DIR)/run.m.in \
		$(DIR)/run_cmd_line.cpp.in \
		$(DIR)/scan.cpp.in \
		$(DIR)/slha_io.hpp.in \
		$(DIR)/slha_io.cpp.in \
		$(DIR)/standard_model_high_scale_spectrum_generator.hpp.in \
		$(DIR)/standard_model_high_scale_spectrum_generator.cpp.in \
		$(DIR)/standard_model_low_scale_spectrum_generator.hpp.in \
		$(DIR)/standard_model_low_scale_spectrum_generator.cpp.in \
		$(DIR)/standard_model_matching.hpp.in \
		$(DIR)/standard_model_matching.cpp.in \
		$(DIR)/standard_model_two_scale_high_scale_initial_guesser.cpp.in \
		$(DIR)/standard_model_two_scale_high_scale_initial_guesser.hpp.in \
		$(DIR)/standard_model_two_scale_low_scale_initial_guesser.cpp.in \
		$(DIR)/standard_model_two_scale_low_scale_initial_guesser.hpp.in \
		$(DIR)/standard_model_two_scale_matching.hpp.in \
		$(DIR)/standard_model_two_scale_matching.cpp.in \
		$(DIR)/susy_scale_constraint.hpp.in \
		$(DIR)/two_scale_convergence_tester.hpp.in \
		$(DIR)/two_scale_convergence_tester.cpp.in \
		$(DIR)/two_scale_high_scale_constraint.hpp.in \
		$(DIR)/two_scale_high_scale_constraint.cpp.in \
		$(DIR)/two_scale_high_scale_initial_guesser.hpp.in \
		$(DIR)/two_scale_high_scale_initial_guesser.cpp.in \
		$(DIR)/two_scale_low_scale_constraint.hpp.in \
		$(DIR)/two_scale_low_scale_constraint.cpp.in \
		$(DIR)/two_scale_low_scale_initial_guesser.hpp.in \
		$(DIR)/two_scale_low_scale_initial_guesser.cpp.in \
		$(DIR)/two_scale_model.hpp.in \
		$(DIR)/two_scale_model.cpp.in \
		$(DIR)/two_scale_susy_scale_constraint.hpp.in \
		$(DIR)/two_scale_susy_scale_constraint.cpp.in

TEMPLATES    := \
		$(BASE_TEMPLATES) \
		$(TWO_SCALE_TEMPLATES)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):
		@true

clean-$(MODNAME):
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
