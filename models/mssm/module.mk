DIR          := models/mssm
MODNAME      := libmssm

LIBMSSM_HDR  := \
		$(DIR)/mssm_solver.h \
		$(DIR)/physpars.h \
		$(DIR)/softpars.h \
		$(DIR)/softsusy.h \
		$(DIR)/susy.h \
		$(DIR)/twoloophiggs.h

LIBMSSM_SRC  := \
		$(DIR)/physpars.cpp \
		$(DIR)/softpars.cpp \
		$(DIR)/softsusy.cpp \
		$(DIR)/susy.cpp \
		$(DIR)/twoloophiggs.f

LIBMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSM_SRC)))

LIBMSSM_DEP  := \
		$(LIBMSSM_OBJ:.o=.d)

LIBMSSM      := $(DIR)/$(MODNAME).a

SOFTSUSY_SRC := \
		$(DIR)/main.cpp

SOFTSUSY_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(SOFTSUSY_SRC)))

SOFTSUSY_DEP := \
		$(SOFTSUSY_OBJ:.o=.d)

SOFTSUSY     := \
		$(DIR)/softsusy.x

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBMSSM) $(SOFTSUSY)

clean-$(MODNAME):
		rm -rf $(LIBMSSM_OBJ)
		rm -rf $(SOFTSUSY_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBMSSM_DEP)
		rm -rf $(LIBMSSM)
		rm -rf $(SOFTSUSY_DEP)
		rm -rf $(SOFTSUSY)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBMSSM): $(LIBMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBMSSM_DEP)
ALLLIB += $(LIBMSSM)

$(SOFTSUSY): $(SOFTSUSY_OBJ) $(ALLLIB)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

ALLDEP += $(SOFTSUSY_DEP)
ALLLIB += $(SOFTSUSY)
