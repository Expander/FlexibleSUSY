#include <complex>
#include <string>
#include <cstdlib>

#include "consts.hpp"
// #include "mssm_parameter_point.hpp"
#include "fmssm_oneloop.hpp"
#include "fmssm_lattice.hpp"
#include "fmssm_lattice_msusy_constraint.hpp"
#include "fmssm_lattice_mz_constraint.hpp"
#include "fmssmn_lattice.hpp"
#include "fmssmn_lattice_mx_constraint.hpp"
// #include "fmssmn_lattice_convergence_tester.hpp"
// #include "two_scale_running_precision.hpp"
#include "lattice_initial_guesser.hpp"
#include "lattice_solver.hpp"
#include "logger.hpp"
// #include "coupling_monitor.hpp"


using namespace std;

namespace flexiblesusy {

#define fortran_fmssm_fmssmn_mc(name)					\
									\
extern "C" void name##_							\
(const Real& scale0, const Real *w, const Real *x, const int& i,	\
 Real *row, Real *rhs);

fortran_fmssm_fmssmn_mc(fmssm_fmssmn_gauge_couplings)
fortran_fmssm_fmssmn_mc(fmssm_fmssmn_yukawas)
fortran_fmssm_fmssmn_mc(fmssm_fmssmn_mu_b)
fortran_fmssm_fmssmn_mc(fmssm_fmssmn_gaugino_masses)
fortran_fmssm_fmssmn_mc(fmssm_fmssmn_higgs_masses)
fortran_fmssm_fmssmn_mc(fmssm_fmssmn_sfermion_masses)
fortran_fmssm_fmssmn_mc(fmssm_fmssmn_trilinears)

class Fmssm_fmssmn_gauge_couplings : public ForeignMatching {
public:
    Fmssm_fmssmn_gauge_couplings() : ForeignMatching(3) {}
    void operator()() {
	for (size_t i = 0; i < 3; i++) {
	    fmssm_fmssmn_gauge_couplings_
	    (f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_yukawas : public ForeignMatching {
public:
    Fmssm_fmssmn_yukawas() : ForeignMatching(54) {}
    void operator()() {
	for (size_t i = 0; i < 54; i++) {
	    fmssm_fmssmn_yukawas_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_mu_b : public ForeignMatching {
public:
    Fmssm_fmssmn_mu_b() : ForeignMatching(4) {}
    void operator()() {
	for (size_t i = 0; i < 4; i++) {
	    fmssm_fmssmn_mu_b_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_gaugino_masses : public ForeignMatching {
public:
    Fmssm_fmssmn_gaugino_masses() : ForeignMatching(6) {}
    void operator()() {
	for (size_t i = 0; i < 6; i++) {
	    fmssm_fmssmn_gaugino_masses_
	    (f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_higgs_masses : public ForeignMatching {
public:
    Fmssm_fmssmn_higgs_masses() : ForeignMatching(2) {}
    void operator()() {
	for (size_t i = 0; i < 2; i++) {
	    fmssm_fmssmn_higgs_masses_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_sfermion_masses : public ForeignMatching {
public:
    Fmssm_fmssmn_sfermion_masses() : ForeignMatching(45) {}
    void operator()() {
	for (size_t i = 0; i < 45; i++) {
	    fmssm_fmssmn_sfermion_masses_
	    (f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_trilinears : public ForeignMatching {
public:
    Fmssm_fmssmn_trilinears() : ForeignMatching(54) {}
    void operator()() {
	for (size_t i = 0; i < 54; i++) {
	    fmssm_fmssmn_trilinears_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_initial_guesser : public Initial_guesser<Lattice> {
public:
    Fmssm_fmssmn_initial_guesser
    (double mxGuess_, double mu_ini, double b_ini,
     Fmssm_mz_constraint& mzc_,
     Fmssm_msusy_constraint& msc_, double mu_match_,
     Fmssmn_mx_constraint& mxc_) :
	mxGuess(mxGuess_),
	mu(mu_ini),
	b(b_ini),
	mzc(mzc_),
	msc(msc_),
	mu_match(mu_match_),
	mxc(mxc_)
    {}

    virtual void operator()() {
	size_t fmssm_height  = f->efts[0].height;
	size_t fmssmn_height = f->efts[1].height;

	Real gY = mzc.gcs.g1 * sqrt(3/5.0);
	Real g2 = mzc.gcs.g2;
	Comp Yt = mzc.ycs.Yu(2,2); // neglecting CKM rotation
	Comp At = mxc.tfc.Au(2,2);
	Real vu = msc.ewsb.vu;
	Real vd = msc.ewsb.vd;
	mxc.mhc.m2Hu_ini = b*vd/vu - sqr(mu) +
	    (sqr(g2)+sqr(gY))*(Pow4(vd)-Pow4(vu))/(4*(sqr(vd)+sqr(vu)));
	mxc.mhc.m2Hd_ini = b*vu/vd - sqr(mu) +
	    (sqr(g2)+sqr(gY))*(Pow4(vu)-Pow4(vd))/(4*(sqr(vd)+sqr(vu)));

	Real m2Q33 = real(mxc.mfc.m2Q(2,2));
	Real m2U33 = real(mxc.mfc.m2U(2,2));

	Real M2stLL =
	    m2Q33+sqr(vu)*norm(Yt)+(3*sqr(g2)-sqr(gY))*(sqr(vd)-sqr(vu))/12;
	Real M2stRR = m2U33+sqr(vu)*norm(Yt)+sqr(gY)*(sqr(vd)-sqr(vu))/3;
	Comp M2stLR = vu*conj(At*Yt) - vd*conj(Yt)*mu;
	Real detM2st = M2stLL*M2stRR - norm(M2stLR);
	Real MS_guess = pow(detM2st, 1/4.0);

	vector<Real> mu_fmssm = { mZ, MS_guess, mu_match };
	vector<Real> mu_fmssmn = { mu_match, mxGuess };
	assert(fmssm_height == mu_fmssm.size());
	assert(fmssmn_height == mu_fmssmn.size());

	const Fmssm<Lattice>& fmssm =
	    *static_cast<const Fmssm<Lattice>*>(f->efts[0].w);
	const Fmssmn<Lattice>& fmssmn =
	    *static_cast<const Fmssmn<Lattice>*>(f->efts[1].w);

	// when a = 0, no Wilson coefficient runs except for g1, g2, g3

	fmssm(0).m2Hu() = fmssmn(0).m2Hu() = mxc.mhc.m2Hu_ini;
	fmssm(0).m2Hd() = fmssmn(0).m2Hd() = mxc.mhc.m2Hd_ini;
	fmssm(0).mu()   = fmssmn(0).mu()   = mu;
	fmssm(0).b()    = fmssmn(0).b()    = b;
	fmssm(0).M1()   = fmssmn(0).M1()   = mxc.mgc.M1;
	fmssm(0).M2()   = fmssmn(0).M2()   = mxc.mgc.M2;
	fmssm(0).M3()   = fmssmn(0).M3()   = mxc.mgc.M3;
	for (size_t i = 0; i < 3; i++)
	    for (size_t j = 0; j < 3; j++) {
		fmssm(0).Yu(i,j)  = fmssmn(0).Yu(i,j)  = mzc.ycs.Yu(i,j);
		fmssm(0).Yd(i,j)  = fmssmn(0).Yd(i,j)  = mzc.ycs.Yd(i,j);
		fmssm(0).Ye(i,j)  = fmssmn(0).Ye(i,j)  = mzc.ycs.Ye(i,j);
		fmssm(0).m2Q(i,j) = fmssmn(0).m2Q(i,j) = mxc.mfc.m2Q(i,j);
		fmssm(0).m2U(i,j) = fmssmn(0).m2U(i,j) = mxc.mfc.m2U(i,j);
		fmssm(0).m2D(i,j) = fmssmn(0).m2D(i,j) = mxc.mfc.m2D(i,j);
		fmssm(0).m2L(i,j) = fmssmn(0).m2L(i,j) = mxc.mfc.m2L(i,j);
		fmssm(0).m2E(i,j) = fmssmn(0).m2E(i,j) = mxc.mfc.m2E(i,j);
		fmssm(0).TAu(i,j) = fmssmn(0).TAu(i,j) = mxc.tfc.Au(i,j)*mzc.ycs.Yu(i,j);
		fmssm(0).TAd(i,j) = fmssmn(0).TAd(i,j) = mxc.tfc.Ad(i,j)*mzc.ycs.Yd(i,j);
		fmssm(0).TAe(i,j) = fmssmn(0).TAe(i,j) = mxc.tfc.Ae(i,j)*mzc.ycs.Ye(i,j);

		fmssmn(0).Yn(i,j)  = mxc.ync.Yn(i,j);
		fmssmn(0).m2N(i,j) = mxc.mfc.m2N(i,j);
		fmssmn(0).TAn(i,j) = mxc.tfc.An(i,j)*mxc.ync.Yn(i,j);
	    }
	for (size_t m = 1; m < fmssm_height; m++)
	    for (size_t i = 0; i < fmssm.width; i++) f->y(0,m,i) = f->y(0,0,i);
	for (size_t m = 1; m < fmssmn_height; m++)
	    for (size_t i = 0; i < fmssmn.width; i++) f->y(1,m,i) = f->y(1,0,i);

	for (size_t m = 0; m < mu_fmssm.size(); m++) {
	    fmssm(m).t()  =  log(mu_fmssm[m]/f->scl0);
	    fmssm(m).g1() = g1L1(mu_fmssm[m]);
	    fmssm(m).g2() = g1L2(mu_fmssm[m]);
	    fmssm(m).g3() = g1L3(mu_fmssm[m]);
	}
	for (size_t m = 0; m < mu_fmssmn.size(); m++) {
	    fmssmn(m).t()  =  log(mu_fmssmn[m]/f->scl0);
	    fmssmn(m).g1() = g1L1(mu_fmssmn[m]);
	    fmssmn(m).g2() = g1L2(mu_fmssmn[m]);
	    fmssmn(m).g3() = g1L3(mu_fmssmn[m]);
	}
    }

private:
    // const QedQcd oneset;       ///< low-energy parameters
    double mxGuess;            ///< guessed GUT scale
    Real mu;
    Real b;
    Fmssm_mz_constraint& mzc;
    Fmssm_msusy_constraint& msc;
    Real mu_match;
    Fmssmn_mx_constraint& mxc;
};

class Fmssmn_cmssmn_constraint : public Fmssmn_mx_constraint {
public:
    Fmssmn_cmssmn_constraint(double m0, double m12, double a0, const CM33& Yn){
	mgc.M1 = mgc.M2 = mgc.M3 = m12;
	mfc.m2Q << sqr(m0), 0, 0,
	    	   0, sqr(m0), 0,
	    	   0, 0, sqr(m0);
	mfc.m2U = mfc.m2D = mfc.m2L = mfc.m2N = mfc.m2E = mfc.m2Q;
	mhc.m2Hu_fin = mhc.m2Hd_fin = sqr(m0);
	// to be set by initial guesser
	// mhc.m2Hu_ini = m2Hu_ini;
	// mhc.m2Hd_ini = m2Hd_ini;
	tfc.Au << a0, a0, a0,
	    	  a0, a0, a0,
	    	  a0, a0, a0;
	tfc.Ad = tfc.An = tfc.Ae = tfc.Au;
	ync.Yn = Yn;
    }
};

}

using namespace flexiblesusy;

int main(int argc, char *argv[])
{
   // Mssm_parameter_point pp;
   Real pp_m0  = 400*GeV;
   Real pp_m12 = 300*GeV;
   Real pp_a0  = pp_m0;
   Real pp_tanBeta = 10;
   Real pp_mxGuess = MX1L();
   Real mu = pp_m0;
   Real b  = sqr(pp_m0);

   Fmssm<Lattice> fmssm;

   Fmssm_mz_constraint fmssm_mz_constraint(pp_tanBeta);
   Fmssm_msusy_constraint fmssm_msusy_constraint(pp_tanBeta);
   // Fmssm_convergence_tester fmssm_convergence_tester(1.0e-4);

   std::vector<Constraint<Lattice>*> fmssm_constraints;
   fmssm_constraints.push_back(&fmssm_mz_constraint);
   fmssm_constraints.push_back(&fmssm_msusy_constraint);

   Fmssmn<Lattice> fmssmn;

   Real matching_scale = 1e12*GeV;
   Fixed_t set_scale(matching_scale);
   SingleSiteInterTheoryConstraint set_matching_scale(&set_scale, 1);
   Match_t match_t;
   Fmssm_fmssmn_gauge_couplings fmssm_fmssmn_gauge_couplings;
   Fmssm_fmssmn_yukawas fmssm_fmssmn_yukawas;
   Fmssm_fmssmn_mu_b fmssm_fmssmn_mu_b;
   Fmssm_fmssmn_gaugino_masses fmssm_fmssmn_gaugino_masses;
   Fmssm_fmssmn_higgs_masses fmssm_fmssmn_higgs_masses;
   Fmssm_fmssmn_sfermion_masses fmssm_fmssmn_sfermion_masses;
   Fmssm_fmssmn_trilinears fmssm_fmssmn_trilinears;
   CompoundMatching<Lattice> fmssm_fmssmn_matching({
	   &set_matching_scale,
	   &match_t,
	   &fmssm_fmssmn_gauge_couplings,
	   &fmssm_fmssmn_yukawas,
	   &fmssm_fmssmn_mu_b,
	   &fmssm_fmssmn_gaugino_masses,
	   &fmssm_fmssmn_higgs_masses,
	   &fmssm_fmssmn_sfermion_masses,
	   &fmssm_fmssmn_trilinears
       });

   CM33 UPMNS;			// tribimaximal mixing
   UPMNS <<  sqrt(2.0/3),  sqrt(1.0/3),	    0.0   ,
       	    -sqrt(1.0/6),  sqrt(1.0/3), sqrt(1.0/2),
       	     sqrt(1.0/6), -sqrt(1.0/3), sqrt(1.0/2);
   CM33 YnDiag;
   YnDiag << 0.2, 0, 0,
       	     0, 0.6, 0,
       	     0, 0, 1.2;
   CM33 Yn;
   Yn = UPMNS * YnDiag;
   Fmssmn_cmssmn_constraint fmssmn_cmssmn_constraint
       (pp_m0, pp_m12, pp_a0, Yn);
   std::vector<Constraint<Lattice>*> fmssmn_constraints;
   fmssmn_constraints.push_back(&fmssmn_cmssmn_constraint);

   // LATTICE: guess of initial profile requires knowledge of entire
   // EFT tower
   Fmssm_fmssmn_initial_guesser initial_guesser(pp_mxGuess, mu, b,
						fmssm_mz_constraint,
						fmssm_msusy_constraint,
						matching_scale,
						fmssmn_cmssmn_constraint);
   // Two_scale_increasing_precision two_scale_increasing_precision(10.0, 1.0e-5);

   RGFlow<Lattice> solver;
   // Q: does each EFT have an associated convergence tester?
   // solver.set_convergence_tester(&mssm_convergence_tester);
   // solver.set_running_precision(&two_scale_increasing_precision);
   solver.set_initial_guesser(&initial_guesser);
   // Q: is fmssm_constraints allowed to be modified after add_model()ed?
   solver.add_model(&fmssm, &fmssm_fmssmn_matching, fmssm_constraints);
   solver.add_model(&fmssmn, fmssmn_constraints);

   if (argc > 1) {
       string argv1(argv[1]);
       if (argv1.find("h") != string::npos) solver.enable_hybrid();
       if (argv1.find("1") != string::npos) solver.disable_multithreading();
   }

   // INFO("Running: " << pp);
   try {
      solver.solve();
   } catch (Error& e) {
      ERROR("no solution found: " << e.what());
      exit(EXIT_FAILURE);
   }

   fmssm.calculate_spectrum();

   INFO("Solution found:\n" << solver);
   // fmssm.print(std::cout);

   stringstream m2Qprint, m2Lprint;
   for (size_t i = 0; i < 3; i++) {
       for (size_t j = 0; j < 3; j++) {
	   m2Qprint << (j ? ", " : "") << Comp(fmssm(0).m2Q(i,j))/sqr(TeV);
	   m2Lprint << (j ? ", " : "") << Comp(fmssm(0).m2L(i,j))/sqr(TeV);
       }
       m2Qprint << '\n';
       m2Lprint << '\n';
   }
   INFO("At scale = " << exp(fmssm(0).t())*solver.scl0/TeV << " TeV:\n\n"
	"mu = " << Comp(fmssm(0).mu())/TeV << " TeV\n\n"
	"m2Q/TeV^2 =\n" << m2Qprint.str() << "\n"
	"m2L/TeV^2 =\n" << m2Lprint.str());

   return 0;
}
