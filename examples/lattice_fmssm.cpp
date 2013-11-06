#include <complex>
#include <string>
#include <cstdlib>

#include "consts.hpp"
// #include "mssm_parameter_point.hpp"
#include "fmssm_oneloop.hpp"
#include "fmssm_lattice.hpp"
#include "lattice_initial_guesser.hpp"
#include "fmssm_lattice_mx_constraint.hpp"
#include "fmssm_lattice_msusy_constraint.hpp"
#include "fmssm_lattice_mz_constraint.hpp"
// #include "fmssm_lattice_convergence_tester.hpp"
// #include "two_scale_running_precision.hpp"
#include "lattice_solver.hpp"
#include "logger.hpp"
// #include "coupling_monitor.hpp"


using namespace std;

namespace flexiblesusy {

class Fmssm_cmssm_constraint : public Fmssm_mx_constraint {
public:
    Fmssm_cmssm_constraint
    (double m0, double m12, double a0) {
	mgc.M1 = mgc.M2 = mgc.M3 = m12;
	mfc.m2Q << sqr(m0), 0, 0,
	    	   0, sqr(m0), 0,
	    	   0, 0, sqr(m0);
	mfc.m2U = mfc.m2D = mfc.m2L = mfc.m2E = mfc.m2Q;
	mhc.m2Hu_fin = mhc.m2Hd_fin = sqr(m0);
	// to be set by initial guesser
	// mhc.m2Hu_ini = m2Hu_ini;
	// mhc.m2Hd_ini = m2Hd_ini;
	tfc.Au << a0, a0, a0,
	    	  a0, a0, a0,
	    	  a0, a0, a0;
	tfc.Ad = tfc.Ae = tfc.Au;
    }
};

class Fmssm_initial_guesser : public Initial_guesser<Lattice> {
public:
    Fmssm_initial_guesser
    (double mxGuess_, double mu_ini, double b_ini,
     Fmssm_mz_constraint& mzc_,
     Fmssm_msusy_constraint& msc_,
     Fmssm_mx_constraint& mxc_) :
	mxGuess(mxGuess_),
	mu(mu_ini),
	b(b_ini),
	mzc(mzc_),
	msc(msc_),
	mxc(mxc_)
    {}

    void operator()() {
	size_t fmssm_height = f->efts[0].height;

	Real gY = mzc.gcs.g1 * sqrt(3/5.0);
	Real g2 = mzc.gcs.g2;
	Comp Yt = mzc.ycs.Yu(2,2); // neglecting CKM rotation
	Comp At = mxc.tfc.Au(2,2);
	Real vu = msc.ewsb.vu;
	Real vd = msc.ewsb.vd;
	// invent Higgs soft masses meeting EWSB conditions for a=0
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

	Real t0 = log(mZ/f->scl0);
	Real t1 = log(MS_guess/f->scl0);
	Real t2 = log(mxGuess/f->scl0);

	const Fmssm<Lattice>& fmssm =
	    *static_cast<const Fmssm<Lattice>*>(f->efts[0].w);

	// when a = 0, no Wilson coefficient runs except for g1, g2, g3

	fmssm(0).m2Hu() = mxc.mhc.m2Hu_ini;
	fmssm(0).m2Hd() = mxc.mhc.m2Hd_ini;
	fmssm(0).mu()   = mu;
	fmssm(0).b()    = b;
	fmssm(0).M1()   = mxc.mgc.M1;
	fmssm(0).M2()   = mxc.mgc.M2;
	fmssm(0).M3()   = mxc.mgc.M3;
	for (size_t i = 0; i < 3; i++)
	    for (size_t j = 0; j < 3; j++) {
		fmssm(0).Yu(i,j)  = mzc.ycs.Yu(i,j);
		fmssm(0).Yd(i,j)  = mzc.ycs.Yd(i,j);
		fmssm(0).Ye(i,j)  = mzc.ycs.Ye(i,j);
		fmssm(0).m2Q(i,j) = mxc.mfc.m2Q(i,j);
		fmssm(0).m2U(i,j) = mxc.mfc.m2U(i,j);
		fmssm(0).m2D(i,j) = mxc.mfc.m2D(i,j);
		fmssm(0).m2L(i,j) = mxc.mfc.m2L(i,j);
		fmssm(0).m2E(i,j) = mxc.mfc.m2E(i,j);
		fmssm(0).TAu(i,j) = mxc.tfc.Au(i,j)*mzc.ycs.Yu(i,j);
		fmssm(0).TAd(i,j) = mxc.tfc.Ad(i,j)*mzc.ycs.Yd(i,j);
		fmssm(0).TAe(i,j) = mxc.tfc.Ae(i,j)*mzc.ycs.Ye(i,j);
	    }
	for (size_t m = 1; m < fmssm_height; m++)
	    for (size_t i = 0; i < fmssm.width; i++) f->y(0,m,i) = f->y(0,0,i);

	fmssm(0).t()  = t0;
	fmssm(0).g1() = g1L1(exp(t0)*f->scl0);
	fmssm(0).g2() = g1L2(exp(t0)*f->scl0);
	fmssm(0).g3() = g1L3(exp(t0)*f->scl0);
	fmssm(1).t()  = t1;
	fmssm(1).g1() = g1L1(exp(t1)*f->scl0);
	fmssm(1).g2() = g1L2(exp(t1)*f->scl0);
	fmssm(1).g3() = g1L3(exp(t1)*f->scl0);
	fmssm(2).t()  = t2;
	fmssm(2).g1() = g1L1(exp(t2)*f->scl0);
	fmssm(2).g2() = g1L2(exp(t2)*f->scl0);
	fmssm(2).g3() = g1L3(exp(t2)*f->scl0);
    }

private:
    // const QedQcd oneset;       ///< low-energy parameters
    double mxGuess;            ///< guessed GUT scale
    Real mu;
    Real b;
    Fmssm_mz_constraint& mzc;
    Fmssm_msusy_constraint& msc;
    Fmssm_mx_constraint& mxc;
};

}

using namespace flexiblesusy;

int main(int argc, char *argv[])
{
   // Mssm_parameter_point pp;
   // FIXME: solver does not work if GeV = 1 instead of 1e-3,
   // this has sth to do with guessing units of vanishing Wilson coefficients
   Real pp_m0  = 400*GeV;
   Real pp_m12 = 300*GeV;
   // METACODE: Fmssm_constraint_on_trilinears divides by 0
   // when a0 = 0, this can be avoided by rewriting the equation
   // without divisions
   // UPDATE: supposed to be fixed in writeRGE.m
   Real pp_a0  = pp_m0;
   Real pp_tanBeta = 10;
   Real pp_mxGuess = MX1L();
   Real mu = pp_m0;		// sign wish enters here
   // one may want to check sign(mu) after solution is found
   Real b  = sqr(pp_m0);

   Fmssm<Lattice> fmssm;

   Fmssm_mz_constraint fmssm_mz_constraint(pp_tanBeta);
   Fmssm_msusy_constraint fmssm_msusy_constraint(pp_tanBeta);
   Fmssm_cmssm_constraint fmssm_cmssm_constraint(pp_m0, pp_m12, pp_a0);
   // Fmssm_convergence_tester fmssm_convergence_tester(1.0e-4);

   // LATTICE: guess of initial profile requires knowledge of entire
   // EFT tower
   Fmssm_initial_guesser initial_guesser(pp_mxGuess, mu, b,
					 fmssm_mz_constraint,
					 fmssm_msusy_constraint,
					 fmssm_cmssm_constraint);
   // Two_scale_increasing_precision two_scale_increasing_precision(10.0, 1.0e-5);

   // LATTICE: there is no distinction between upward and downward
   // constraints
   // std::vector<Constraint<Lattice>*> fmssm_upward_constraints;
   // fmssm_upward_constraints.push_back(&fmssm_mz_constraint);
   // fmssm_upward_constraints.push_back(&fmssm_sugra_constraint);

   std::vector<Constraint<Lattice>*> fmssm_constraints;
   // lowest scale first
   fmssm_constraints.push_back(&fmssm_mz_constraint);
   fmssm_constraints.push_back(&fmssm_msusy_constraint);
   fmssm_constraints.push_back(&fmssm_cmssm_constraint);

   RGFlow<Lattice> solver;
   // Q: does each EFT have an associated convergence tester?
   // solver.set_convergence_tester(&mssm_convergence_tester);
   // solver.set_running_precision(&two_scale_increasing_precision);
   solver.set_initial_guesser(&initial_guesser);
   // Q: is fmssm_constraints allowed to be modified after add_model()ed?
   solver.add_model(&fmssm, fmssm_constraints);

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
