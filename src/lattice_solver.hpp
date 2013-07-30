#ifndef lattice_solver_hpp
#define lattice_solver_hpp


#include <functional>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include "mathdefs.hpp"

#include "rg_flow.hpp"
#include "error.hpp"

namespace flexiblesusy {

#if 1
template <class T> class Constraint;
template <class T> class Matching;
template <class T> class Convergence_tester;
template <class T> class Initial_guesser;
#endif

class Lattice;
class Lattice_constraint;
// class SingleSiteConstraint;
// class InterTheoryConstraint;
typedef Constraint<Lattice> SingleSiteConstraint;
typedef Matching<Lattice> InterTheoryConstraint;


typedef std::vector<Real> RVec;


template<class T>
class band_matrix {
public:
    band_matrix(size_t N, size_t KL, size_t KU, size_t LD) :
	n(N), kl(KL), ku(KU), ld(LD) {
	assert(ld >= kl+ku+1);
	AB = new T[ld*n];
    }
    ~band_matrix() { delete[] AB; }
    T *pointer() { return AB; }
    void clear() { for (size_t i = 0; i < ld*n; i++) AB[i] = 0; }
    T& operator()(size_t i, size_t j) {
	assert(i < n && j < n && i+ku >= j && i <= j+kl);
	size_t idx = kl+ku+i-j + ld*j;
	assert(idx < ld*n);
	return AB[idx];
    }

private:
    size_t n, kl, ku, ld;
    T *AB;
};

			      
struct Wilson {
    Wilson(size_t nvars) : width(nvars) {}
    virtual ~Wilson() {}
    size_t width;		// 1 + number of Wilson coefficients
    // derivative wrt t == x[0]
    virtual Real  dx(const Real *x, size_t i) const = 0;
    // d dx[i] / d x[j] == d dx[i] / d y[j] / unit[j]
    virtual void ddx(const Real *x, size_t i, Real *ddx) const = 0;
};

struct ParWilson {
    ParWilson(size_t width_) : width(width_) {}
    virtual ~ParWilson() {}
    size_t width;		// 1 + number of Wilson coefficients
    // derivative wrt t == x[0]
    virtual Real  dx(const Real a, const Real *x, size_t i) const = 0;
    // d dx[i] / d x[j] == d dx[i] / d y[j] / unit[j]
    virtual void ddx(const Real a, const Real *x, size_t i, Real *ddx) const=0;
};

class Lattice_model: public ParWilson {
public:
    Lattice_model(size_t width) : ParWilson(width) {}
    virtual ~Lattice_model() {}
    virtual void init(RGFlow<Lattice> *flow, size_t theory)
    { f = flow; T = theory; }
    virtual void calculate_spectrum() = 0;
    virtual std::string name() const { return "unnamed"; }
    // virtual int run_to(double, double eps = -1.0) = 0;
    virtual void print(std::ostream& out) const { out << "Model: " << name(); }
    friend std::ostream& operator<<(std::ostream& out, const Lattice_model& model) {
	model.print(out);
	return out;
    }

protected:
    RGFlow<Lattice> *f;
    size_t T;
};

template<>
class RGFlow<Lattice> {
public:
    enum Inner_status { JUMPED=-2, ENDLESS=-1, CONVERGED=0 };

    class MemoryError : public Error {
    public:
	MemoryError(const std::string& message_) : message(message_) {}
	virtual ~MemoryError() {}
	virtual std::string what() const { return message; }
    private:
	std::string message;
    };

    class SetupError : public Error {
    public:
	SetupError(const std::string& message_) : message(message_) {}
	virtual ~SetupError() {}
	virtual std::string what() const { return message; }
    private:
	std::string message;
    };

    class NonInvertibleMatrixError : public Error {
    public:
	NonInvertibleMatrixError(const std::string& message_) :
	    message(message_) {}
	virtual ~NonInvertibleMatrixError() {}
	virtual std::string what() const { return message; }
    private:
	std::string message;
    };

    class DivergenceError : public Error {
    public:
	DivergenceError(const std::string& message_) : message(message_) {}
	virtual ~DivergenceError() {}
	virtual std::string what() const { return message; }
    private:
	std::string message;
    };

    class NoConvergenceError : public Error {
    public:
	NoConvergenceError(size_t number_of_iterations_)
	    : number_of_iterations(number_of_iterations_) {}
	virtual ~NoConvergenceError() {}
	virtual std::string what() const {
	    std::stringstream message;
	    message << "RGFlow<Lattice>::NoConvergenceError: no convergence"
		    << " after " << number_of_iterations << " iterations";
	    return message.str();
	}
    private:
	size_t number_of_iterations;
    };

    class NonPerturbativeRunningError : public Error {
    public:
	NonPerturbativeRunningError(Lattice_model* model_, double scale_)
	    : model(model_)
	    , scale(scale_)
	    {}
	virtual ~NonPerturbativeRunningError() {}
	virtual std::string what() const;
    private:
	Lattice_model* model;
	double scale;
    };

    struct Translator {
	template<class T>
	struct Var {
	    Var(std::function<T()> get, std::function<void(T)> set) :
		get_(get), set_(set) {}
	    operator T() { return get_(); }
	    T operator=(T value) { set_(value); return get_(); }
	    T operator=(Var<T>& var) { return operator=(T(var)); }
	    std::function<T()> get_;
	    std::function<void(T)> set_;
	};

	Translator(RGFlow *flow, size_t theory, size_t site) :
	    f(flow), T(theory), m(site) {}
	Real  u(size_t i) { return f->efts[T].units[i]; }
	Real& y(size_t i) { return f->y(T, m, i); }
	RGFlow *f;
	size_t T;
	size_t m;
    };

    struct EFT {
	EFT(Lattice_model *model, RGFlow *flow) :
	    w(model), units(w->width, 1), f(flow)
	    {}
	Lattice_model *w;
	RVec units;		// x normalizations
	// std::vector<ptrdiff_t> r_offset;
	RGFlow *f;
	size_t height;
	size_t T;
	size_t offset;		// location within y_
    };

    union EqRow {
	EqRow *next;		// next free element
	struct {
	    size_t T, m;	// first site
	    size_t n;		// number of occupied sites
	    size_t realRow;	// row within A_
	} rowSpec;
    };

    // RGFlow(std::vector<EFT *> efts, std::vector<Constraint *> constraints,
    // 	   Real scale0, int verbose, std::ostream *lout);
    RGFlow();
    ~RGFlow();

    /// add a model and constraints
    void add_model(Lattice_model* model,
		   const std::vector<SingleSiteConstraint*>& constraints);
    /// add a model, constraints and matching condition
    void add_model(Lattice_model*,
		   InterTheoryConstraint *m = NULL,
		   const std::vector<SingleSiteConstraint*>& constraints = std::vector<SingleSiteConstraint*>());
    /// add a model and up- and downwards constraints
    void add_model(Lattice_model*,
		   const std::vector<SingleSiteConstraint*>&,
		   const std::vector<SingleSiteConstraint*>&);
    /// add a model, up- and downward constraints and matching condition
    void add_model(Lattice_model*,
		   InterTheoryConstraint *m,
		   const std::vector<SingleSiteConstraint*>& upwards_constraints,
		   const std::vector<SingleSiteConstraint*>& downwards_constraints);
    void set_initial_guesser(Initial_guesser<Lattice>*);

    void enable_hybrid() { hybrid = true; }
    void disable_multithreading() { multithreading = false; };
    void solve();
    friend std::ostream& operator<<(std::ostream &out, const RGFlow& flow);

    struct EFTspec : public EFT {
	EFTspec(Lattice_model *model,
		const std::vector<SingleSiteConstraint*>& cs,
		InterTheoryConstraint* m,
		RGFlow *flow);

	// Lattice_model *model;
	std::vector<SingleSiteConstraint*> constraints;
	InterTheoryConstraint *matching;
    };
    std::vector<EFTspec> efts;
    Initial_guesser<Lattice> *init_profile;
    Real tiny_dy;
    Real huge_dy;
    Real a;			// continuation parameter btw 0 & 1
    size_t max_a_steps;
    size_t max_iter;
    bool units_set;
    bool hybrid;
    Real scl0;			// t == log(mu/scl0)
    RVec y_;			// entire RG flow
    band_matrix<Real> *A_;
    RVec z;			// RHS of linear equations
    std::vector<EqRow> row_pool;
    EqRow *free_row_list_head;
    int N, KL, KU, LDA;		// for lapack
    std::vector<int> IPIV;	// ditto
    int verb;
    // std::ostream *log;
    size_t site_offset(size_t T, size_t m) const {
	assert(m < efts[T].height);
	return efts[T].offset + m*efts[T].w->width;
    }
    Real& A(size_t r, size_t T, size_t m, size_t i)
    { return (*A_)(r, site_offset(T,m) + i); }
    Real& y(size_t T, size_t m, size_t i) { return y_[site_offset(T,m) + i]; }
    Real  y(size_t T, size_t m, size_t i) const
    { return y_[site_offset(T,m) + i]; }
    Real  x(size_t T, size_t m, size_t i) const
    { return y(T,m,i) * efts[T].units[i]; }
    std::vector<Lattice_constraint*> constraints;
    std::vector<size_t> teqidx;
    std::vector<size_t> rgeidx;
    std::unordered_set<Lattice_constraint*> elementary_constraints;
    bool multithreading;
    // as of v1.52 boost::thread gives valgrind the impression that
    // 8 bytes are not freed at exit
    // see http://www.cplusplus.com/forum/unices/83480/
    boost::thread_group *threads;
    boost::barrier *threads_begin, *threads_end;
    bool keep_threads;
    void set_units();
    void apply_constraints();
    void apply_constraints_thread(Lattice_constraint *c);
    void create_threads();
    void join_threads();
    Real maxdiff(const RVec& y0, const RVec& y1);
    void init_lattice();
    void increase_a();
    void increase_density();
    void rk_stage();
    Inner_status iterate();
    std::vector<std::vector<size_t>> refine_lattice();
    void enable_Runge_Kutta();
    void disable_Runge_Kutta();
    void resample(const std::vector<std::vector<size_t>>& site_maps);
    EqRow *ralloc(size_t T, size_t m, size_t span);
    void rfree(EqRow *r);
    void init_free_row_list();
    void sort_rows();
};

}

#endif // lattice_solver_hpp
