// special case for Higgs -> Z Z
// TODO: implement higher order corrections

struct my_f_params {double mHOS; double mVOS; double GammaV;};

double hVV_4body(double *q2, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  struct my_f_params * fp = (struct my_f_params *)params;
  const double mHOS = fp->mHOS;
  if (q2[1] > Sqr(mHOS - std::sqrt(q2[0]))) return 0.;
  const double mVOS = fp->mVOS;
  const double GammaV = fp->GammaV;
  const double kl = KallenLambda(1., q2[0]/Sqr(mHOS), q2[1]/Sqr(mHOS));
  return
     mVOS*GammaV/(Sqr(q2[0] - Sqr(mVOS)) + Sqr(mVOS*GammaV))
     * mVOS*GammaV/(Sqr(q2[1] - Sqr(mVOS)) + Sqr(mVOS*GammaV))
     * std::sqrt(kl)*(kl + 12.*q2[0]*q2[1]/Power4(mHOS));
}

template <>
double CLASSNAME::get_partial_width<H,Z,Z>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<Z>::type const& indexOut1,
   typename field_indices<Z>::type const& indexOut2
   ) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   // There might be large differences between mZ from mass block
   // and one from slha input, especially in the decoupling limit
   // so we use the latter one. There might be a problem with
   // models where Z mixes with something else.
   // const double mZOS = context.physical_mass<Z>(indexOut1);
   const double mZOS = qedqcd.displayPoleMZ();
   const double x = Sqr(mZOS/mHOS);
   double res = 0;

   // mH < mZ
   if ((x > 1.0 && flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays) != 0) ||
       (4.*x > 1.0 && flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays) == 2)) {

      // integrand
      constexpr double GammaZ = 2.4952;
      struct my_f_params params = {mHOS, mZOS, GammaZ};
      gsl_monte_function G = {&hVV_4body, 2, &params};

      // setup integration
      gsl_rng_env_setup ();
      double xl[2] = {0, 0};
      double xu[2] = {Sqr(mHOS), Sqr(mHOS)};
      // this gives relative error < 0.05%
      constexpr size_t calls = 1'000'000;
      double err;
      const gsl_rng_type *T = gsl_rng_default;
      gsl_rng *r = gsl_rng_alloc (T);
      gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
      gsl_monte_miser_integrate (&G, xl, xu, 2, calls, r, s,
                                 &res, &err);

      // clean-up
      gsl_monte_miser_free (s);
      gsl_rng_free (r);

      // prefactor
      const auto indices = concatenate(indexIn, indexOut2, indexOut1);
      const auto ghZZ =
         Vertex<H, Z, Z>::evaluate(indices, context).value();
      const double normalization = 1./(64.*Cube(Pi)) * std::norm(ghZZ) * Cube(mHOS)/Power4(mZOS)/2.;
      res *= normalization;
      err *= normalization;
      const double rel_err = err/res;
      if (rel_err > 1e-3) {
         std::ostringstream oss;
         oss <<
            "Warning: result of H->ZZ is "
            << res << " +/- " << err << " GeV (" << 1e+2*rel_err
            << "%). Relative error > 0.1%.\n";
         WARNING (oss.str());
      }
   // mZ < mH < 2*mZ
   // three-body decay
   }
   else if (4 * x > 1.0 && flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays) != 0) {

      if (check_3body_Vff_decay<BSMForZdecay,Z>(context, mHOS, indexOut1)) {
         const std::string index_as_string = indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn.at(0)) + ")";
         WARNING("Warning in H" + index_as_string + "->ZZ decays: Single off-shell decays H->Zff' assume no possible BSM particles in the final state. Turning off.");
         return 0.;
      }

      const double sw2 = Sqr(std::sin(context.model.ThetaW()));
      const double deltaV = 7.0/12.0 - 10.0/9.0*sw2 + 40.0/27.0*Sqr(sw2);

      res = 3./(512.*Power3(Pi)) * 1./mHOS * deltaV * RT(x)/x;

      const auto indices = concatenate(indexIn, indexOut2, indexOut1);
      const auto ghZZ =
         Vertex<H, Z, Z>::evaluate(indices, context).value();

      const double g2 = context.model.get_g2();

      res *= std::norm(ghZZ*g2)/(1-sw2);
   // mH > 2mZ
   // two-body decay
   }
   else if (4.*x < 1.0) {

      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * std::sqrt(KallenLambda(1., x, x));

      // phase space symmetry factor
      const double ps_symmetry = 1./2.;

      // matrix element squared
      const auto mat_elem = calculate_amplitude<H, Z, Z>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      res = flux * ps * ps_symmetry * mat_elem_sq;
   }

   return res;
}
