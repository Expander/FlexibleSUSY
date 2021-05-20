template <>
double CLASSNAME::get_partial_width<H, G, G>(
      const context_base& context,
      const typename field_indices<H>::type& in_idx,
      const typename field_indices<G>::type& out1_idx,
      const typename field_indices<G>::type& out2_idx) const
{
   const auto amp = calculate_amplitude<H, G, G>(context, in_idx, out1_idx, out2_idx);
   const double mH = context.physical_mass<H>(in_idx);
   constexpr double ps {1./(8.*Pi)};
   constexpr double ps_symmetry {1./2.};
   constexpr double color_fact = squared_color_generator<H, G, G>();
   const double flux = 0.5/mH;

   // full LO (SM+BSM) result
   double result = flux * color_fact * ps * ps_symmetry * amp.square();

   const double mt {context.mass<uq>({2})};
   const double tau = Sqr(mH/(2.*mt));
   // the analytic form o corrections is valid for small tau
   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) && tau < 0.7) {
      auto qedqcd_ = qedqcd;
      qedqcd_.to(mH);
      // 5-flavour SM alpha_s
      const double alpha_s_5f = qedqcd_.displayAlpha(softsusy::ALPHAS);

      const auto indices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, in_idx);
      const auto HGGVertex = Vertex<bar<uq>::type, uq, H>::evaluate(indices, context);
      std::complex<double> const HGGVertexSVal = 0.5*(HGGVertex.left() + HGGVertex.right());

      // eq. 2.46 of 0503172
      const std::complex<double> A12_H = 2.*(tau + (tau -1)*f(tau))/Sqr(tau);

      // LO width comming only from the top-loop
      // agrees up to a full double precision with automatically generated one
      // when compensated for different alpha_s
      // heavy top limit (tau -> 0)
      //    3./4.*A12 -> 1
      //    HGGVertexSVal*std::sqrt(tau) is mt independent
      const double Gamma_SM_LO_S =
         mH/(18.*Cube(Pi))*std::norm(alpha_s_5f * HGGVertexSVal*std::sqrt(tau) * 3./4.*A12_H);

      const double mu = mH;
      const double LH = std::log(Sqr(mu/mH));
      const double Lt = std::log(Sqr(mu/mt));

      // eq. 5 of 0708.0916
      constexpr int Nf = 5;
      const double hnlo0 =
         95./4. - 7./6.*Nf + (33 - 2*Nf)/6.*LH;
      const double hnlo1 =
         5803./540. + 77./30.*LH - 14./15.*Lt + Nf*(-29./60. - 7./45.*LH);
      const double hnlo2 =
         1029839./189000. + 16973./12600.*LH - 1543./1575.*Lt + Nf*(-89533./378000. - 1543./18900.*LH);
      const double hnlo3 =
         9075763./2976750. + 1243./1575*LH - 452./525.*Lt + Nf*(-3763./28350. - 226./4725.*LH);
      const double hnlo4 =
         50854463./27783000. + 27677./55125.*LH - 442832./606375.*Lt + Nf*(-10426231./127338750. - 55354./1819125.*LH);
      const double hnlo5 =
         252432553361./218513295000. + 730612./2149875.*LH - 2922448./4729725.*Lt + Nf*(-403722799./7449316875. - 1461224./70945875.*LH);
      const double deltaNLO_S {
         hnlo0 + tau*(hnlo1 + tau*(hnlo2 + tau*(hnlo3 + tau*(hnlo4 + tau*hnlo5))))
      };

      // eq. 9 of 0708.0916
      const double hnnlo0 =
         149533./288. - 121./16.*Sqr(Pi) - 495./8.*zeta3 + 3301./16.*LH + 363./16.*Sqr(LH) + 19./8.*Lt
         + Nf*(-4157./72. + 11./12.*Sqr(Pi) + 5./4.*zeta3 - 95./4.*LH - 11./4.*Sqr(LH) + 2./3.*Lt)
         + Sqr(Nf)*(127./108. - Sqr(Pi)/36. + 7./12.*LH + Sqr(LH)/12);
      const double hnnlo1 =
         104358341./1555200. - 847./240*Sqr(Pi) + 7560817./69120.*zeta3 + LH*(203257./2160. - 77./15.*Lt) + 847./80.*Sqr(LH)
         - 24751./1080.*Lt - 77./180.*Sqr(Lt) + Nf*(-9124273./388800. + 77./180.*Sqr(Pi) + 7./12.*zeta3 + LH*(-67717./6480.
                  + 14./45*Lt) - 77./60.*Sqr(LH) + 586./405.*Lt + 7./90.*Sqr(Lt)) + Sqr(Nf)*(5597./12960. - 7./540.*Sqr(Pi) + 29./120.*LH
                  + 7./180.*Sqr(LH));
      const double hnnlo2 =
         -1279790053883./12192768000. - 186703./100800.*Sqr(Pi) + 39540255113./232243200.*zeta3 + LH*(9158957./189000.
               - 16973./3150.*Lt) + 186703./33600.*Sqr(LH) - 10980293./453600.*Lt + 20059./37800.*Sqr(Lt) + Nf*(-64661429393./5715360000.
                  - 16973./25200.*Sqr(LH) + 16973./75600.*Sqr(Pi) + 1543./5040.*zeta3 + LH*(- 10306537./1944000. + 1543./4725.*Lt)
                  +8973773./6804000.*Lt + 1543./18900.*Sqr(Lt)) + Sqr(Nf)*(3829289./19440000. - 1543./226800.*Sqr(Pi) + 89533./756000.*LH + 1543./75600.*Sqr(LH));

      const double deltaNNLO_S {
         hnnlo0 + tau*(hnnlo1 + tau*hnnlo2)
      };
      // eq. 4.20 from Adam's thesis + conversion to MSbar top mass
      const double deltaNNNLO_S {467.683620788 - 8/3.*(19/8. + 2/3.*Nf) + (122.440972222 - 2*(19/8. + 2/3.*Nf))*Lt + 10.9409722222*Sqr(Lt)};

      const double alpha_s_red = alpha_s_5f/Pi;

      double scalar_corr = 0;
      const double Gamma0 = std::norm(3./4*A12_H);
      switch (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections))) {
         case 4:
         case 3:
            scalar_corr += deltaNNNLO_S*alpha_s_red;
         case 2:
            scalar_corr += deltaNNLO_S;
            scalar_corr *= alpha_s_red;
         case 1:
            scalar_corr += deltaNLO_S;
            scalar_corr *= alpha_s_red/Gamma0;
            //convert LO from 6 to 5 flavour scheme
            scalar_corr += 1. - Sqr(get_alphas(context)/alpha_s_5f);
            scalar_corr *= Gamma_SM_LO_S;
            break;
         default:
            WARNING("Unknow correcion in Phi->gg");
      }

      result += scalar_corr;

      if(info::is_CP_violating_Higgs_sector) {

         std::complex<double> const HGGVertexValP = 0.5*(HGGVertex.right()-HGGVertex.left());

         const std::complex<double> A12_A = 2.*f(tau)/tau;
         // LO width comming only from the top-loop
         // agrees up to a full double precision with automatically generated one
         const double Gamma_SM_LO_P = mH/(18.*Power3(Pi))*std::norm(alpha_s_5f * HGGVertexValP*sqrt(tau) * 3./4*A12_A);

         const double deltaNLO_P {
            97./4. - 7./6.*Nf + (33.-2*Nf)/6*LH
         };

         const double mtpole {qedqcd.displayPoleMt()};
         const double Lt = std::log(Sqr(mu/mtpole));
         const double deltaNNLO_P {
            237311./864. - 529./24.*zeta2 - 445./8.*zeta3 + 5.*Lt
         };

         double pseudoscalar_corr = 0.0;
         switch (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections))) {
            case 4:
            case 3:
            case 2:
               pseudoscalar_corr += deltaNNLO_P*alpha_s_red;
            case 1:
               pseudoscalar_corr += deltaNLO_P;
               pseudoscalar_corr *= alpha_s_red/std::norm(0.5*A12_A);
               pseudoscalar_corr += 1. - Sqr(get_alphas(context)/alpha_s_5f);
               pseudoscalar_corr *= Gamma_SM_LO_P;
               break;
            default:
               WARNING("Unknow correcion in Phi->gg");
         }
         result += pseudoscalar_corr;
      }
   }

   return result;
}


