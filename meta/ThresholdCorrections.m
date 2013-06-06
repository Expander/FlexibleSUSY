
BeginPackage["ThresholdCorrections`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Constraint`"}];

SetDRbarCouplings::usage="";
CalculateDRbarColorCoupling::usage="";
CalculateDRbarElectromagneticCoupling::usage="";

Begin["Private`"];

DRbarConversion[SARAH`U[1]] := 0;

DRbarConversion[SARAH`SU[n_Integer]] := n/6;

DRbarConversion[group_] := Null;

CalculateDRbarHyperchargeCoupling[] :=
    CalculateDRbarCoupling[SARAH`hypercharge];

CalculateDRbarColorCoupling[] :=
    CalculateDRbarCoupling[TreeMasses`FindColorGaugeGroup[]];

CalculateDRbarElectromagneticCoupling[] :=
    CalculateDRbarCoupling[{SARAH`electricCharge, FlexibleSUSY`electricCharge, SARAH`U[1]}];

CalculateDRbarLeftCoupling[] :=
    CalculateDRbarCoupling[SARAH`left];

CalculateDRbarCoupling[{coupling_, name_, group_}] :=
    Module[{susyParticles, prefactor, i, result = 0, particle, dynkin,
            dim, dimStart, nc, casimir},
           susyParticles = Select[TreeMasses`GetSusyParticles[], (!TreeMasses`IsGhost[#])&];
           For[i = 1, i <= Length[susyParticles], i++,
               particle = susyParticles[[i]];
               dim = GetDimension[particle];
               dimStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
               Which[IsMajoranaFermion[particle], prefactor = 2/3,
                     IsDiracFermion[   particle], prefactor = 4/3,
                     IsRealScalar[     particle], prefactor = 1/6,
                     IsComplexScalar[  particle], prefactor = 1/3,
                     IsVector[         particle], prefactor = 1/3,
                     True                       , prefactor = 1];
               If[coupling === SARAH`electricCharge,
                  casimir = SA`Casimir[particle, SARAH`color];
                  nc = If[NumericQ[casimir] && casimir =!= 0, 3, 1];
                  dynkin = SARAH`getElectricCharge[particle]^2 nc;
                  ,
                  dynkin = SA`Dynkin[particle, Position[SARAH`Gauge, name][[1, 1]]];
                 ];
               (* some dynkin indices are not defined in SARAH *)
               If[!NumericQ[dynkin], dynkin = 0];
               If[dim == 1,
                  result -= prefactor dynkin Log[particle/Global`currentScale];,
                  result -= Sum[prefactor dynkin Log[particle[i]/Global`currentScale], {i,dimStart,dim}];
                 ];
              ];
           Simplify[coupling + (result + DRbarConversion[group]) CConversion`oneOver16PiSqr (coupling)^3]
          ];

SetDRbarCouplings[] :=
    Module[{result,
            couplingG1, couplingG2, couplingG3, couplingEm,
            drBarCouplingG3, drBarCouplingEm, drBarCouplingSi,
            drBarCouplingG1CVariable, drBarCouplingG2CVariable,
            drBarCouplingG3CVariable,
            drBarCouplingEmCVariable, drBarCouplingSiCVariable},
           (* get couplings *)
           couplingG1 = FindHyperchargeGaugeCoupling[];
           couplingG2 = FindLeftGaugeCoupling[];
           couplingG3 = FindColorGaugeCoupling[];
           couplingEm = SARAH`electricCharge;
           (* calculate drBar couplings *)
           drBarCouplingG3 = CalculateDRbarColorCoupling[] /.
              Rule[couplingG3, SARAH`SM[couplingG3]];
           drBarCouplingEm = CalculateDRbarElectromagneticCoupling[] /.
              Rule[couplingEm, SARAH`SM[couplingEm]];
           drBarCouplingSi = Global`sinThetaW /.
              Rule[Global`sinThetaW, SARAH`SM[Global`sinThetaW]];
           (* create C variables *)
           drBarCouplingG3CVariable = ToValidCSymbolString[couplingG3] <> "_drbar";
           drBarCouplingEmCVariable = ToValidCSymbolString[couplingEm] <> "_drbar";
           drBarCouplingSiCVariable = "sinThetaW_drbar";
           (* create local const refs *)
           result = Constraint`CreateLocalConstRefs[drBarCouplingG3 + drBarCouplingEm] <> "\n";
           (* calculate g3 *)
           result = result <> "const double " <> drBarCouplingG3CVariable <> " = " <>
                    RValueToCFormString[drBarCouplingG3] <> ";\n\n";
           (* calculate e *)
           result = result <> "const double " <> drBarCouplingEmCVariable <> " = " <>
                    RValueToCFormString[drBarCouplingEm] <> ";\n\n";
           (* calculate sin(thetaW) *)
           result = result <> "const double " <> drBarCouplingSiCVariable <> " = " <>
                    RValueToCFormString[drBarCouplingSi] <> ";\n\n";
           drBarCouplingG1CVariable = drBarCouplingEmCVariable <> " * " <>
                                      RValueToCFormString[1/Parameters`GetGUTNormalization[couplingG1]] <>
                                      " / cos(asin(" <> drBarCouplingSiCVariable <> "))";
           drBarCouplingG2CVariable = drBarCouplingEmCVariable <> "/" <> drBarCouplingSiCVariable;
           (* write results back to the model *)
           result = result <>
                    Constraint`SetParameter[couplingG1, drBarCouplingG1CVariable, "model"] <>
                    Constraint`SetParameter[couplingG2, drBarCouplingG2CVariable, "model"] <>
                    Constraint`SetParameter[couplingG3, drBarCouplingG3CVariable, "model"];
           Return[result];
          ];

End[];

EndPackage[];
