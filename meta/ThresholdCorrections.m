
BeginPackage["ThresholdCorrections`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Constraint`"}];

SetDRbarGaugeCouplings::usage="";
SetDRbarYukawaCouplings::usage="";
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
           susyParticles = TreeMasses`GetSusyParticles[];
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
                  result -= prefactor dynkin Global`FiniteLog[Global`Sqr[FlexibleSUSY`M[particle]/Global`currentScale]] / 2;,
                  result -= Sum[prefactor dynkin Global`FiniteLog[Global`Sqr[FlexibleSUSY`M[particle][i]/Global`currentScale]] / 2,
                                {i,dimStart,dim}];
                 ];
              ];
           Return[coupling + (result + DRbarConversion[group]) CConversion`oneOver16PiSqr (coupling)^3];
          ];

SetDRbarGaugeCouplings[] :=
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
                    CConversion`RValueToCFormString[drBarCouplingG3] <> ";\n\n";
           (* calculate e *)
           result = result <> "const double " <> drBarCouplingEmCVariable <> " = " <>
                    CConversion`RValueToCFormString[drBarCouplingEm] <> ";\n\n";
           (* calculate sin(thetaW) *)
           (* result = result <> "const double " <> drBarCouplingSiCVariable <> " = " <> *)
           (*          CConversion`RValueToCFormString[drBarCouplingSi] <> ";\n\n"; *)
           drBarCouplingG1CVariable = drBarCouplingEmCVariable <> " * " <>
                                      RValueToCFormString[1/Parameters`GetGUTNormalization[couplingG1]] <>
                                      " / Sqrt(1.0 - Sqr(" <> drBarCouplingSiCVariable <> "))";
           drBarCouplingG2CVariable = drBarCouplingEmCVariable <> "/" <> drBarCouplingSiCVariable;
           (* write results back to the model *)
           result = result <>
                    Constraint`SetParameter[couplingG1, drBarCouplingG1CVariable, "model"] <>
                    Constraint`SetParameter[couplingG2, drBarCouplingG2CVariable, "model"] <>
                    Constraint`SetParameter[couplingG3, drBarCouplingG3CVariable, "model"];
           Return[result];
          ];

GetPrefactor[expr_Plus, yukawas_List] := 1;

GetPrefactor[expr_Integer, yukawas_List] := 1;

GetPrefactor[expr_Symbol, yukawas_List] := 1;

ContainsNoYukawa[expr_, yukawas_List] :=
    And @@ (FreeQ[expr, #]& /@ yukawas);

GetPrefactor[expr_Times, yukawas_List] :=
    Module[{factors, prefactors},
           factors = List @@ expr;
           prefactors = Select[factors, ContainsNoYukawa[#, yukawas]&];
           Times @@ prefactors
          ];

StripMatrixIndices[sym_Symbol] := sym;

StripMatrixIndices[sym_[_Integer, _Integer]] := sym;

ToMatrixSymbol[{}] := Null;

ToMatrixSymbol[list_List] :=
    Module[{dim, symbol, matrix, i, k, diag},
           dim = Length[list];
           symbol = StripMatrixIndices[list[[1,1]]];
           matrix = Table[symbol[i,k], {i,1,dim}, {k,1,dim}];
           diag = DiagonalMatrix[Table[symbol[i,i], {i,1,dim}]];
           Which[matrix === list, symbol,
                 Transpose[matrix] === list, Transpose[symbol],
                 diag === list, FlexibleSUSY`Diag[symbol],
                 True, Null
                ]
          ];

InvertRelation[Transpose[sym_], expr_] := {sym, SARAH`Tp[expr]};
InvertRelation[ConjugateTranspose[sym_], expr_] := {sym, SARAH`Adj[expr]};
InvertRelation[FlexibleSUSY`Diag[sym_], expr_] := {sym, FlexibleSUSY`Diag[expr]};
InvertRelation[sym_, expr_] := {sym, expr};

InvertMassRelation[fermion_, yukawas_List] :=
    Module[{massMatrix, polynom, prefactor, matrixSymbol},
           massMatrix = SARAH`MassMatrix[fermion];
           polynom = Factor[massMatrix /. List -> Plus];
           prefactor = GetPrefactor[polynom, yukawas];
           matrixSymbol = ToMatrixSymbol[massMatrix / prefactor];
           If[matrixSymbol === Null,
              Print["Error: could not convert expression to matrix symbol: ",
                    massMatrix / prefactor];
              Return[{Null, fermion}];
             ];
           InvertRelation[matrixSymbol, fermion / prefactor]
          ];

SetDRbarYukawaCouplings[] :=
    Module[{result, yTop, top, yBot, bot, yTau, tau, yukawas},
           yukawas = {SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa};
           {yTop, top} = InvertMassRelation[SARAH`TopQuark   , yukawas];
           {yBot, bot} = InvertMassRelation[SARAH`BottomQuark, yukawas];
           {yTau, tau} = InvertMassRelation[SARAH`Electron   , yukawas];
           top = top /. SARAH`TopQuark    -> Global`topDRbar;
           bot = bot /. SARAH`BottomQuark -> Global`bottomDRbar;
           tau = tau /. SARAH`Electron    -> Global`electronDRbar;
           result = Constraint`CreateLocalConstRefs[top + bot + tau] <>
                    Constraint`SetParameter[yTop, RValueToCFormString[top], "model"] <>
                    Constraint`SetParameter[yBot, RValueToCFormString[bot], "model"] <>
                    Constraint`SetParameter[yTau, RValueToCFormString[tau], "model"];
           Return[result];
          ];

End[];

EndPackage[];
