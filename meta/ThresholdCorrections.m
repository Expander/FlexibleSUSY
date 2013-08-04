
BeginPackage["ThresholdCorrections`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Constraint`"}];

CalculateDeltaAlphaEm::usage="";
CalculateDeltaAlphaS::usage="";
SetDRbarYukawaCouplings::usage="";

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
                  result -= prefactor dynkin Global`FiniteLog[Abs[FlexibleSUSY`M[particle]/Global`currentScale]];,
                  result -= Sum[prefactor dynkin Global`FiniteLog[Abs[FlexibleSUSY`M[particle][i]/Global`currentScale]],
                                {i,dimStart,dim}];
                 ];
              ];
           Return[result + DRbarConversion[group]];
          ];

CalculateDeltaAlphaEm[] :=
    Module[{result, deltaSusy, deltaSM, prefactor},
           prefactor = Global`alphaEm / (2 Pi);
           deltaSM = 1/3 - 16/9 Global`FiniteLog[Abs[FlexibleSUSY`M[SARAH`TopQuark][3]/Global`currentScale]];
           deltaSusy = CalculateDRbarElectromagneticCoupling[];
           result = Parameters`CreateLocalConstRefs[deltaSusy + deltaSM] <> "\n" <>
                    "const double delta_alpha_em_SM = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSM] <> ";\n\n" <>
                    "const double delta_alpha_em = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSusy] <> ";\n\n" <>
                    "return delta_alpha_em + delta_alpha_em_SM;\n";
           Return[result];
          ];

CalculateDeltaAlphaS[] :=
    Module[{result, deltaSusy, deltaSM, prefactor},
           prefactor = Global`alphaS / (2 Pi);
           deltaSM = - 2/3 Global`FiniteLog[Abs[FlexibleSUSY`M[SARAH`TopQuark][3]/Global`currentScale]];
           deltaSusy = CalculateDRbarColorCoupling[];
           result = Parameters`CreateLocalConstRefs[deltaSusy + deltaSM] <> "\n" <>
                    "const double delta_alpha_s_SM = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSM] <> ";\n\n" <>
                    "const double delta_alpha_s = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSusy] <> ";\n\n" <>
                    "return delta_alpha_s + delta_alpha_s_SM;\n";
           Return[result];
          ];

GetPrefactor[expr_Plus, _] := 1;

GetPrefactor[expr_Integer, _] := 1;

GetPrefactor[expr_Symbol, _] := 1;

ContainsNoYukawa[expr_, yukawa_] := FreeQ[expr, yukawa];

GetPrefactor[expr_Times, yukawa_] :=
    Module[{factors, prefactors},
           factors = List @@ expr;
           prefactors = Select[factors, ContainsNoYukawa[#, yukawa]&];
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

(* Solve the equation #1 == #2 for #3 *)
InvertRelation[Transpose[sym_], expr_, sym_] :=
    {sym, SARAH`Tp[expr]};

InvertRelation[ConjugateTranspose[sym_], expr_, sym_] :=
    {sym, SARAH`Adj[expr]};

InvertRelation[FlexibleSUSY`Diag[sym_], expr_, sym_] :=
    {sym, FlexibleSUSY`Diag[expr]};

InvertRelation[sym_, expr_, sym_] :=
    {sym, expr};

InvertMassRelation[fermion_, yukawa_] :=
    Module[{massMatrix, polynom, prefactor, matrixSymbol},
           massMatrix = SARAH`MassMatrix[fermion];
           polynom = Factor[massMatrix /. List -> Plus];
           prefactor = GetPrefactor[polynom, yukawa];
           matrixSymbol = ToMatrixSymbol[massMatrix / prefactor];
           If[matrixSymbol === Null,
              Print["Error: could not convert expression to matrix symbol: ",
                    massMatrix / prefactor];
              Quit[1];
              Return[{Null, fermion}];
             ];
           InvertRelation[matrixSymbol, fermion / prefactor, yukawa]
          ];

SetDRbarYukawaCouplings[] :=
    Module[{result, yTop, top, yBot, bot, yTau, tau},
           {yTop, top} = InvertMassRelation[SARAH`TopQuark   , SARAH`UpYukawa];
           {yBot, bot} = InvertMassRelation[SARAH`BottomQuark, SARAH`DownYukawa];
           {yTau, tau} = InvertMassRelation[SARAH`Electron   , SARAH`ElectronYukawa];
           top = top /. SARAH`TopQuark    -> Global`topDRbar;
           bot = bot /. SARAH`BottomQuark -> Global`bottomDRbar;
           tau = tau /. SARAH`Electron    -> Global`electronDRbar;
           result = Parameters`CreateLocalConstRefs[top + bot + tau] <>
                    "new_Yu = " <> RValueToCFormString[top] <> ";\n" <>
                    "new_Yd = " <> RValueToCFormString[bot] <> ";\n" <>
                    "new_Ye = " <> RValueToCFormString[tau] <> ";\n";
           Return[result];
          ];

End[];

EndPackage[];
