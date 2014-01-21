
BeginPackage["ThresholdCorrections`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Constraint`"}];

CalculateGaugeCouplings::usage="";
CalculateDeltaAlphaEm::usage="";
CalculateDeltaAlphaS::usage="";
SetDRbarYukawaCouplings::usage="";

Begin["`Private`"];

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
                  result -= Sum[prefactor dynkin Global`FiniteLog[Abs[FlexibleSUSY`M[particle][i-1]/Global`currentScale]],
                                {i,dimStart,dim}];
                 ];
              ];
           Return[result + DRbarConversion[group]];
          ];

CalculateDeltaAlphaEm[] :=
    Module[{result, deltaSusy, deltaSM, prefactor},
           prefactor = Global`alphaEm / (2 Pi);
           deltaSM = 1/3 - 16/9 Global`FiniteLog[Abs[FlexibleSUSY`M[SARAH`TopQuark][2]/Global`currentScale]];
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
           deltaSM = - 2/3 Global`FiniteLog[Abs[FlexibleSUSY`M[SARAH`TopQuark][2]/Global`currentScale]];
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

GetPrefactor[expr_Times, yukawa_] :=
    Module[{factors, prefactors},
           factors = List @@ expr;
           prefactors = Select[factors, FreeQ[#, yukawa]&];
           Times @@ prefactors
          ];

ExtractSymbols[sym_?NumberQ] := {};

ExtractSymbols[sym_Symbol] := {sym};

ExtractSymbols[HoldPattern[SARAH`sum[_,_,_,expr_]]] := ExtractSymbols[expr];

ExtractSymbols[expr_Plus] := Flatten[ExtractSymbols /@ (List @@ expr)];

ExtractSymbols[expr_Times] := Flatten[ExtractSymbols /@ (List @@ expr)];

ExtractSymbols[sym_[_,_]] := {sym};

ToMatrixExpression[{}] := Null;

ToMatrixExpression[list_List] :=
    Module[{dim, symbol, matrix, i, k, diag, expression = Null,
            expandedList, permutations},
           dim = Length[list];
           symbol = ExtractSymbols[list[[1,1]]];
           (* expand sum[] *)
           expandedList = Expand[list //. SARAH`sum[idx_,start_,stop_,expr_] :> Sum[expr,{idx,start,stop}]];
           If[Length[symbol] == 1,
              symbol = symbol[[1]];
              matrix = Table[symbol[i,k], {i,1,dim}, {k,1,dim}];
              diag = DiagonalMatrix[Table[symbol[i,i], {i,1,dim}]];
              Which[matrix === expandedList, expression = symbol;,
                    Transpose[matrix] === expandedList, expression = SARAH`Tp[symbol];,
                    diag === expandedList, expression = FlexibleSUSY`Diag[symbol];
                   ];
              ,
              (* create all permutations of matrix products *)
              permutations = Permutations[symbol];
              For[i = 1, i <= Length[permutations], i++,
                  matrix = Table[#[i,k],{i,1,dim},{k,1,dim}]& /@ permutations[[i]];
                  matrix = Expand[Dot @@ matrix];
                  If[matrix === expandedList,
                     (* found a combination which yields our list *)
                     expression = SARAH`MatMul @@ permutations[[i]];
                     Break[];
                    ];
                 ];
             ];
           Return[expression];
          ];

(* Solve the equation #1 == #2 for #3 *)
InvertRelation[Transpose[sym_], expr_, sym_] :=
    {sym, SARAH`Tp[expr]};

InvertRelation[SARAH`Tp[sym_], expr_, sym_] :=
    {sym, SARAH`Tp[expr]};

InvertRelation[ConjugateTranspose[sym_], expr_, sym_] :=
    {sym, SARAH`Adj[expr]};

InvertRelation[SARAH`Adj[sym_], expr_, sym_] :=
    {sym, SARAH`Adj[expr]};

InvertRelation[FlexibleSUSY`Diag[sym_], expr_, sym_] :=
    {sym, FlexibleSUSY`Diag[expr]};

InvertRelation[sym_, expr_, sym_] :=
    {sym, expr};

(* remove matrices from the left *)
InvertRelation[SARAH`MatMul[SARAH`Adj[U_],X___,sym_,V___], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[X,sym,V], SARAH`MatMul[U,expr], sym];

InvertRelation[SARAH`MatMul[U_,X___,sym_,V___], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[X,sym,V], SARAH`MatMul[SARAH`Adj[U],expr], sym];

(* remove matrices from the right *)
InvertRelation[SARAH`MatMul[sym_,V___,SARAH`Adj[U_]], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[sym,V], SARAH`MatMul[expr,U], sym];

InvertRelation[SARAH`MatMul[sym_,U__], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[sym], SARAH`MatMul[expr,SARAH`Adj[U]], sym];

InvertRelation[SARAH`MatMul[sym_], expr_, sym_] :=
    InvertRelation[sym, expr, sym];

InvertRelation[sym_, expr_, other_] :=
    Block[{},
          Print["Error: InvertRelation: don't know how to solve equation: ",
                sym, " == ", expr, " for ", other];
          Quit[1];
         ];

InvertMassRelation[fermion_, yukawa_] :=
    Module[{massMatrix, polynom, prefactor, matrixExpression},
           massMatrix = SARAH`MassMatrix[fermion];
           polynom = Factor[massMatrix /. List -> Plus];
           prefactor = GetPrefactor[polynom, yukawa];
           matrixExpression = ToMatrixExpression[massMatrix / prefactor];
           If[matrixExpression === Null,
              Print["Error: could not convert list to matrix expression: ",
                    massMatrix / prefactor];
              Quit[1];
             ];
           InvertRelation[matrixExpression, fermion / prefactor, yukawa]
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

CalculateGaugeCouplings[] :=
    Module[{subst, weinbergAngle, g1Def, g2Def, g3Def, result},
           subst = { SARAH`Mass[SARAH`VectorW] -> FlexibleSUSY`MWDRbar,
                     SARAH`Mass[SARAH`VectorZ] -> FlexibleSUSY`MZDRbar,
                     SARAH`electricCharge      -> FlexibleSUSY`EDRbar };
           weinbergAngle = Parameters`FindSymbolDef[SARAH`Weinberg] /. subst;
           g1Def = (Parameters`FindSymbolDef[SARAH`hyperchargeCoupling]
                    / Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]) /. subst;
           g2Def = (Parameters`FindSymbolDef[SARAH`leftCoupling]
                    / Parameters`GetGUTNormalization[SARAH`leftCoupling]) /. subst;
           g3Def = (Parameters`FindSymbolDef[SARAH`strongCoupling]
                    / Parameters`GetGUTNormalization[SARAH`strongCoupling]) /. subst;
           result = Parameters`CreateLocalConstRefs[{weinbergAngle, g1Def, g2Def, g3Def}] <>
                    "const double " <> CConversion`ToValidCSymbolString[SARAH`Weinberg] <>
                    " = " <> CConversion`RValueToCFormString[weinbergAngle] <> ";\n" <>
                    "new_g1 = " <> CConversion`RValueToCFormString[g1Def] <> ";\n" <>
                    "new_g2 = " <> CConversion`RValueToCFormString[g2Def] <> ";\n" <>
                    "new_g3 = " <> CConversion`RValueToCFormString[g3Def] <> ";\n";
           Return[result];
          ];

End[];

EndPackage[];
