
BeginPackage["ThresholdCorrections`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Constraint`", "Vertices`"}];

CalculateGaugeCouplings::usage="";
CalculateDeltaAlphaEm::usage="";
CalculateDeltaAlphaS::usage="";
CalculateThetaW::usage="";
SetDRbarYukawaCouplingTop::usage="";
SetDRbarYukawaCouplingBottom::usage="";
SetDRbarYukawaCouplingElectron::usage="";

Begin["`Private`"];

DRbarConversion[SARAH`U[1]] := 0;

DRbarConversion[SARAH`SU[n_Integer]] := n/6;

DRbarConversion[group_] := Null;

CalculateColorCoupling[scheme_] :=
    CalculateCoupling[TreeMasses`FindColorGaugeGroup[], scheme];

CalculateElectromagneticCoupling[scheme_] :=
    CalculateCoupling[{SARAH`electricCharge, FlexibleSUSY`electricCharge, SARAH`U[1]}, scheme];

CalculateCoupling::UnknownRenormalizationScheme = "Unknown\
 renormalization scheme `1`.";

(* Calculate threshold correction for a gauge coupling from SM
   (MS-bar) to a given renormalization scheme in the given model. *)
CalculateCoupling[{coupling_, name_, group_}, scheme_] :=
    Module[{susyParticles, prefactor, i, result = 0, particle, dynkin,
            dim, dimStart, nc, casimir, conversion},
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
           conversion = Switch[scheme,
                               FlexibleSUSY`DRbar, DRbarConversion[group],
                               FlexibleSUSY`MSbar, 0,
                               _, Message[CalculateCoupling::UnknownRenormalizationScheme, scheme]; 0
                              ];
           Return[result + conversion];
          ];

CalculateDeltaAlphaEm[renormalizationScheme_] :=
    Module[{result, deltaSusy, deltaSM, prefactor, topQuark},
           topQuark = TreeMasses`GetThirdGenerationMass[SARAH`TopQuark];
           prefactor = Global`alphaEm / (2 Pi);
           deltaSM = 1/3 - 16/9 Global`FiniteLog[Abs[topQuark/Global`currentScale]];
           deltaSusy = CalculateElectromagneticCoupling[renormalizationScheme];
           result = Parameters`CreateLocalConstRefs[deltaSusy + deltaSM] <> "\n" <>
                    "const double delta_alpha_em_SM = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSM] <> ";\n\n" <>
                    "const double delta_alpha_em = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSusy] <> ";\n\n" <>
                    "return delta_alpha_em + delta_alpha_em_SM;\n";
           Return[result];
          ];

CalculateDeltaAlphaS[renormalizationScheme_] :=
    Module[{result, deltaSusy, deltaSM, prefactor, topQuark},
           topQuark = TreeMasses`GetThirdGenerationMass[SARAH`TopQuark];
           prefactor = Global`alphaS / (2 Pi);
           deltaSM = - 2/3 Global`FiniteLog[Abs[topQuark/Global`currentScale]];
           deltaSusy = CalculateColorCoupling[renormalizationScheme];
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

ToMatrixExpression[expr_ /; Head[expr] =!= List] := expr;

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

InvertRelation[sym_[i1_,i2_], expr_, sym_] :=
    {sym[i1,i2], expr};

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
    Module[{massMatrix, polynom, prefactor, matrixExpression, dim},
           If[TreeMasses`IsUnmixed[fermion],
              massMatrix = TreeMasses`GetMassOfUnmixedParticle[fermion];
              massMatrix = TreeMasses`ReplaceDependencies[massMatrix];
              massMatrix = Vertices`StripGroupStructure[massMatrix, {SARAH`ct1, SARAH`ct2}];
              ,
              massMatrix = SARAH`MassMatrix[fermion];
             ];
           dim = Length[massMatrix];
           If[massMatrix === Table[0, {i,1,dim}, {k,1,dim}],
              Return[{yukawa,FlexibleSUSY`ZEROMATRIX[dim,dim]}];
             ];
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

SetDRbarYukawaCouplingTop[settings_] :=
    SetDRbarYukawaCouplingFermion[SARAH`TopQuark, SARAH`UpYukawa, Global`topDRbar, settings];

SetDRbarYukawaCouplingBottom[settings_] :=
    SetDRbarYukawaCouplingFermion[SARAH`BottomQuark, SARAH`DownYukawa, Global`bottomDRbar, settings];

SetDRbarYukawaCouplingElectron[settings_] :=
    SetDRbarYukawaCouplingFermion[SARAH`Electron, SARAH`ElectronYukawa, Global`electronDRbar, settings];

SetDRbarYukawaCouplingFermion[fermion_, yukawa_, mass_, settings_] :=
    Module[{y, f},
           f = Cases[settings, {yukawa, s_} :> s];
           If[f === {}, Return[""];];
           f = f[[1]];
           If[f === Automatic,
              {y, f} = InvertMassRelation[fermion, yukawa];
              f = f /. fermion -> mass;
              ,
              y = yukawa;
             ];
           Parameters`CreateLocalConstRefs[f] <>
           Parameters`SetParameter[yukawa, f, "MODEL"]
          ];

MultiplyBy[factor_ /; factor == 1] := "";

MultiplyBy[factor_] :=
    " * " <> CConversion`RValueToCFormString[factor];

GetParameter[par_[idx1_,idx2_], factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[par] <>
    "(" <> CConversion`RValueToCFormString[idx1] <> "," <>
    CConversion`RValueToCFormString[idx2] <> ")" <>
    MultiplyBy[factor];

GetParameter[FlexibleSUSY`M[par_], factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[FlexibleSUSY`M[par]] <> "()" <>
    MultiplyBy[factor];

GetParameter[par_[idx_], factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[par] <>
    "(" <> CConversion`RValueToCFormString[idx] <> ")" <>
    MultiplyBy[factor];

GetParameter[par_, factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[par] <> "()" <>
    MultiplyBy[factor];

CalculateThetaWFromFermiConstant[] :=
    Module[{g1Str, g2Str, g3Str, vuStr, vdStr,
            mTop, mBot, mHiggs,
            mtStr, mbStr,
            mnStr, mcStr, znStr, umStr, upStr,
            mhStr, zhStr,
            mseStr, zseStr,
            msvStr, zsvStr,
            zStr, wStr, ymStr
           },
           mTop    = TreeMasses`GetThirdGenerationMass[SARAH`TopQuark];
           mBot    = TreeMasses`GetThirdGenerationMass[SARAH`BottomQuark];
           mHiggs  = TreeMasses`GetLightestMass[SARAH`HiggsBoson];
           mtStr   = GetParameter[mTop];
           mbStr   = GetParameter[mBot];
           g1Str   = GetParameter[SARAH`hyperchargeCoupling, Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]];
           g2Str   = GetParameter[SARAH`leftCoupling       , Parameters`GetGUTNormalization[SARAH`leftCoupling]];
           g3Str   = GetParameter[SARAH`strongCoupling     , Parameters`GetGUTNormalization[SARAH`strongCoupling]];
           vuStr   = GetParameter[SARAH`VEVSM2];
           vdStr   = GetParameter[SARAH`VEVSM1];
           mnStr   = GetParameter[FlexibleSUSY`M[Parameters`GetParticleFromDescription["Neutralinos"]]];
           mcStr   = GetParameter[FlexibleSUSY`M[Parameters`GetParticleFromDescription["Charginos"]]];
           znStr   = GetParameter[SARAH`NeutralinoMM];
           umStr   = GetParameter[SARAH`CharginoMinusMM];
           upStr   = GetParameter[SARAH`CharginoPlusMM];
           mhStr   = GetParameter[mHiggs];
           zhStr   = GetParameter[SARAH`HiggsMixingMatrix[0,1]];
           mseStr  = GetParameter[FlexibleSUSY`M[SARAH`Selectron]];
           zseStr  = GetParameter[SARAH`SleptonMM];
           msvStr  = GetParameter[FlexibleSUSY`M[SARAH`Sneutrino]];
           zsvStr  = GetParameter[SARAH`SneutrinoMM];
           wStr    = CConversion`RValueToCFormString[SARAH`VectorW];
           zStr    = CConversion`RValueToCFormString[SARAH`VectorZ];
           ymStr   = GetParameter[SARAH`ElectronYukawa[1,1]];
    "\
using namespace weinberg_angle;

const double scale         = MODEL->get_scale();
const double gfermi        = softsusy::GMU; // @todo
const double mw_pole       = oneset.displayPoleMW();
const double mz_pole       = oneset.displayPoleMZ();
const double mt_pole       = oneset.displayPoleMt();
const double mt_drbar      = " <> mtStr <> ";
const double mb_drbar      = " <> mbStr <> ";
const double mh_drbar      = " <> mhStr <> ";
const double alphaDrbar    = ALPHA_EM_DRBAR;
const double gY            = " <> g1Str <> ";
const double g2            = " <> g2Str <> ";
const double g3            = " <> g3Str <> ";
const double ymu           = " <> ymStr <> ";
const double alpha         = " <> zhStr <> ";
const double tanBeta       = " <> vuStr <> " / " <> vdStr <> ";
double mselL               = 0.;
double msmuL               = 0.;
double msnue               = 0.;
double msnumu              = 0.;
const Eigen::ArrayXd mneut(" <> mnStr <> ");
const Eigen::MatrixXcd zn (" <> znStr <> ");
const Eigen::ArrayXd mch  (" <> mcStr <> ");
const Eigen::MatrixXcd um (" <> umStr <> ");
const Eigen::MatrixXcd up (" <> upStr <> ");
const auto MSe(" <> mseStr <> ");
const auto ZE(" <> zseStr <> ");
const auto MSv(" <> msvStr <> ");
const auto ZV(" <> zsvStr <> ");

for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
   mselL += AbsSqr(ZE(i,0))*MSe(i);
   msmuL += AbsSqr(ZE(i,1))*MSe(i);
}

for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
   msnue  += AbsSqr(ZV(i,0))*MSv(i);
   msnumu += AbsSqr(ZV(i,1))*MSv(i);
}

const double pizztMZ = Re(MODEL->self_energy_" <> zStr <> "(mz_pole));
const double piwwt0  = Re(MODEL->self_energy_" <> wStr <> "(0));
const double piwwtMW = Re(MODEL->self_energy_" <> wStr <> "(mw_pole));

Weinberg_angle::Self_energy_data se_data;
se_data.scale    = scale;
se_data.mt_pole  = mt_pole;
se_data.mt_drbar = mt_drbar;
se_data.mb_drbar = mb_drbar;
se_data.gY       = gY;
se_data.g2       = g2;

const double pizztMZ_corrected =
   Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole, se_data);

const double piwwtMW_corrected =
   Weinberg_angle::replace_mtop_in_self_energy_w(piwwtMW, mw_pole, se_data);

const double piwwt0_corrected =
   Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0., se_data);

Weinberg_angle::Data data;
data.scale               = scale;
data.alpha_em_drbar      = alphaDrbar;
data.fermi_contant       = gfermi;
data.self_energy_z_at_mz = pizztMZ_corrected;
data.self_energy_w_at_mw = piwwtMW_corrected;
data.self_energy_w_at_0  = piwwt0_corrected;
data.mw_pole             = mw_pole;
data.mz_pole             = mz_pole;
data.mt_pole             = mt_pole;
data.mh_drbar            = mh_drbar;
data.hmix_12             = alpha;
data.msel_drbar          = mselL;
data.msmul_drbar         = msmuL;
data.msve_drbar          = msnue;
data.msvm_drbar          = msnumu;
data.mn_drbar            = mneut;
data.mc_drbar            = mch;
data.zn                  = zn;
data.um                  = um;
data.up                  = up;
data.gY                  = gY;
data.g2                  = g2;
data.g3                  = g3;
data.tan_beta            = tanBeta;
data.ymu                 = ymu;

Weinberg_angle weinberg;
weinberg.set_data(data);

THETAW = ArcSin(weinberg.calculate());"
          ];

CalculateThetaWFromMW[] :=
    Module[{subst, weinbergAngle, result},
           subst = { SARAH`Mass[SARAH`VectorW] -> FlexibleSUSY`MWDRbar,
                     SARAH`Mass[SARAH`VectorZ] -> FlexibleSUSY`MZDRbar,
                     SARAH`electricCharge      -> FlexibleSUSY`EDRbar };
           weinbergAngle = Parameters`FindSymbolDef[SARAH`Weinberg] /. subst;
           result = Parameters`CreateLocalConstRefs[{weinbergAngle}] <>
                    "THETAW = " <>
                    CConversion`RValueToCFormString[weinbergAngle] <> ";\n";
           Return[result];
          ];

CalculateThetaW[input_] :=
    Switch[input,
           FlexibleSUSY`FSFermiConstant, CalculateThetaWFromFermiConstant[],
           FlexibleSUSY`FSMassW, CalculateThetaWFromMW[],
           _,
           Print["Error: CalculateThetaW[", input, "]: unknown input ", input];
           Print["   Using default input: ", FlexibleSUSY`FSMassW];
           CalculateThetaWFromMW[]
          ];

CalculateGaugeCouplings[] :=
    Module[{subst, g1Def, g2Def, g3Def, result},
           subst = { SARAH`Mass[SARAH`VectorW] -> FlexibleSUSY`MWDRbar,
                     SARAH`Mass[SARAH`VectorZ] -> FlexibleSUSY`MZDRbar,
                     SARAH`electricCharge      -> FlexibleSUSY`EDRbar,
                     SARAH`Weinberg            -> FlexibleSUSY`ThetaW };
           g1Def = (Parameters`FindSymbolDef[SARAH`hyperchargeCoupling]
                    / Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]) /. subst;
           g2Def = (Parameters`FindSymbolDef[SARAH`leftCoupling]
                    / Parameters`GetGUTNormalization[SARAH`leftCoupling]) /. subst;
           g3Def = (Parameters`FindSymbolDef[SARAH`strongCoupling]
                    / Parameters`GetGUTNormalization[SARAH`strongCoupling]) /. subst;
           result = Parameters`CreateLocalConstRefs[{g1Def, g2Def, g3Def}] <>
                    "new_" <> CConversion`ToValidCSymbolString[SARAH`hyperchargeCoupling] <>
                    " = " <> CConversion`RValueToCFormString[g1Def] <> ";\n" <>
                    "new_" <> CConversion`ToValidCSymbolString[SARAH`leftCoupling] <>
                    " = " <> CConversion`RValueToCFormString[g2Def] <> ";\n" <>
                    "new_" <> CConversion`ToValidCSymbolString[SARAH`strongCoupling] <>
                    " = " <> CConversion`RValueToCFormString[g3Def] <> ";\n";
           Return[result];
          ];

End[];

EndPackage[];
