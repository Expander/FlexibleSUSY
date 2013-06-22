
BeginPackage["CConversion`", {"SARAH`", "TextFormatting`"}];

MatrixType::usage="";
VectorType::usage="";
ScalarType::usage="";

UNITMATRIX::usage="";
oneOver16PiSqr::usage="";
twoLoop::usage="";
AbsSqr::usage="";

ToValidCSymbol::usage="creates a valid C variable name from a symbol";

ToValidCSymbolString::usage="returns the result of ToValidCSymbol[] as
a string";

RValueToCFormString::usage="converts a Mathematica expression to a C
expression.";

GetHead::usage="returns the head of a symbol
GetHead[s]        ->  s
GetHead[s[a]]     ->  s
GetHead[s[a][b]]  ->  s";

CreateUnitMatrix::usage="creates a unit matrix for a given parameter
type";

GetCParameterType::usage="returns the C/C++ data type of a given
Mathematica symbol";

CreateGetterPrototype::usage="creates C/C++ getter prototype";

CreateInlineSetter::usage="creates C/C++ inline setter"

CreateInlineGetter::usage="creates C/C++ inline getter"

CreateGetterReturnType::usage="creates C/C++ getter return type";

CreateDefaultConstructor::usage="creates C/C++ default constructor for
a given parameter type";

CreateDefaultDefinition::usage="creates C/C++ variable definition
using the default constructor";

ExpandSums::usage="expands expressions that contain sum symbols of the
form sum[index,1,3,expression]"

MakeUnique::usage="create a unique symbol from a string";

Begin["Private`"];

(* This rule is essential for the ExpandSums[] function.
 * It prevents the following bug:
 *
 *   In[]:= expr = Power[ThetaStep[i1,5],2];
 *   In[]:= DeleteCases[expr, ThetaStep[_,_]]
 *   Out[]= 2
 *)
SARAH`ThetaStep /: Power[SARAH`ThetaStep[a_,b_],_] := SARAH`ThetaStep[a,b];

CreateGetterReturnType[type_] :=
    Print["Error: unknown type: " <> ToString[type]];

CreateGetterReturnType[CConversion`ScalarType[type_]] :=
    ToString[type];

CreateGetterReturnType[CConversion`VectorType[type_, entries_]] :=
    "const " <> ToString[type] <> "&";

CreateGetterReturnType[CConversion`MatrixType[type_, rows_, cols_]] :=
    "const " <> ToString[type] <> "&";

CreateSetterInputType[type_] :=
    CreateGetterReturnType[type];

(* Creates a C++ setter *)
CreateInlineSetter[parameter_String, type_String] :=
    "void set_" <> parameter <> "(" <> type <>
    " " <> parameter <> "_) { " <> parameter <> " = " <>
    parameter <> "_; }\n";

CreateInlineSetter[parameter_String, type_] :=
    CreateInlineSetter[parameter, CreateSetterInputType[type]];

(* Creates a C++ inline getter *)
CreateInlineGetter[parameter_String, type_String] :=
    type <> " get_" <> parameter <>
    "() const { return " <> parameter <> "; }\n";

CreateInlineGetter[parameter_, type_] :=
    CreateInlineGetter[parameter, CreateGetterReturnType[type]];

(* Creates C++ getter prototype *)
CreateGetterPrototype[parameter_String, type_String] :=
    type <> " get_" <> parameter <> "() const;\n";

CreateGetterPrototype[parameter_, type_] :=
    CreateGetterPrototype[parameter, CreateGetterReturnType[type]];

(* create default constructor initialization list *)
CreateDefaultConstructor[parameter_, type_] :=
    Print["Error: unknown parameter type: " <> ToString[type]];

CreateDefaultConstructor[parameter_, CConversion`ScalarType[type_]] :=
    parameter <> "(0)";

CreateDefaultConstructor[parameter_, CConversion`VectorType[type_, entries_]] :=
    parameter <> "(" <> type <> "::Zero())";

CreateDefaultConstructor[parameter_, CConversion`MatrixType[type_, rows_, cols_]] :=
    parameter <> "(" <> type <> "::Zero())";

CreateDefaultDefinition[parameter_, type_] :=
    Print["Error: unknown parameter type: " <> ToString[type]];

CreateDefaultDefinition[parameter_, CConversion`ScalarType[type_]] :=
    type <> " " <> parameter <> " = 0";

CreateDefaultDefinition[parameter_, CConversion`VectorType[type_, _]] :=
    type <> " " <> parameter;

CreateDefaultDefinition[parameter_, CConversion`MatrixType[type_, _, _]] :=
    type <> " " <> parameter;

CreateDefaultDefinition[parameter_, type_] :=
    Print["Error: unknown parameter type: " <> ToString[type]];

CreateDefaultDefinition[parameter_, CConversion`ScalarType[type_]] :=
    type <> " " <> parameter <> " = 0";

CreateDefaultDefinition[parameter_, CConversion`VectorType[type_, entries_]] :=
    type <> " " <> parameter <> "(" <> ToString[entries] <> ")";

CreateDefaultDefinition[parameter_, CConversion`MatrixType[type_, rows_, cols_]] :=
    type <> " " <> parameter <> "(" <> ToString[rows] <> "," <> ToString[cols] <> ")";

GetCParameterType[parameterType_] :=
    ToString[parameterType[[1]]];

(* create unitary matrix *)
CreateUnitMatrix[type_] :=
    Block[{},
          Print["Error: CreateUnitMatrix: can't create unity matrix for type: ", type];
          Quit[];
         ];

CreateUnitMatrix[CConversion`ScalarType[t_]] := 1;

CreateUnitMatrix[CConversion`MatrixType[t_, rows_, rows_]] := CConversion`UNITMATRIX[rows];

MakeUnique[name_String] :=
    Module[{appendix = ""},
           While[NameQ[Evaluate[name <> appendix]] &&
                 (MemberQ[Attributes[Evaluate[name <> appendix]], Protected] ||
                  Context[Evaluate[name <> appendix]] == "SARAH`"),
                 appendix = appendix <> "x";
             ];
           Clear[Evaluate[name <> appendix]];
           Symbol[Evaluate[name <> appendix]]
          ];

MakeUniqueStr[name_String] :=
    ToString[MakeUnique[name]];

ConvertGreekLetters[text_] :=
   Symbol[StringReplace[ToString[text], {
       (* replace greek symbol by uniqe greek symbol string only if
          the symbol appears alone *)
       StartOfString ~~ "\[Alpha]" ~~ EndOfString        -> MakeUniqueStr["Alpha"],
       StartOfString ~~ "\[Beta]" ~~ EndOfString         -> MakeUniqueStr["Beta"],
       StartOfString ~~ "\[Gamma]" ~~ EndOfString        -> MakeUniqueStr["Gamma"],
       StartOfString ~~ "\[Delta]" ~~ EndOfString        -> MakeUniqueStr["Delta"],
       StartOfString ~~ "\[Epsilon]" ~~ EndOfString      -> MakeUniqueStr["Epsilon"],
       StartOfString ~~ "\[CurlyEpsilon]" ~~ EndOfString -> MakeUniqueStr["CurlyEpsilon"],
       StartOfString ~~ "\[Zeta]" ~~ EndOfString         -> MakeUniqueStr["Zeta"],
       StartOfString ~~ "\[Eta]" ~~ EndOfString          -> MakeUniqueStr["Eta"],
       StartOfString ~~ "\[Theta]" ~~ EndOfString        -> MakeUniqueStr["Theta"],
       StartOfString ~~ "\[CurlyTheta]" ~~ EndOfString   -> MakeUniqueStr["CurlyTheta"],
       StartOfString ~~ "\[Iota]" ~~ EndOfString         -> MakeUniqueStr["Iota"],
       StartOfString ~~ "\[Kappa]" ~~ EndOfString        -> MakeUniqueStr["Kappa"],
       StartOfString ~~ "\[CurlyKappa]" ~~ EndOfString   -> MakeUniqueStr["CurlyKappa"],
       StartOfString ~~ "\[Lambda]" ~~ EndOfString       -> MakeUniqueStr["Lambda"],
       StartOfString ~~ "\[Mu]" ~~ EndOfString           -> MakeUniqueStr["Mu"],
       StartOfString ~~ "\[Nu]" ~~ EndOfString           -> MakeUniqueStr["Nu"],
       StartOfString ~~ "\[Xi]" ~~ EndOfString           -> MakeUniqueStr["Xi"],
       StartOfString ~~ "\[Omicron]" ~~ EndOfString      -> MakeUniqueStr["Omicron"],
       StartOfString ~~ "\[Pi]" ~~ EndOfString           -> MakeUniqueStr["Pi"],
       StartOfString ~~ "\[CurlyPi]" ~~ EndOfString      -> MakeUniqueStr["CurlyPi"],
       StartOfString ~~ "\[Rho]" ~~ EndOfString          -> MakeUniqueStr["Rho"],
       StartOfString ~~ "\[CurlyRho]" ~~ EndOfString     -> MakeUniqueStr["CurlyRho"],
       StartOfString ~~ "\[Sigma]" ~~ EndOfString        -> MakeUniqueStr["Sigma"],
       StartOfString ~~ "\[FinalSigma]" ~~ EndOfString   -> MakeUniqueStr["FinalSigma"],
       StartOfString ~~ "\[Tau]" ~~ EndOfString          -> MakeUniqueStr["Tau"],
       StartOfString ~~ "\[Upsilon]" ~~ EndOfString      -> MakeUniqueStr["Upsilon"],
       StartOfString ~~ "\[Phi]" ~~ EndOfString          -> MakeUniqueStr["Phi"],
       StartOfString ~~ "\[CurlyPhi]" ~~ EndOfString     -> MakeUniqueStr["CurlyPhi"],
       StartOfString ~~ "\[Chi]" ~~ EndOfString          -> MakeUniqueStr["Chi"],
       StartOfString ~~ "\[Psi]" ~~ EndOfString          -> MakeUniqueStr["Psi"],
       StartOfString ~~ "\[Omega]" ~~ EndOfString        -> MakeUniqueStr["Omega"],
       StartOfString ~~ "\[Digamma]" ~~ EndOfString      -> MakeUniqueStr["Digamma"],
       StartOfString ~~ "\[Koppa]" ~~ EndOfString        -> MakeUniqueStr["Koppa"],
       StartOfString ~~ "\[Stigma]" ~~ EndOfString       -> MakeUniqueStr["Stigma"],
       StartOfString ~~ "\[Sampi]" ~~ EndOfString        -> MakeUniqueStr["Sampi"],
       (* otherwise, i.e. if the symbol appears in combination with others,
          we can replace it by it's proper name *)
       "\[Alpha]"        -> "Alpha",
       "\[Beta]"         -> "Beta",
       "\[Gamma]"        -> "Gamma",
       "\[Delta]"        -> "Delta",
       "\[Epsilon]"      -> "Epsilon",
       "\[CurlyEpsilon]" -> "CurlyEpsilon",
       "\[Zeta]"         -> "Zeta",
       "\[Eta]"          -> "Eta",
       "\[Theta]"        -> "Theta",
       "\[CurlyTheta]"   -> "CurlyTheta",
       "\[Iota]"         -> "Iota",
       "\[Kappa]"        -> "Kappa",
       "\[CurlyKappa]"   -> "CurlyKappa",
       "\[Lambda]"       -> "Lambda",
       "\[Mu]"           -> "Mu",
       "\[Nu]"           -> "Nu",
       "\[Xi]"           -> "Xi",
       "\[Omicron]"      -> "Omicron",
       "\[Pi]"           -> "Pi",
       "\[CurlyPi]"      -> "CurlyPi",
       "\[Rho]"          -> "Rho",
       "\[CurlyRho]"     -> "CurlyRho",
       "\[Sigma]"        -> "Sigma",
       "\[FinalSigma]"   -> "FinalSigma",
       "\[Tau]"          -> "Tau",
       "\[Upsilon]"      -> "Upsilon",
       "\[Phi]"          -> "Phi",
       "\[CurlyPhi]"     -> "CurlyPhi",
       "\[Chi]"          -> "Chi",
       "\[Psi]"          -> "Psi",
       "\[Omega]"        -> "Omega",
       "\[Digamma]"      -> "Digamma",
       "\[Koppa]"        -> "Koppa",
       "\[Stigma]"       -> "Stigma",
       "\[Sampi]"        -> "Sampi"
    }]];

ToValidCSymbol[symbol_Symbol] := ConvertGreekLetters[symbol];

ToValidCSymbol[symbol_Integer] := symbol;

ToValidCSymbol[symbol_Real] := symbol;

ToValidCSymbol[symbol_[i1,i2]] := ToValidCSymbol[symbol];

ToValidCSymbol[symbol_ /; Length[symbol] > 0] :=
    Module[{result = "", i},
           For[i = 0, i <= Length[symbol], i++,
               result = result <> ToString[ToValidCSymbol[symbol[[i]]]];
              ];
           Return[Symbol[result]];
          ];

(* creates a valid C parameter name string by converting the symbol to
   a valid C variable name and removing matrix indices *)
ToValidCSymbolString[symbol_] :=
    ToString[ToValidCSymbol[symbol]];

Format[Susyno`LieGroups`conj[x_],CForm] :=
    If[SARAH`getDimParameters[x] === {},
       Format["Conj(" <> ToString[CForm[x]] <> ")", OutputForm],
       Format[ToString[CForm[x]] <> ".conjugate()", OutputForm]
      ];

Format[SARAH`Conj[x_],CForm]            :=
    If[SARAH`getDimParameters[x] === {} || SARAH`getDimParameters[x] === {0},
       Format["Conj(" <> ToString[CForm[x]] <> ")", OutputForm],
       Format[ToString[CForm[x]] <> ".conjugate()", OutputForm]
      ];

Format[SARAH`Adj[x_Symbol],CForm]       := Format[ToString[CForm[x]] <> ".adjoint()"  , OutputForm];

Format[SARAH`Adj[x_],CForm]             :=
    Format["(" <> ToString[CForm[x]] <> ").adjoint()"  , OutputForm];

Format[SARAH`Tp[x_Symbol],CForm]        := Format[ToString[CForm[x]] <> ".transpose()", OutputForm];

Format[SARAH`Tp[x_],CForm]              :=
    Format["(" <> ToString[CForm[x]] <> ").transpose()", OutputForm];

Format[SARAH`trace[HoldPattern[x_Symbol]],CForm] :=
    Format[ToString[CForm[HoldForm[x]]] <> ".trace()", OutputForm];

Format[SARAH`trace[HoldPattern[x_]],CForm] :=
    Format["(" <> ToString[CForm[HoldForm[x]]] <> ").trace()", OutputForm];

(* Converts an expression to CForm and expands SARAH symbols
 *
 *   MatMul[A]      ->   A
 *   MatMul[A,B]    ->   A * B
 *   MatMul[A,B,C]  ->   A * B * C
 *   MatMul[A,B,A]  ->   A * B * A
 *   trace[A,B]     ->   trace[A * B]
 *   trace[A,B,A]   ->   trace[A * B * A]
 *
 * etc.
 *)
RValueToCFormString[expr_] :=
    Module[{times, result},
           result = expr //. {
                    SARAH`A0[b_]             :> SARAH`A0[FlexibleSUSY`Mass[b]]        /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`B0[a__ , b_, c___] :> SARAH`B0[a , FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`B1[a__ , b_, c___] :> SARAH`B1[a , FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`B00[a__, b_, c___] :> SARAH`B00[a, FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`B22[a__, b_, c___] :> SARAH`B22[a, FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`F0[a__ , b_, c___] :> SARAH`F0[a , FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`G0[a__ , b_, c___] :> SARAH`G0[a , FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass,
                    SARAH`H0[a__ , b_, c___] :> SARAH`H0[a , FlexibleSUSY`Mass[b], c] /; Head[b] =!= FlexibleSUSY`Mass } /.
                    SARAH`Mass2[a_?NumberQ]  :> Global`Sqr[a] /.
                    SARAH`Mass2[a_]          :> Global`Sqr[FlexibleSUSY`Mass[a]] /.
                    FlexibleSUSY`Mass[a_?NumberQ]   :> a /.
                    FlexibleSUSY`Mass[bar[a_]]      :> FlexibleSUSY`Mass[a] /.
                    FlexibleSUSY`Mass[a_[idx_]]     :> ToValidCSymbol[FlexibleSUSY`Mass[a]][idx] /.
                    FlexibleSUSY`Mass[a_]           :> ToValidCSymbol[FlexibleSUSY`Mass[a]] /.
                    Susyno`LieGroups`conj    -> SARAH`Conj /.
                    SARAH`Conj[a_] a_        :> AbsSqr[a] /.
                    a_[SARAH`i1,SARAH`i2]    :> a /.
                    SARAH`Delta[a_,a_]       -> 1 /.
                    Power[a_,2]              :> Global`Sqr[a] /.
                    Power[E,a_]              :> exp[a];
           result = Apply[Function[code, Hold[CForm[code]], HoldAll],
                          Hold[#] &[result /. { SARAH`MatMul[a__] :> times @@ SARAH`MatMul[a],
                                                SARAH`trace[a__]  :> SARAH`trace[times[a]] }]
                          /. times -> Times
                         ];
           ToString[HoldForm @@ result]
          ];

(* returns the head of a symbol
 * GetHead[s]        ->  s
 * GetHead[s[a]]     ->  s
 * GetHead[s[a][b]]  ->  s
 *)
GetHead[sym_] :=
    Module[{result},
           result = sym;
           If[IntegerQ[result] || RealQ[result], Return[result]];
           While[Head[result] =!= Symbol, result = Head[result]];
           Return[result];
          ];

(* this variable is increased during each call of
   CreateUniqueCVariable[] *)
sumVariableCounter = 0;

CreateUniqueCVariable[] :=
    Module[{variable},
           variable = "tmp_" <> ToString[sumVariableCounter];
           sumVariableCounter++;
           Return[variable];
          ];

ExpandSums[sum[index_, start_, stop_, expr_], variable_String, type_String:"Complex", initialValue_String:""] :=
    Module[{result, tmpSum, idxStr, startStr, stopStr},
           idxStr   = ToValidCSymbolString[index];
           startStr = ToValidCSymbolString[start];
           stopStr  = ToValidCSymbolString[stop];
           tmpSum   = CreateUniqueCVariable[];
           result = type <> " " <> tmpSum <> initialValue <> ";\n" <>
                    "for (unsigned " <> idxStr <> " = " <>
                    startStr <> "; " <> idxStr <> " <= " <> stopStr <>
                    "; ++" <> idxStr <> ") {\n" <>
                    IndentText[ExpandSums[expr,tmpSum,type,initialValue]] <> "}\n" <>
                    variable <> " += " <> tmpSum <> ";\n";
           Return[result];
          ];

ExpandSums[expr_Plus, variable_String, type_String:"Complex", initialValue_String:""] :=
    Module[{summands},
           summands = List @@ expr;
           StringJoin[ExpandSums[#,variable,type,initialValue]& /@ summands]
          ];

ToCondition[SARAH`ThetaStep[i1_,i2_]] := ToString[i1] <> " <= " <> ToString[i2];

StripThetaStep[expr_] :=
      Module[{thetas, strippedExpr, i, condition = ""},
           thetas = Cases[ThetaMark expr, SARAH`ThetaStep[_,_], Infinity];
           strippedExpr = DeleteCases[ThetaMark expr, SARAH`ThetaStep[_,_], Infinity] /. ThetaMark -> 1;
           (* create condition *)
           For[i = 1, i <= Length[thetas], i++,
               If[i > 1, condition = condition <> " && ";];
               condition = condition <> ToCondition[thetas[[i]]];
              ];
           Return[{strippedExpr, condition}];
          ];

ExpandSums[expr_Times /; !FreeQ[expr,SARAH`ThetaStep], variable_String,
           type_String:"Complex", initialValue_String:""] :=
    Module[{strippedExpr, condition, result, expandedExpr},
           expandedExpr = Expand[expr];
           If[expandedExpr === expr,
              {strippedExpr, condition} = StripThetaStep[expr];
              result = "if (" <> condition <> ") {\n" <>
                       IndentText[ExpandSums[strippedExpr, variable, type, initialValue]] <>
                       "}\n";
              ,
              result = ExpandSums[expandedExpr, variable, type, initialValue];
             ];
           Return[result];
          ];

ExpandSums[expr_Times, variable_String, type_String:"Complex", initialValue_String:""] :=
    Module[{factors, sums, rest, expandedSums, sumProduct, result = "", i},
           factors = List @@ expr;
           sums = Select[factors, (!FreeQ[#,sum[__]])&];
           rest = Complement[factors, sums];
           expandedSums = ({#, CreateUniqueCVariable[]})& /@ sums;
           expandedSums = ({ExpandSums[#[[1]], #[[2]], type, initialValue], #[[2]]})& /@ expandedSums;
           (* add for loops *)
           result = StringJoin[(type <> " " <> #[[2]] <> initialValue <> ";\n" <> #[[1]])& /@ expandedSums];
           result = result <> variable <> " += ";
           (* multiply the rest *)
           For[i = 1, i <= Length[rest], i++,
               If[i > 1, result = result <> " * ";];
               If[Head[rest[[i]]] === Plus,
                  result = result <> "(" <> RValueToCFormString[rest[[i]]] <> ")";
                  ,
                  result = result <> RValueToCFormString[rest[[i]]];
                 ];
              ];
           (* multiply the sums *)
           For[i = 1, i <= Length[expandedSums], i++,
               If[Length[rest] > 0 || i > 1, result = result <> " * ";];
               result = result <> expandedSums[[i,2]];
              ];
           result = result <>";\n";
           Return[result];
          ];

ExpandSums[expr_?((!FreeQ[#,SARAH`ThetaStep])&), variable_String, type_String:"Complex", initialValue_String:""] :=
    Module[{strippedExpr, i, condition = "", result},
           {strippedExpr, condition} = StripThetaStep[expr];
           result = "if (" <> condition <> ") {\n" <>
                    IndentText[variable <> " += " <>
                               RValueToCFormString[strippedExpr] <> ";\n"] <>
                    "}\n";
           Return[result];
          ];

ExpandSums[expr_, variable_String, type_String:"Complex", initialValue_String:""] :=
    variable <> " += " <> RValueToCFormString[expr] <> ";\n";

End[];

EndPackage[];
