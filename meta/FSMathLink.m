BeginPackage["FSMathLink`", {"CConversion`", "Parameters`", "TreeMasses`", "Utils`"}];

GetNumberOfInputParameterRules::usage = "";
GetNumberOfSpectrumEntries::usage = "";
PutInputParameters::usage = "";
SetInputParametersFromArguments::usage = "";
SetInputParameterDefaultArguments::usage = "";
SetInputParameterArgumentTypes::usage = "";
SetInputParameterArgumentCTypes::usage = "";
SetInputParameterArguments::usage = "";
PutSpectrum::usage = "";

Begin["`Private`"];

GetNumberOfInputParameterRules[inputPars_List] :=
    Length[inputPars];

PutInputParameter[{par_, _}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameters[inputPars_List, linkName_String] :=
    StringJoin[PutInputParameter[#, linkName]& /@ inputPars];

CreateComponent[parStr_String, CConversion`integerScalarCType] := parStr;
CreateComponent[parStr_String, CConversion`realScalarCType] := parStr;
CreateComponent[parStr_String, CConversion`complexScalarCType] := "std::complex<double>(Re" <> parStr <> ", Im" <> parStr <> ")";

SetInputParameterFromArguments[{par_, CConversion`ScalarType[st_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "INPUTPARAMETER(" <> parStr <> ") = " <> CreateComponent[parStr, st] <> ";\n"
          ];

SetInputParameterFromArguments[{par_, (CConversion`ArrayType | CConversion`VectorType)[st_,dim_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[("INPUTPARAMETER(" <> parStr <> "(" <> ToString[#-1] <> ")) = " <>
                       CreateComponent[parStr <> "_" <> ToString[#], st] <> ";\n")& /@ Table[i, {i,1,dim}]]
          ];

SetInputParameterFromArguments[{par_, CConversion`MatrixType[st_,dim1_,dim2_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[Flatten[Outer[("INPUTPARAMETER(" <> parStr <> "(" <> ToString[#1-1] <> "," <> ToString[#2-1] <> ")) = " <>
                                     CreateComponent[parStr <> "_" <> ToString[#1] <> ToString[#2], st] <> ";\n")&,
                                    Table[i, {i, 1, dim1}],
                                    Table[j, {j, 1, dim2}]], 1]]
          ];

SetInputParameterFromArguments[{par_, CConversion`TensorType[st_,dim1_,dim2_,dim3_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[Flatten[Outer[("INPUTPARAMETER(" <> parStr <> "(" <> ToString[#1-1] <> "," <> ToString[#2-1] <> "," <> ToString[#3-1] <> ")) = " <>
                                     CreateComponent[parStr <> "_" <> ToString[#1] <> ToString[#2] <> ToString[#3], st] <> ";\n")&,
                                    Table[i, {i, 1, dim1}],
                                    Table[j, {j, 1, dim2}],
                                    Table[k, {k, 1, dim3}]], 2]]
          ];

SetInputParameterFromArguments[{par_, CConversion`TensorType[st_,dim1_,dim2_,dim3_,dim4_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[Flatten[Outer[("INPUTPARAMETER(" <> parStr <> "(" <> ToString[#1-1] <> "," <> ToString[#2-1] <> "," <> ToString[#3-1] <> "," <> ToString[#4-1] <> ")) = " <>
                                     CreateComponent[parStr <> "_" <> ToString[#1] <> ToString[#2] <> ToString[#3] <> ToString[#4], st] <> ";\n")&,
                                    Table[i, {i, 1, dim1}],
                                    Table[j, {j, 1, dim2}],
                                    Table[k, {k, 1, dim3}],
                                    Table[l, {l, 1, dim4}]], 3]]
          ];

SetInputParametersFromArguments[inputPars_List] :=
    StringJoin[SetInputParameterFromArguments /@ inputPars];

SetInputParameterDefaultArgument[{par_, CConversion`ScalarType[_]}] :=
    CConversion`ToValidCSymbol[par] -> 0;

SetInputParameterDefaultArgument[{par_, (CConversion`ArrayType | CConversion`VectorType)[_,dim_]}] :=
    CConversion`ToValidCSymbol[par] -> Table[0, {dim}];

SetInputParameterDefaultArgument[{par_, CConversion`MatrixType[_,dim1_,dim2_]}] :=
    CConversion`ToValidCSymbol[par] -> Table[0, {dim1}, {dim2}];

SetInputParameterDefaultArgument[{par_, CConversion`TensorType[_,dims__]}] :=
    CConversion`ToValidCSymbol[par] -> Table[0, Evaluate[Sequence @@ ({#}& /@ {dims})]];

SetInputParameterDefaultArguments[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[SetInputParameterDefaultArgument[#]]& /@ inputPars, ",\n"];

ConcatIndices[idx__] := StringJoin[ToString /@ {idx}];

ParAndType[par_, CConversion`realScalarCType   ] := {HoldForm[OptionValue[par]], Real, "double " <> CConversion`ToValidCSymbolString[par]};
ParAndType[par_, CConversion`complexScalarCType] := Sequence[{HoldForm[Re[OptionValue[par]]], Real, "double Re" <> CConversion`ToValidCSymbolString[par]},
                                                             {HoldForm[Im[OptionValue[par]]], Real, "double Im" <> CConversion`ToValidCSymbolString[par]}];
ParAndType[par_, CConversion`integerScalarCType] := {HoldForm[OptionValue[par]], Integer, "int " <> CConversion`ToValidCSymbolString[par]};

ParAndType[par_, CConversion`realScalarCType, idx__   ] := {HoldForm[OptionValue[par][[idx]]], Real,
                                                            "double " <> CConversion`ToValidCSymbolString[par] <> "_" <> ConcatIndices[idx]};
ParAndType[par_, CConversion`complexScalarCType, idx__] := Sequence[{HoldForm[Re[OptionValue[par][[idx]]]], Real,
                                                                     "double Re" <> CConversion`ToValidCSymbolString[par] <> "_" <> ConcatIndices[idx]},
                                                                    {HoldForm[Im[OptionValue[par][[idx]]]], Real,
                                                                     "double Im" <> CConversion`ToValidCSymbolString[par] <> "_" <> ConcatIndices[idx]}];
ParAndType[par_, CConversion`integerScalarCType, idx__] := {HoldForm[OptionValue[par][[idx]]], Integer,
                                                            "int " <> CConversion`ToValidCSymbolString[par] <> "_" <> ConcatIndices[idx]};

SetInputParameterArgumentsAndType[{par_, CConversion`ScalarType[st_]}] :=
    {ParAndType[CConversion`ToValidCSymbol[par], st]};

SetInputParameterArgumentsAndType[{par_, (CConversion`ArrayType | CConversion`VectorType)[st_, dim_]}] :=
    Table[ParAndType[CConversion`ToValidCSymbol[par], st, i], {i, 1, dim}];

SetInputParameterArgumentsAndType[{par_, CConversion`MatrixType[st_, dim1_, dim2_]}] :=
    Flatten[Outer[ParAndType[CConversion`ToValidCSymbol[par], st, #1, #2] &,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}]], 1];

SetInputParameterArgumentsAndType[{par_, CConversion`TensorType[st_, dim1_, dim2_, dim3_]}] :=
    Flatten[Outer[ParAndType[CConversion`ToValidCSymbol[par], st, #1, #2, #] &,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}],
                  Table[k, {k, 1, dim3}]], 2];

SetInputParameterArgumentsAndType[{par_, CConversion`TensorType[st_, dim1_, dim2_, dim3_, dim4_]}] :=
    Flatten[Outer[ParAndType[CConversion`ToValidCSymbol[par], st, #1, #2, #3, #4] &,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}],
                  Table[k, {k, 1, dim3}],
                  Table[l, {l, 1, dim4}]], 3];

SetInputParameterArgumentsAndTypes[inputPars_List] :=
    Join @@ SetInputParameterArgumentsAndType /@ inputPars;

SetInputParameterArguments[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[#[[1]]]& /@ SetInputParameterArgumentsAndTypes[inputPars], ",\n"];

SetInputParameterArgumentTypes[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[#[[2]]]& /@ SetInputParameterArgumentsAndTypes[inputPars], ",\n"];

SetInputParameterArgumentCTypes[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[#[[3]]]& /@ SetInputParameterArgumentsAndTypes[inputPars], ",\n"];

GetNumberOfSpectrumEntries[pars_List] :=
    Length[pars];

(* returns all heads of a nested expression of the form f[g[h[x]]] -> {f,g,h} *)
GetHeads[h_[p___]] := Join[{h}, GetHeads[p]];
GetHeads[p___] := {};

HeadStr[par_] :=
    Module[{heads = ToString /@ GetHeads[par]},
           If[heads === {}, "",
              ", {\"" <> Utils`StringJoinWithSeparator[heads, "\", \""] <> "\"}"
             ]
          ];

ToVaildOutputParStr[FlexibleSUSY`Pole[par_]] := ToVaildOutputParStr[par]; (* Pole[x] is not a valid parameter name *)
ToVaildOutputParStr[p:FlexibleSUSY`M[_]] := CConversion`ToValidCSymbolString[p];
ToVaildOutputParStr[FlexibleSUSY`SCALE] := "scale";
ToVaildOutputParStr[par_] := CConversion`ToValidCSymbolString[par];

WithoutHeads[_[par_]] := WithoutHeads[par];
WithoutHeads[par_] := par;

GetMacroFor[FlexibleSUSY`Pole[par_]] := "PHYSICALPARAMETER";
GetMacroFor[par_] := "MODELPARAMETER";

PutParameter[par_, _, link_String] :=
    Module[{parStr = ToVaildOutputParStr[par],
            parWithoutHeads = CConversion`ToValidCSymbolString[WithoutHeads[par]]},
           "MLPutRuleTo(" <> link <> ", " <> GetMacroFor[par] <> "(" <> parStr <>
           "), \"" <> parWithoutHeads <> "\"" <> HeadStr[par] <> ");\n"
          ];

PutParameter[par_, link_String] :=
    PutParameter[par, Parameters`GetType[par /. FlexibleSUSY`Pole -> Identity], link];

PutSpectrum[pars_List, link_String] :=
    StringJoin[PutParameter[#,link]& /@ pars];

End[];

EndPackage[];
