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

ScalarTypeName[st_] :=
    Switch[st,
           CConversion`complexScalarCType, "Complex",
           CConversion`realScalarCType, "Real",
           CConversion`integerScalarCType, "Integer"
          ];

PutInputParameter[{par_, CConversion`ScalarType[st_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameter[{par_, CConversion`ArrayType[st_,dim_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenArray(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <> ToString[dim] <> ");\n"
          ];

PutInputParameter[{par_, CConversion`VectorType[st_,dim_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenVector(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <> ToString[dim] <> ");\n"
          ];

PutInputParameter[{par_, CConversion`MatrixType[st_,dim1_,dim2_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenMatrix(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <>
           ToString[dim1] <> ", " <> ToString[dim2] <> ");\n"
          ];

PutInputParameter[{par_, CConversion`TensorType[st_,dims__]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenTensor(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <>
           Utils`StringJoinWithSeparator[ToString /@ {dims}, ", "] <> ");\n"
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

HeadOpt[FlexibleSUSY`Pole[par_]] := "Head";
HeadOpt[par_] := "";

HeadStr[FlexibleSUSY`Pole[par_]] := ", \"Pole\"";
HeadStr[par_] := "";

ToVaildOutputParStr[FlexibleSUSY`Pole[par_]] := CConversion`ToValidCSymbolString[par];
ToVaildOutputParStr[par_] := CConversion`ToValidCSymbolString[par];

WrapMLAround[par_, CConversion`ScalarType[st_], str_String] :=
    "MLPutRuleTo" <> HeadOpt[par] <> ScalarTypeName[st] <> "(" <> str <> HeadStr[par] <> ");\n";

WrapMLAround[par_, CConversion`ArrayType[st_,__], str_String] :=
    "MLPutRuleTo" <> HeadOpt[par] <> ScalarTypeName[st] <> "EigenArray(" <> str <> HeadStr[par] <> ");\n";

WrapMLAround[par_, CConversion`VectorType[st_,__], str_String] :=
    "MLPutRuleTo" <> HeadOpt[par] <> ScalarTypeName[st] <> "EigenVector(" <> str <> HeadStr[par] <> ");\n";

WrapMLAround[par_, CConversion`MatrixType[st_,__], str_String] :=
    "MLPutRuleTo" <> HeadOpt[par] <> ScalarTypeName[st] <> "EigenMatrix(" <> str <> HeadStr[par] <> ");\n";

WrapMLAround[par_, CConversion`TensorType[st_,__], str_String] :=
    "MLPutRuleTo" <> HeadOpt[par] <> ScalarTypeName[st] <> "EigenTensor(" <> str <> HeadStr[par] <> ");\n";

GetMacroFor[FlexibleSUSY`Pole[par_]] := "PHYSICALPARAMETER";
GetMacroFor[par_] := "MODELPARAMETER";

PutParameter[par_, t:CConversion`ScalarType[st_], link_String] :=
    Module[{parStr = ToVaildOutputParStr[par]},
           WrapMLAround[par, t, link <> ", " <> GetMacroFor[par] <> "(" <> parStr <> "), \"" <> parStr <> "\""]
          ];

PutParameter[par_, t:(CConversion`ArrayType | CConversion`VectorType)[st_,dim_], link_String] :=
    Module[{parStr = ToVaildOutputParStr[par]},
           WrapMLAround[par, t, link <> ", " <> GetMacroFor[par] <> "(" <> parStr <> "), \"" <> parStr <> "\", " <> ToString[dim]]
          ];

PutParameter[par_, t:CConversion`MatrixType[st_,dim1_,dim2_], link_String] :=
    Module[{parStr = ToVaildOutputParStr[par]},
           WrapMLAround[par, t, link <> ", " <> GetMacroFor[par] <> "(" <> parStr <> "), \"" <> parStr <> "\", " <> ToString[dim1] <> ", " <> ToString[dim2]]
          ];

PutParameter[par_, t:CConversion`TensorType[st_,dims__], link_String] :=
    Module[{parStr = ToVaildOutputParStr[par]},
           WrapMLAround[par, t, link <> ", " <> GetMacroFor[par] <> "(" <> parStr <> "), \"" <> parStr <> "\", " <> Utils`StringJoinWithSeparator[ToString /@ {dims}, ", "]]
          ];

PutParameter[par_, link_String] :=
    PutParameter[par, Parameters`GetType[par /. FlexibleSUSY`Pole -> Identity], link];

PutSpectrum[pars_List, link_String] :=
    StringJoin[PutParameter[#,link]& /@ pars];

End[];

EndPackage[];
