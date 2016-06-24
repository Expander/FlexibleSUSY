BeginPackage["FSMathLink`", {"CConversion`", "Parameters`", "TreeMasses`", "Utils`"}];

GetNumberOfInputParameterRules::usage = "";
PutInputParameters::usage = "";
SetInputParameterDefaultArguments::usage = "";
SetInputParameterArgumentTypes::usage = "";

Begin["`Private`"];

GetNumberOfInputParameterRules[inputPars_List] :=
    Length[inputPars];

ScalarTypeName[st_] :=
    Switch[st,
           CConversion`complexScalarCType, "Complex",
           CConversion`realScalarCType, "Real",
           CConversion`integerScalarCType, "Integer"
          ];

PutInputParameter[{par_, _, CConversion`ScalarType[st_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameter[{par_, _, CConversion`ArrayType[st_,dim_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenArray(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <> ToString[dim] <> ");\n"
          ];

PutInputParameter[{par_, _, CConversion`VectorType[st_,dim_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenVector(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <> ToString[dim] <> ");\n"
          ];

PutInputParameter[{par_, _, CConversion`MatrixType[st_,dim1_,dim2_]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenMatrix(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <>
           ToString[dim1] <> ", " <> ToString[dim2] <> ");\n"
          ];

PutInputParameter[{par_, _, CConversion`TensorType[st_,dims__]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo" <> ScalarTypeName[st] <> "EigenTensor(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\", " <>
           Utils`StringJoinWithSeparator[ToString /@ {dims}, ", "] <> ");\n"
          ];

PutInputParameters[inputPars_List, linkName_String] :=
    StringJoin[PutInputParameter[#, linkName]& /@ inputPars];

SetInputParameterDefaultArgument[{par_, _, CConversion`ScalarType[_]}] :=
    CConversion`ToValidCSymbol[par] -> 0;

SetInputParameterDefaultArgument[{par_, _, (CConversion`ArrayType | CConversion`VectorType)[_,dim_]}] :=
    CConversion`ToValidCSymbol[par] -> Table[0, {dim}];

SetInputParameterDefaultArgument[{par_, _, CConversion`MatrixType[_,dim1_,dim2_]}] :=
    CConversion`ToValidCSymbol[par] -> Table[0, {dim1}, {dim2}];

SetInputParameterDefaultArgument[{par_, _, CConversion`TensorType[_,dims__]}] :=
    CConversion`ToValidCSymbol[par] -> Table[0, Evaluate[Sequence @@ ({#}& /@ {dims})]];

SetInputParameterDefaultArguments[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[SetInputParameterDefaultArgument[#]]& /@ inputPars, ",\n"];

ParAndType[par_, CConversion`realScalarCType] := {par, Real};
ParAndType[par_, CConversion`complexScalarCType] := Sequence[{Re[par], Real}, {Im[par], Real}];
ParAndType[par_, CConversion`integerScalarCType] := {par, Integer};

SetInputParameterArgumentsAndType[{par_, _, CConversion`ScalarType[st_]}] := {ParAndType[par, st]};

SetInputParameterArgumentsAndType[{par_, _, (CConversion`ArrayType | CConversion`VectorType)[st_, dim_]}] :=
    Table[ParAndType[par[i], st], {i, 1, dim}];

SetInputParameterArgumentsAndType[{par_, _, CConversion`MatrixType[st_, dim1_, dim2_]}] :=
    Flatten[Outer[ParAndType[par[#1, #2], st]&,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}]], 1];

SetInputParameterArgumentsAndType[{par_, _, CConversion`TensorType[st_, dim1_, dim2_, dim3_]}] :=
    Flatten[Outer[ParAndType[par[#1, #2, #3], st]&,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}],
                  Table[k, {k, 1, dim3}]], 2];

SetInputParameterArgumentsAndType[{par_, _, CConversion`TensorType[st_, dim1_, dim2_, dim3_, dim4_]}] :=
    Flatten[Outer[ParAndType[par[#1, #2, #3, #4], st]&,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}],
                  Table[k, {k, 1, dim3}],
                  Table[l, {l, 1, dim4}]], 3];

SetInputParameterArgumentsAndTypes[inputPars_List] :=
    Join @@ SetInputParameterArgumentsAndType /@ inputPars;

SetInputParameterArgumentTypes[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[#[[2]]]& /@ SetInputParameterArgumentsAndTypes[inputPars], ",\n"];

End[];

EndPackage[];
