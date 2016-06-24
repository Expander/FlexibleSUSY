BeginPackage["FSMathLink`", {"CConversion`", "Parameters`", "TreeMasses`", "Utils`"}];

GetNumberOfInputParameterRules::usage = "";
PutInputParameters::usage = "";

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

End[];

EndPackage[];
