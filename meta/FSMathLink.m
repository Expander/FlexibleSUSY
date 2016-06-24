BeginPackage["FSMathLink`", {"CConversion`", "Parameters`", "TreeMasses`", "Utils`"}];

GetNumberOfInputParameterRules::usage = "";
PutInputParameters::usage = "";

Begin["`Private`"];

GetNumberOfInputParameterRules[inputPars_List] :=
    Length[inputPars];

PutInputParameter[{par_, _, CConversion`ScalarType[CConversion`complexScalarCType]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleToComplex(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameter[{par_, _, CConversion`ScalarType[CConversion`integerScalarCType]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleToInteger(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameter[{par_, _, CConversion`ScalarType[CConversion`realScalarCType]}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleToReal(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameters[inputPars_List, linkName_String] :=
    StringJoin[PutInputParameter[#, linkName]& /@ inputPars];

End[];

EndPackage[];
