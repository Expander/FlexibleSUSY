
BeginPackage["Tadpoles`", {"SARAH`", "TextFormatting`", "CConversion`"}];

CreateTadpoleEqPrototype::usage="creates C function prototype for a
given tadpole equation";

CreateTadpoleEqFunction::usage="creates C function definition for a
given tadpole equation";

FillArrayWithTadpoles::usage="fills array of doubles with the values
of the tadpoles";

FillInitialGuessArray::usage="fills a C array with initial values for the
EWSB eqs. solver";

Begin["Private`"];

CreateTadpoleEqPrototype[vev_Symbol] :=
    Module[{result = ""},
           result = "double get_tadpole_" <> ToValidCSymbolString[vev] <>
                    "() const;\n";
           Return[result];
          ];

CreateTadpoleEqFunction[vev_Symbol, equation_] :=
    Module[{result = "", body = "", variableName = ""},
           variableName = "tadpole_" <> ToValidCSymbolString[vev];
           result = "double CLASSNAME::get_" <> variableName <>
                    "() const\n{\n";
           body = "double " <> variableName <> " = " <>
                  RValueToCFormString[equation] <> ";\n";
           body = body <> "\nreturn " <> variableName <> ";\n";
           body = IndentText[WrapLines[body]];
           Return[result <> body <> "}\n\n"];
          ];

FillArrayWithTadpoles[tadpoleEquations_List, parametersFixedByEWSB_List,
                      gslIntputVector_String:"x", gslOutputVector_String:"tadpole"] :=
    Module[{i, result = "", vev},
           If[Length[tadpoleEquations] =!= Length[parametersFixedByEWSB],
              Print["Error: number of tadpole equations (",Length[tadpoleEquations],
                    ") is not equal to the number of fixed parameters (",
                    Length[parametersFixedByEWSB],")"];
              Return[""];
             ];
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               result = result <> "model->set_" <>
                        ToValidCSymbolString[parametersFixedByEWSB[[i]]] <>
                        "(gsl_vector_get(" <> gslIntputVector <> ", " <>
                        ToString[i-1] <> "));\n";
              ];
           result = result <> "\n";
           For[i = 1, i <= Length[tadpoleEquations], i++,
               vev = tadpoleEquations[[i,1]];
               result = result <> gslOutputVector <> "[" <> ToString[i-1] <>
                        "] = " <> "model->get_tadpole_" <>
                        ToValidCSymbolString[vev] <> "();\n";
              ];
           Return[result];
          ];

FillInitialGuessArray[parametersFixedByEWSB_List, arrayName_String:"x_init"] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               result = result <> arrayName <> "[" <> ToString[i-1] <> "] = " <>
                        ToValidCSymbolString[parametersFixedByEWSB[[i]]] <>
                        ";\n";
              ];
           Return[result];
          ];

End[];

EndPackage[];
