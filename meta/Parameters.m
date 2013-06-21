
BeginPackage["Parameters`", {"SARAH`", "CConversion`"}];

CreateSetAssignment::usage="";
CreateDisplayAssignment::usage="";

ApplyGUTNormalization::usage="Returns a list of repacement rules for
gauge couplings, which replace non-normalized gauge couplings
(e.g. gY) by normalized ones (e.g. g1).";

GetGUTNormalization::usage="Returns GUT normalization of the given
coupling";

Begin["Private`"];

GetGUTNormalization[coupling_Symbol] :=
    Module[{pos, norm},
           pos = Position[SARAH`Gauge, coupling];
           If[pos =!= {},
              norm = SARAH`GUTren[pos[[1,1]]];
              If[NumericQ[norm],
                 Return[norm];
                ];
             ];
           Return[1];
          ];

ApplyGUTNormalization[] :=
    Module[{i, rules = {}, coupling},
           For[i = 1, i <= Length[SARAH`Gauge], i++,
               If[NumericQ[SARAH`GUTren[i]],
                  coupling = SARAH`Gauge[[i,4]];
                  AppendTo[rules, Rule[coupling, coupling SARAH`GUTren[i]]];
                 ];
              ];
           Return[rules];
          ];

CreateSetAssignment[name_, startIndex_, parameterType_] :=
    Block[{},
          Print["Error: CreateSetAssignment: unknown parameter type: ", ToString[parameterType]];
          Quit[];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType["double"]] :=
    Module[{ass = ""},
           ass = name <> " = v(" <> ToString[startIndex] <> ");\n";
           Return[{ass, 1}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType["Complex"]] :=
    Module[{ass = ""},
           ass = name <> " = Complex(v(" <> ToString[startIndex] <>
                 ", v(" <> ToString[startIndex + 1] <> "));\n";
           Return[{ass, 2}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[type_, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> ") = v(" <> ToString[startIndex + count] <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, rows * cols}];
          ];

CreateDisplayAssignment[name_, startIndex_, parameterType_] :=
    Block[{},
          Print["Error: CreateDisplayAssignment: unknown parameter type: ",
                ToString[parameterType]];
          Quit[];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`ScalarType["double"]] :=
    Module[{ass = ""},
           ass = "pars(" <> ToString[startIndex] <> ") = "
                 <> name <> ";\n";
           Return[{ass, 1}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`ScalarType["Complex"]] :=
    Module[{ass = ""},
           ass = "pars(" <> ToString[startIndex] <> ") = Re(" <> name <> ");\n" <>
                 "pars(" <> ToString[startIndex + 1] <> ") = Im(" <> name <> ");\n";
           Return[{ass, 2}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`MatrixType[type_, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> "pars(" <> ToString[startIndex + count] <> ") = "
                          <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                          <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, rows * cols}];
          ];

End[];

EndPackage[];
