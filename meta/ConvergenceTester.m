
BeginPackage["ConvergenceTester`", {"CConversion`", "TextFormatting`", "TreeMasses`"}];

CreateCompareFunction::usage="";

Begin["Private`"];

CountNumberOfMasses[particles_List] :=
    Plus @@ (TreeMasses`GetDimension /@ particles);

CountNumberOfMasses[particle_] :=
    TreeMasses`GetDimension[particle];

CalcDifference[particle_, offset_Integer, diff_String] :=
    Module[{result, body, dim, esStr, comma},
           dim = TreeMasses`GetDimension[particle];
           esStr = ToValidCSymbolString[particle];
           If[dim == 1,
              comma = "OLD1(" <> esStr <> "),NEW1(" <> esStr <> ")";
              body = diff <> "(" <> ToString[1 + offset] <> ") = " <>
                     "std::fabs(1.0 - std::min(" <> comma <> ")/std::max(" <> comma <> "));";
              result = body <> "\n";
              ,
              comma = "OLD(" <> esStr <> ",i),NEW(" <> esStr <> ",i)";
              result = "for (unsigned i = 1; i <= " <> ToString[dim] <> "; ++i) {\n";
              body = diff <> "(i + " <> ToString[offset] <> ") = " <>
                     "std::fabs(1.0 - std::min(" <> comma <> ")/std::max(" <> comma <> "));";
              result = result <> IndentText[body] <> "\n}\n";
             ];
           Return[result];
          ];

CreateCompareFunction[particles_List] :=
    Module[{result, numberOfMasses, i, offset = 0, massiveParticles},
           massiveParticles = Select[particles, (!TreeMasses`IsMassless[#])&];
           numberOfMasses = CountNumberOfMasses[massiveParticles];
           result = "DoubleVector diff(" <> ToString[numberOfMasses] <> ");\n\n";
           For[i = 1, i <= Length[massiveParticles], i++,
               result = result <> CalcDifference[massiveParticles[[i]], offset, "diff"];
               offset += CountNumberOfMasses[massiveParticles[[i]]];
              ];
           If[offset != numberOfMasses,
              Print["Error: something is wrong with the counting of masses:"];
              Print["  numberOfMasses = ", numberOfMasses, ", offset = ", offset];
             ];
           result = result <> "\nreturn diff.max();\n";
           Return[result];
          ];

End[];

EndPackage[];
