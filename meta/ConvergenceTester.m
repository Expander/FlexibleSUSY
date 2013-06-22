
BeginPackage["ConvergenceTester`", {"CConversion`", "TextFormatting`", "TreeMasses`"}];

CreateCompareFunction::usage="";

Begin["Private`"];

CountNumberOfMasses[particles_List] :=
    Plus @@ (TreeMasses`GetDimension /@ particles);

CountNumberOfMasses[particle_] :=
    TreeMasses`GetDimension[particle];

CalcDifference[particle_, offset_Integer, diff_String] :=
    Module[{result, body, dim, dimStart, esStr, comma},
           dim = TreeMasses`GetDimension[particle];
           esStr = ToValidCSymbolString[FlexibleSUSY`Mass[particle]];
           If[dim == 1,
              comma = "OLD1(" <> esStr <> "),NEW1(" <> esStr <> ")";
              body = diff <> "(" <> ToString[1 + offset] <> ") = " <>
                     "std::fabs(1.0 - std::min(" <> comma <> ")/std::max(" <> comma <> "));";
              result = body <> "\n";
              ,
              comma = "OLD(" <> esStr <> ",i),NEW(" <> esStr <> ",i)";
              dimStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
              result = "for (unsigned i = " <> ToString[dimStart] <>
                       "; i <= " <> ToString[dim] <> "; ++i) {\n";
              body = diff <> "(i + " <> ToString[offset] <> ") = " <>
                     "std::fabs(1.0 - std::min(" <> comma <> ")/std::max(" <> comma <> "));";
              result = result <> IndentText[body] <> "\n}\n";
             ];
           Return[result];
          ];

CreateCompareFunction[particles_List] :=
    Module[{result, numberOfMasses, i, offset = 0, massiveSusyParticles},
           massiveSusyParticles = Select[particles /. FlexibleSUSY`Mass -> Identity,
                                         (!TreeMasses`IsMassless[#] && !SARAH`SMQ[#])&];
           numberOfMasses = CountNumberOfMasses[massiveSusyParticles];
           result = "DoubleVector diff(" <> ToString[numberOfMasses] <> ");\n\n";
           For[i = 1, i <= Length[massiveSusyParticles], i++,
               result = result <> CalcDifference[massiveSusyParticles[[i]], offset, "diff"];
               offset += CountNumberOfMasses[massiveSusyParticles[[i]]];
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
