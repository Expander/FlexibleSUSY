
BeginPackage["ConvergenceTester`", {"CConversion`", "TextFormatting`", "TreeMasses`"}];

CreateCompareFunction::usage="";

Begin["`Private`"];

CountNumberOfParameters[FlexibleSUSY`M[particle_]] :=
    TreeMasses`GetDimension[particle];

CountNumberOfParameters[particle_] :=
    TreeMasses`GetDimension[particle];

CountNumberOfParameters[particles_List] :=
    Plus @@ (CountNumberOfParameters /@ particles);

CalcDifference[FlexibleSUSY`M[particle_], offset_Integer, diff_String] :=
    Module[{result, body, dim, dimStart, esStr},
           dim = TreeMasses`GetDimension[particle];
           esStr = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[dim == 1,
              result = diff <> "[" <> ToString[offset] <> "] = " <>
                       "MaxRelDiff(OLD1(" <> esStr <> "),NEW1(" <> esStr <> "));\n";
              ,
              dimStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle] - 1;
              result = "for (unsigned i = " <> ToString[dimStart] <>
                       "; i < " <> ToString[dim] <> "; ++i) {\n";
              body = diff <> "[i + " <> ToString[offset] <> "] = " <>
                     "MaxRelDiff(OLD(" <> esStr <> ",i),NEW(" <> esStr <> ",i));";
              result = result <> IndentText[body] <> "\n}\n";
             ];
           Return[result];
          ];

CreateCompareFunction[crit_ /; crit === Automatic] :=
    Module[{particles},
           If[SARAH`SupersymmetricModel,
              particles = TreeMasses`GetSusyParticles[];,
              particles = TreeMasses`GetParticles[];
             ];
           particles = Select[particles, (!TreeMasses`IsMassless[#] && !IsGhost[#] && !IsGoldstone[#])&];
           particles = FlexibleSUSY`M /@ particles;
           CreateCompareFunction[particles]
          ];

CreateCompareFunction[parameters_List] :=
    Module[{result, numberOfParameters, i, offset = 0},
           numberOfParameters = CountNumberOfParameters[parameters];
           If[numberOfParameters == 0,
              Print["Error: no parameters specified for the convergence test!"];
              Return["return 0.;"];
             ];
           result = "double diff[" <> ToString[numberOfParameters] <> "] = { 0 };\n\n";
           For[i = 1, i <= Length[parameters], i++,
               result = result <> CalcDifference[parameters[[i]], offset, "diff"];
               offset += CountNumberOfParameters[parameters[[i]]];
              ];
           If[offset != numberOfParameters,
              Print["Error: something is wrong with the counting of masses:"];
              Print["  numberOfParameters = ", numberOfParameters, ", offset = ", offset];
             ];
           result = result <>
                    "\nreturn *std::max_element(diff, diff + " <>
                    ToString[numberOfParameters] <> ");\n";
           Return[result];
          ];

End[];

EndPackage[];
