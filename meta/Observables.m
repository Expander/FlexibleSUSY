BeginPackage["Observables`", {"FlexibleSUSY`", "Utils`"}];

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { aMuonGM2Calc, aMuonGM2CalcUncertainty };
End[];

CalculateObservables::usage="";

Begin["`Private`"];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    structName <> ".AMUGM2CALC = gm2calc_calculate_amu();\n";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty, structName_String] :=
    structName <> ".AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty();\n";

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}];
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

End[];

EndPackage[];
