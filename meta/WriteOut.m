
BeginPackage["WriteOut`", {"CConversion`"}];

PrintParameters::usage="Creates parameter printout statements";

Begin["Private`"];

PrintParameter[Null, streamName_String] := "";

PrintParameter[parameter_, streamName_String] :=
    Module[{parameterName},
           parameterName = ToValidCSymbolString[parameter /. a_[i1,i2] :> a];
           Return[streamName <> " << \"" <> parameterName <> " = \" << " <>
                  parameterName <> " << '\\n';\n"];
          ];

PrintParameters[parameters_List, streamName_String] :=
    Module[{result = ""},
           (result = result <> PrintParameter[#,streamName])& /@ parameters;
           Return[result];
          ];

End[];

EndPackage[];
