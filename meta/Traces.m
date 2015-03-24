
BeginPackage["Traces`", {"SARAH`", "CConversion`", "Parameters`"}];

CreateTraceRules::usage="takes a list of traces and returns a list of
rules to replace the traces by their C/C++ variables.";

ConvertSARAHTraces::usage="takes SARAH's `TraceAbbr' and returns a
two-component list, where the first entry is string of C/C++ variable
definitions that hold the trace values.  The second entry is a list of
rules to replace the traces by their C/C++ variables.";

CreateTraceDefs::usage="";
CreateTraceCalculation::usage="";
CreateSARAHTraceDefs::usage="";
CreateSARAHTraceCalculation::usage="";
CreateLocalCopiesOfTraces::usage="";
CreateLocalCopiesOfSARAHTraces::usage="";
FindSARAHTraces::usage="";

SARAHTrace;

Begin["`Private`"];

GetSARAHTraceName[Traces`SARAHTrace[name_, expr_]] := name;
GetSARAHTraceExpr[Traces`SARAHTrace[name_, expr_]] := expr;

ConvertToScalar[expr_] :=
    expr /. {
        SARAH`MatMul[a_, b_] :> SARAH`ScalarProd[a, b],
        SARAH`MatMul[a_, b__] :> SARAH`ScalarProd[a, SARAH`MatMul[b]]
            };

FindSARAHTraces[expr_, sarahTraces_List] :=
    Module[{traceSymbols},
           traceSymbols = DeleteDuplicates[Flatten[GetSARAHTraceName /@ sarahTraces]];
           Select[traceSymbols, (!FreeQ[expr,#])&]
          ];

ConvertSARAHTraces[abbrs_] :=
    Module[{i, j, name, expr, traces = {}},
           For[i = 1, i <= Length[abbrs], i++,
               For[j = 1, j <= Length[abbrs[[i]]], j++,
                   {name, expr} = abbrs[[i,j]];
                   (* replace MatMul by ScalarProd *)
                   expr = ConvertToScalar[expr];
                   AppendTo[traces, Traces`SARAHTrace[name, expr]];
                  ];
              ];
           Return[traces];
          ];

GetTraceType[trace_] :=
    CConversion`ScalarType @
    If[Parameters`IsRealExpression[trace],
       CConversion`realScalarCType,
       CConversion`complexScalarCType
      ];

GetTraceCType[trace_] :=
    CConversion`CreateCType[GetTraceType[trace]];

FindMultipleTraces[list_List] :=
    Module[{traces},
           traces = Flatten[Cases[list, trace[__], Infinity]];
           traces = (#[[1]])& /@ Select[Tally[traces], (#[[2]] > 1)&];
           Return[traces];
          ];

FindAllTraces[list_List] :=
    DeleteDuplicates[Flatten[Cases[list, trace[__], Infinity]]];

CreateTraceRules[traces_List] :=
    (Rule[#, ToValidCSymbol[#]])& /@ FindAllTraces[traces];

CreateLocalCopiesOfSARAHTraces[expr_, sarahTraces_List, structName_String] :=
    Module[{defs = "", traces},
           traces = FindSARAHTraces[expr, sarahTraces];
           (defs = defs <> "const " <> GetTraceCType[#] <>
            " " <> ToValidCSymbolString[#] <>
            " = " <> structName <> "." <> ToValidCSymbolString[#] <>
            ";\n")& /@ traces;
           Return[defs];
          ];

CreateLocalCopiesOfTraces[list_List, structName_String] :=
    Module[{defs = "", traces},
           traces = FindAllTraces[list];
           (defs = defs <> "const " <> GetTraceCType[#] <>
            " " <> ToValidCSymbolString[#] <> " = " <>
            structName <> "." <> ToValidCSymbolString[#] <> ";\n")& /@ traces;
           Return[defs];
          ];

CreateTraceDefs[list_List] :=
    Module[{defs = "", traces},
           traces = FindAllTraces[list];
           (defs = defs <> "" <> GetTraceCType[#] <> " " <>
            ToValidCSymbolString[#] <> ";\n")& /@ traces;
           Return[defs];
          ];

CreateSARAHTraceDefs[list_List] :=
    Module[{defs = ""},
           (defs = defs <> "" <> GetTraceCType[GetSARAHTraceName[#]] <> " " <>
            ToValidCSymbolString[GetSARAHTraceName[#]] <> ";\n")& /@ list;
           Return[defs];
          ];

CreateCastedTraceExprStr[expr_] :=
    CConversion`CastTo[RValueToCFormString[expr], GetTraceType[expr]];

CreateTraceCalculation[list_List, structName_String] :=
    Module[{defs = "", traces},
           traces = FindAllTraces[list];
           (defs = defs <> structName <> "." <> ToValidCSymbolString[#] <>
            " = " <> CreateCastedTraceExprStr[#] <> ";\n")& /@ traces;
           Return[defs];
          ];

CreateSARAHTraceCalculation[list_List, structName_String] :=
    Module[{defs = ""},
           (defs = defs <> structName <> "." <> ToValidCSymbolString[GetSARAHTraceName[#]] <>
            " = " <> CreateCastedTraceExprStr[GetSARAHTraceExpr[#]] <> ";\n")& /@ list;
           defs = Parameters`CreateLocalConstRefsForInputParameters[list] <>
                  "\n" <> defs;
           Return[defs];
          ];

End[];

EndPackage[];
