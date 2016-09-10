BeginPackage["ReadSLHA`"];

ReadSLHAFile::usage="reads SLHA file and returns list of input and
 output parameters";

ReadSLHAStream::usage="reads stream with SLHA input and returns list
 of input and output parameters";

ReadSLHAString::usage="reads string with SLHA input and returns list
 of input and output parameters";

Begin["`Private`"];

IsDataLine[str_String] :=
    !StringMatchQ[str, StartOfString ~~ "#" ~~ ___];

BlockStarts[str_String, blockName_String] :=
    StringMatchQ[str, StartOfString ~~ "BLOCK" ~~ Whitespace ~~ blockName ~~ ___,
                 IgnoreCase -> True];

BlockStarts[str_String] :=
    StringMatchQ[str, StartOfString ~~ "BLOCK" ~~ Whitespace ~~ ___,
                 IgnoreCase -> True];

floatRegex = "[+\-]?(?:[0-9]*)(?:\\.[0-9]*)?(?:[eE][+\-]?[0-9]+)?";

ReadBlock[stream_, blockName_String] :=
    Module[{line = "", block = "", inBlock = True},
           SetStreamPosition[stream, 0];
           While[(line = Read[stream, String]) =!= EndOfFile,
                 If[BlockStarts[line, blockName],
                    inBlock = True,
                    If[BlockStarts[line],
                       inBlock = False
                      ];
                   ];
                 If[!inBlock || !IsDataLine[line], Continue[]];
                 block = block <> line <> "\n";
                ];
           block
          ];

ComposeFixedRegEx[idx__] :=
    StringExpression[StartOfString, Whitespace,
                     Sequence @@ (StringExpression[ToString[#], Whitespace]& /@ {idx})];

ToWolframExpr[str_String] :=
    ToExpression[StringReplace[str, {"E"|"e" -> "*10^", "=" -> ""}]];

ReadIndexFromBlock[block_String, idx___] :=
    Module[{value = 0, values, line = "", val,
            stream = StringToStream[block], patt = ComposeFixedRegEx[idx]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 If[BlockStarts[line], Continue[]];
                 values = StringCases[line, patt ~~ v:RegularExpression[floatRegex] :> v];
                 If[values =!= {},
                    val = ToWolframExpr[First[values]];
                    If[NumberQ[val],
                       value = val;
                      ];
                   ];
                ];
           value
          ];

IsIndex[i_?IntegerQ] := True;
IsIndex[_] := False;
IsIndex[indices_List] := And @@ (IsIndex /@ indices);
IsIndex[indices__] := IsIndex[{indices}];

WithinRange[i_Integer, dim_Integer] := i >= 1 && i <= dim;
WithinRange[idx_List, dim_List] := Length[idx] == Length[dim] && (And @@ MapThread[WithinRange, {idx, dim}]);

ReadMatrixFromBlock[block_String, dims__] :=
    Module[{matrix = Array[0&, {dims}],
            ndims = Length[{dims}],
            line = "", values, indices, value,
            stream = StringToStream[block]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 If[BlockStarts[line], Continue[]];
                 values = StringSplit[line];
                 If[Length[values] > ndims,
                    indices = ToWolframExpr /@ Take[values, ndims];
                    value = ToWolframExpr[values[[ndims + 1]]];
                    If[IsIndex[indices] && NumberQ[value] && WithinRange[indices, {dims}],
                       matrix[[Sequence @@ indices]] = value;
                      ];
                   ];
                ];
           matrix
          ];

ReadParameter[stream_, par_, {} | {0|1} | 0 | 1, {block_, idx___}] :=
    ReadIndexFromBlock[ReadBlock[stream, ToString[block]], idx];

ReadParameter[stream_, par_, {dims__}, block_] :=
    ReadMatrixFromBlock[ReadBlock[stream, ToString[block]], dims];

ReadParameter[stream_, {par_, type_, block_}] :=
    par -> ReadParameter[stream, par, type, block];

(*
 The elements of the parameters list are 3-component lists, where
 #1 is the parameter name
 #2 is the parameter type
 #3 is the block / entry in the SLHA input

 parameters = {
    {Qin  , {0}   , {EXTPAR, 0}},
    {m0   , {0}   , {MINPAR, 1}},
    {CpHPP, {0}   , {EFFHIGGSCOUPLINGS, 25, 22, 22}},
    {k    , {3}   , KappaIn},
    {Yu   , {3, 3}, YuIN}
 };
 *)
ReadSLHAStream[stream_, parameters_List] :=
    ReadParameter[stream, #]& /@ parameters;

ReadSLHAFile[fileName_String, parameters_List] :=
    Module[{file, values},
           file = OpenRead[fileName];
           values = ReadSLHAStream[file, parameters];
           Close[file];
           values
          ];

ReadSLHAString[slha_String, parameters_List] :=
    ReadSLHAStream[StringToStream[slha], parameters];

End[];

EndPackage[];