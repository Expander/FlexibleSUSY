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

ReadIndexFromBlock[block_String, idx_] :=
    Module[{value = 0, values, line = "", val,
            stream = StringToStream[block], idxStr = ToString[idx]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 values = StringCases[line, StartOfString ~~ Whitespace ~~ idxStr ~~ Whitespace ~~ v:NumberString :> v];
                 If[values =!= {},
                    val = ToExpression[First[values]];
                    If[NumberQ[val],
                       value = val;
                      ];
                   ];
                ];
           value
          ];

ReadVectorFromBlock[block_String, dim_] :=
    Module[{vector = Table[0, {i,1,dim}], line = "", value, index, values, 
            stream = StringToStream[block]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 values = StringCases[line, StartOfString ~~ Whitespace ~~ i:DigitCharacter.. ~~ Whitespace ~~ v:NumberString :> {i,v}];
                 If[values =!= {},
                    {index, value} = First[values];
                    index = ToExpression[index];
                    value = ToExpression[value];
                    If[IntegerQ[index] && NumberQ[value],
                       vector[[index]] = value;
                      ];
                   ];
                ];
           vector
          ];

ReadMatrixFromBlock[block_String, dim1_, dim2_] :=
    Module[{matrix = Table[0, {i,1,dim1}, {j,1,dim2}], line = "",
            value, index1, index2, values, 
            stream = StringToStream[block]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 values = StringCases[line, StartOfString ~~ Whitespace ~~ i:DigitCharacter.. ~~ Whitespace ~~ j:DigitCharacter.. ~~ Whitespace ~~ v:NumberString :> {i,j,v}];
                 If[values =!= {},
                    {index1, index2, value} = First[values];
                    index1 = ToExpression[index1];
                    index2 = ToExpression[index2];
                    value  = ToExpression[value];
                    If[IntegerQ[index1] && IntegerQ[index2] && NumberQ[value],
                       matrix[[index1, index2]] = value;
                      ];
                   ];
                ];
           matrix
          ];

ReadTensorFromBlock[block_String, dim1_, dim2_, dim3_] :=
    Module[{tensor = Table[0, {i,1,dim1}, {j,1,dim2}, {k,1,dim3}], line = "",
            value, index1, index2, index3, values,
            stream = StringToStream[block]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 values = StringCases[line, StartOfString ~~ Whitespace ~~ i:DigitCharacter.. ~~ Whitespace ~~ j:DigitCharacter.. ~~ Whitespace ~~ k:DigitCharacter.. ~~ Whitespace ~~ v:NumberString :> {i,j,k,v}];
                 If[values =!= {},
                    {index1, index2, index3, value} = First[values];
                    index1 = ToExpression[index1];
                    index2 = ToExpression[index2];
                    index3 = ToExpression[index3];
                    value  = ToExpression[value];
                    If[IntegerQ[index1] && IntegerQ[index2] && IntegerQ[index3] && NumberQ[value],
                       tensor[[index1, index2, index3]] = value;
                      ];
                   ];
                ];
           tensor
          ];

ReadTensorFromBlock[block_String, dim1_, dim2_, dim3_, dim4_] :=
    Module[{tensor = Table[0, {i,1,dim1}, {j,1,dim2}, {k,1,dim3}, {l,1,dim4}], line = "",
            value, index1, index2, index3, index4, values,
            stream = StringToStream[block]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 values = StringCases[line, StartOfString ~~ Whitespace ~~ i:DigitCharacter.. ~~ Whitespace ~~ j:DigitCharacter.. ~~ Whitespace ~~ k:DigitCharacter.. ~~ Whitespace ~~ l:DigitCharacter.. ~~ Whitespace ~~ v:NumberString :> {i,j,k,l,v}];
                 If[values =!= {},
                    {index1, index2, index3, index4, value} = First[values];
                    index1 = ToExpression[index1];
                    index2 = ToExpression[index2];
                    index3 = ToExpression[index3];
                    index4 = ToExpression[index4];
                    value  = ToExpression[value];
                    If[IntegerQ[index1] && IntegerQ[index2] && IntegerQ[index3] && IntegerQ[index4] && NumberQ[value],
                       tensor[[index1, index2, index3, index4]] = value;
                      ];
                   ];
                ];
           tensor
          ];

ReadParameter[stream_, par_, {} | {0|1} | 0 | 1, {block_, idx_}] :=
    ReadIndexFromBlock[ReadBlock[stream, ToString[block]], idx];

ReadParameter[stream_, par_, {dim_}, block_] :=
    ReadVectorFromBlock[ReadBlock[stream, ToString[block]], dim];

ReadParameter[stream_, par_, {dim1_, dim2_}, block_] :=
    ReadMatrixFromBlock[ReadBlock[stream, ToString[block]], dim1, dim2];

ReadParameter[stream_, par_, {dims__}, block_] :=
    ReadTensorFromBlock[ReadBlock[stream, ToString[block]], dims];

ReadParameter[stream_, {par_, type_, block_}] :=
    par -> ReadParameter[stream, par, type, block];

(*
 The elements of the parameters list are 3-component lists, where
 #1 is the parameter name
 #2 is the parameter type
 #3 is the block / entry in the SLHA input

 parameters = {
    {Qin, {0}   , {EXTPAR, 0}},
    {m0 , {0}   , {MINPAR, 1}},
    {k  , {3}   , KappaIn},
    {Yu , {3, 3}, YuIN}
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
