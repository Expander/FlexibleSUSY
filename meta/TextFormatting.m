
BeginPackage["TextFormatting`"];

WrapLines::usage="breaks text lines";
IndentText::usage="indents text by a given number of spaces";

Begin["`Private`"];

GetBestSplitPoint[line_String, maxWidth_:79] :=
    Module[{chars, lastSplitPoint, i, char, nextChar, numberOfQuotes = 0,
            splitChars = StringSplit["(){}[]*+,;<> ", ""]},
           If[StringLength[line] <= maxWidth, Return[StringLength[line]]];
           chars = StringSplit[line, ""];
           lastSplitPoint = StringLength[line];
           For[i = 1, i < Length[chars] && i <= maxWidth, i++,
               char = chars[[i]];
               If[char == "\"", numberOfQuotes++;];
               nextChar = chars[[i+1]];
               If[(MemberQ[splitChars, char] || MemberQ[splitChars, nextChar])
                  && EvenQ[numberOfQuotes],
                  lastSplitPoint = i;
                 ];
              ];
           Return[lastSplitPoint];
          ];

SplitLine[line_String, maxWidth_:79] :=
    Module[{bestSplitPoint, first, rest, result = {}},
           bestSplitPoint = GetBestSplitPoint[line, maxWidth];
           first = StringTake[line, bestSplitPoint];
           AppendTo[result, first];
           rest = StringDrop[line, bestSplitPoint];
           While[rest != "",
                 bestSplitPoint = GetBestSplitPoint[rest, maxWidth];
                 first = StringTake[rest, bestSplitPoint];
                 AppendTo[result, first];
                 rest = StringDrop[rest, bestSplitPoint];
                ];
           Return[result];
          ];

WrapLines[text_String, maxWidth_:79, offset_:"   "] :=
    Module[{result = "", lines, line, i, k, indent, splitLines},
           lines = StringSplit[text, "\n"];
           For[i = 1, i <= Length[lines], i++,
               line = lines[[i]];
               (* get intentation of line *)
               indent = StringReplace[line, x:(StartOfString ~~ Whitespace...) ~~ ___ -> x];
               splitLines = SplitLine[line, maxWidth - StringLength[indent]
                                      - StringLength[offset]];
               (* remove strings that consist of whitespace only *)
               splitLines = Cases[StringReplace[#,StartOfString ~~ Whitespace.. ~~ EndOfString -> EmptyString]& /@ splitLines, _String];
               For[k = 1, k <= Length[splitLines], k++,
                   (* strip whitespace *)
                   splitLines[[k]] = StringTrim[splitLines[[k]]];
                   result = result <> indent;
                   If[k > 1, result = result <> offset];
                   result = result <> splitLines[[k]] <> "\n";
                  ];
              ];
           Return[result];
          ];

IndentText[text_String, spaces_Integer:3] :=
    Module[{i, whiteSpace = ""},
           If[text == "", Return[text];];
           For[i = 0, i < spaces, i++, whiteSpace = whiteSpace <> " "];
           Return[StringReplace[whiteSpace <> text,
                                { x:("\n" ~~ ("\n"..)) :> x <> whiteSpace,
                                  "\n" ~~ EndOfString -> "\n",
                                  "\n" -> "\n" <> whiteSpace } ]];
          ];

End[];

EndPackage[];
