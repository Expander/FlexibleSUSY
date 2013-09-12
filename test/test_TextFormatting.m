Needs["TestSuite`", "TestSuite.m"];
Needs["TextFormatting`", "TextFormatting.m"];

Print["testing GetBestSplitPoint[] ..."];

TestEquality[TextFormatting`Private`GetBestSplitPoint["",0], 0];
TestEquality[TextFormatting`Private`GetBestSplitPoint["a",1], 1];
TestEquality[TextFormatting`Private`GetBestSplitPoint["a b",1], 1];
TestEquality[TextFormatting`Private`GetBestSplitPoint["a b",2], 2];
TestEquality[TextFormatting`Private`GetBestSplitPoint["a b c",4], 4];
TestEquality[TextFormatting`Private`GetBestSplitPoint["a b c",10], StringLength["a b c"]];

(* test that we do not split a C string *)
TestEquality[TextFormatting`Private`GetBestSplitPoint["\"a,b\"\"cd\"",2], StringLength["\"a,b\"\"cd\""]];
TestEquality[TextFormatting`Private`GetBestSplitPoint["1 + \"a,b\"",7], 4];

Print["testing SplitLine[] ..."];

TestEquality[TextFormatting`Private`SplitLine["",0], {""}];
TestEquality[TextFormatting`Private`SplitLine["",1], {""}];
TestEquality[TextFormatting`Private`SplitLine["a b",1], {"a", " ", "b"}];
TestEquality[TextFormatting`Private`SplitLine["a b",2], {"a ", "b"}];
TestEquality[TextFormatting`Private`SplitLine["a b c",3], {"a b", " c"}];
TestEquality[TextFormatting`Private`SplitLine["a b c",4], {"a b ", "c"}];
(* no break possible *)
TestEquality[TextFormatting`Private`SplitLine["abcdefghij",5], {"abcdefghij"}];

Print["testing WrapLines[] ..."];

TestEquality[WrapLines["a"  ,1,""], "a\n"];
TestEquality[WrapLines["a\n",1,""], "a\n"];
TestEquality[WrapLines["a  ",1,""], "a\n"];
TestEquality[WrapLines["a b",2,""], "a\nb\n"];
TestEquality[WrapLines["a b c",1,""], "a\nb\nc\n"];
TestEquality[WrapLines["abc def",3,""], "abc\ndef\n"];
TestEquality[WrapLines["abc def",4,""], "abc\ndef\n"];

(* test indentation *)
TestEquality[WrapLines[" abc,def",5,""], " abc\n ,def\n"];
TestEquality[WrapLines[" abc,def",6,""], " abc,\n def\n"];
TestEquality[WrapLines[" abc def",5,""], " abc\n def\n"];

(* test indentation + offset *)
TestEquality[WrapLines[" abc,def",5,"x"], " abc\n x,\n xdef\n"];
TestEquality[WrapLines[" abc,def",6,"x"], " abc\n x,def\n"];
(* no break possible because offset is too long *)
TestEquality[WrapLines[" abc,def",5,"xxx"], " abc,def\n"];

Print["testing IndentText[] ..."];

TestEquality[IndentText["abc",1], " abc"];
TestEquality[IndentText["abc\n",1], " abc\n"];
TestEquality[IndentText["abc\ndef",1], " abc\n def"];
TestEquality[IndentText["abc\n\ndef",1], " abc\n\n def"];
TestEquality[IndentText["abc\n\n\ndef",1], " abc\n\n\n def"];
TestEquality[IndentText["abc\n\n\n\ndef",1], " abc\n\n\n\n def"];

PrintTestSummary[];
