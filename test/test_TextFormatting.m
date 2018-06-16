(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

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

TestEquality[WrapLines["a"  ,1,""], "a"];
TestEquality[WrapLines["a\n",1,""], "a\n"];
TestEquality[WrapLines["a  ",1,""], "a"];
TestEquality[WrapLines["a b",2,""], "a\nb"];
TestEquality[WrapLines["a b c",1,""], "a\nb\nc"];
TestEquality[WrapLines["abc def",3,""], "abc\ndef"];
TestEquality[WrapLines["abc def",4,""], "abc\ndef"];

(* test indentation *)
TestEquality[WrapLines[" abc,def",5,""], " abc,\n def"];
TestEquality[WrapLines[" abc,def",6,""], " abc,\n def"];
TestEquality[WrapLines[" abc def",5,""], " abc\n def"];

(* test indentation + offset *)
TestEquality[WrapLines[" abc,def",5,"x"], " abc,\n  def"];
TestEquality[WrapLines[" abc,def",6,"x"], " abc,\n  def"];
TestEquality[WrapLines[" abc,def",5,"xxx"], " abc,\n    def"];

Print["testing IndentText[] ..."];

TestEquality[IndentText["abc",1], " abc"];
TestEquality[IndentText["abc\n",1], " abc\n"];
TestEquality[IndentText["abc\ndef",1], " abc\n def"];
TestEquality[IndentText["abc\n\ndef",1], " abc\n\n def"];
TestEquality[IndentText["abc\n\n\ndef",1], " abc\n\n\n def"];
TestEquality[IndentText["abc\n\n\n\ndef",1], " abc\n\n\n\n def"];

PrintTestSummary[];
