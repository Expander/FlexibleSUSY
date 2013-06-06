
BeginPackage["TestSuite`"];

TestEquality::usage="tests equality of two expressions";
TestCPPCode::usage="tests a C/C++ code snippet for an expected
result";
PrintTestSummary::usage="prints test summary";

Begin["Private`"];

numberOfFailedTests := 0;
numberOfPassedTests := 0;

TestEquality[val_, expr_, msg_:""] := 
    If[val =!= expr,
       numberOfFailedTests++;
       Print["Error: expressions are not equal: ",
             InputForm[val], " =!= ", InputForm[expr]];,
       numberOfPassedTests++;
      ];

TestCPPCode[{preface_String, expr_String}, value_String, type_String, expected_String] :=
    Module[{code = expr, output},
           code = code <> "\n" <>
                  type <> " result__ = " <> value <> ";\n" <>
                  "std::cout << result__ << std::endl;";
           output = RunCPPProgram[{preface, code}];
           TestEquality[output, expected];
          ];

PrintTestSummary[] :=
    Block[{},
          Print["Test summary"];
          Print["============"];
          Print["Number of passed tests: ", numberOfPassedTests];
          Print["Number of failed tests: ", numberOfFailedTests];
          Print["Total number of tests: ", numberOfPassedTests + numberOfFailedTests];
         ];

RunCPPProgram[{preface_String, expr_String}, fileName_String:"tmp.cpp"] :=
    Module[{code, output = "", errorCode},
           code = "#include <iostream>\n" <>
                  preface <> "\n" <>
                  "int main() {\n" <>
                  expr <>
                  "\nreturn 0;\n}\n";
           Export[fileName, code, "String"];
           errorCode = Run["g++ -o a.out " <> fileName];
           If[errorCode != 0,
              Print["Error: could not compile the following: ", code];
              Return[""];
             ];
           Run["./a.out > a.out.log"];
           If[errorCode != 0, Return[""]];
           If[MemberQ[FileNames[], "a.out.log"],
              output = Import["a.out.log"];,
              Print["Error: output file \"a.out.log\" not found"];
             ];
           DeleteFile[{"a.out", "a.out.log", fileName}];
           Return[output];
          ];

End[];

EndPackage[];
