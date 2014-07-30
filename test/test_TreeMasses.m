Needs["TestSuite`", "TestSuite.m"];
Needs["TreeMasses`", "TreeMasses.m"];

Print["testing IsSymmetric[] ..."];

TestEquality[TreeMasses`Private`IsSymmetric[{}], True];
TestEquality[TreeMasses`Private`IsSymmetric[{a}], False];
TestEquality[TreeMasses`Private`IsSymmetric[{{a}}], True];
TestEquality[TreeMasses`Private`IsSymmetric[{{a,b},{b,a}}], True];
TestEquality[TreeMasses`Private`IsSymmetric[{{a,b},{c,a}}], False];

PrintTestSummary[];
