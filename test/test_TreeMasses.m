Needs["TestSuite`", "TestSuite.m"];
Needs["TreeMasses`", "TreeMasses.m"];

Print["testing IsSymmetric[] ..."];

TestEquality[Private`IsSymmetric[{}], True];
TestEquality[Private`IsSymmetric[{a}], False];
TestEquality[Private`IsSymmetric[{{a}}], True];
TestEquality[Private`IsSymmetric[{{a,b},{b,a}}], True];
TestEquality[Private`IsSymmetric[{{a,b},{c,a}}], False];

PrintTestSummary[];
