Needs["TestSuite`", "TestSuite.m"];
Needs["TreeMasses`", "TreeMasses.m"];

Print["testing IsSymmetric[] ..."];

TestEquality[TreeMasses`IsSymmetric[{}], True];
TestEquality[TreeMasses`IsSymmetric[{a}], False];
TestEquality[TreeMasses`IsSymmetric[{{a}}], True];
TestEquality[TreeMasses`IsSymmetric[{{a,b},{b,a}}], True];
TestEquality[TreeMasses`IsSymmetric[{{a,b},{c,a}}], False];

PrintTestSummary[];
