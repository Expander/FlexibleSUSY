Needs["TestSuite`", "TestSuite.m"];
Get["model_files/HSSUSY/FlexibleSUSY.m.in"];

limitMQMUM3Deg = \
    Limit[
        Limit[
            lambda2LPhiHSSAlphaTAlphaSFull,
            msu2[3, 3] -> M3Input^2
        ],
        msq2[3, 3] -> M3Input^2
    ];

limitMQMUM3QDeg = Limit[limitMQMUM3Deg, M3Input -> SCALE];

limitMQM3Deg = \
    Limit[
        lambda2LPhiHSSAlphaTAlphaSFull,
        msq2[3, 3] -> M3Input^2
    ];

limitMUM3Deg = \
    Limit[
        lambda2LPhiHSSAlphaTAlphaSFull,
        msu2[3, 3] -> M3Input^2
    ];

TestEquality[Simplify[limitMQMUM3Deg  - lambda2LPhiHSSAlphaTAlphaSMQMUM3Degenerate], 0];
TestEquality[Simplify[limitMQMUM3QDeg - lambda2LPhiHSSAlphaTAlphaSDegenerate]      , 0];
TestEquality[Simplify[limitMQM3Deg    - lambda2LPhiHSSAlphaTAlphaSMQM3Degenerate]  , 0];
TestEquality[Simplify[limitMUM3Deg    - lambda2LPhiHSSAlphaTAlphaSMUM3Degenerate]  , 0];

PrintTestSummary[];
