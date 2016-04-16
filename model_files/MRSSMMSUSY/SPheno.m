SA`TadpoleEquationsWithField[EWSB] = Drop[SA`TadpoleEquationsWithField[EWSB],-3];
TadpoleEquations[EWSB] = Drop[TadpoleEquations[EWSB],-3];

MINPAR = {
    {1, Ms},
    {3, TanBeta}
};

RealParameters = { TanBeta, Ms, vS, vT };

ParametersToSolveTadpoles = { mHd2, mHu2, vS, vT };
AssumptionsTadpoleEquations = { conj[x_] :> x };

UseParameterAsGUTscale = { Ms };

RenormalizationScaleFirstGuess = Ms^2;
RenormalizationScale = Ms^2;

BoundarySUSYScale = {};

BoundaryHighScale = {
    {mq2, DIAGONAL Ms^2},
    {ml2, DIAGONAL Ms^2},
    {md2, DIAGONAL Ms^2},
    {mu2, DIAGONAL Ms^2},
    {me2, DIAGONAL Ms^2},
    {mS2, Ms^2},
    {mT2, Ms^2},
    {moc2, Ms^2},
    {mRd2, Ms^2},
    {mRu2, Ms^2},
    {\[Mu], 0},
    {B[\[Mu]], Ms^2/(TanBeta + 1/TanBeta)},
    {LamSD, LHInput[LamSD]},
    {LamSU, LHInput[LamSU]},
    {LamTD, LHInput[LamTD]},
    {LamTU, LHInput[LamTU]},
    {MDBS, Ms},
    {MDGoc, Ms},
    {MDWBT, Ms},
    {MuD, Ms},
    {MuU, Ms},
    {B[MuD], 0},
    {B[MuU], 0},
    {vd, Sqrt[vu^2 + vd^2] Cos[ArcTan[TanBeta]]},
    {vu, Sqrt[vu^2 + vd^2] Sin[ArcTan[TanBeta]]}
};

BoundaryLowScaleInput = {
    {vd, Sqrt[4 mz2/(g1^2+g2^2)]*Cos[ArcTan[vu/vd]]},
    {vu, Sqrt[4 mz2/(g1^2+g2^2)]*Sin[ArcTan[vu/vd]]}
};

ListDecayParticles = Automatic;
ListDecayParticles3B = Automatic;

DefaultInputValues = {
    Ms -> 1000,
    TanBeta -> 5
};
