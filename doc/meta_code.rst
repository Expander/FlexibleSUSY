Meta code
=========

.. contents:: Table of Contents
..    :depth: 2

Re-usable meta code packages
----------------------------

FlexibleSUSY contains several packages with general algorithms needed
in Quantum Field Theory calculations and implementations of analytic
expressions from the literature.

Passarino-Veltman loop functions
````````````````````````````````

``meta/LoopFunctions.m`` contains an implementation of
Passarino-Veltman 1-loop functions :math:`A_0`, :math:`B_0`,
:math:`B_1`, :math:`B_{00}`, :math:`B_{11}`, etc.

Renormalization group equation integrator
`````````````````````````````````````````

``meta/RGIntegrator.m`` contains a routine, which perturbatively
integrates a system of coupled renormalization group equations.

Beta functions
``````````````

``meta/ThreeLoopMSSM.m`` contains a routine which returns the beta
functions of the MSSM up to the 3-loop level.  The beta functions are
stored in the files ``meta/MSSM/beta_*.m``, which have been obtained
from http://www.liv.ac.uk/~dij/betas/allgennb.log [hep-ph:0308231]_.

``meta/ThreeLoopSM.m`` contains a routine which returns the beta
functions of the Standard Model up to the 3-loop level.  The beta
functions are stored in the files ``meta/SM/beta_*.m``, which have
been obtained from [1303.4364]_, [1504.05200]_.

Loop corrections to masses
``````````````````````````

Standard Model
''''''''''''''

``meta/SM/Mh2_effpot.m`` contains the QCD contributions to the 4-loop
effective Higgs potential in the Standard Model from [1508.00912]_

``meta/ThreeLoopQCD.m`` contains a routine, which returns the ratio of
the :math:`\overline{\text{MS}}` top mass over the top pole mass,
[hep-ph:9912391]_, Eq. (10).  The expression contains the full
renormalization scale dependence, which has been taken from
[hep-ph:9911434]_.

``meta/TwoLoopQCD.m`` contains routines, which return the ratio of the
top pole mass over the running top mass up to the 2-loop level in the
:math:`\overline{\text{MS}}` and :math:`\overline{\text{DR}}` schemes
[hep-ph:0210258]_, [hep-ph:9803493]_.

MSSM
''''

``meta/TwoLoopMSSM.m`` contains routines, which return the
analytic 2-loop corrections to the Higgs masses in the CP-conserving
MSSM [hep-ph:0105096]_.


Threshold corrections
`````````````````````

Standard Model
''''''''''''''

``meta/SM/mf_3loop_qcd.m`` contains the 3-loop relation between a quark pole
mass and the corresponding running :math:`\overline{\text{MS}}` mass
from [hep-ph:9912391]_, [hep-ph:9911434]_.

2HDM
''''

``meta/THDM/Thresholds_1L_full.m`` contains the implementation of the
complete analytic 1-loop threshold corrections of the THDM and the
THDM + Higgsinos + gauginos to the MSSM [0901.2065]_.

MSSM
''''

``meta/MSSM/tquark_2loop_strong.m`` contains the analytic
expression for the 2-loop relation :math:`O(\alpha_s^2)` between the top
quark pole mass and the :math:`\overline{\text{DR}}` top mass in the
MSSM [hep-ph:0210258]_, [hep-ph:0507139]_.

``meta/MSSM/bquark_2loop_sqcd_decoupling.m`` contains the analytic
expression for the 2-loop relation :math:`O(\alpha_s^2)` between the
:math:`\overline{\text{MS}}` bottom quark mass in the Standard Model
(without the top quark) and the :math:`\overline{\text{DR}}` bottom
mass in the MSSM [0707.0650]_.

``meta/MSSM/dmtauas2.m`` contains the analytic expression for the
2-loop relation between the tau lepton pole mass and the
:math:`\overline{\text{DR}}` tau mass in the MSSM.

``meta/MSSM/das2.m`` contains the analytic expression for the 2-loop
relation between the :math:`\overline{\text{MS}}` :math:`\alpha_s` in
the Standard Model (without the top quark) and the
:math:`\overline{\text{DR}}` value in the MSSM [hep-ph:0509048]_,
[0810.5101]_, [1009.5455]_.

.. [hep-ph:9803493] `Nucl.Phys. B539 (1999) 671-690 <https://inspirehep.net/record/468752>`_ [`arXiv:hep-ph/9803493 <https://arxiv.org/abs/hep-ph/9803493>`_]
.. [hep-ph:9911434] `Nucl.Phys. B573 (2000) 617-651 <https://inspirehep.net/record/510551>`_ [`arXiv:hep-ph/9911434 <https://arxiv.org/abs/hep-ph/9911434>`_]
.. [hep-ph:9912391] `Phys.Lett. B482 (2000) 99-108 <https://inspirehep.net/record/522686>`_ [`arXiv:hep-ph/9912391 <https://arxiv.org/abs/hep-ph/9912391>`_]
.. [hep-ph:0105096] `Nucl.Phys. B611 (2001) 403-422 <https://inspirehep.net/record/556417>`_ [`arXiv:hep-ph/0105096 <https://arxiv.org/abs/hep-ph/0105096>`_]
.. [hep-ph:0210258] `Eur.Phys.J. C29 (2003) 87-101 <https://inspirehep.net/record/600038>`_ [`arXiv:hep-ph/0210258 <https://arxiv.org/abs/hep-ph/0210258>`_]
.. [hep-ph:0308231] `Phys.Lett. B579 (2004) 180-188 <https://inspirehep.net/record/626390>`_ [`arXiv:hep-ph/0308231 <https://arxiv.org/abs/hep-ph/0308231>`_]
.. [hep-ph:0507139] `Phys.Atom.Nucl. 71 (2008) 343-350 <https://inspirehep.net/record/687205>`_ [`arXiv:hep-ph/0507139 <https://arxiv.org/abs/hep-ph/0507139>`_]
.. [hep-ph:0509048] `Phys.Rev. D72 (2005) 095009 <https://inspirehep.net/record/691479>`_ [`arXiv:hep-ph/0509048 <https://arxiv.org/abs/hep-ph/0509048>`_]
.. [0707.0650] `Int.J.Mod.Phys. A22 (2007) 5245-5277 <https://inspirehep.net/record/755029>`_ [`arXiv:0707.0650 <https://arxiv.org/abs/0707.0650>`_]
.. [0810.5101] `JHEP 0902 (2009) 037 <https://inspirehep.net/record/800842>`_ [`arXiv:0810.5101 <https://arxiv.org/abs/0810.5101>`_]
.. [0901.2065] `Phys.Rev. D84 (2011) 034030 <https://inspirehep.net/record/811006>`_ [`arXiv:0901.2065 <https://arxiv.org/abs/0901.2065>`_]
.. [1009.5455] `C10-06-06.1 <https://inspirehep.net/record/871111>`_ [`arXiv:1009.5455 <https://arxiv.org/abs/1009.5455>`_]
.. [1303.4364] `Nucl.Phys. B875 (2013) 552-565 <https://inspirehep.net/record/1224266>`_ [`arXiv:1303.4364 <https://arxiv.org/abs/1303.4364>`_]
.. [1504.05200] `JHEP 1507 (2015) 159 <https://inspirehep.net/record/1362483>`_ [`arXiv:1504.05200 <https://arxiv.org/abs/1504.05200>`_]
.. [1508.00912] `Phys.Rev. D92 (2015) no.5, 054029 <https://inspirehep.net/record/1386688>`_ [`arXiv:1508.00912 <https://arxiv.org/abs/1508.00912>`_]
