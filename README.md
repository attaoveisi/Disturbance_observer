# Disturbance_observer
In this repo, disturbance rejection control (DRC) based on unknown input observation (UIO), and disturbance-observer based control (DOBC) methods are revisited for a class of MIMO systems with mismatch disturbance conditions. In both of these methods, the estimated disturbance is considered to be in the feedback channel. The disturbance term could represent either unknown mismatched signals penetrating the states, or unknown dynamics not captured in the modeling process, or physical parameter variations not accounted for in the mathematical model of the plant. Unlike the high-gain approaches and variable structure methods, a systematic synthesis of the state/disturbance observer-based controller is carried out. For this purpose, first, using a series of singular value decompositions, the linearized plant is transformed into disturbance-free and disturbance-dependent subsystems. Then, functional state reconstruction based on generalized detectability concept is proposed for the disturbance-free part. Then, a DRC based on quadratic stability theorem is employed to guarantee the performance of the closed-loop system. An important contribution offered in this article is the independence of the estimated disturbance from the control input which seem to be missing in the literature for disturbance decoupling problems. In the second method, DOBC is reconsidered with the aim of achieving a high level of robustness against modeling uncertainties and matched/mismatched disturbances, while at the same time retaining performance. Accordingly, unlike the first method, DRC, full information state observation is developed independent of the disturbance estimation. An advantage of such a combination is that disturbance estimation does not involve output derivatives. Finally, the case of systems with matched disturbances is presented as a corollary of the main results.

Keywords: Unknown input observation, Disturbance rejection control, Disturbance observer-based control, noise decoupling.


References
[1]  Guan, M. Saif, A novel approach to the design of unknown input observers, IEEE T. Automat. Contr. 36 (1991) 632-635.

[2]  M. Hou, P.C. Müller, Disturbance decoupled observer design: a unified viewpoint, IEEE T. Automat. Contr. 39 (1994) 1338-1341.

[3]  J.E. Kurek, The state vector reconstruction for linear systems with unknown inputs, IEEE T. Automat. Contr. 28 (1983) 1120-1122.

[4]  M. Darouach, M. Zasadzinski, S.J. Xu, Full-order observers for linear systems with unknown inputs, IEEE T. Automat. Contr. 39 (1994) 606-609.

[5]  S. Mondal, G. Chakraborty, K. Bhattacharyya, LMI approach to robust unknown input observer design for continuous systems with noise and uncertainties, Int. J. Control Autom. 8 (2010) 210-219.

[6]  A. Oveisi, T. Nestorović, Robust observer-based adaptive fuzzy sliding mode controller, Mech. Syst. Signal Pr. 76-77 (2016) 58-71.

[7]  H.J. Sussmann, P.V. Kokotovic, The peaking phenomenon and the global stabilization of nonlinear systems, IEEE T. Automat. Contr. 36 (1991) 424-440.

[8]  H. Shim H, N.H. Jo, An almost necessary and sufficient condition for robust stability of closed-loop systems with disturbance observer, Automatica 45 (2009) 296-299.

[9]  W-H. Chen, J. Yang, L. Guo, S. Li, Disturbance-observer-based control and related methods-an overview, IEEE T. Ind. Electron. 63 (2016) 1083-1095.

[10]  Z. Gao, Y. Huang, J. Han, An alternative paradigm for control system design, In: Proceedings of the 40th IEEE Conference on Decision and Control, 2001.

[11]  Q-C. Zhong, A. Kuperman, R-K. Stobart, Design of UDE based controllers from their two-degree-of-freedom nature, Int. J. Robust Nonlin. 17 (2011) 1994-2008.

[12]  J. Yang, W.-H. Chen, S. Li, X. Chen, Static disturbance-to-output decoupling for nonlinear systems with arbitrary disturbance relative degree, Int. J. Robust Nonlin. 23 (2013) 562–577.

[13]  J. Yang, S. Li, W.-H. Chen, Nonlinear disturbance observer-based control for multi-input multi-output nonlinear systems subject to mismatching condition, Int. J. Control. 85 (2012) 1071–1082.

[14]  R. Rajamani, Observer for Lipschitz nonlinear systems. IEEE T. Automat. Contr., 43 (1998) 397-401.

[15]  H. Trinh, T.D. Tran, T. Fernando, Disturbance decoupled observers for systems with unknown inputs, IEEE T. Automat. Contr. 53 (2008) 2397-2402.

[16]  C.J. Kempf, S. Kobayashi, Disturbance observer and feedforward design for a high-speed direct-drive position table, IEEE T. Contr. Syst. T. 7 (1999) 513-526.

[17]  W-H. Chen, Disturbance observer based control for nonlinear systems, IEEE-ASME T. Mech. 9 (2004) 706-710.

[18]  C.K. Thum, C. Du, F.L. Lewis, B.M. Chen, E.H. Ong, H-infinity disturbance observer design for high precision track following in hard disk drives, IET Control Theory A. 3 (2009) 1591-1598.

[19]  X.S. Chen, J. Yang, S.H. Li, Q. Li, Disturbance observer based multi-variable control of ball mill grinding circuits. J. Process Contr. 19 (2009) 1205-1213.

[20]  J. Yang, A. Zolotas, W-H. Chen, K. Michail, S. Li, Robust control of nonlinear MAGLEV suspension system with mismatched uncertainties via DOBC approach, ISA T. 50 (2011) 389-396.

[21]  B. Barmish, G. Leitmann, On ultimate boundedness control of uncertain systems in the absence of matching assumptions, IEEE T. Automat. Contr. 27 (1982) 153-158.

[22]  L. Guo, W-H. Chen, Disturbance attenuation and rejection for systems with nonlinearity via DOBC approach, Int. J. Robust Nonlin. 15 (2005) 109-125.

[23]  X.J. Wei, L. Guo, Composite disturbance-observer-based control and terminal sliding model control for non-linear systems with disturbances, Int. J. Control 82 (2009) 1082-1098. 

[24]  X.J. Wei, L. Guo, Composite disturbance-observer-based control and H-infinity control for complex continuous models, Int. J. Robust Nonlin. 20 (2010) 106-118.

[25]  W-H. Chen, D.J. Ballance, P.J. Gawthrop, J.J Gribble, J. O’Reilly, Nonlinear PID predictive controller, IEE P-Contr. Theor. Ap. 146 (1999) 603-611.

[26]  W-H. Chen, Nonlinear disturbance observer-enhanced dynamic inversion control of missiles, J. Guid. Control Dynam. 26 (2003) 161-166.

[27]  J.L. Chang, Robust output feedback disturbance rejection control by simultaneously estimating state and disturbance, J. Contr. Sci. Eng. 2011 (2011) 1-13.

[28]  D. Koenig, S. Mammar, Design of a class of reduced order unknown inputs nonlinear observer for fault diagnosis, In: Proceedings of the American Control Conference, Arlington, 2001.

[29]  J.A. Moreno, Simultaneous observation of linear systems: a state-space interpretation, IEEE T. Automat. Contr. 50 (2005) 1021-1025.

[30]  G. Li, G. Herrman, D.P. Stoten, J. Tu, M.C. Turner, A novel robust disturbance rejection anti-windup framework, Int. J. Control 84 (2011) 123-137. 

[31]  C.-H. Lien, H∞ non-fragile observer-based controls of dynamical systems via LMI optimization approach, Chaos Soliton. Fract. 34 (2007) 428-436.

[32]  S. Boyd, L. El Ghaoui, E. Feron, V. Balakrishnan, Linear matrix inequalities in system and control theory, Vol 15 of studies in applied mathematics, SIAM, Philadelphia, PA, 1994.

[33]  Z. Gao, H. Wang, Descriptor observer approaches for multivariable systems with measurement noises and application in fault detection and diagnosis, Syst. Control Lett. 55 (2006) 304-313.

[34]  M. Aldeen, R. Sharma, Estimation of states, faults and unknown disturbances in non-linear systems, Int. J. Control 81 (2008) 1195-1201.

[35]  G. R. Duan, “Solution to matrix equation AV+BW = VF and their applications to eigenstructure assignment in linear systems, IEEE T. Automat. Contr. 38 (1993) 276-280.

[36]  M.R. Garey, D.S. Johnson, Computers and intractability: A guide to NP-completeness, W. H. Freeman and company, New York, 1983.

[37]  C.H. Papadimitriou, K. Steiglitz, Combinatorial global optimization: algorithms and complexity, Prentice-Hall, Englewood Cliffs, New Jersey, 1982.

[38]  J.G. VanAntwerp, R.D. Braatz, A tutorial on linear and bilinear matrix inequalities, J. Process Contr. 10 (2000) 363-385.

[39]  E.L. Lawler, D.E. Wood, Branch-and-Bound methods: A survey, Oper. Res. 14 (1966) 699-719.

[40]  J. Löfberg, YALMIP: A toolbox for modeling and optimization in MATLAB. In: Proceedings of the CACSD Conference, Taipei, Taiwan, 2004.

[41]  J.J. Rubio, F. Meléndez, M. Figueroa, An observer with controller to detect and reject disturbances, Int. J. Control 87 (2014) 524-536.

[42]  J.J. Rubio, M. Figueroa, J.H.P. Cruz, F.J. Bejarano, Geometric approach and structure at infinity controls for the disturbance rejection, IET Control Theory A. 6 (2012) 2528-2537.

[43]  R. Sharma, M. Aldeen, Estimation of unknown disturbances in nonlinear systems, In: Control, Bath, UK, 2004.

[44]  A. Oveisi, T. Nestorović, Robust nonfragile observer-based  controller, J. Vib. Control (2016) 1-17.

[45]  M.L.J. Hautus, Strong detectability and observers, Linear Alg. Appl. 50 (1983) 353-368.

[46]  J. Moreno, Quasi unknown input observers for linear systems, In: Proceedings of IEEE International Conference on Control Applications, 2001.

[47]  F. Alizadeh, Optimization over the positive semi-definite cone: interior-point methods and combinatorial algorithms, Advances in Optimization and Parallel Computing, Elsevier Science, 1992.

[48]  L. Vandenberghe, S. Boyd, Semidefinite programming, SIAM Rev. 38 (1996) 49-95.

[49]  W-H. Chen, P.J. Gawthrop, J. O’Reilly, A nonlinear disturbance observer for robotic manipulators, IEEE T. Ind. Electron. 47 (2000) 932-938.

[50]  A. Oveisi, T. Nestorović, Mu-Synthesis based active robust vibration control of an MRI inlet, FACTA Univ. Ser. Mech. Eng. 14 (2016) 37-53.

[51]  J.P. Noël, G. Kerschen, Nonlinear system identification in structural dynamics: 10 more years of progress, Mech. Syst. Signal Process. 83 (2017) 2-35.

[52]  A. Oveisi, T. Nestorović, Transient response of an active nonlinear sandwich piezolaminated plate, Commun. Nonlinear Sci. Numer. Simul. 45 (2017) 158-175.

[53]  J.P. Noël, J. Schoukens, Grey-box state-space identification of nonlinear mechanical vibrations, Int. J. Control. 1–22 (2017).
