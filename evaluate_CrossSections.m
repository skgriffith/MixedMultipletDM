(* ::Package:: *)

(* ::Input::Initialization:: *)
t=AbsoluteTime[];


(* ::Section::Initialization:: *)
(*Define Particles*)


(* ::Input::Initialization:: *)
(*Representation sizes*)
(*Requires multiplets with the representation sizes Majorana=(1+Dirac)*)
(*Must have Majorana representation >=3*)
Mrep=13;
Drep=12;
Dbarrep=Drep;

(*Size of Higgs coupling ("small", "mid", or "large")*)
yHval="large";


(* ::Section::Initialization:: *)
(*Preliminary*)


(* ::Input::Initialization:: *)
(*Make sure package mixed_multiplet.m is in the same folder as this notebook*)
SetDirectory[NotebookDirectory[]];
Get["./mixed_multiplet.m"];
Get["./overlap_integrals.m"];

(*coupling strength to the higgs, if yH=0, set to something extremely small but nonzero*)
yHAbs=Which[yHval=="large",1,yHval=="mid",0.005,yHval=="small",10^-9];
\[Alpha]Hactual=yHAbs^2/(4*\[Pi]);

(*Averaging factor for cross sections*)
Avg=1/(g[Mrep]+g[Drep])^2;

(*Replacement of couplings*)
\[Alpha]2MZ=\[Alpha]2MZactual;
\[Alpha]1MZ=\[Alpha]1MZactual;\[Alpha]H=\[Alpha]Hactual;


(* ::Subsection::Initialization::Closed:: *)
(*Load Pre-Computed Overlap Integrals*)


(* ::Input::Initialization:: *)
(*all of the relevant most bound states are n=1, \[ScriptL]=0 so we can use the pre-computed overlap integrals*)
Clear[Jintcalc,Tintcalc,Fintcalc,nbs,\[ScriptL]bs];
nbs=1;
\[ScriptL]bs=0;
Jintcalc=JintCalc1100;
Tintcalc=TintCalc1100;
Fintcalc=FintCalc100;
Print["Overlap integrals loaded"]


(* ::Section::Initialization::Closed:: *)
(*Annihilation of Scattering States*)


(* ::Subsection::Initialization::Closed:: *)
(*Tree Level Cross Section (averaged)*)


(* ::Input::Initialization:: *)
\[Sigma]vAnnTreelevel=Module[{\[Sigma]vAnnMix1,\[Sigma]vAnnMix2,\[Sigma]vAnnMM,\[Sigma]vAnnDD,\[Sigma]vAnnDbarDbar,\[Sigma]vAnnMDbar,\[Sigma]vAnnMD},If[Mrep<5,
\[Sigma]vAnnMix1=Sum[Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin],{Rep,{1,3}},{Spin,{0,1}}];
\[Sigma]vAnnMix2=Sum[Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin],{Rep,{1,3}},{Spin,{0,1}}];
\[Sigma]vAnnMM=Sum[\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin],{Rep,{5}},{Spin,{0,1}}];
\[Sigma]vAnnDD=Sum[\[Sigma]vAnn[Rep,0,Drep,0,0,Spin],{Rep,{1,3}},{Spin,{0,1}}];
\[Sigma]vAnnDbarDbar=Sum[\[Sigma]vAnn[Rep,0,0,Dbarrep,0,Spin],{Rep,{1,3}},{Spin,{0,1}}];
\[Sigma]vAnnMDbar=Sum[\[Sigma]vAnn[Rep,Mrep,0,Dbarrep,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
\[Sigma]vAnnMD=Sum[\[Sigma]vAnn[Rep,Mrep,Drep,0,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
Avg*2*(\[Sigma]vAnnMix1+\[Sigma]vAnnMix2+\[Sigma]vAnnMM+\[Sigma]vAnnDD+\[Sigma]vAnnDbarDbar+\[Sigma]vAnnMDbar+\[Sigma]vAnnMD)
,
\[Sigma]vAnnMix1=Sum[Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnMix2=Sum[Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnDD=Sum[\[Sigma]vAnn[Rep,0,Drep,0,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnDbarDbar=Sum[\[Sigma]vAnn[Rep,0,0,Dbarrep,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnMDbar=Sum[\[Sigma]vAnn[Rep,Mrep,0,Dbarrep,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
\[Sigma]vAnnMD=Sum[\[Sigma]vAnn[Rep,Mrep,Drep,0,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
Avg*2*(\[Sigma]vAnnMix1+\[Sigma]vAnnMix2+\[Sigma]vAnnDD+\[Sigma]vAnnDbarDbar+\[Sigma]vAnnMDbar+\[Sigma]vAnnMD)]];
Print["Tree-level done"]


(* ::Subsection::Initialization::Closed:: *)
(*Sommerfeld Enhanced*)


(* ::Input::Initialization:: *)
\[Sigma]vAnnSomm=Module[{\[Sigma]vAnnMix1Somm,\[Sigma]vAnnMix2Somm,\[Sigma]vAnnMMSomm,\[Sigma]vAnnDDSomm,\[Sigma]vAnnDbarDbarSomm,\[Sigma]vAnnMDbarSomm,\[Sigma]vAnnMDSomm},If[Mrep<5,
\[Sigma]vAnnMix1Somm=Sum[Somm[Chop[\[Alpha]eff[Mrep,Drep,Dbarrep,Rep,1]],v]*(Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin]),{Rep,{1,3}},{Spin,{0,1}}];
\[Sigma]vAnnMix2Somm=Sum[Somm[Chop[\[Alpha]eff[Mrep,Drep,Dbarrep,Rep,2]],v]*(Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin]),{Rep,{1,3}},{Spin,{0,1}}];
\[Sigma]vAnnMMSomm=Sum[Somm[Chop[\[Alpha]eff[Mrep,0,0,Rep,0]],v]*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin],{Rep,{5}},{Spin,{0,1}}];
\[Sigma]vAnnDDSomm=Sum[Somm[Chop[\[Alpha]eff[0,Drep,0,Rep,0]],v]*\[Sigma]vAnn[Rep,0,Drep,0,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnDbarDbarSomm=Sum[Somm[Chop[\[Alpha]eff[0,0,Dbarrep,Rep,0]],v]*\[Sigma]vAnn[Rep,0,0,Dbarrep,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnMDbarSomm=Sum[Somm[Chop[\[Alpha]eff[Mrep,0,Dbarrep,Rep,0]],v]*\[Sigma]vAnn[Rep,Mrep,0,Dbarrep,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
\[Sigma]vAnnMDSomm=Sum[Somm[Chop[\[Alpha]eff[Mrep,Drep,0,Rep,0]],v]*\[Sigma]vAnn[Rep,Mrep,Drep,0,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
Avg*2*(\[Sigma]vAnnMix1Somm+\[Sigma]vAnnMix2Somm+\[Sigma]vAnnMMSomm+\[Sigma]vAnnDDSomm+\[Sigma]vAnnDbarDbarSomm+\[Sigma]vAnnMDbarSomm+\[Sigma]vAnnMDSomm)
,
\[Sigma]vAnnMix1Somm=Sum[Somm[Chop[\[Alpha]eff[Mrep,Drep,Dbarrep,Rep,1]],v]*(Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[1]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin]),{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnMix2Somm=Sum[Somm[Chop[\[Alpha]eff[Mrep,Drep,Dbarrep,Rep,2]],v]*(Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]MM]^2*\[Sigma]vAnn[Rep,Mrep,0,0,0,Spin]+Coefficient[mixStates[Mrep,Drep,Dbarrep,Rep,0,Spin][[2]],\[Phi]DDbar]^2*\[Sigma]vAnn[Rep,0,Drep,Dbarrep,0,Spin]),{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnDDSomm=Sum[Somm[Chop[\[Alpha]eff[0,Drep,0,Rep,0]],v]*\[Sigma]vAnn[Rep,0,Drep,0,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnDbarDbarSomm=Sum[Somm[Chop[\[Alpha]eff[0,0,Dbarrep,Rep,0]],v]*\[Sigma]vAnn[Rep,0,0,Dbarrep,0,Spin],{Rep,{1,3,5}},{Spin,{0,1}}];
\[Sigma]vAnnMDbarSomm=Sum[Somm[Chop[\[Alpha]eff[Mrep,0,Dbarrep,Rep,0]],v]*\[Sigma]vAnn[Rep,Mrep,0,Dbarrep,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
\[Sigma]vAnnMDSomm=Sum[Somm[Chop[\[Alpha]eff[Mrep,Drep,0,Rep,0]],v]*\[Sigma]vAnn[Rep,Mrep,Drep,0,0,Spin],{Rep,{2,4}},{Spin,{0,1}}];
Avg*2*(\[Sigma]vAnnMix1Somm+\[Sigma]vAnnMix2Somm+\[Sigma]vAnnDDSomm+\[Sigma]vAnnDbarDbarSomm+\[Sigma]vAnnMDbarSomm+\[Sigma]vAnnMDSomm)]];
Print["Sommerfeld done"]


(* ::Section::Initialization::Closed:: *)
(*Bound States*)


(* ::Subsection::Initialization::Closed:: *)
(*Binding energies of bound states*)


(* ::Input::Initialization:: *)
BSenergies=Module[{intlist,halflist,Mix1List,Mix2List,DDList,MDList,MDbarList},
intlist={1,3,5};
halflist={2,4};
Mix1List={};
Mix2List={};
DDList={};
MDList={};
MDbarList={};

For[index=1,index<=Length[intlist],index++,If[(\[Lambda][Mrep,Drep,Drep,intlist[[index]],1])>0,AppendTo[Mix1List,{StringJoin["Mix1R",ToString[intlist[[index]]]],Ebs[Mrep,Drep,Drep,intlist[[index]],1,1]/MDM}]]&&If[(\[Lambda][Mrep,Drep,Drep,intlist[[index]],2])>0,AppendTo[Mix2List,{StringJoin["Mix2R",ToString[intlist[[index]]]],
Ebs[Mrep,Drep,Drep,intlist[[index]],2,1]/MDM}]]&&If[(\[Lambda][0,Drep,0,intlist[[index]],0])>0,AppendTo[DDList,{StringJoin["DDR",ToString[intlist[[index]]]],
Ebs[0,Drep,0,intlist[[index]],0,1]/MDM}]]];

For[index=1,index<=Length[halflist],index++,If[(\[Lambda][Mrep,Drep,0,intlist[[index]],0])>0,AppendTo[MDList,{StringJoin["MDR",ToString[halflist[[index]]]],
Ebs[Mrep,Drep,0,intlist[[index]],0,1]/MDM}]]&&If[(\[Lambda][Mrep,0,Drep,intlist[[index]],0])>0,AppendTo[MDbarList,{StringJoin["MDbarR",ToString[halflist[[index]]]],
Ebs[Mrep,0,Drep,intlist[[index]],0,1]/MDM}]]];

BSenergies=Join[Mix1List,Mix2List,DDList,MDList,MDbarList];
BSenergies = Reverse[Sort[BSenergies,#1[[2]]<#2[[2]]&]];
(*Select only bound states with E_BS>=E_BS,max/4 and eliminate DDR5 since it doesnt annihilate*)
BSenergies=Select[BSenergies,#1[[2]]>=BSenergies[[1]][[2]]/4 &];
BSenergies=Select[BSenergies,#1[[1]]!="DDR5" &];
BSenergies];


(* ::Subsection::Initialization::Closed:: *)
(*Individual Bound State Annihilation Cross Sections*)


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState1 1-plet, n=1, \[ScriptL]=0, m=0, spin-0*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 1plet, n=1, \[ScriptL]=0, m=0, spin-0*)
(*W emission is from 3-plet to 1-plet, B emission is from 1-plet to 1-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 2-plet to 1-plet l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix11s0,\[Sigma]vBSBMix11s0,\[Sigma]vBSHMix11s0,\[Sigma]vBSMix11s0}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iH1,\[Lambda]iH2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,F1,F2,\[Sigma]vBSWMix11s0,\[Sigma]vBSBMix11s0,\[Sigma]vBSHMix11s0,\[Sigma]vBSMix11s0},If[MemberQ[BSenergies[[All,1]],"Mix1R1"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=0;
RBSbs=1;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=1;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iH1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iH2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
F1=If[\[Lambda]iH1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
F2=If[\[Lambda]iH2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vWcalcDDbar=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcDDbar=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSWMix11s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix11s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix11s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix11s0,\[Sigma]vBSBMix11s0,\[Sigma]vBSHMix11s0,\[Sigma]vBSWMix11s0+\[Sigma]vBSBMix11s0+\[Sigma]vBSHMix11s0}]];
Print["Mix11s0 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState1 1-plet, n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 1plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 3-plet to 1-plet, B emission is from 1-plet to 1-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 2-plet to 1-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix11s1,\[Sigma]vBSBMix11s1,\[Sigma]vBSHMix11s1,\[Sigma]vBSMix11s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iH1,\[Lambda]iH2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,F1,F2,\[Sigma]vBSWMix11s1,\[Sigma]vBSBMix11s1,\[Sigma]vBSHMix11s1,\[Sigma]vBSMix11s1},If[MemberQ[BSenergies[[All,1]],"Mix1R1"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=1;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=1;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iH1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iH2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
F1=If[\[Lambda]iH1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
F2=If[\[Lambda]iH2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vWcalcDDbar=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcDDbar=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSWMix11s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix11s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix11s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix11s1,\[Sigma]vBSBMix11s1,\[Sigma]vBSHMix11s1,\[Sigma]vBSWMix11s1+\[Sigma]vBSBMix11s1+\[Sigma]vBSHMix11s1}]];
Print["Mix11s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState1 3-plet, n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 3-plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 5-plet or 1-let to 3-plet, B emission is from 3-plet to 3-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 4-plet or 2-plet to 3-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix13s1,\[Sigma]vBSBMix13s1,\[Sigma]vBSHMix13s1,\[Sigma]vBSMix13s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWquint,\[Lambda]iWsing,\[Lambda]iB,\[Lambda]iHquad1,\[Lambda]iHquad2,\[Lambda]iHdoubl1,\[Lambda]iHdoubl2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,Fquad1,Fquad2,Fdoubl1,Fdoubl2,\[Sigma]vBSWMix13s1,\[Sigma]vBSBMix13s1,\[Sigma]vBSHMix13s1,\[Sigma]vBSMix13s1},If[MemberQ[BSenergies[[All,1]],"Mix1R3"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=3;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=1;
\[Lambda]iWquint=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWsing=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHquad1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iHquad2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]iHdoubl1=\[Lambda][Mrep,Drep,0,RBSbs-1,0];
\[Lambda]iHdoubl2=\[Lambda][Mrep,0,Dbarrep,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Fquad1=If[\[Lambda]iHquad1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fquad2=If[\[Lambda]iHquad2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fdoubl1=If[\[Lambda]iHdoubl1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fdoubl2=If[\[Lambda]iHdoubl2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=(If[\[Lambda]iWquint<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWquint/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWsing<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWsing/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vWcalcDDbar=(If[\[Lambda]iWquint<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWquint/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWsing<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsing/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fdoubl1,Fdoubl2,0,RBSbs-1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHdoubl1/.\[Lambda]i2->\[Lambda]iHdoubl2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vHcalcDDbar=(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fdoubl1,Fdoubl2,0,RBSbs-1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHdoubl1/.\[Lambda]i2->\[Lambda]iHdoubl2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSWMix13s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix13s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix13s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix13s1,\[Sigma]vBSBMix13s1,\[Sigma]vBSHMix13s1,\[Sigma]vBSWMix13s1+\[Sigma]vBSBMix13s1+\[Sigma]vBSHMix13s1}]];
Print["Mix13s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState1 5-plet, n=1, \[ScriptL]=0, m=0, spin-0*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 3-plet, n=1, \[ScriptL]=0, m=0, spin-0*)
(*W emission is from 7-plet or 3-let to 5-plet, B emission is from 5-plet to 5-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 6-plet or 4-plet to 5-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix15s0,\[Sigma]vBSBMix15s0,\[Sigma]vBSHMix15s0,\[Sigma]vBSMix15s0}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWsept,\[Lambda]iWtrip,\[Lambda]iB,\[Lambda]iHsix1,\[Lambda]iHsix2,\[Lambda]iHquad1,\[Lambda]iHquad2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,Fsix1,Fsix2,Fquad1,Fquad2,\[Sigma]vBSWMix15s0,\[Sigma]vBSBMix15s0,\[Sigma]vBSHMix15s0,\[Sigma]vBSMix15s0},If[MemberQ[BSenergies[[All,1]],"Mix1R5"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=0;
RBSbs=5;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=1;
\[Lambda]iWsept=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWtrip=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHsix1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iHsix2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]iHquad1=\[Lambda][Mrep,Drep,0,RBSbs-1,0];
\[Lambda]iHquad2=\[Lambda][Mrep,0,Dbarrep,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Fsix1=If[\[Lambda]iHsix1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fsix2=If[\[Lambda]iHsix2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fquad1=If[\[Lambda]iHquad1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fquad2=If[\[Lambda]iHquad2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=(If[\[Lambda]iWsept<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWsept/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWtrip<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWtrip/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vWcalcDDbar=(If[\[Lambda]iWsept<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsept/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWtrip<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWtrip/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=(\[Sigma]vH[Fsix1,Fsix2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHsix1/.\[Lambda]i2->\[Lambda]iHsix2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs-1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vHcalcDDbar=(\[Sigma]vH[Fsix1,Fsix2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHsix1/.\[Lambda]i2->\[Lambda]iHsix2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs-1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSWMix15s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix15s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix15s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix15s0,\[Sigma]vBSBMix15s0,\[Sigma]vBSHMix15s0,\[Sigma]vBSWMix15s0+\[Sigma]vBSBMix15s0+\[Sigma]vBSHMix15s0}]];
Print["Mix15s0 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState2 1-plet, n=1, \[ScriptL]=0, m=0, spin-0*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 1plet, n=1, \[ScriptL]=0, m=0, spin-0*)
(*W emission is from 3-plet to 1-plet, B emission is from 1-plet to 1-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 2-plet to 1-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix21s0,\[Sigma]vBSBMix21s0,\[Sigma]vBSHMix21s0,\[Sigma]vBSMix21s0}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iH1,\[Lambda]iH2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,F1,F2,\[Sigma]vBSWMix21s0,\[Sigma]vBSBMix21s0,\[Sigma]vBSHMix21s0,\[Sigma]vBSMix21s0},If[MemberQ[BSenergies[[All,1]],"Mix2R1"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=0;
RBSbs=1;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=2;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iH1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iH2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
F1=If[\[Lambda]iH1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
F2=If[\[Lambda]iH2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vWcalcDDbar=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcDDbar=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSWMix21s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix21s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix21s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix21s0,\[Sigma]vBSBMix21s0,\[Sigma]vBSHMix21s0,\[Sigma]vBSWMix21s0+\[Sigma]vBSBMix21s0+\[Sigma]vBSHMix21s0}]];
Print["Mix21s0 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState2 1-plet, n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 1plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 3-plet to 1-plet, B emission is from 1-plet to 1-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 2-plet to 1-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix21s1,\[Sigma]vBSBMix21s1,\[Sigma]vBSHMix21s1,\[Sigma]vBSMix21s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iH1,\[Lambda]iH2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,F1,F2,\[Sigma]vBSWMix21s1,\[Sigma]vBSBMix21s1,\[Sigma]vBSHMix21s1,\[Sigma]vBSMix21s1},If[MemberQ[BSenergies[[All,1]],"Mix2R1"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=1;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=2;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iH1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iH2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
F1=If[\[Lambda]iH1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
F2=If[\[Lambda]iH2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vWcalcDDbar=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcDDbar=\[Sigma]vH[F1,F2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iH1/.\[Lambda]i2->\[Lambda]iH2/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSWMix21s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix21s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix21s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix21s1,\[Sigma]vBSBMix21s1,\[Sigma]vBSHMix21s1,\[Sigma]vBSWMix21s1+\[Sigma]vBSBMix21s1+\[Sigma]vBSHMix21s1}]];
Print["Mix21s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState2 3-plet, n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 3-plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 5-plet or 1-let to 3-plet, B emission is from 3-plet to 3-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 4-plet or 2-plet to 3-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix23s1,\[Sigma]vBSBMix23s1,\[Sigma]vBSHMix23s1,\[Sigma]vBSMix23s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWquint,\[Lambda]iWsing,\[Lambda]iB,\[Lambda]iHquad1,\[Lambda]iHquad2,\[Lambda]iHdoubl1,\[Lambda]iHdoubl2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,Fquad1,Fquad2,Fdoubl1,Fdoubl2,\[Sigma]vBSWMix23s1,\[Sigma]vBSBMix23s1,\[Sigma]vBSHMix23s1,\[Sigma]vBSMix23s1},
If[MemberQ[BSenergies[[All,1]],"Mix2R3"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=3;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=2;
\[Lambda]iWquint=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWsing=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHquad1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iHquad2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]iHdoubl1=\[Lambda][Mrep,Drep,0,RBSbs-1,0];
\[Lambda]iHdoubl2=\[Lambda][Mrep,0,Dbarrep,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Fquad1=If[\[Lambda]iHquad1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fquad2=If[\[Lambda]iHquad2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fdoubl1=If[\[Lambda]iHdoubl1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fdoubl2=If[\[Lambda]iHdoubl2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=(If[\[Lambda]iWquint<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWquint/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWsing<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWsing/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vWcalcDDbar=(If[\[Lambda]iWquint<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWquint/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWsing<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsing/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0||\[Lambda]fbs<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fdoubl1,Fdoubl2,0,RBSbs-1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHdoubl1/.\[Lambda]i2->\[Lambda]iHdoubl2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vHcalcDDbar=(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fdoubl1,Fdoubl2,0,RBSbs-1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHdoubl1/.\[Lambda]i2->\[Lambda]iHdoubl2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSWMix23s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix23s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix23s1=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix23s1,\[Sigma]vBSBMix23s1,\[Sigma]vBSHMix23s1,\[Sigma]vBSWMix23s1+\[Sigma]vBSBMix23s1+\[Sigma]vBSHMix23s1}]];
Print["Mix23s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*MixState2 5-plet, n=1, \[ScriptL]=0, m=0, spin-0*)


(* ::Input::Initialization:: *)
(*Cross section for MM/DDbar-bound in the 3-plet, n=1, \[ScriptL]=0, m=0, spin-0*)
(*W emission is from 7-plet or 3-let to 5-plet, B emission is from 5-plet to 5-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 6-plet or 4-plet to 5-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWMix25s0,\[Sigma]vBSBMix25s0,\[Sigma]vBSHMix25s0,\[Sigma]vBSMix25s0}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWsept,\[Lambda]iWtrip,\[Lambda]iB,\[Lambda]iHsix1,\[Lambda]iHsix2,\[Lambda]iHquad1,\[Lambda]iHquad2,\[Lambda]fbs,\[Sigma]vWcalcMM,\[Sigma]vWcalcDDbar,\[Sigma]vBcalcMM,\[Sigma]vBcalcDDbar,\[Sigma]vHcalcMM,\[Sigma]vHcalcDDbar,Fsix1,Fsix2,Fquad1,Fquad2,\[Sigma]vBSWMix25s0,\[Sigma]vBSBMix25s0,\[Sigma]vBSHMix25s0,\[Sigma]vBSMix25s0},If[MemberQ[BSenergies[[All,1]],"Mix2R5"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=0;
RBSbs=5;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=Dbarrep;
Mixbs=2;
\[Lambda]iWsept=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWtrip=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHsix1=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iHsix2=\[Lambda][Mrep,0,Dbarrep,RBSbs+1,0];
\[Lambda]iHquad1=\[Lambda][Mrep,Drep,0,RBSbs-1,0];
\[Lambda]iHquad2=\[Lambda][Mrep,0,Dbarrep,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Fsix1=If[\[Lambda]iHsix1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fsix2=If[\[Lambda]iHsix2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fquad1=If[\[Lambda]iHquad1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fquad2=If[\[Lambda]iHquad2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
\[Sigma]vWcalcMM=(If[\[Lambda]iWsept<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWsept/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWtrip<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iWtrip/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vWcalcDDbar=(If[\[Lambda]iWsept<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsept/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWtrip<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWtrip/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBcalcMM=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,0,0,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBcalcDDbar=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,0,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vHcalcMM=(\[Sigma]vH[Fsix1,Fsix2,0,RBSbs+1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHsix1/.\[Lambda]i2->\[Lambda]iHsix2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs-1,RBSbs,RMbs,0,0,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vHcalcDDbar=(\[Sigma]vH[Fsix1,Fsix2,0,RBSbs+1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHsix1/.\[Lambda]i2->\[Lambda]iHsix2/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fquad1,Fquad2,0,RBSbs-1,RBSbs,0,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHquad1/.\[Lambda]i2->\[Lambda]iHquad2/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSWMix25s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vWcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vWcalcDDbar);
\[Sigma]vBSBMix25s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vBcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vBcalcDDbar);
\[Sigma]vBSHMix25s0=Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]MM]^2*(\[Sigma]vHcalcMM)+Coefficient[mixStates[RMbs,RDbs,RDbarbs,RBSbs,\[ScriptL]bs,Sbs][[Mixbs]],\[Phi]DDbar]^2*(\[Sigma]vHcalcDDbar);
{\[Sigma]vBSWMix25s0,\[Sigma]vBSBMix25s0,\[Sigma]vBSHMix25s0,\[Sigma]vBSWMix25s0+\[Sigma]vBSBMix25s0+\[Sigma]vBSHMix25s0}]];
Print["Mix25s0 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*DbarM 2-plet , n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for DbarM-bound in the 2plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 4-plet to 2-plet, B emission is from 2-plet to 2-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 3-plet or 1-plet to 2-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSDbarM2s1,\[Sigma]vBSWDbarM2s1,\[Sigma]vBSBDbarM2s1,\[Sigma]vBSHDbarM2s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iHtrip1,\[Lambda]iHtrip2,\[Lambda]iHtrip3,\[Lambda]iHsing1,\[Lambda]iHsing2,\[Lambda]iHsing3,\[Lambda]fbs,Ftrip1,Ftrip2,Ftrip3,Fsing1,Fsing2,Fsing3,\[Sigma]vWcalc,\[Sigma]vBcalc,\[Sigma]vHcalc,\[Sigma]vBSWDbarM2s1,\[Sigma]vBSBDbarM2s1,\[Sigma]vBSHDbarM2s1,\[Sigma]vBSDbarM2s1},If[MemberQ[BSenergies[[All,1]],"MDbarR2"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=2;
RMbs=Mrep;
RDbs=0;
RDbarbs=Dbarrep;
Mixbs=0;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHtrip1=\[Lambda][0,Drep,Dbarrep,RBSbs+1,0];
\[Lambda]iHtrip2=\[Lambda][0,0,Dbarrep,RBSbs+1,0];
\[Lambda]iHtrip3=\[Lambda][Mrep,0,0,RBSbs+1,0];
\[Lambda]iHsing1=\[Lambda][0,Drep,Dbarrep,RBSbs-1,0];
\[Lambda]iHsing2=\[Lambda][0,0,Dbarrep,RBSbs-1,0];
\[Lambda]iHsing3=\[Lambda][Mrep,0,0,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Ftrip1=If[\[Lambda]iHtrip1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Ftrip2=If[\[Lambda]iHtrip2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Ftrip3=If[\[Lambda]iHtrip3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
Fsing1=If[\[Lambda]iHsing1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fsing2=If[\[Lambda]iHsing2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fsing3=If[\[Lambda]iHsing3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
\[Sigma]vBSWDbarM2s1=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSBDbarM2s1=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSHDbarM2s1=(\[Sigma]vH[Ftrip1,Ftrip2,Ftrip3,RBSbs+1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHtrip1/.\[Lambda]i2->\[Lambda]iHtrip2/.\[Lambda]i3->\[Lambda]iHtrip3/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fsing1,Fsing2,Fsing3,RBSbs-1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHsing1/.\[Lambda]i2->\[Lambda]iHsing2/.\[Lambda]i3->\[Lambda]iHsing3/.\[Lambda]f->\[Lambda]fbs);
{\[Sigma]vBSWDbarM2s1,\[Sigma]vBSBDbarM2s1,\[Sigma]vBSHDbarM2s1,\[Sigma]vBSWDbarM2s1+\[Sigma]vBSBDbarM2s1+\[Sigma]vBSHDbarM2s1}]];
Print["DbarM2s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*DbarM 4-plet , n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for DbarM-bound in the 4-plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 6-plet or 2-plet to 4-plet, B emission is from 4-plet to 4-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 5-plet or 3-plet to 2-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWDbarM4s1,\[Sigma]vBSBDbarM4s1,\[Sigma]vBSHDbarM4s1,\[Sigma]vBSDbarM4s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWsix,\[Lambda]iWquad,\[Lambda]iB,\[Lambda]iHquint1,\[Lambda]iHquint2,\[Lambda]iHquint3,\[Lambda]iHtrip1,\[Lambda]iHtrip2,\[Lambda]iHtrip3,\[Lambda]fbs,Fquint1,Fquint2,Fquint3,Ftrip1,Ftrip2,Ftrip3,\[Sigma]vWcalc,\[Sigma]vBcalc,\[Sigma]vHcalc,\[Sigma]vBSWDbarM4s1,\[Sigma]vBSBDbarM4s1,\[Sigma]vBSHDbarM4s1,\[Sigma]vBSDbarM4s1},If[MemberQ[BSenergies[[All,1]],"MDbarR4"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=4;
RMbs=Mrep;
RDbs=0;
RDbarbs=Dbarrep;
Mixbs=0;
\[Lambda]iWsix=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWquad=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHquint1=\[Lambda][0,Drep,Dbarrep,RBSbs+1,0];
\[Lambda]iHquint2=\[Lambda][0,0,Dbarrep,RBSbs+1,0];
\[Lambda]iHquint3=\[Lambda][Mrep,0,0,RBSbs+1,0];
\[Lambda]iHtrip1=\[Lambda][0,Drep,Dbarrep,RBSbs-1,0];
\[Lambda]iHtrip2=\[Lambda][0,0,Dbarrep,RBSbs-1,0];
\[Lambda]iHtrip3=\[Lambda][Mrep,0,0,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Fquint1=If[\[Lambda]iHquint1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fquint2=If[\[Lambda]iHquint2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fquint3=If[\[Lambda]iHquint3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
Ftrip1=If[\[Lambda]iHtrip1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Ftrip2=If[\[Lambda]iHtrip2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Ftrip3=If[\[Lambda]iHtrip3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
\[Sigma]vBSWDbarM4s1=(If[\[Lambda]iWsix<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsix/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWquad<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWquad/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSBDbarM4s1=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSHDbarM4s1=(\[Sigma]vH[Fquint1,Fquint2,Fquint3,RBSbs+1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHquint1/.\[Lambda]i2->\[Lambda]iHquint2/.\[Lambda]i3->\[Lambda]iHquint3/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Ftrip1,Ftrip2,Ftrip3,RBSbs-1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHtrip1/.\[Lambda]i2->\[Lambda]iHtrip2/.\[Lambda]i3->\[Lambda]iHtrip3/.\[Lambda]f->\[Lambda]fbs);
{\[Sigma]vBSWDbarM4s1,\[Sigma]vBSBDbarM4s1,\[Sigma]vBSHDbarM4s1,\[Sigma]vBSWDbarM4s1+\[Sigma]vBSBDbarM4s1+\[Sigma]vBSHDbarM4s1}]];
Print["DbarM2s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*DM 2-plet , n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for DbarM-bound in the 2plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 4-plet to 2-plet, B emission is from 2-plet to 2-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 3-plet or 1-plet to 2-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWDM2s1,\[Sigma]vBSBDM2s1,\[Sigma]vBSHDM2s1,\[Sigma]vBSDM2s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iHtrip1,\[Lambda]iHtrip2,\[Lambda]iHtrip3,\[Lambda]iHsing1,\[Lambda]iHsing2,\[Lambda]iHsing3,\[Lambda]fbs,Ftrip1,Ftrip2,Ftrip3,Fsing1,Fsing2,Fsing3,\[Sigma]vWcalc,\[Sigma]vBcalc,\[Sigma]vHcalc,\[Sigma]vBSWDM2s1,\[Sigma]vBSBDM2s1,\[Sigma]vBSHDM2s1,\[Sigma]vBSDM2s1},If[MemberQ[BSenergies[[All,1]],"MDR2"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=2;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=0;
Mixbs=0;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHtrip1=\[Lambda][0,Drep,Dbarrep,RBSbs+1,0];
\[Lambda]iHtrip2=\[Lambda][0,Drep,0,RBSbs+1,0];
\[Lambda]iHtrip3=\[Lambda][Mrep,0,0,RBSbs+1,0];
\[Lambda]iHsing1=\[Lambda][0,Drep,Dbarrep,RBSbs-1,0];
\[Lambda]iHsing2=\[Lambda][0,Drep,0,RBSbs-1,0];
\[Lambda]iHsing3=\[Lambda][Mrep,0,0,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Ftrip1=If[\[Lambda]iHtrip1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Ftrip2=If[\[Lambda]iHtrip2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Ftrip3=If[\[Lambda]iHtrip3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
Fsing1=If[\[Lambda]iHsing1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fsing2=If[\[Lambda]iHsing2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fsing3=If[\[Lambda]iHsing3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
\[Sigma]vBSWDM2s1=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSBDM2s1=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSHDM2s1=(\[Sigma]vH[Ftrip1,Ftrip2,Ftrip3,RBSbs+1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHtrip1/.\[Lambda]i2->\[Lambda]iHtrip2/.\[Lambda]i3->\[Lambda]iHtrip3/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Fsing1,Fsing2,Fsing3,RBSbs-1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHsing1/.\[Lambda]i2->\[Lambda]iHsing2/.\[Lambda]i3->\[Lambda]iHsing3/.\[Lambda]f->\[Lambda]fbs);
{\[Sigma]vBSWDM2s1,\[Sigma]vBSBDM2s1,\[Sigma]vBSHDM2s1,\[Sigma]vBSWDM2s1+\[Sigma]vBSBDM2s1+\[Sigma]vBSHDM2s1}]];
Print["DM2s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*DM 4-plet , n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for DbarM-bound in the 4-plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 6-plet or 2-plet to 4-plet, B emission is from 4-plet to 4-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 5-plet or 3-plet to 2-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWDM4s1,\[Sigma]vBSBDM4s1,\[Sigma]vBSHDM4s1,\[Sigma]vBSDM4s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWsix,\[Lambda]iWquad,\[Lambda]iB,\[Lambda]iHquint1,\[Lambda]iHquint2,\[Lambda]iHquint3,\[Lambda]iHtrip1,\[Lambda]iHtrip2,\[Lambda]iHtrip3,\[Lambda]fbs,Fquint1,Fquint2,Fquint3,Ftrip1,Ftrip2,Ftrip3,\[Sigma]vWcalc,\[Sigma]vBcalc,\[Sigma]vHcalc,\[Sigma]vBSWDM4s1,\[Sigma]vBSBDM4s1,\[Sigma]vBSHDM4s1,\[Sigma]vBSDM4s1},If[MemberQ[BSenergies[[All,1]],"MDR4"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=4;
RMbs=Mrep;
RDbs=Drep;
RDbarbs=0;
Mixbs=0;
\[Lambda]iWsix=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWquad=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHquint1=\[Lambda][0,Drep,Dbarrep,RBSbs+1,0];
\[Lambda]iHquint2=\[Lambda][0,Drep,0,RBSbs+1,0];
\[Lambda]iHquint3=\[Lambda][Mrep,0,0,RBSbs+1,0];
\[Lambda]iHtrip1=\[Lambda][0,Drep,Dbarrep,RBSbs-1,0];
\[Lambda]iHtrip2=\[Lambda][0,Drep,0,RBSbs-1,0];
\[Lambda]iHtrip3=\[Lambda][Mrep,0,0,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
Fquint1=If[\[Lambda]iHquint1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Fquint2=If[\[Lambda]iHquint2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Fquint3=If[\[Lambda]iHquint3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
Ftrip1=If[\[Lambda]iHtrip1<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i1];
Ftrip2=If[\[Lambda]iHtrip2<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i2];
Ftrip3=If[\[Lambda]iHtrip3<=0,0,Fintcalc/.\[Lambda]i->\[Lambda]i3];
\[Sigma]vBSWDM4s1=(If[\[Lambda]iWsix<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsix/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWquad<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWquad/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSBDM4s1=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSHDM4s1=(\[Sigma]vH[Fquint1,Fquint2,Fquint3,RBSbs+1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHquint1/.\[Lambda]i2->\[Lambda]iHquint2/.\[Lambda]i3->\[Lambda]iHquint3/.\[Lambda]f->\[Lambda]fbs)+(\[Sigma]vH[Ftrip1,Ftrip2,Ftrip3,RBSbs-1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]/.\[Lambda]i1->\[Lambda]iHtrip1/.\[Lambda]i2->\[Lambda]iHtrip2/.\[Lambda]i3->\[Lambda]iHtrip3/.\[Lambda]f->\[Lambda]fbs);
{\[Sigma]vBSWDM4s1,\[Sigma]vBSBDM4s1,\[Sigma]vBSHDM4s1,\[Sigma]vBSWDM4s1+\[Sigma]vBSBDM4s1+\[Sigma]vBSHDM4s1}]];
Print["DM4s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*DD 1-plet, n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for DD-bound in the 1-plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 3-plet to 1-plet, B emission is from 1-plet to 1-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 2-plet to 1-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWDD1s1,\[Sigma]vBSBDD1s1,\[Sigma]vBSHDD1s1,\[Sigma]vBSDD1s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iW,\[Lambda]iB,\[Lambda]iH,\[Lambda]fbs,\[Sigma]vWcalc,\[Sigma]vBcalc,\[Sigma]vHcalc,\[Sigma]vBSWDD1s1,\[Sigma]vBSBDD1s1,\[Sigma]vBSHDD1s1,\[Sigma]vBSDD1s1},If[MemberQ[BSenergies[[All,1]],"DDR1"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=1;
RMbs=0;
RDbs=Drep;
RDbarbs=0;
Mixbs=0;
\[Lambda]iW=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iH=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Sigma]vBSWDD1s1=If[\[Lambda]iW<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iW/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSBDD1s1=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSHDD1s1=If[\[Lambda]iH<=0,0,\[Sigma]vH[Fintcalc,0,0,RBSbs+1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iH/.\[Lambda]f->\[Lambda]fbs;
{\[Sigma]vBSWDD1s1,\[Sigma]vBSBDD1s1,\[Sigma]vBSHDD1s1,\[Sigma]vBSWDD1s1+\[Sigma]vBSBDD1s1+\[Sigma]vBSHDD1s1}]];
Print["DD1s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*DD 3-plet , n=1, \[ScriptL]=0, m=0, spin-1*)


(* ::Input::Initialization:: *)
(*Cross section for DD-bound in the 3-plet, n=1, \[ScriptL]=0, m=0, spin-1*)
(*W emission is from 1-plet or 5-plet to 3-plet, B emission is from 3-plet to 3-plet, both l=1->\[ScriptL]=0*)
(*H emission is from 2-plet or 4-plet to 3-plet, l=0->\[ScriptL]=0*)
{\[Sigma]vBSWDD3s1,\[Sigma]vBSBDD3s1,\[Sigma]vBSHDD3s1,\[Sigma]vBSDD3s1}=Module[{Sbs,RBSbs,RMbs,RDbs,RDbarbs,Mixbs,\[Lambda]iWquint,\[Lambda]iWsing,\[Lambda]iB,\[Lambda]iHquad,\[Lambda]iHdoubl,\[Lambda]fbs,\[Sigma]vWcalc,\[Sigma]vBcalc,\[Sigma]vHcalc,\[Sigma]vBSWDD3s1,\[Sigma]vBSBDD3s1,\[Sigma]vBSHDD3s1,\[Sigma]vBSDD3s1},
If[MemberQ[BSenergies[[All,1]],"DDR3"]==False,{0,0,0,0},
(*Define the bound state*)
Sbs=1;
RBSbs=3;
RMbs=0;
RDbs=Drep;
RDbarbs=0;
Mixbs=0;
\[Lambda]iWquint=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs+2,Mixbs];
\[Lambda]iWsing=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs-2,Mixbs];
\[Lambda]iB=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Lambda]iHquad=\[Lambda][Mrep,Drep,0,RBSbs+1,0];
\[Lambda]iHdoubl=\[Lambda][Mrep,Drep,0,RBSbs-1,0];
\[Lambda]fbs=\[Lambda][RMbs,RDbs,RDbarbs,RBSbs,Mixbs];
\[Sigma]vBSWDD3s1=(If[\[Lambda]iWquint<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs+2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWquint/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iWsing<=0,0,\[Sigma]vW[Jintcalc,Tintcalc,RBSbs-2,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iWsing/.\[Lambda]f->\[Lambda]fbs);
\[Sigma]vBSBDD3s1=If[\[Lambda]iB<=0,0,\[Sigma]vB[Jintcalc,RBSbs,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iB/.\[Lambda]f->\[Lambda]fbs;
\[Sigma]vBSHDD3s1=(If[\[Lambda]iHquad<=0,0,\[Sigma]vH[Fintcalc,0,0,RBSbs+1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iHquad/.\[Lambda]f->\[Lambda]fbs)+(If[\[Lambda]iHdoubl<=0,0,\[Sigma]vH[Fintcalc,0,0,RBSbs-1,RBSbs,RMbs,RDbs,RDbarbs,Sbs]]/.\[Lambda]i->\[Lambda]iHdoubl/.\[Lambda]f->\[Lambda]fbs);
{\[Sigma]vBSWDD3s1,\[Sigma]vBSBDD3s1,\[Sigma]vBSHDD3s1,\[Sigma]vBSWDD3s1+\[Sigma]vBSBDD3s1+\[Sigma]vBSHDD3s1}]];
Print["DD3s1 done"]


(* ::Subsubsection::Initialization::Closed:: *)
(*Totals by gauge boson*)


(* ::Input::Initialization:: *)
\[Sigma]vBSW=\[Sigma]vBSWMix11s0+\[Sigma]vBSWMix11s1+\[Sigma]vBSWMix13s1+\[Sigma]vBSWMix15s0+\[Sigma]vBSWMix21s0+\[Sigma]vBSWMix21s1+\[Sigma]vBSWMix23s1+\[Sigma]vBSWMix25s0+\[Sigma]vBSWDbarM2s1+\[Sigma]vBSWDbarM4s1+\[Sigma]vBSWDM2s1+\[Sigma]vBSWDM4s1+\[Sigma]vBSWDD1s1+\[Sigma]vBSWDD3s1;
\[Sigma]vBSB=\[Sigma]vBSBMix11s0+\[Sigma]vBSBMix11s1+\[Sigma]vBSBMix13s1+\[Sigma]vBSBMix15s0+\[Sigma]vBSBMix21s0+\[Sigma]vBSBMix21s1+\[Sigma]vBSBMix23s1+\[Sigma]vBSBMix25s0+\[Sigma]vBSBDbarM2s1+\[Sigma]vBSBDbarM4s1+\[Sigma]vBSBDM2s1+\[Sigma]vBSBDM4s1+\[Sigma]vBSBDD1s1+\[Sigma]vBSBDD3s1;
\[Sigma]vBSH=\[Sigma]vBSHMix11s0+\[Sigma]vBSHMix11s1+\[Sigma]vBSHMix13s1+\[Sigma]vBSHMix15s0+\[Sigma]vBSHMix21s0+\[Sigma]vBSHMix21s1+\[Sigma]vBSHMix23s1+\[Sigma]vBSHMix25s0+\[Sigma]vBSHDbarM2s1+\[Sigma]vBSHDbarM4s1+\[Sigma]vBSHDM2s1+\[Sigma]vBSHDM4s1+\[Sigma]vBSHDD1s1+\[Sigma]vBSHDD3s1;
Print["Gauge boson totals done"]


(* ::Section::Initialization::Closed:: *)
(*Export cross sections*)


(* ::Input::Initialization:: *)
time=AbsoluteTime[]-t;
filename=StringJoin[ToString[Mrep],"M",ToString[Drep],"D_CrossSections_",yHval,"aH.m"];
DeleteFile[filename];
Save[filename,time,yHval,yHAbs,BSenergies,\[Sigma]vAnnTreelevel,\[Sigma]vAnnSomm,\[Sigma]vBSMix11s0,\[Sigma]vBSMix11s1,\[Sigma]vBSMix13s1,\[Sigma]vBSMix15s0,\[Sigma]vBSMix21s0,\[Sigma]vBSMix21s1,\[Sigma]vBSMix23s1,\[Sigma]vBSMix25s0,\[Sigma]vBSDbarM2s1,\[Sigma]vBSDbarM4s1,\[Sigma]vBSDM2s1,\[Sigma]vBSDM4s1,\[Sigma]vBSDD1s1,\[Sigma]vBSDD3s1,\[Sigma]vBSW,\[Sigma]vBSB,\[Sigma]vBSH];
Print["Saved"]
