(* ::Package:: *)

(* ::Input::Initialization:: *)
Off[ClebschGordan::phy]


(* ::Input::Initialization:: *)
(*Standard model parameters*)
\[Alpha]2MZactual=1/29.6;
\[Alpha]1MZactual=1/98.4;


(* ::Input::Initialization:: *)
(*SU(2) Generators up to the 13-plet*)
Do[
(*Define spin j*)
dim=2 j+1;

(*Define m values:from j to-j*)
mValues=Table[j-i,{i,0,dim-1}];

(*Initialize J+and J-as zero matrices*)
Jplus=ConstantArray[0,{dim,dim}];
Jminus=ConstantArray[0,{dim,dim}];

(*Populate J+and J-matrices*)
Do[If[i>1,mDown=mValues[[i]];
coef=Sqrt[(j-mDown) (j+mDown+1)];
Jplus[[i-1,i]]=coef;];
If[i<dim,mUp=mValues[[i]];
coef=Sqrt[(j+mUp) (j-mUp+1)];
Jminus[[i+1,i]]=coef;];,{i,1,dim}];

(*Define J1,J2,J3*)
J[dim][1]=(Jplus+Jminus)/2;
J[dim][2]=(Jplus-Jminus)/(2 I);
J[dim][3]=DiagonalMatrix[mValues];
,{j,{1/2,1,3/2,2,5/2,3,7/2,4,9/2,5,11/2,6}}]


(* ::Input::Initialization:: *)
(*Particle degrees of freedom*)
g[R_Integer]:=Which[OddQ[R],2*R,EvenQ[R],4*R];

(*Weak hypercharge*)
Y[R_Integer]:=Which[OddQ[R],0,EvenQ[R],1/2];
Ybar[R_Integer]:=Which[OddQ[R],0,EvenQ[R],-1/2];

(*Isospin for a representation of size R and opposite*)
i[R_Integer]:=(R-1)/2;
R[i_Integer]:=2i+1;
i3[R_Integer,n_Integer]:=If[n>R,Abort[],i[R]-(n-1)];

(*Quadratic casimir*)
C2R[R_Integer]:=i[R]*(i[R]+1);

(*Gauge factors for combining 2 multiplets into a representation R*)
F[R1_Integer,R2_Integer,R_Integer]:=-(1/2)*(C2R[R1]+C2R[R2]-C2R[R]);

(*Generalization of Levi Civita symbol for SU(2) tensor formalism*)
eps[d_]:=Module[{epsMat},
epsMat=Normal[SparseArray[{},{d,d}]];
For[q=0,q<d,q++,Which[EvenQ[q],epsMat[[q+1,d-q]]=1,OddQ[q],epsMat[[q+1,d-q]]=-1]];
epsMat];


(* ::Input::Initialization:: *)
(*Potentials for all possible interactions, V=-Subscript[\[Alpha], eff]/r. Temporary since we need to consider the eigenvalues for the mixed states*)
\[Alpha]efftemp[RM_Integer,RD_Integer,RDbar_Integer,R_Integer]:=
Which[
RM==0&&RD!=0&&RDbar!=0,(\[Alpha]1MZ/4-\[Alpha]2MZ*F[RD,RDbar,R]),
RM==0&&RD!=0&&RDbar==0,-(\[Alpha]1MZ/4+\[Alpha]2MZ*F[RD,RD,R]),
RM==0&&RD==0&&RDbar!=0,-(\[Alpha]1MZ/4+\[Alpha]2MZ*F[RDbar,RDbar,R]),
RM!=0&&RD==0&&RDbar==0,-\[Alpha]2MZ*F[RM,RM,R],
RM!=0&&RD!=0&&RDbar==0,-(\[Alpha]2MZ*F[RM,RD,R]+\[Alpha]H),
RM!=0&&RD==0&&RDbar!=0,-(\[Alpha]2MZ*F[RM,RDbar,R]-\[Alpha]H),
RM!=0&&RD!=0&&RDbar!=0,-\[Alpha]H];

(*Mixing matrix for the MM<->DDbar states for large \[Alpha]H*)
MixingMatrix={{\[Alpha]effMM,\[Alpha]effMix},{\[Alpha]effMix,\[Alpha]effDDbar}};

(*Final Subscript[\[Alpha], eff], same as \[Alpha]efftemp for non-mixed states (set Mix=0). For mixed, Mix=1,2 denotes which eigenvector we are choosing*)
\[Alpha]eff[RM_Integer,RD_Integer,RDbar_Integer,R_Integer,Mix_Integer]:=Which[
Mix==0,\[Alpha]efftemp[RM,RD,RDbar,R],
Mix!=0,Eigenvalues[MixingMatrix][[Mix]]/.\[Alpha]effMM->\[Alpha]efftemp[RM,0,0,R]/.\[Alpha]effMix->\[Alpha]efftemp[RM,RD,RDbar,R]/.\[Alpha]effDDbar->\[Alpha]efftemp[0,RD,RDbar,R]];

(*Mixed MM DDbar states corresponding to the diagonalized mixing matrix*)
mixStates[RM_Integer,RD_Integer,RDbar_Integer,R_Integer,\[ScriptL]_Integer,S_Integer]:=Module[{Evecs,Pmatrix},
Evecs=Eigenvectors[MixingMatrix/.\[Alpha]effMM->\[Alpha]efftemp[RM,0,0,R]/.\[Alpha]effMix->\[Alpha]efftemp[RM,RD,RDbar,R]/.\[Alpha]effDDbar->\[Alpha]efftemp[0,RD,RDbar,R]];
Pmatrix=Transpose[{Normalize[Evecs[[1]]],Normalize[Evecs[[2]]]}];
Inverse[Pmatrix] . {\[Phi]MM,\[Phi]DDbar}];

(*\[Lambda] defined by \[Lambda]=\[Alpha]/\[Alpha]2Mz to simplify cross section calculations*)
\[Lambda][RM_Integer,RD_Integer,RDbar_Integer,R_Integer,Mix_Integer]:=\[Alpha]eff[RM,RD,RDbar,R,Mix]/\[Alpha]2MZ;

(*Bound state energy*)
Ebs[RM_Integer,RD_Integer,RDbar_Integer,R_Integer,Mix_Integer,n_Integer]:=MDM*(\[Alpha]2MZ*\[Lambda][RM,RD,RDbar,R,Mix])^2/(4*n^2);
Ebs[\[Alpha]_,n_Integer]:=MDM/(4*n^2) \[Alpha]^2;


(* ::Input::Initialization:: *)
(*J-integral coefficient for BS formation through W emission*)
JCoeffWtemp[RIS_Integer,MIS_,RBS_Integer,MBS_,R1_Integer,R2_Integer,b_Integer]:=Module[{CGi,CGf},
CGi=Table[ClebschGordan[{i[R1],i3[R1,n]},{i[R2],i3[R2,m]},{i[RIS],MIS}],{n,1,R1},{m,1,R2}];
CGf=Table[ClebschGordan[{i[R1],i3[R1,n]},{i[R2],i3[R2,m]},{i[RBS],MBS}],{n,1,R1},{m,1,R2}];
1/2*Tr[CGi . Transpose[eps[R2]] . (Transpose[CGf] . J[R1][b] . eps[R1]-Transpose[J[R2][b]] . Transpose[CGf] . eps[R1])]];
JCoeffW[RIS_Integer,MIS_,RBS_Integer,MBS_,RM_Integer,RD_Integer,RDbar_Integer,b_Integer]:=Which[
RM==0&&RD!=0&&RDbar!=0,JCoeffWtemp[RIS,MIS,RBS,MBS,RD,RDbar,b],
RM==0&&RD!=0&&RDbar==0,JCoeffWtemp[RIS,MIS,RBS,MBS,RD,RD,b],
RM==0&&RD==0&&RDbar!=0,JCoeffWtemp[RIS,MIS,RBS,MBS,RDbar,RDbar,b],
RM!=0&&RD==0&&RDbar==0,JCoeffWtemp[RIS,MIS,RBS,MBS,RM,RM,b],
RM!=0&&RD!=0&&RDbar==0,JCoeffWtemp[RIS,MIS,RBS,MBS,RM,RD,b],
RM!=0&&RD==0&&RDbar!=0,JCoeffWtemp[RIS,MIS,RBS,MBS,RM,RDbar,b]];


(* ::Input::Initialization:: *)
(*T-integral coefficient for BS formation through W emission*)
TCoeffWtemp[RIS_Integer,MIS_,RBS_Integer,MBS_,R1_Integer,R2_Integer,b_Integer]:=Module[{CGi,CGf},
CGi=Table[ClebschGordan[{i[R1],i3[R1,n]},{i[R2],i3[R2,m]},{i[RIS],MIS}],{n,1,R1},{m,1,R2}];
CGf=Table[ClebschGordan[{i[R1],i3[R1,n]},{i[R2],i3[R2,m]},{i[RBS],MBS}],{n,1,R1},{m,1,R2}];
I*Sum[LeviCivitaTensor[3][[b,c,d]]*Tr[CGi . Transpose[eps[R2]] . Transpose[J[R2][d]] . Transpose[CGf] . J[R1][c] . eps[R1]],{c,1,3},{d,1,3}]];
TCoeffW[RIS_Integer,MIS_,RBS_Integer,MBS_,RM_Integer,RD_Integer,RDbar_Integer,b_Integer]:=Which[
RM==0&&RD!=0&&RDbar!=0,TCoeffWtemp[RIS,MIS,RBS,MBS,RD,RDbar,b],
RM==0&&RD!=0&&RDbar==0,TCoeffWtemp[RIS,MIS,RBS,MBS,RD,RD,b],
RM==0&&RD==0&&RDbar!=0,TCoeffWtemp[RIS,MIS,RBS,MBS,RDbar,RDbar,b],
RM!=0&&RD==0&&RDbar==0,TCoeffWtemp[RIS,MIS,RBS,MBS,RM,RM,b],
RM!=0&&RD!=0&&RDbar==0,TCoeffWtemp[RIS,MIS,RBS,MBS,RM,RD,b],
RM!=0&&RD==0&&RDbar!=0,TCoeffWtemp[RIS,MIS,RBS,MBS,RM,RDbar,b]];


(* ::Input::Initialization:: *)
(*J-integral coefficient for BS formation through B emission*)
JCoeffBtemp[RIS_Integer,MIS_,RBS_Integer,MBS_,R1_Integer,R2_Integer]:=Module[{CGi,CGf},
CGi=Table[ClebschGordan[{i[R1],i3[R1,n]},{i[R2],i3[R2,m]},{i[RIS],MIS}],{n,1,R1},{m,1,R2}];
CGf=Table[ClebschGordan[{i[R1],i3[R1,n]},{i[R2],i3[R2,m]},{i[RBS],MBS}],{n,1,R1},{m,1,R2}];
(1/2)*(Y1-Y2)*Tr[CGi . Transpose[eps[R2]] . Transpose[CGf] . eps[R1]]];
JCoeffB[RIS_Integer,MIS_,RBS_Integer,MBS_,RM_Integer,RD_Integer,RDbar_Integer]:=Which[
RM==0&&RD!=0&&RDbar!=0,JCoeffBtemp[RIS,MIS,RBS,MBS,RD,RDbar]/.{Y1->Y[RD],Y2->Ybar[RDbar]},
RM==0&&RD!=0&&RDbar==0,0,
RM==0&&RD==0&&RDbar!=0,0,
RM!=0&&RD==0&&RDbar==0,0,
RM!=0&&RD!=0&&RDbar==0,JCoeffBtemp[RIS,MIS,RBS,MBS,RM,RD]/.{Y1->Y[RM],Y2->Y[RD]},
RM!=0&&RD==0&&RDbar!=0,JCoeffBtemp[RIS,MIS,RBS,MBS,RM,RDbar]/.{Y1->Y[RM],Y2->Ybar[RDbar]}];


(* ::Input::Initialization:: *)
(*F-integral coefficients for BS formation through H emission*)
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,DM,MM]:=Sum[ClebschGordan[{i[Drep],n},{i[Mrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Mrep],n+a},{i[Mrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Drep],i[Drep]},{m,-i[Mrep],i[Mrep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,MD,DD]:=Sum[ClebschGordan[{i[Mrep],n},{i[Drep],m},{i[RIS],MIS}]*ClebschGordan[{i[Drep],n+a},{i[Drep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Mrep],i[Mrep]},{m,-i[Drep],i[Drep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,MDbar,DbarDbar]:=Sum[ClebschGordan[{i[Mrep],n},{i[Dbarrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Dbarrep],n+a},{i[Dbarrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Mrep],i[Mrep]},{m,-i[Dbarrep],i[Dbarrep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,MDbar,DDbar]:=Sum[ClebschGordan[{i[Mrep],n},{i[Dbarrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Drep],n+a},{i[Dbarrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Mrep],i[Mrep]},{m,-i[Dbarrep],i[Dbarrep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,DbarD,MD]:=Sum[ClebschGordan[{i[Dbarrep],n},{i[Drep],m},{i[RIS],MIS}]*ClebschGordan[{i[Mrep],n+a},{i[Drep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Dbarrep],i[Dbarrep]},{m,-i[Drep],i[Drep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,DD,MD]:=Sum[ClebschGordan[{i[Drep],n},{i[Drep],m},{i[RIS],MIS}]*ClebschGordan[{i[Mrep],n+a},{i[Drep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Drep],i[Drep]},{m,-i[Drep],i[Drep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,MM,DM]:=Sum[ClebschGordan[{i[Mrep],n},{i[Mrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Drep],n+a},{i[Mrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Mrep],i[Mrep]},{m,-i[Mrep],i[Mrep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,DDbar,MDbar]:=Sum[ClebschGordan[{i[Drep],n},{i[Dbarrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Mrep],n+a},{i[Dbarrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Drep],i[Drep]},{m,-i[Dbarrep],i[Dbarrep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,DbarDbar,MDbar]:=Sum[ClebschGordan[{i[Dbarrep],n},{i[Dbarrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Mrep],n+a},{i[Dbarrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Dbarrep],i[Dbarrep]},{m,-i[Dbarrep],i[Dbarrep]}];
FCoeffH[RIS_Integer,MIS_,RBS_Integer,MBS_,MM,DbarM]:=Sum[ClebschGordan[{i[Mrep],n},{i[Mrep],m},{i[RIS],MIS}]*ClebschGordan[{i[Dbarrep],n+a},{i[Mrep],m},{i[RBS],MBS}],{a,-1/2,1/2},{n,-i[Mrep],i[Mrep]},{m,-i[Mrep],i[Mrep]}];


(* ::Input::Initialization:: *)
(*Sommerfeld factor, massless limit*)
Somm[\[Alpha]_,v_]:=If[\[Alpha]==0,1,(2\[Pi]*\[Alpha]/v)/(1-E^(-2\[Pi]*\[Alpha]/v))];


(* ::Input::Initialization:: *)
(*Bound state wavefunction Rnl, and scattering state wavefunction Rpl*)
Rnl[r_,\[Alpha]_,n_Integer,\[ScriptL]_Integer]:=((\[Alpha]*MDM)/n)^(3/2) Sqrt[(n-\[ScriptL]-1)!/(2n (n+\[ScriptL])!)] E^-(r*\[Alpha]*MDM/(2n)) ((\[Alpha]*MDM*r)/(n))^\[ScriptL] LaguerreL[n-\[ScriptL]-1,2\[ScriptL]+1,(\[Alpha]*MDM*r)/(n)];
Rpl[r_,p_,\[Alpha]_,v_,\[ScriptL]_Integer]:=Sqrt[4\[Pi](2\[ScriptL]+1)Somm[\[Alpha],v]]/Gamma[2\[ScriptL]+2] E^(-I*p*r) (2p*r)^\[ScriptL] Hypergeometric1F1[\[ScriptL]+1+I*\[Alpha]/v,2\[ScriptL]+2,2I*p*r]\!\(
\*SubsuperscriptBox[\(\[Product]\), \(k = 1\), \(\[ScriptL]\)]\((\[ScriptL] - k + 1 - I*\[Alpha]/v)\)\);


(* ::Input::Initialization:: *)
(*Basis for the vector spherical harmonics*)
e[-1]=1/Sqrt[2] {1,-I,0};
e[0]={0,0,1};
e[1]=-(1/Sqrt[2]){1,I,0};


(* ::Input::Initialization:: *)
(*J_nlm overlap integral, \[ScriptL] is for final state, l is for initial*)
Jang[r_,l_Integer,n_Integer,\[ScriptL]_Integer,m_Integer,\[Lambda]i_,\[Lambda]f_]:=Rpl[r,MDM*v/2,\[Alpha]2MZ*\[Lambda]i,v,l]*Which[l==\[ScriptL]+1,-Sqrt[((\[ScriptL]+1)/(2\[ScriptL]+1))]*(D[Rnl[r,\[Alpha]2MZ*\[Lambda]f,n,\[ScriptL]],r]\[Conjugate]-\[ScriptL]/r Rnl[r,\[Alpha]2MZ*\[Lambda]f,n,\[ScriptL]]\[Conjugate])*Sum[ClebschGordan[{l,0},{1,j},{\[ScriptL],m}]*e[j],{j,-1,1}],l==\[ScriptL]-1,Sqrt[\[ScriptL]/(2\[ScriptL]+1)]*(D[Rnl[r,\[Alpha]2MZ*\[Lambda]f,n,\[ScriptL]],r]\[Conjugate]+(\[ScriptL]+1)/r Rnl[r,\[Alpha]2MZ*\[Lambda]f,n,\[ScriptL]]\[Conjugate])*Sum[ClebschGordan[{l,0},{1,j},{\[ScriptL],m}]*e[j],{j,-1,1}],l!=\[ScriptL]+1&&l!=\[ScriptL]-1,0];
Jint[l_Integer,n_Integer,\[ScriptL]_Integer,m_Integer,\[Lambda]i_,\[Lambda]f_]:=-Integrate[r^2*Jang[r,l,n,\[ScriptL],m,\[Lambda]i,\[Lambda]f],{r,0,Infinity}];


(* ::Input::Initialization:: *)
(*T_nlm overlap integral, \[ScriptL] is for final state, l is for initial*)
Tang[r_,l_Integer,n_Integer,\[ScriptL]_Integer,m_Integer,\[Lambda]i_,\[Lambda]f_]:=(\[Alpha]2MZ*MDM)/2*Integrate[Sin[\[Theta]]*Conjugate[Rnl[r,\[Alpha]2MZ*\[Lambda]f,n,\[ScriptL]]*SphericalHarmonicY[\[ScriptL],m,\[Theta],\[Phi]]]*((Sin[\[Theta]]*E^(I*\[Phi]))/Sqrt[2]
*e[-1]-(Sin[\[Theta]]*E^(-I*\[Phi]))/Sqrt[2]*e[1]+Cos[\[Theta]]*e[0])*Rpl[r,MDM*v/2,\[Alpha]2MZ*\[Lambda]i,v,l]*SphericalHarmonicY[l,0,\[Theta],\[Phi]],{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
Tint[l_Integer,n_Integer,\[ScriptL]_Integer,m_Integer,\[Lambda]i_,\[Lambda]f_]:=-Integrate[r^2*Tang[r,l,n,\[ScriptL],m,\[Lambda]i,\[Lambda]f],{r,0,Infinity}];


(* ::Input::Initialization:: *)
(*F_nlm overlap integral, \[ScriptL] is for final state, l is for initial*)
Fint[n_Integer,\[ScriptL]_Integer,m_,\[Lambda]i_,\[Lambda]f_]:=(MDM/2)*4*\[Pi]*Integrate[r^2*Conjugate[Rnl[r,\[Alpha]2MZ*\[Lambda]f,n,\[ScriptL]]]*Rpl[r,MDM*v/2,\[Alpha]2MZ*\[Lambda]i,v,\[ScriptL]],{r,0,Infinity}];


(* ::Input::Initialization:: *)
(*Cross section to form a bound state through W emission averaged over initial states, summed over spins*)
(*Valid for \[ScriptL]=0,n=0*)
\[Sigma]vW[J_,T_,RIS_Integer,RBS_Integer,RM_Integer,RD_Integer,RDbar_Integer,S_Integer]:=Avg*(2*S+1)*
If[(RM==0&&RD!=0&&RDbar==0)||(RM==0&&RD==0&&RDbar!=0)||(RM!=0&&RD==0&&RDbar==0),Which[(RM==0&&RD!=0&&RDbar==0),(1+(-1)^(S+2*i[RD]-i[RBS]))^2/4,(RM==0&&RD==0&&RDbar!=0),(1+(-1)^(S+2*i[RDbar]-i[RBS]))^2/4,(RM!=0&&RD==0&&RDbar==0),(1+(-1)^(S+2*i[RM]-i[RBS]))^2/4],1]*Sum[(8*\[Alpha]2MZ)/MDM^2*(MDM*v^2)/2*(1+(\[Alpha]2MZ^2*\[Lambda]f^2)/v^2)*2/3*Abs[(JCoeffW[RIS,MIS,RBS,MBS,RM,RD,RDbar,b]*J-TCoeffW[RIS,MIS,RBS,MBS,RM,RD,RDbar,b]*T) . (JCoeffW[RIS,MIS,RBS,MBS,RM,RD,RDbar,b]*J-TCoeffW[RIS,MIS,RBS,MBS,RM,RD,RDbar,b]*T)\[Conjugate]],{b,1,3},{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}];


(* ::Input::Initialization:: *)
(*Cross section to form a bound state through B emission averaged over initial states, summed over spins*)
(*Valid for \[ScriptL]=0,n=0*)
\[Sigma]vB[J_,RIS_Integer,RBS_Integer,RM_Integer,RD_Integer,RDbar_Integer,S_Integer]:=Avg*(2*S+1)*
Sum[(8*\[Alpha]1MZ)/MDM^2*(MDM*v^2)/2*(1+(\[Alpha]2MZ^2*\[Lambda]f^2)/v^2)*2/3 *Abs[(JCoeffB[RIS,MIS,RBS,MBS,RM,RD,RDbar]*J) . (JCoeffB[RIS,MIS,RBS,MBS,RM,RD,RDbar]*J)\[Conjugate] ],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}];


(* ::Input::Initialization:: *)
(*Cross section to form a bound state through H emission averaged over initial states, summed over spins*)
(*For states without only one or two overlap integrals set others to 0*)
(*Valid for \[ScriptL]=0,n=0*)
\[Sigma]vH[F1_,F2_,F3_,RIS_Integer,RBS_Integer,RM_Integer,RD_Integer,RDbar_Integer,S_Integer]:=Avg*(2*S+1)*(4*\[Alpha]H)/MDM^2*(MDM*v^2)/2*(1+(\[Alpha]2MZ^2*\[Lambda]f^2)/v^2)*
Which[
RM==0&&RD!=0&&RDbar!=0,Sum[Abs[F1*F1\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,MDbar,DDbar]*FCoeffH[RIS,MIS,RBS,MBS,MDbar,DDbar]\[Conjugate]]+Abs[F2*F2\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,MDbar,DDbar]*FCoeffH[RIS,MIS,RBS,MBS,MDbar,DDbar]\[Conjugate]],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}],
RM==0&&RD!=0&&RDbar==0,(1+(-1)^(S+2*i[RD]-i[RBS]))^2*Abs[F1*F1\[Conjugate]]*Sum[Abs[FCoeffH[RIS,MIS,RBS,MBS,MD,DD]*FCoeffH[RIS,MIS,RBS,MBS,MD,DD]\[Conjugate]],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}],
RM==0&&RD==0&&RDbar!=0,(1+(-1)^(S+2*i[RD]-i[RBS]))^2*Abs[F1*F1\[Conjugate]]*Sum[Abs[FCoeffH[RIS,MIS,RBS,MBS,MDbar,DbarDbar]*FCoeffH[RIS,MIS,RBS,MBS,MDbar,DbarDbar]\[Conjugate]],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}],
RM!=0&&RD==0&&RDbar==0,(1+(-1)^(S+2*i[RM]-i[RBS]))^2*Sum[Abs[F1*F1\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,DM,MM]*FCoeffH[RIS,MIS,RBS,MBS,DM,MM]\[Conjugate]]+Abs[F2*F2\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,DM,MM]*FCoeffH[RIS,MIS,RBS,MBS,DM,MM]\[Conjugate]],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}],
RM!=0&&RD!=0&&RDbar==0,Sum[2*Abs[F1*F1\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,DbarD,MD]*FCoeffH[RIS,MIS,RBS,MBS,DbarD,MD]\[Conjugate]]+(1+(-1)^(S+2*i[RD]-i[RIS]))^2*(1+(-1)^(i[RM]+i[RD]-(i[RIS]-i[RBS])))^2*Abs[F2*F2\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,DD,MD]*FCoeffH[RIS,MIS,RBS,MBS,DD,MD]\[Conjugate]]+(1+(-1)^(S+2*i[RM]-i[RIS]))^2*(1-(-1)^(i[RM]+i[RD]-(i[RIS]-i[RBS])))^2*Abs[F3*F3\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,MM,DM]*FCoeffH[RIS,MIS,RBS,MBS,MM,DM]\[Conjugate]],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}],
RM!=0&&RD==0&&RDbar!=0,Sum[2*Abs[F1*F1\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,DDbar,MDbar]*FCoeffH[RIS,MIS,RBS,MBS,DDbar,MDbar]\[Conjugate]]+(1+(-1)^(S+2*i[RDbar]-i[RIS]))^2*(1+(-1)^(i[RM]+i[RD]-(i[RIS]-i[RBS])))^2*Abs[F2*F2\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,DbarDbar,MDbar]*FCoeffH[RIS,MIS,RBS,MBS,DbarDbar,MDbar]\[Conjugate]]+(1+(-1)^(S+2*i[RM]-i[RIS]))^2*(1-(-1)^(i[RM]+i[RD]-(i[RIS]-i[RBS])))^2*Abs[F3*F3\[Conjugate]]*Abs[FCoeffH[RIS,MIS,RBS,MBS,MM,DbarM]*FCoeffH[RIS,MIS,RBS,MBS,MM,DbarM]\[Conjugate]],{MIS,-i[RIS],i[RIS]},{MBS,-i[RBS],i[RBS]}]];


(* ::Input::Initialization:: *)
(*Annihilation to SM cross sections*)
(*Only consider annihilation of states with \[ScriptL]=0*)
(*Ignoring contribution from B boson*)
\[Sigma]vAnn[RBS_Integer,RM_Integer,RD_Integer,RDbar_Integer,\[ScriptL]_Integer,S_Integer]:=If[\[ScriptL]==0,
Which[
RM==0&&RD!=0&&RDbar!=0,If[(RBS==1||RBS==3||RBS==5)&&S==0,(\[Pi]*\[Alpha]2MZ^2)/MDM^2*Sum[Abs[Sum[ClebschGordan[{i[RD],i3[RD,m]},{i[RD],i3[RD,k]},{i[RBS],MBS}]*(eps[RD] . J[RD][b] . J[RD][c]+eps[RD] . J[RD][c] . J[RD][b])[[m,k]],{m,1,RD},{k,1,RD}]]^2,{MBS,-i[RBS],i[RBS]},{b,1,3},{c,1,3}],0]+If[(RBS==1||RBS==3)&&S==1,(\[Alpha]H^2*4*\[Pi])/MDM^2,0]+If[RBS==3&&S==1,(2*\[Pi]*\[Alpha]2MZ^2)/MDM^2*25/2*Sum[Abs[Sum[ClebschGordan[{i[RD],i3[RD,m]},{i[RD],i3[RD,k]},{i[RBS],MBS}]*(eps[RD] . J[RD][b])[[m,k]],{m,1,RD},{k,1,RD}]*Sum[ClebschGordan[{i[2],i3[2,n]},{i[2],i3[2,s]},{i[RBS],MFS}]*(eps[2] . J[2][b])[[n,s]],{n,1,2},{s,1,2}]]^2,{MBS,-i[RBS],i[RBS]},{MFS,-i[RBS],i[RBS]},{b,1,3}],0],
RM==0&&RD!=0&&RDbar==0,If[(RBS==1||RBS==3)&&S==1,(4*\[Pi]*\[Alpha]H^2)/MDM^2,0],
RM==0&&RD==0&&RDbar!=0,If[(RBS==1||RBS==3)&&S==1,(4*\[Pi]*\[Alpha]H^2)/MDM^2,0],
RM!=0&&RD==0&&RDbar==0,
If[(RBS==1||RBS==3||RBS==5)&&S==0,(\[Pi]*\[Alpha]2MZ^2)/(2*MDM^2)*Sum[Abs[Sum[ClebschGordan[{i[RM],i3[RM,m]},{i[RM],i3[RM,k]},{i[RBS],MBS}]*(eps[RM] . J[RM][b] . J[RM][c]+eps[RM] . J[RM][c] . J[RM][b])[[m,k]],{m,1,RM},{k,1,RM}]]^2,{MBS,-i[RBS],i[RBS]},{b,1,3},{c,1,3}],0]+If[(RBS==1||RBS==3)&&S==1,(8*\[Alpha]H^2*\[Pi])/MDM^2,0]+If[RBS==3&&S==1,(\[Pi]*\[Alpha]2MZ^2)/MDM^2*25/2*Sum[Abs[Sum[ClebschGordan[{i[RM],i3[RM,m]},{i[RM],i3[RM,k]},{i[RBS],MBS}]*(eps[RM] . J[RM][b])[[m,k]],{m,1,RM},{k,1,RM}]*Sum[ClebschGordan[{i[2],i3[2,n]},{i[2],i3[2,s]},{i[RBS],MFS}]*(eps[2] . J[2][b])[[n,s]],{n,1,2},{s,1,2}]]^2,{MBS,-i[RBS],i[RBS]},{MFS,-i[RBS],i[RBS]},{b,1,3}],0],
RM!=0&&RD!=0&&RDbar==0,
If[(RBS==2||RBS==4)&&S==1,(2*\[Pi]*\[Alpha]2MZ*\[Alpha]H)/MDM^2*Abs[Sum[ClebschGordan[{i[RM],i3[RM,k]},{i[RD],i3[RM,m]+c},{i[RBS],MBS}]*(eps[RM] . J[RM][b])[[k,m]]+ClebschGordan[{i[RM],i3[RD,k]+c},{i[RD],i3[RD,m]},{i[RBS],MBS}]*J[RD][b][[k,m]],{k,1,RD},{m,1,RD},{b,1,3},{c,-1/2,1/2},{MBS,-i[RBS],i[RBS]}]]^2,0],
RM!=0&&RD==0&&RDbar!=0,
If[(RBS==2||RBS==4)&&S==1,(2*\[Pi]*\[Alpha]2MZ*\[Alpha]H)/MDM^2*Abs[Sum[ClebschGordan[{i[RM],i3[RM,k]},{i[RDbar],i3[RM,m]+c},{i[RBS],MBS}]*(eps[RM] . J[RM][b])[[k,m]]+ClebschGordan[{i[RM],i3[RDbar,k]+c},{i[RDbar],i3[RDbar,m]},{i[RBS],MBS}]*J[RDbar][b][[k,m]],{k,1,RDbar},{m,1,RDbar},{b,1,3},{c,-1/2,1/2},{MBS,-i[RBS],i[RBS]}]]^2,0]
],0];


(* ::Input::Initialization:: *)
(*Thermal averaging for cross sections*)
\[Sigma]vAvg[\[Sigma]v_,z_]:=(
\[Beta]min=10^-3.;(* naive \[Beta] at which thermal masses shut off the effect *)
\[Beta]max=5/Sqrt[z];
(4 z^(3/2))/Sqrt[\[Pi]] NIntegrate[\[Beta]^2  Exp[-z \[Beta]^2]\[Sigma]v/.v->2\[Beta],{\[Beta], \[Beta]min,Min[\[Beta]max,1]}]);


(* ::Input::Initialization:: *)
(*Bound state annihilation rate*)
(*First factor averages over the degrees of freedom of the bound state*)
\[CapitalGamma]ann[RBS_Integer,RM_Integer,RD_Integer,RDbar_Integer,n_Integer,\[ScriptL]_Integer,S_Integer]:=1/(RBS*(2*S+1))*\[Sigma]vAnn[RBS,RM,RD,RDbar,\[ScriptL],S]*((MDM^3*\[Alpha]eff[RM,RD,RDbar,RBS,0]^3)/(2^3*\[Pi]*n^3));


(* ::Input::Initialization:: *)
(*Branching ratio of BS to SM particles, factor of 4 in numerator is DOF^2 of the constituent particles*)
BR[\[Sigma]vAvg_,RBS_Integer,Ebs_,\[CapitalGamma]ann_,S_Integer,zf_]:=If[\[CapitalGamma]ann==0,0,(1+(\[Sigma]vAvg*MDM^3*E^(-zf*Ebs/MDM))/(2*RBS*(2*S+1)*(4*\[Pi]*zf)^(3/2)*\[CapitalGamma]ann))^-1];
BRpureMaj[\[Sigma]vAvg_,RBS_Integer,Ebs_,\[CapitalGamma]ann_,S_Integer,zf_]:=If[\[CapitalGamma]ann==0,0,(1+(\[Sigma]vAvg*MDM^3*E^(-zf*Ebs/MDM))/(2*RBS*(2*S+1)*(4*\[Pi]*zf)^(3/2)*\[CapitalGamma]ann))^-1];
BRmix[Mix_,\[Sigma]vAvg_,RBS_Integer,RM_Integer,RD_Integer,RDbar_Integer,n_Integer,\[ScriptL]_Integer,S_Integer,zf_]:=If[Mix==1||Mix==2,Coefficient[mixStates[RM,RD,RDbar,RBS,\[ScriptL],S][[Mix]],\[Phi]MM]^2*2*BR[\[Sigma]vAvg,RBS,Ebs[\[Alpha]eff[RM,RD,RDbar,RBS,Mix],n],\[CapitalGamma]ann[RBS,RM,0,0,n,\[ScriptL],S],S,zf] + Coefficient[mixStates[RM,RD,RDbar,RBS,\[ScriptL],S][[Mix]],\[Phi]DDbar]^2*BR[\[Sigma]vAvg,RBS,Ebs[\[Alpha]eff[RM,RD,RDbar,RBS,Mix],n],\[CapitalGamma]ann[RBS,0,RD,RDbar,n,\[ScriptL],S],S,zf],Abort[]];


(* ::Input::Initialization:: *)
(*Tree level annihilation rates to SM*)
\[Sigma]vTreeAnn[RM_Integer,RD_Integer]:=Which[RM!=0&&RD==0,(16*\[Pi]^2)/(256*\[Pi]*MDM^2*2*RM)*\[Alpha]2MZ^2*(2*RM^4+17*RM^2-19),
RM==0&&RD!=0,(16*\[Pi]^2)/(256*\[Pi]*MDM^2*4*RD)*(\[Alpha]2MZ^2*(2*RD^4+17*RD^2-19)+4*(1/2)^2*\[Alpha]1MZ^2*(41+8*(1/2)^2)+16*\[Alpha]2MZ*\[Alpha]1MZ*(1/2)^2*(RD^2-1))];
\[Sigma]vTreeEff[RM_Integer,RD_Integer]:=Which[RM!=0&&RD==0,1 /RM^2*\[Sigma]vTreeAnn[RM,0],RM==0&&RD!=0,1 /RD^2*\[Sigma]vTreeAnn[0,RD],RM!=0&&RD!=0,1 /(RM+RD)^2*(\[Sigma]vTreeAnn[RM,0]+\[Sigma]vTreeAnn[0,RD])];
