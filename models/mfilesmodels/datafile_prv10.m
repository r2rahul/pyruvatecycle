function [output] = datafile_prv10(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pyruvate cycle Model
% 
% [output] = datafile_prv10(statevector) => time=0
% [output] = datafile_prv10(time,statevector)
% 
% output: structure containing information about the values of the variables
%         and reaction rates for given time and state information.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameterValuesNew_ALL = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,
	time_vector = [0];
	statevector = varargin{1};
elseif nargin == 2,
	time_vector = varargin{1};
	statevector = varargin{2};
elseif nargin == 3,
	time_vector = varargin{1};
	statevector = varargin{2};
	parameterValuesNew_ALL = varargin{3};
else
	error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableValuesTotal = [];
reactionValuesTotal = [];
if length(time_vector) == 1,
	statevector = statevector(:)';
	[variableValuesTotal, reactionValuesTotal] = getValues(time_vector,statevector,parameterValuesNew_ALL);
else
	for k = 1:length(time_vector),
		if ~isempty(parameterValuesNew_ALL),
			[variableValues, reactionValues] = getValues(time_vector(k),statevector(k,:),parameterValuesNew_ALL(k,:));
		else
			[variableValues, reactionValues] = getValues(time_vector(k),statevector(k,:),[]);
		end
		variableValuesTotal = [variableValuesTotal; variableValues];
		reactionValuesTotal = [reactionValuesTotal; reactionValues];
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT OUTPUT VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output.time = time_vector;
output.states = {'glc_c','F6P_c','FBP_c','GAP_c','DPG_c','PEP_c','LAC_c','PYR_c','MAL_c','CIT_c','ICIT_c','AKG_c','OAA_c','NADPH_c','PYR_m','ACO_m','CIT_m','ICIT_m','AKG_m','SCOA_m','SUC_m','FUM_m','MAL_m','OAA_m',};
output.statevalues = statevector(:,1:24);
output.variables = {'NADP_c','ATP','ADP',};
output.variablevalues = variableValuesTotal;
output.reactions = {'J0entry_c','J1gk_c','J2pfk_c','J3fba_c','J4gapd_c','J5pgp_c','J6pk_c','J7ldh_c','J9pyr_s','J10cit_s','J11icit_s','J12akhmal_s','J13malh_s','J14nadph_c','J15citl_c','J16mdh_c','J17acon_c','J18isod_c','J19me_c','J20pdh_m','J21pc_m','J22cs_m','J23ac_m','J24icd_m','J25akg_m','J26sco_m','J27sdh_m','J28fum_m','J29mdh_m','J30me_m','akgflow','lacsink',};
output.reactionvalues = reactionValuesTotal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear global variables used for delay handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global delaybase_tpc26a34*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getValues FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variableValues, reactionValues] = getValues(time_local,statevector,parameterValuesNew)
global time
variableValues = [];
reactionValues = [];

time = time_local;

% STATES
glc_c = statevector(1);
F6P_c = statevector(2);
FBP_c = statevector(3);
GAP_c = statevector(4);
DPG_c = statevector(5);
PEP_c = statevector(6);
LAC_c = statevector(7);
PYR_c = statevector(8);
MAL_c = statevector(9);
CIT_c = statevector(10);
ICIT_c = statevector(11);
AKG_c = statevector(12);
OAA_c = statevector(13);
NADPH_c = statevector(14);
PYR_m = statevector(15);
ACO_m = statevector(16);
CIT_m = statevector(17);
ICIT_m = statevector(18);
AKG_m = statevector(19);
SCOA_m = statevector(20);
SUC_m = statevector(21);
FUM_m = statevector(22);
MAL_m = statevector(23);
OAA_m = statevector(24);
% PARAMETERS
if isempty(parameterValuesNew),
	inglc = 0.003;
	Vr = 20;
	NADPtot = 0.0005;
	v0_1 = 2.7271e-06;
	v0_2 = 101.054;
	v1_1 = 0.00013424;
	v1_2 = 0.00030118;
	v1_3 = 1.34;
	v2_1 = 3.16667e-05;
	v2_2 = 0.0089;
	v2_3 = 0.9;
	v3_1 = 7.97833e-05;
	v3_2 = 4e-06;
	v4_1 = 0.001;
	v4_2 = 0.00032;
	v4_3 = 1.5;
	v5_1 = 3.33e-05;
	v5_2 = 8e-06;
	v6_1 = 5.33e-05;
	v6_2 = 0.00015;
	v6_3 = 2.9;
	v7_1 = 29.0969;
	v7_eq = 21.121;
	v7_2 = 448491;
	v7_3 = 449.094;
	v9_1 = 3.7674e-08;
	v9_2 = 0.004;
	v9_3 = 49.2637;
	v9_4 = 187.379;
	v10_1 = 32514;
	v10_2 = 84267;
	v10_3 = 0.00013;
	v10_4 = 0.00044;
	v10_5 = 0.00033;
	v10_6 = 4.18e-05;
	v12_1 = 5811.9;
	v12_2 = 6739.9;
	v12_3 = 0.0003;
	v12_4 = 0.0007;
	v12_5 = 0.0014;
	v12_6 = 0.00017;
	v13_1 = 0.8512;
	v13_2 = 0.2721;
	v13_3 = 1388.9;
	v13_4 = 1111.1;
	v14_1 = 1e-06;
	v14_2 = 1e-05;
	v15_1 = 0.8333;
	v15_2 = 0.00035714;
	v16_1 = 0.0036;
	v16_2 = 50000;
	v16_3 = 1111.1;
	v16_4 = 16667;
	v17_1 = 0.0158;
	v17_2 = 0.1104;
	v17_3 = 9090.9;
	v17_4 = 2000;
	v18_1 = 1.1521e+06;
	v18_2 = 5.4825e+06;
	v18_3 = 2.30415e+10;
	v18_4 = 142857;
	v18_5 = 161290;
	v18_6 = 1.09649e+11;
	v18_7 = 416667;
	v18_8 = 263158;
	v19_1 = 913070;
	v19_eq = 1000;
	v19_2 = 0.00012;
	v19_3 = 1.39e-06;
	v19_4 = 0.0048;
	v19_5 = 5.3e-06;
	v20_1 = 1.41667e-07;
	v20_2 = 3.5e-05;
	v20_3 = 2e-05;
	v21_1 = 1.7163e-07;
	v21_2 = 27;
	v21_3 = 1.2944e-05;
	v21_4 = 0.00018413;
	v21_5 = 2.22e-05;
	v21_6 = 1.1333e-06;
	v21_7 = 0.0926;
	v21_8 = 0.48787;
	v21_9 = 0.036228;
	v21_10 = 1.0444;
	v22_1 = 50.85;
	v22_2 = 2.5e+10;
	v22_3 = 295000;
	v22_4 = 120000;
	v23_1 = 0.0158;
	v23_2 = 0.1104;
	v23_3 = 9090.9;
	v23_4 = 2000;
	v24_1 = 0.11263;
	v24_2 = 82222;
	v24_3 = 2777.8;
	v24_4 = 0.63889;
	v24_5 = 0.21667;
	v24_6 = 1.7778;
	v25_1 = 0.0311;
	v25_2 = 1.4546e+09;
	v25_3 = 1.4546e+09;
	v25_4 = 24691;
	v26_1 = 3.2512e-05;
	v26_2 = 6.4;
	v26_3 = 4.1876;
	v26_4 = 37.5348;
	v26_5 = 1478.2;
	v26_6 = 9509.6;
	v26_7 = 236.711;
	v27_1 = 2.941;
	v27_2 = 11.5344;
	v27_3 = 449.661;
	v27_4 = 199130;
	v27_5 = 1.7921e+08;
	v28_1 = 13.9024;
	v28_2 = 13.9024;
	v28_3 = 4e+06;
	v28_4 = 880000;
	v29_1 = 0.0016;
	v29_2 = 2.5563e-10;
	v29_3 = 33263;
	v29_4 = 14.0754;
	v29_5 = 2.7429e+06;
	v30_1 = 157.812;
	v30_eq = 1000;
	v30_2 = 500;
	v30_3 = 227.273;
	v15_3 = 8333.3;
	v15_4 = 5000;
	flow = 0;
	sink = 0;
else
	inglc = parameterValuesNew(1);
	Vr = parameterValuesNew(2);
	NADPtot = parameterValuesNew(3);
	v0_1 = parameterValuesNew(4);
	v0_2 = parameterValuesNew(5);
	v1_1 = parameterValuesNew(6);
	v1_2 = parameterValuesNew(7);
	v1_3 = parameterValuesNew(8);
	v2_1 = parameterValuesNew(9);
	v2_2 = parameterValuesNew(10);
	v2_3 = parameterValuesNew(11);
	v3_1 = parameterValuesNew(12);
	v3_2 = parameterValuesNew(13);
	v4_1 = parameterValuesNew(14);
	v4_2 = parameterValuesNew(15);
	v4_3 = parameterValuesNew(16);
	v5_1 = parameterValuesNew(17);
	v5_2 = parameterValuesNew(18);
	v6_1 = parameterValuesNew(19);
	v6_2 = parameterValuesNew(20);
	v6_3 = parameterValuesNew(21);
	v7_1 = parameterValuesNew(22);
	v7_eq = parameterValuesNew(23);
	v7_2 = parameterValuesNew(24);
	v7_3 = parameterValuesNew(25);
	v9_1 = parameterValuesNew(26);
	v9_2 = parameterValuesNew(27);
	v9_3 = parameterValuesNew(28);
	v9_4 = parameterValuesNew(29);
	v10_1 = parameterValuesNew(30);
	v10_2 = parameterValuesNew(31);
	v10_3 = parameterValuesNew(32);
	v10_4 = parameterValuesNew(33);
	v10_5 = parameterValuesNew(34);
	v10_6 = parameterValuesNew(35);
	v12_1 = parameterValuesNew(36);
	v12_2 = parameterValuesNew(37);
	v12_3 = parameterValuesNew(38);
	v12_4 = parameterValuesNew(39);
	v12_5 = parameterValuesNew(40);
	v12_6 = parameterValuesNew(41);
	v13_1 = parameterValuesNew(42);
	v13_2 = parameterValuesNew(43);
	v13_3 = parameterValuesNew(44);
	v13_4 = parameterValuesNew(45);
	v14_1 = parameterValuesNew(46);
	v14_2 = parameterValuesNew(47);
	v15_1 = parameterValuesNew(48);
	v15_2 = parameterValuesNew(49);
	v16_1 = parameterValuesNew(50);
	v16_2 = parameterValuesNew(51);
	v16_3 = parameterValuesNew(52);
	v16_4 = parameterValuesNew(53);
	v17_1 = parameterValuesNew(54);
	v17_2 = parameterValuesNew(55);
	v17_3 = parameterValuesNew(56);
	v17_4 = parameterValuesNew(57);
	v18_1 = parameterValuesNew(58);
	v18_2 = parameterValuesNew(59);
	v18_3 = parameterValuesNew(60);
	v18_4 = parameterValuesNew(61);
	v18_5 = parameterValuesNew(62);
	v18_6 = parameterValuesNew(63);
	v18_7 = parameterValuesNew(64);
	v18_8 = parameterValuesNew(65);
	v19_1 = parameterValuesNew(66);
	v19_eq = parameterValuesNew(67);
	v19_2 = parameterValuesNew(68);
	v19_3 = parameterValuesNew(69);
	v19_4 = parameterValuesNew(70);
	v19_5 = parameterValuesNew(71);
	v20_1 = parameterValuesNew(72);
	v20_2 = parameterValuesNew(73);
	v20_3 = parameterValuesNew(74);
	v21_1 = parameterValuesNew(75);
	v21_2 = parameterValuesNew(76);
	v21_3 = parameterValuesNew(77);
	v21_4 = parameterValuesNew(78);
	v21_5 = parameterValuesNew(79);
	v21_6 = parameterValuesNew(80);
	v21_7 = parameterValuesNew(81);
	v21_8 = parameterValuesNew(82);
	v21_9 = parameterValuesNew(83);
	v21_10 = parameterValuesNew(84);
	v22_1 = parameterValuesNew(85);
	v22_2 = parameterValuesNew(86);
	v22_3 = parameterValuesNew(87);
	v22_4 = parameterValuesNew(88);
	v23_1 = parameterValuesNew(89);
	v23_2 = parameterValuesNew(90);
	v23_3 = parameterValuesNew(91);
	v23_4 = parameterValuesNew(92);
	v24_1 = parameterValuesNew(93);
	v24_2 = parameterValuesNew(94);
	v24_3 = parameterValuesNew(95);
	v24_4 = parameterValuesNew(96);
	v24_5 = parameterValuesNew(97);
	v24_6 = parameterValuesNew(98);
	v25_1 = parameterValuesNew(99);
	v25_2 = parameterValuesNew(100);
	v25_3 = parameterValuesNew(101);
	v25_4 = parameterValuesNew(102);
	v26_1 = parameterValuesNew(103);
	v26_2 = parameterValuesNew(104);
	v26_3 = parameterValuesNew(105);
	v26_4 = parameterValuesNew(106);
	v26_5 = parameterValuesNew(107);
	v26_6 = parameterValuesNew(108);
	v26_7 = parameterValuesNew(109);
	v27_1 = parameterValuesNew(110);
	v27_2 = parameterValuesNew(111);
	v27_3 = parameterValuesNew(112);
	v27_4 = parameterValuesNew(113);
	v27_5 = parameterValuesNew(114);
	v28_1 = parameterValuesNew(115);
	v28_2 = parameterValuesNew(116);
	v28_3 = parameterValuesNew(117);
	v28_4 = parameterValuesNew(118);
	v29_1 = parameterValuesNew(119);
	v29_2 = parameterValuesNew(120);
	v29_3 = parameterValuesNew(121);
	v29_4 = parameterValuesNew(122);
	v29_5 = parameterValuesNew(123);
	v30_1 = parameterValuesNew(124);
	v30_eq = parameterValuesNew(125);
	v30_2 = parameterValuesNew(126);
	v30_3 = parameterValuesNew(127);
	v15_3 = parameterValuesNew(128);
	v15_4 = parameterValuesNew(129);
	flow = parameterValuesNew(130);
	sink = parameterValuesNew(131);
end

% VARIABLES
NADP_c = NADPtot-NADPH_c;
atpeq=(0.44444*(inglc-1e-3)+3e-3);
adpeq=(-0.066667*(inglc-1e-3)+1.2e-3);
ATP=piecewiseATP(inglc, atpeq);
ADP=piecewiseADP(inglc, adpeq);
% REACTIONS
J0entry_c = v0_1*(inglc-glc_c)/(1+(inglc+glc_c)*v0_2+(v0_2)^2*inglc*glc_c);
J1gk_c = v1_1*glc_c^v1_3/(v1_2^v1_3+glc_c^v1_3);
J2pfk_c = v2_1*F6P_c^v2_3/(v2_2^v2_3+F6P_c^v2_3);
J3fba_c = v3_1*FBP_c/(v3_2+FBP_c);
J4gapd_c = v4_1*GAP_c^v4_3/(v4_2^v4_3+GAP_c^v4_3);
J5pgp_c = v5_1*DPG_c/(v5_2+DPG_c);
J6pk_c = v6_1*PEP_c^v6_3/(v6_2^v6_3+PEP_c^v6_3);
J7ldh_c = v7_1*(PYR_c-(LAC_c/v7_eq))/(1+v7_2*PYR_c+v7_3*LAC_c);
J9pyr_s = ((v9_1*PYR_m-PYR_c*v9_2)/(1+v9_3*PYR_m+PYR_c*v9_4));
J10cit_s = ((CIT_c*MAL_m*v10_1-MAL_c*CIT_m*v10_2)/(1+CIT_c/v10_3+MAL_m/v10_4+MAL_c/v10_5+CIT_m/v10_6+CIT_c*MAL_m/(v10_3*v10_4)+MAL_c*CIT_m/(v10_5*v10_6)+MAL_m*CIT_m/(v10_4*v10_6)+CIT_c*MAL_c/(v10_3*v10_5)));
J11icit_s = ((ICIT_c*MAL_m*v10_1-MAL_c*ICIT_m*v10_2)/(1+ICIT_c/v10_3+MAL_m/v10_4+MAL_c/v10_5+ICIT_m/v10_6+ICIT_c*MAL_m/(v10_3*v10_4)+MAL_c*ICIT_m/(v10_5*v10_6)+MAL_m*ICIT_m/(v10_4*v10_6)+ICIT_c*MAL_c/(v10_3*v10_5)));
J12akhmal_s = ((AKG_m*MAL_c*v12_1-MAL_m*AKG_c*v12_2)/(1+AKG_m/v12_3+MAL_c/v12_4+MAL_m/v12_5+AKG_c/v12_6+AKG_m*MAL_c/(v12_3*v12_4)+MAL_m*AKG_c/(v12_5*v12_6)+MAL_c*AKG_c/(v12_4*v12_6)+AKG_m*MAL_m/(v12_3*v12_5)));
J13malh_s = ((MAL_m*v13_1-MAL_c*v13_2)/(1+v13_3*MAL_m+MAL_c*v13_4));
J14nadph_c = v14_1*NADPH_c/(v14_2+NADPH_c);
J15citl_c = v15_1*(CIT_c-v15_2*OAA_c)/(1+v15_3*CIT_c+v15_4*OAA_c);
J16mdh_c = v16_1*(MAL_c-OAA_c*v16_2)/(1+v16_3*MAL_c+v16_4*OAA_c);
J17acon_c = ((v17_1*CIT_c-v17_2*ICIT_c)/(v17_3*ICIT_c+v17_4*CIT_c+1));
J18isod_c = ((v18_1*ICIT_c*NADP_c)/(v18_3*ICIT_c*NADP_c+v18_4*NADP_c+v18_5*ICIT_c+1)-(v18_2*AKG_c*NADPH_c)/(v18_6*AKG_c*NADPH_c+v18_7*NADPH_c+v18_8*AKG_c+1));
J19me_c = (v19_1)*((MAL_c*NADP_c-(PYR_c*NADPH_c/v19_eq))/(1+MAL_c/v19_2+NADP_c/v19_3+PYR_c/v19_4+NADPH_c/v19_5+MAL_c*NADP_c/(v19_2*v19_3)+PYR_c*NADPH_c/(v19_4*v19_5)+MAL_c*NADPH_c/(v19_2*v19_5)+PYR_c*NADP_c/(v19_4*v19_3)));
J20pdh_m = (v20_1*PYR_m)/(v20_2*(1+ACO_m/v20_3)+PYR_m);
J21pc_m = (v21_1*(ADP*OAA_m-v21_2*ATP*PYR_m)/(v21_3*PYR_m+v21_4*OAA_m+v21_5*ATP+v21_6*ADP+v21_7*ATP*PYR_m+v21_8*ATP*OAA_m+v21_9*ADP*PYR_m+v21_10*ADP*OAA_m));
J22cs_m = (v22_1*ACO_m*OAA_m)/(v22_2*ACO_m*OAA_m+v22_3*OAA_m+v22_4*ACO_m+1);
J23ac_m = ((v23_1*CIT_m-v23_2*ICIT_m)/(v23_3*ICIT_m+v23_4*CIT_m+1));
J24icd_m = (v24_1*(ICIT_m*ICIT_m+v24_2*ICIT_m)/(v24_3*ICIT_m*ICIT_m+v24_4*ICIT_m+v24_5*ADP+v24_6*ADP*ICIT_m+1));
J25akg_m = (v25_1*AKG_m/(1+v25_2*AKG_m+v25_3*SCOA_m+v25_4*AKG_m*SCOA_m));
J26sco_m = (v26_1*(v26_3*SCOA_m-SUC_m)*((v26_2*SUC_m+1))/(1+v26_4*SCOA_m+v26_5*SUC_m+v26_6*SUC_m*SUC_m+v26_7*SCOA_m*SUC_m));
J27sdh_m = (v27_1*(SUC_m-FUM_m*v27_2)/(1+v27_3*SUC_m+v27_4*FUM_m+v27_5*SUC_m*FUM_m));
J28fum_m = ((v28_1*FUM_m-v28_2*MAL_m)/(v28_3*MAL_m+v28_4*FUM_m+1));
J29mdh_m = ((v29_1*MAL_m-v29_2*OAA_m)/(1+MAL_m*v29_3+v29_4*OAA_m+v29_5*MAL_m*OAA_m));
J30me_m = v30_1*(MAL_m-(PYR_m/v30_eq))/(1+v30_2*MAL_m+v30_3*PYR_m);
akgflow = AKG_c*flow;
lacsink = sink*LAC_c;
% OUTPUT
variableValues(1) = NADP_c;
variableValues(2) = ATP;
variableValues(3) = ADP;
reactionValues(1) = J0entry_c;
reactionValues(2) = J1gk_c;
reactionValues(3) = J2pfk_c;
reactionValues(4) = J3fba_c;
reactionValues(5) = J4gapd_c;
reactionValues(6) = J5pgp_c;
reactionValues(7) = J6pk_c;
reactionValues(8) = J7ldh_c;
reactionValues(9) = J9pyr_s;
reactionValues(10) = J10cit_s;
reactionValues(11) = J11icit_s;
reactionValues(12) = J12akhmal_s;
reactionValues(13) = J13malh_s;
reactionValues(14) = J14nadph_c;
reactionValues(15) = J15citl_c;
reactionValues(16) = J16mdh_c;
reactionValues(17) = J17acon_c;
reactionValues(18) = J18isod_c;
reactionValues(19) = J19me_c;
reactionValues(20) = J20pdh_m;
reactionValues(21) = J21pc_m;
reactionValues(22) = J22cs_m;
reactionValues(23) = J23ac_m;
reactionValues(24) = J24icd_m;
reactionValues(25) = J25akg_m;
reactionValues(26) = J26sco_m;
reactionValues(27) = J27sdh_m;
reactionValues(28) = J28fum_m;
reactionValues(29) = J29mdh_m;
reactionValues(30) = J30me_m;
reactionValues(31) = akgflow;
reactionValues(32) = lacsink;
return
