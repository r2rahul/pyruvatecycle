%%
clc;
%order of metabolite
%pyr,cit,akg,suc,fum,mal,lac
%%****Pillai LG*************
plg=[9.6e-05,3.2e-05,5.71e-06,7.54e-05,2.63e-05,3.89e-05,0.000107];
%****Pillai HG*************
phg=[0.000496,0.00019086,0.00024457,0.00018629,0.00012,0.00026171,0.000408];
%****Pillai Knockdown LG*************
kplg=[8.8e-05,5.7143e-06,1.1429e-06,3.4286e-05,1.2571e-05,1.3714e-05,1.8667e-05];
%****Pillai Knockdown HG*************
kphg=[0.00017067,9.7143e-05,2.4e-05,5.9429e-05,3.8857e-05,8.8e-05,0.00024];
%%****Ronnebaum Knockdown LG*************
%reduction in ICDc Vmax 0.609 param37 index 38
krlg=[6.4810e-04,3.92E-004,3.36E-005,4.80E-005,1.60E-005,1.20E-004,4.96E-004];
%****Ronnebaum Knockdown HG*************
krhg=[2.7747e-03,2.42E-003,6.00E-004,1.07E-004,7.20E-005,1.31E-003,2.08E-003];
xobs1=1e-2*[3.44E-004,2.08E-004,3.52E-005,6.40E-005,1.60E-005,8.80E-005;2.69E-003,1.98E-003,7.62E-004,1.28E-004,1.22E-004,1.50E-003];
er=[5.950,0.711,0.150,0.187,0.15,0.187,0.374,10;35.70,2.36,1.31,0.299,0.299,2.10,1.23,12; ...
    9.520,1.16,0.112,0.187,0.15,0.187,0.486,9;32.10,3.17,0.785,0.187,0.075,1.38,2.77,12];
xobs=1e-2*[3.44E-004,2.08E-004,3.52E-005,6.40E-005,1.60E-005,8.80E-005,3.36E-004,0.0000248*1e2;2.69E-003,1.98E-003,7.62E-004,1.28E-004,1.22E-004,1.50E-003,7.90E-004,0.00004848*1e2;...
    6.4810e-04,3.92E-004,3.36E-005,4.80E-005,1.60E-005,1.20E-004,4.96E-004,0.00002608*1e2;2.7747e-03,2.42E-003,6.00E-004,1.07E-004,7.20E-005,1.31E-003,2.08E-003,0.0000328*1e2];
sc=abs(xobs(1,1:7)./(1e-1*plg));
%load complete_data.mat
%load ARNT1.mat
%%
NADPtot =5e-4;
%error data
ronnebaumICDcF7A_control_3_error = [0.187,0.15,0.150,0.187,0.374,0.711,5.95];
ronnebaumICDcF7A_control_12_error = [0.299,0.299,1.31,2.10,1.23,2.36,35.7];
ronnebaumICDcF7A_si_3_error = [0.187,0.15,0.112,0.187,0.486,1.16,9.52];
ronnebaumICDcF7A_si_12_error = [0.187,0.075,0.785,1.38,2.77,3.17,32.1];
format shortG;
%+++++++++CIC data++++++++++++++++++++++++++
cwtot=[1.08800000000000e-05, 5.44000000000000e-05];
cktot=[1.088e-05, 3.392e-05];
cwcyt=[1.024e-05,4.24e-05];
ckcyt=[8.32e-06,2.4e-05];
cwmit=[1.408e-05, 4.512e-05];
ckmit=[1.072e-05,3.408e-05];
fluxname={'J0entry','J1gk','J2pfk','J3fba','J4gapd','J5pgp','J6pk','J7ldh','LACsink','J9pyr_s','J10cit_s','J11icit_s','J12akhmal_s',...
    'J13malh_s','J14nadph','J15citl_c','J16mdh_c','J17acon_c','J18isod_c','J19me_c','J20pdh','J21pc','J22cs','J23ac','J24icd',...
    'J25akg','J26sco','J27sdh','J28fum','J29mdh_m','J30me_m','Akgsink'}; %

%==================================================================
%% Start the simulation Experiment
points=[2.5e-3,2.8e-3,3e-3,(4:16)*1e-3,16.7e-3,(17:22)*1e-3]; % Glucose input 4e-6 22e-6 4e-6 22e-6
temp=ics();
numcases=6; % Number of different Knockdowns
% for storing all the steady states
X = repmat(0, [length(points),length(temp),numcases]);
J = repmat(0, [length(points),length(fluxname),numcases]);
%===============Begining of Simulation=============================
for curcase = 1:numcases
    va=prv10('parametervalues');  % re-load model parameters for each case
    indFix=[4     5    46    47];
    fixVal=[8.16980648140150e-06,96.1495175230328,0.00610074140147127,0.672751542742064];
    va(indFix)=fixVal;
    indFix2=[76,115,116,75,42,82,43,85,86,36,89,31,84,72,73,74,55,30,66,69,111,61,64,78,37,87,48,27];
    fixVal2=[0.658423582094916,40.0785499825684,1381.85063062883,4.40120829671354e-08,84.7669746283260,45.2364109830938,0.225353473680178,5047.00324480593,885015928625.771,560434.876696313,0.0160755731174134,842.670000000000,1.53331642398236,1.97875531301938e-06,0.00312887994502887,1.68969778357061e-06,8.48149874150939,425172.536953913,9133.86531089249,4.48999225337987e-05,2.05645551605442,12265237.6651391,12843170.6389473,0.0117243568664210,49855.7047804443,579409.968993638,0.00913865562640276,0.0906028347792317];
    va(indFix2)=fixVal2;
    %va(130) = 1.46e-9;
    %va(131) = 8.89e-08;
    switch curcase % adjust params for each case
        case 1 % control
            display('*** Generating Continuation Diagram for Control Case ***')
            
        case 2 % knockdown of ICDc by Ad-siICDc
            
            display('*** Generating Continuation Diagram for Ad-siICDc Case ***')
            reduce = .609; % compute how much enzyme activity (Vmax) is reduced
            display(['ICDc Vmax reduced by ' num2str(100*(1-reduce)) '%']);
            va(58)=.609*va(58);
            va(59)=.609*va(59);
            
        case 3
            
            disp('Citrate Isocitrate Knockdown')
            va(30)=(.37)*va(30);
            va(31)=(.37)*va(31);
            
        case 4
            disp('Citrate Lyase Knockdown')
            va(48)=(1-.75)*va(48);
        case 5
            disp('Malic Enzyme Knock Down')
            va(66)=(1-.75)*va(66);
        case 6
            disp('Pyruvate Carboxylase Knock down')
            va(75)=(1-0.6)*va(75);
    end
    
    count=0;
    for point = points
        count=count+1;
        
        va(1) = point; % Input glucose level
        
        
        x0 = ics(); % Defining initial conditions
        
        % Run the simulation to equilibrium
        tic
        display(['Starting Simulation... ',int2str(count), '/', int2str(length(points))])
        [x,t] = v10solver(x0,va,1e5,2);
        disp('Simulation Complete.')
        toc
        
        % save simulation experiment data
        X(count,:,curcase) = x;
        flux=fsolprv10(0,x,va,1);
        J(count,:,curcase) = flux;
        %=========================================================================================================================
    end
end
if true
    %pyr,cit,akg,suc,fum,mal,lac
    h1=figure(1);
    [amt,r]=min(abs(points-3e-3));  % takes point closest to 3 mM
    simval1=X(r,[15,17,19,21,22,23],1);
    simval2=X(r,[8,10,12,9,7],1);
    LGsimval=[simval1(1)+simval2(1) simval1(2)+simval2(2) simval1(3)+simval2(3) simval1(4) simval1(5) simval1(6)+simval2(4),simval2(5)];
    lgxobs=xobs(1,1:7);
    knock1=X(r,[15,17,19,21,22,23],2);
    knock2=X(r,[8,10,12,9,7],2);
    LGknock=[knock1(1)+knock2(1),knock1(2)+knock2(2),knock1(3)+knock2(3),knock1(4),knock1(5),knock1(6)+knock2(4),knock2(5)]; %
    lgknockobs=xobs(3,1:7);
    LGT=[LGknock(:),lgknockobs(:)];
    lgbar=[LGsimval(:),lgxobs(:)];
    lger=min(xobs(:))*er(1,1:7);
    lgerbar=[NaN*ones(length(lger),1),lger(:)];
    bb=bar(lgbar);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('3M Glucose - siControl Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 8],'XTick',1:7, 'box', 'off');
    set(gca,'XTickLabel',{'Pyruvate','Citrate','AKG','Succinate','Fumarate','Malate','Lactate'},'FontSize',10)
    ylabel('Metabolite Concentration in M')
    legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    numgroups = size(lgbar, 1);
    numbars = size(lgbar, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, lgbar(:,i), lgerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h1,'-dsvg','AnalysisFigure');
    print(h1,'-dpdf','afig.pdf');
    %% Low Glucose Knock Down
    h2=figure(2);
    bb=bar(LGT);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('3mM Glucose siICDc Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 8],'XTick',[1:7]);
    set(gca,'XTickLabel',{'Pyruvate','Citrate','AKG','Succinate','Fumarate','Malate','Lactate'},'FontSize',12)
    ylabel('Metabolite Concentration in M')
    legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    temp=er(2,1:7);
    temp2=min(xobs(:)).*[zeros(length(lger),1),temp(:)];
    temp3=.1*temp2(1,:);
    Klgerbar=[temp3;temp2(2:end,:)];
    numgroups = size(LGT, 1);
    numbars = size(LGT, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, LGT(:,i), Klgerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h2,'-dsvg','AnalysisFigure');
    print(h2,'-dtiff','LG.tiff');
    
    nadlgbar=[X(r,14,1),xobs(1,8);0,0];
    nadlgerbar=min(xobs(:)).*[0,er(1,8);0,0];
    h3=figure(3);
    bb=bar(nadlgbar);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    set(gca,'FontSize',12);
    title('3mM Glucose - siControl Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 2],'XTick',1);
    set(gca,'XTickLabel',{'NADPH'},'FontSize',12)
    ylabel('Metabolite Concentration in M')
     legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    numgroups = size(nadlgbar, 1);
    numbars = size(nadlgbar, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, nadlgbar(:,i), nadlgerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h3,'-dsvg','AnalysisFigure');
    print(h3,'-dtiff','KLG.tiff');
    knadlgbar=[X(r,14,2),xobs(3,8);0,0];
    knadlgerbar=min(xobs(:)).*[0,er(3,8);0,0];
    h4=figure(4);
    bb=bar(knadlgbar);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('3mM Glucose - siICDc Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 2],'XTick',1);
    set(gca,'FontSize',12);
    set(gca,'XTickLabel',{'NADPH'})
    ylabel('Metabolite Concentration in M')
     legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    numgroups = size(knadlgbar, 1);
    numbars = size(knadlgbar, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, knadlgbar(:,i), knadlgerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h4,'-dsvg','AnalysisFigure');
        
    %************************************************************
    %*******************High Glucose*****************************
    %************************************************************
    h=figure();
    [amt,r]=min(abs(points-12e-3));  % takes point closest to 12 mM
    xr=X(r,:,1);
    simval1 =X(r,[15,17,19,21,22,23],1);
    simval2=X(r,[8,10,12,9,7],1);
    HGsimval=[simval1(1)+simval2(1) simval1(2)+simval2(2) simval1(3)+simval2(3) simval1(4) simval1(5) simval1(6)+simval2(4),simval2(5)]; %
    hgxobs=xobs(2,1:7);
    
    %**********************HG Knock Down***********
    knock1=X(r,[15,17,19,21,22,23],2);
    knock2=X(r,[8,10,12,9,7],2);
    HGknock=[knock1(1)+knock2(1) knock1(2)+knock2(2) knock1(3)+knock2(3) knock1(4) knock1(5) knock1(6)+knock2(4),knock2(5)]; %
    hgknockobs=xobs(4,1:7);
    HGbar=[HGsimval(:),hgxobs(:)];
    % Ploting Each Value
    bb=bar(HGbar);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('12mM Glucose - siControl Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 8],'XTick',[1:7]);
    set(gca,'XTickLabel',{'Pyruvate','Citrate','AKG','Succinate','Fumarate','Malate','Lactate'},'FontSize',10)
    ylabel('Metabolite Concentration in M')
     legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    temp=er(3,1:7);
    HGerbar=min(xobs(:)).*[zeros(length(lger),1),temp(:)];
    numgroups = size(HGbar, 1);
    numbars = size(HGbar, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, HGbar(:,i), HGerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h,'-dsvg','AnalysisFigure');
    h=figure();
    HGT=[HGknock(:),hgknockobs(:)];
    bb=bar(HGT);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('12mM Glucose siICDc Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 8],'XTick',[1:7]);
    set(gca,'XTickLabel',{'Pyruvate','Citrate','AKG','Succinate','Fumarate','Malate','Lactate'},'FontSize',10)
    ylabel('Metabolite Concentration in M')
     legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    temp=er(4,1:7);
    KHGerbar=min(xobs(:)).*[zeros(length(lger),1),temp(:)];
    numgroups = size(HGT, 1);
    numbars = size(HGT, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, HGT(:,i), KHGerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h,'-dsvg','AnalysisFigure');
    print(h,'-dpdf','afig.pdf');
    nadhgbar=[X(r,14,1),xobs(2,8);0,0];
    nadhgerbar=min(xobs(:)).*[0,er(2,8);0,0];
    h3=figure(3);
    bb=bar(nadhgbar);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('12mM Glucose - siControl Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 2],'XTick',1);
    set(gca,'XTickLabel',{'NADPH'},'FontSize',10)
    ylabel('Metabolite Concentration in M')
     legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    numgroups = size(nadhgbar, 1);
    numbars = size(nadhgbar, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, nadhgbar(:,i), nadhgerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h,'-dsvg','AnalysisFigure');
        
    %%High Glucose Knock Down SUBplot
    knadhgbar=[X(r,14,2),xobs(4,8);0,0];
    knadhgerbar=min(xobs(:)).*[0,er(4,8);0,0];
    
    h=figure();
    bb=bar(knadhgbar);
    set(bb,'BarWidth',1); % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    title('12mM Glucose - siICDc Experimental and Simulation Values');
    %xlabel('Metabolite Name')
    set(gca,'XLim',[0 2],'XTick',1);
    set(gca,'XTickLabel',{'NADPH'},'FontSize',10)
    ylabel('Metabolite Concentration in M')
     legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1)
    hold on
    numgroups = size(knadhgbar, 1);
    numbars = size(knadhgbar, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
        errorbar(x, knadhgbar(:,i), knadhgerbar(:,i), 'k', 'linestyle', 'none');
    end
    print(h,'-dsvg','AnalysisFigure');
    insetplotting
end


if true
    disp('***********CIC Analysis*****************')
    ind1=find(points==2.8e-3);
    ind2=find(points==16.7e-3);
    sciclg1=X(ind1,[15,17,19,21,22,23],1);
    sciclg2=X(ind1,[8,10,12,9,7,14],1);
    scichg1=X(ind2,[15,17,19,21,22,23],1);
    scichg2=X(ind2,[8,10,12,9,7,14],1);
    kciclg1=X(ind1,[15,17,19,21,22,23],4);
    kciclg2=X(ind1,[8,10,12,9,7,14],4);
    kcichg1=X(ind2,[15,17,19,21,22,23],4);
    kcichg2=X(ind2,[8,10,12,9,7,14],4);
    simtot=sciclg1(2)+sciclg2(2);
    ktot=kciclg1(2)+kciclg2(2);
    cic.tod=(simtot-ktot)/simtot;
    cic.mitd=(sciclg1(2)-kciclg1(2))/sciclg1(2);
    cic.cyto=(sciclg2(2)-kciclg2(2))/sciclg2(2);
    disp('The total Citrate level decrease by 37% from wild type to knock down')
    disp(cic.tod)
    disp('The cytoslolic Citrate level decrease by 54% from wild type to knockdown')
    disp(cic.cyto)
    disp('The mitochondrial citrate should not have any change ')
    disp(cic.mitd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tv=[cwtot(1),cwcyt(1),cwmit(1)];
    tv1=[cwtot(2),cwcyt(2),cwmit(2)];
    tv2=[cktot(1),ckcyt(1),ckmit(1)];
    tv3=[cktot(2),ckcyt(2),ckmit(2)];
    swlg=[sciclg1(2)+sciclg2(2),sciclg2(2),sciclg1(2)];
    swhg=[scichg1(2)+scichg2(2),scichg2(2),scichg1(2)];
    kwlg=[kciclg1(2)+kciclg2(2),kciclg2(2),kciclg1(2)];
    kwhg=[kcichg1(2)+kcichg2(2),kcichg2(2),kcichg1(2)];
    
    
    disp('OUTPUT Wild type')
    met={'CIT_c+CIT_m','CIT_c+Null','CIT_m+Null'};
    fprintf(1,'%s\t%s\t\t%s\t\t%s\t\t%s\n','Metabolites','ObsLg','SimLg','ObsHg','SimHg');
    for j=1:3
        tm=char(met(j));
        fprintf(1,'%s\t%6.2e\t%6.2e\t%6.2e\t%6.2e\n',tm,tv(j),swlg(j),tv1(j),swhg(j));
    end
    
    disp('OUTPUT Knock down type')
    met={'CIT_c+CIT_m','CIT_c+Null','CIT_m+Null'};
    fprintf(1,'%s\t%s\t\t%s\t\t%s\t\t%s\n','Metabolites','ObsLg','SimLg','ObsHg','SimHg');
    for j=1:3
        tm=char(met(j));
        fprintf(1,'%s\t%6.2e\t%6.2e\t%6.2e\t%6.2e\n',tm,tv2(j),kwlg(j),tv3(j),kwhg(j));
    end
    
end
if true
    disp('***********Citrate lyase Analysis*****************')
    sciclg1=X(1,18,1);
    sciclg2=X(1,8,1);
    scichg1=X(2,18,1);
    scichg2=X(2,8,1);
    kciclg1=X(1,18,5);
    kciclg2=X(1,8,5);
    kcichg1=X(2,18,5);
    kcichg2=X(2,8,5);
    simtot=sciclg1+sciclg2;
    ktot=kciclg1+kciclg2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('OUTPUT Wild type')
    disp(simtot)
    disp(ktot)
    disp(X(1,:,5))
end

