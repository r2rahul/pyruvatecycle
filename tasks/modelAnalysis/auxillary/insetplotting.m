%cmap=[0.6,0.55686,0.76471;0.9451,0.63922,0.25098];
%colormap(cmap)
%%
if true
    %pyr,cit,akg,suc,fum,mal,lac
    %h1=figure(1); clf;
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
    temp=er(2,1:7);
    temp2=min(xobs(:)).*[zeros(length(lger),1),temp(:)];
    temp3=.1*temp2(1,:);
    Klgerbar=[temp3;temp2(2:end,:)];
    Klgerbar(:,1)=NaN*ones(length(lger),1);
    nadlgbar=[X(r,14,1),xobs(1,8);0,0];
    nadlgerbar=min(xobs(:)).*[NaN,er(1,8);0,0];
    knadlgbar=[X(r,14,2),xobs(3,8);0,0];
    knadlgerbar=min(xobs(:)).*[NaN,er(3,8);0,0];
    %h1 = subplot(2,2,1);
    inset_plot(lgbar,lgerbar,nadlgbar,nadlgerbar,'3mMControl')
    %h1 = subplot(2,2,2);
    inset_plot(LGT,Klgerbar,knadlgbar,knadlgerbar,'3mMKnock')
    %%
    %************************************************************
    %*******************High Glucose*****************************
    %************************************************************
    %h=figure(); clf;
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
    temp=er(3,1:7);
    HGerbar=min(xobs(:)).*[NaN*ones(length(lger),1),temp(:)];
    nadhgbar=[X(r,14,1),xobs(2,8);0,0];
    nadhgerbar=min(xobs(:)).*[0,er(2,8);0,0];
    temp=er(4,1:7);
    KHGerbar=min(xobs(:)).*[NaN*ones(length(lger),1),temp(:)];
    knadhgbar=[X(r,14,2),xobs(4,8);0,0];
    knadhgerbar=min(xobs(:)).*[NaN,er(4,8);0,0];
    % Ploting Each Value
    %h1 = subplot(2,2,3);
    inset_plot(HGbar,HGerbar,nadhgbar,nadhgerbar,'12mMControl')
    %h1 = subplot(2,2,4);
    inset_plot(HGT,KHGerbar,knadhgbar,knadhgerbar,'12mMKnock')
    
end