function [nc_rE,nc_gADP]=NullClLIF(rE,par,mgsyn)
% Copyright (C) 2002,2003 Daniel Durstewitz

% computes the nullclines of the LIF-system with gADP
% rE = vector of firing rates
% par = LIF parameter configuration
% mgsyn = average background [gAMPA gNMDA gGABA]

% Scale factor to prevent too large exponentials in the integrals - adjust
% this if you get warnings:
%fac2=1e20;
fac2=1e50;

fac=1.0;
if par(4,1)>0.3
    fac=0.3/par(4,1);
    par([1 2 4 6 9 12 19],1)=fac*par([1 2 4 6 9 12 19],1);
    mgsyn=fac*mgsyn;
end;

lb0=0.0001;
ub0=0.008;
for i=1:length(rE)
    Tisi=1000.0./rE(i);
    ub=ub0;
    lb=lb0;
    if (nc_rEgADP(0,Tisi,par,mgsyn,fac2)<0) nc_rE(i)=0;
    else
        sdiff=sign(nc_rEgADP(lb0,Tisi,par,mgsyn,fac2))*sign(nc_rEgADP(ub0,Tisi,par,mgsyn,fac2));
        while sdiff>0 & ub<5*ub0
            lb=0;
            ub=ub+ub0;
            sdiff=sign(nc_rEgADP(lb,Tisi,par,mgsyn,fac2))*sign(nc_rEgADP(ub,Tisi,par,mgsyn,fac2));
        end;
        if sdiff<=0
            [nc_rE(i),fval,efl]=fzero(@nc_rEgADP,[lb ub],[],Tisi,par,mgsyn,fac2);
        else nc_rE(i)=NaN; end;
    end;
end;

ACa=par(16,1);
tc1Ca=par(17,1);
tc2Ca=par(18,1);
gADP=par(19,1);
tcfADP=par(20,1);
slopefADP=par(21,1);
halffADP=par(22,1);
for i=1:length(rE)
    Tisi=1000.0./rE(i);
    a=exp(-Tisi/tcfADP);
    b=quadl(@Gsglintgr,0,Tisi,[],[],Tisi,par);
    g0=a*b/(1-a);
    nc_gADP(i)=gADP*quadl(@Gdblintgr,0,Tisi,[],[],Tisi,par,g0)/Tisi;
    
    %avCa=ACa*tc1Ca*tc2Ca/Tisi;
    %nc_gADP2(i)=gADP/(1+exp(slopefADP*(halffADP-avCa)));
end;

nc_rE=nc_rE./fac;
nc_gADP=nc_gADP./fac;


function z=nc_rEgADP(mgADP,Tisi,par,mgsyn,fac2)
Cm=par(1);
gL=par(2);
EL=par(3);
gAHP=par(4);
tcAHP=par(5);
gAMPA=par(6);
tc1AMPA=par(7);
tc2AMPA=par(8);
gNMDA=par(9);
tc1NMDA=par(10);
tc2NMDA=par(11);
gGABA=par(12);
tc1GABA=par(13);
tc2GABA=par(14);
Vth=par(15);
ACa=par(16);
tc1Ca=par(17);
tc2Ca=par(18);
gADP=par(19);
tcfADP=par(20);
slopefADP=par(21);
halffADP=par(22);
Uampa=par(23);
trecAMPA=par(24);
tfacAMPA=par(25);
Unmda=par(26);
trecNMDA=par(27);
tfacNMDA=par(28);
Ugaba=par(29);
trecGABA=par(30);
tfacGABA=par(31);

mgAMPA=mgsyn(1);
mgNMDA=mgsyn(2);
mgGABA=mgsyn(3);

wAMPA=EPSPn(1e3/Tisi,Uampa,trecAMPA/1000,tfacAMPA/1000);
wNMDA=EPSPn(1e3/Tisi,Unmda,trecNMDA/1000,tfacNMDA/1000);
wGABA=EPSPn(1e3/Tisi,Ugaba,trecGABA/1000,tfacGABA/1000);

a0=gL*Tisi;
a1=gAMPA*tc1AMPA*tc2AMPA*wAMPA;
a2=gAHP*tcAHP*(1-exp(-Tisi/tcAHP));
a3=-0.29*gNMDA*tc1NMDA*tc2NMDA*wNMDA;
a4=mgADP*Tisi;
a5=gGABA*tc1GABA*tc2GABA*wGABA;
a6=mgGABA*Tisi;
a7=mgAMPA*Tisi;
a8=-0.29*mgNMDA*Tisi;
%A=exp(-(a0+a1+a2+a3+a4+a5+a6+a7+a8)/Cm);
A=exp(-(a0+a1+a2+a3+a4+a5+a6+a7+a8))^(1/Cm);
B=quadl(@Bintgr,0,Tisi,[],[],mgADP,Tisi,par,[wAMPA wNMDA wGABA],mgsyn,fac2);
z=Vth-A*(EL+fac2*B);

function B=Bintgr(t,mgADP,Tisi,par,wALL,mgsyn,fac2)
Cm=par(1);
gL=par(2);
EL=par(3);
gAHP=par(4);
tcAHP=par(5);
gAMPA=par(6);
tc1AMPA=par(7);
tc2AMPA=par(8);
gNMDA=par(9);
tc1NMDA=par(10);
tc2NMDA=par(11);
gGABA=par(12);
tc1GABA=par(13);
tc2GABA=par(14);
Vth=par(15);
ECl=EL-5;

mgAMPA=mgsyn(1);
mgNMDA=mgsyn(2);
mgGABA=mgsyn(3);

b0=gL*EL;
b1=gAHP*exp(-t/tcAHP)*EL;
bpN=gNMDA*wALL(2)*tc1NMDA*tc2NMDA/(tc2NMDA-tc1NMDA);
b2=24*bpN*(exp(-t/tc2NMDA)/(1-exp(-Tisi/tc2NMDA))-exp(-t/tc1NMDA)/(1-exp(-Tisi/tc1NMDA)));
b3=35*mgADP;
bpG=gGABA*wALL(3)*tc1GABA*tc2GABA/(tc2GABA-tc1GABA);
b4=bpG*(exp(-t/tc2GABA)/(1-exp(-Tisi/tc2GABA))-exp(-t/tc1GABA)/(1-exp(-Tisi/tc1GABA)))*ECl;
b5=ECl*mgGABA;
b6=24*mgNMDA;
b=(b0+b1+b2+b3+b4+b5+b6)/Cm;

a0=gL*t;
apA=gAMPA*wALL(1)*tc1AMPA*tc2AMPA/(tc2AMPA-tc1AMPA);
a1=apA*(tc2AMPA*(1-exp(-t/tc2AMPA))/(1-exp(-Tisi/tc2AMPA))-tc1AMPA*(1-exp(-t/tc1AMPA))/(1-exp(-Tisi/tc1AMPA)));
a2=gAHP*tcAHP*(1-exp(-t/tcAHP));
apN=-0.29*gNMDA*wALL(2)*tc1NMDA*tc2NMDA/(tc2NMDA-tc1NMDA);
a3=apN*(tc2NMDA*(1-exp(-t/tc2NMDA))/(1-exp(-Tisi/tc2NMDA))-tc1NMDA*(1-exp(-t/tc1NMDA))/(1-exp(-Tisi/tc1NMDA)));
a4=mgADP*t;
apG=gGABA*wALL(3)*tc1GABA*tc2GABA/(tc2GABA-tc1GABA);
a5=apG*(tc2GABA*(1-exp(-t/tc2GABA))/(1-exp(-Tisi/tc2GABA))-tc1GABA*(1-exp(-t/tc1GABA))/(1-exp(-Tisi/tc1GABA)));
a6=mgGABA*t;
a7=mgAMPA*t;
a8=-0.29*mgNMDA*t;

%B=b.*exp((a0+a1+a2+a3+a4+a5+a6+a7+a8)/Cm);
B=(b.*(exp(a0+a1+a2+a3+a4+a5+a6+a7+a8).^(1/Cm)))/fac2;

function Gs=Gsglintgr(t,Tisi,par)
% could theoretically be reduced to a single integral ...
ACa=par(16);
tc1Ca=par(17);
tc2Ca=par(18);
gADP=par(19);
tcfADP=par(20);
slopefADP=par(21);
halffADP=par(22);

a=exp(t/tcfADP);
bp=ACa*tc1Ca*tc2Ca/(tc2Ca-tc1Ca);
b0=bp*(exp(-t/tc2Ca)/(1-exp(-Tisi/tc2Ca))-exp(-t/tc1Ca)/(1-exp(-Tisi/tc1Ca)));
Gs=a./((1+exp(slopefADP*(halffADP-b0)))*tcfADP);

function Gd=Gdblintgr(t,Tisi,par,g0)
tcfADP=par(20);

A=exp(-t/tcfADP);
for i=1:length(t)
    if (t(i)==0) B(i)=0;
    else B(i)=quadl(@Gsglintgr,0,t(i),[],[],Tisi,par);
    end;
end;
Gd=A.*(g0+B);

function a=EPSPn(r,U,trec,tfac)
q=exp(-1./(r*trec));
ust=U./(1-(1-U)*exp(-1./(r*tfac)));
Rst=(1-q)./(1-(1-ust).*q);
a=Rst.*ust;
