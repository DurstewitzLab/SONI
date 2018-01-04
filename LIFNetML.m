function [tall,Vall,T,R,Gav,tCaiAll,CaiAll]=LIFNetML(par,tiv,Iaff,TList,gsynAvg,NoisePar,s,Flags)
% Copyright (C) 2002,2003 Daniel Durstewitz

% inputs: par - vector of neuron parameters
%         tiv - vector with start & stop time of sim.
%         Iaff - current stimuli applied
%         TList - times & durations of current stimuli
%         gsynAvg - constant syn. conductances
%         NoisePar - vector specifying params. for noisy background syn. input
%         s - random seed
%         Flags - flags determining various sim. options
% outputs: tall - times of Vm, m-gate etc. recordings
%          Vall - matrix containing various sim. variables (specified on line 172 below)
%          T - vector with all spike times
%          R - instantaneous firing rates computed from spike times
%          Gav - gADP conductance averaged over interspike-intervals
%          tCaiAll, CaiAll - similar to tall/Vall for [Ca2+]i

rand('state',s);
randn('state',s);

% transfer parameters to meaningful par. names
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
tcADP=par(20);
slope_ADP=par(21);
half_ADP=par(22);
utilE=par(23);
trecE=par(24);
tfacE=par(25);
utilI=par(29);
trecI=par(30);
tfacI=par(31);

% transfer parameters for noisy syn. inp.
Tnoise=NoisePar(1);
RnoiseE=NoisePar(2);
RnoiseI=NoisePar(3);
gAMPAnoise=NoisePar(4);
gNMDAnoise=NoisePar(5);
gGABAnoise=NoisePar(6);
% compute Poisson-distributions in advance
err=1e-4;
alpha=RnoiseE*Tnoise/1000;
pcumE(1)=exp(-alpha);
ps=pcumE(1);
i=0;
while (1-max(pcumE))>err
    i=i+1;
    ps=ps*alpha/i;
    pcumE(i+1)=pcumE(i)+ps;
end;
pcumE(i+2)=1;
alpha=RnoiseI*Tnoise/1000;
pcumI(1)=exp(-alpha);
ps=pcumI(1);
i=0;
while (1-max(pcumI))>err
    i=i+1;
    ps=ps*alpha/i;
    pcumI(i+1)=pcumI(i)+ps;
end;
pcumI(i+2)=1;

% initialize variables
pdep=[utilE trecE tfacE utilE trecE tfacE utilI trecI tfacI];
Rprev=ones(1,3);
uprev=zeros(1,3);
m_lim=1/(1+exp(slope_ADP*half_ADP));
V0=[EL m_lim 0 0 0];
gsyn0=zeros(1,8);
tall=[];
Vall=[];
Gall=[];
CaiAll=[];
tCaiAll=[];
t0=tiv(1);
t1=tiv(2);
t(1)=t0;
Nsp=0;
ST=-100000.0;
Tnext=Tnoise;
TL=[TList(1,:) (TList(1,:)+TList(2,:)) t1];
Iext=0;

% Flag decides on whether to use implicit or explicit method for ODE-solver
if Flags(2)==1
    opt=odeset('Events',@evLIF,'Jacobian',@JacLIF,'RelTol',1e-5,'AbsTol', ...
                (1e-6)*[1 1e-2 1e-3 1e-3 1e-10]);
else
    opt=odeset('Events',@evLIF,'RelTol',1e-5,'AbsTol', ...
                (1e-6)*[1 1e-2 1e-3 1e-3 1e-10]);
end;

% MAIN LOOP ########################
kk=1;
i=0;
while t(length(t))<t1
    
    % call solver and stop at each time of synaptic or external event, or spike
    k=find(TL>t0);
    Tstop=min(TL(k));
    Tstop=min(Tstop,Tnext);
    if (Nsp>0) ST=Tsp(Nsp); end;
    if Flags(2)==1
        [t,V,te,Ve,ie]=ode23s(@LIF,[t0 Tstop],V0,opt,t0,ST,gsyn0,Iext,par,gsynAvg,Flags);
    else
        [t,V,te,Ve,ie]=ode23(@LIF,[t0 Tstop],V0,opt,t0,ST,gsyn0,Iext,par,gsynAvg,Flags);
    end;

    % apply stimulus currents
    [F,idx]=ismember(t(length(t)),TL);
    if F>0
        k=length(TList(1,:));
        if idx<=2*k
            if idx<=k
                Iext=Iaff(idx);
            else
                Iext=0;
            end;
        end;
    end;

    % just for recording [Ca]i
    tCa=(t0:0.2:t(length(t)))';
    Cai=gsyn(gsyn0(7:8),tCa-t0,[tc1Ca tc2Ca]);
    CaiAll=[CaiAll;Cai];
    tCaiAll=[tCaiAll;tCa];
    
    % just to monitor simulation time
    if (t(length(t))>=kk*5000) & (t0<kk*5000)
        t(length(t))
        kk=kk+1;
    end;

    % decay of synaptic currents & [Ca]i over integration interval dt
    dt=t(length(t))-t0;
    gsyn0=gsyn0.*exp(-dt./par([8 7 11 10 14 13 18 17])');

    % in case of spike: apply recurrent synaptic inputs etc.
    if ~isempty(ie) 
        Nsp=Nsp+1;
        Tsp(Nsp)=t(length(t));  % save spike time
        V(length(t),1)=EL;  % reset Vm
        [u,R]=syndep(uprev,Rprev,t(length(t))-ST,pdep); % compute synaptic depression/facilitation 
        uprev=u;
        Rprev=R;
        X=[u.*R 1].*[gAMPA gNMDA gGABA ACa];
        gsyn0([1 3 5 7])=gsyn0([1 3 5 7])+X;    % update syn. inputs
        gsyn0([2 4 6 8])=gsyn0([2 4 6 8])+X;
    end;

    % just for recording sim. variables
    t0=t(length(t));
    V0=V(length(t),:)';
    tall=[tall;t];
    Vall=[Vall;V(:,[2 4 5])];
    
    % apply Poisson background syn. input
    if abs(t0-Tnext)<1e-8
        Tnext=Tnext+Tnoise;
        q=rand;
        Ne=min(find(pcumE>=q))-1;
        q=rand;
        Ni=min(find(pcumI>=q))-1;
        X=[gAMPAnoise gNMDAnoise gGABAnoise].*[Ne Ne Ni];
        gsyn0([1 3 5])=gsyn0([1 3 5])+X;
        gsyn0([2 4 6])=gsyn0([2 4 6])+X;
    end;
    
end;
% END OF MAIN LOOP #################

% compute instantaneous firing rates & average gADP
if Nsp>0
    R=1000./diff(Tsp);
    T=Tsp(1:length(Tsp)-1);
    for i=1:length(R)
        k=find(tall>=Tsp(i) & tall<=Tsp(i+1));
        tx=tall(k);
        gx=Vall(k,1);
        dt=((tx(2:length(tx)-1)-tx(1:length(tx)-2))+(tx(3:length(tx))-tx(2:length(tx)-1)))/2;
        Gav(i)=gADP*dt'*gx(2:length(gx)-1)/sum(dt);
    end;
else R=[]; T=[]; Gav=[];
end;


% system of differential equations
function dV=LIF(t,V,t0,spt_last,gsyn0,Iext,par,gsynAvg,Flags)
Cm=par(1);
gL=par(2);
EL=par(3);
gAHP=par(4);
tcAHP=par(5);
gAMPA=par(6);
tc1AMPA=par(7);
tc2AMPA=par(8);
gfAMPA=gsyn0(1:2);
gNMDA=par(9);
tc1NMDA=par(10);
tc2NMDA=par(11);
gfNMDA=gsyn0(3:4);
gGABA=par(12);
tc1GABA=par(13);
tc2GABA=par(14);
gfGABA=gsyn0(5:6);
Vth=par(15);
ACa=par(16);
tc1Ca=par(17);
tc2Ca=par(18);
gfCa=gsyn0(7:8);
gADP=par(19);
tcADP=par(20);
slope_ADP=par(21);
half_ADP=par(22);

% LIF neuron (eqn. 1-5)
Iahp=gAHP*exp(-(t-spt_last)/tcAHP)*(EL-V(1));
Iampa=(gsyn(gfAMPA,t-t0,[tc1AMPA tc2AMPA])+gsynAvg(1))*(0-V(1));
if Flags(1)==1
    a1=par(32:34);
    Inmda=a1(1)/(1+a1(2)*exp(-a1(3)*V(1)))*(gsyn(gfNMDA,t-t0,[tc1NMDA tc2NMDA])+gsynAvg(2))*(0-V(1));
else
    Inmda=(gsyn(gfNMDA,t-t0,[tc1NMDA tc2NMDA])+gsynAvg(2))*(24+0.29*V(1));
end;
Igaba=(gsyn(gfGABA,t-t0,[tc1GABA tc2GABA])+gsynAvg(3))*((EL-5)-V(1));
Iadp=gADP*V(2)*(35-V(1));
dV(1)=(gL*(EL-V(1))+Iahp+Iampa+Inmda+Igaba+Iadp+Iext)/Cm;

Cai=gsyn(gfCa,t-t0,[tc1Ca tc2Ca]);
m_lim=1/(1+exp(slope_ADP*(half_ADP-Cai)));  % slope>0
dV(2)=(m_lim-V(2))/tcADP;

% Ca-var.-signals (eqn. 7-8)
dV(3)=(Cai-V(3))/10000;
dV(4)=((Cai-V(3))^2-V(4))/10;
dV(5)=(dV(3)^2-V(5))/10;

dV=dV';

% compute Jacobian matrix in case implicit solver is used
function J=JacLIF(t,V,t0,spt_last,gsyn0,Iext,par,gsynAvg,Flags)
Cm=par(1);
gL=par(2);
EL=par(3);
gAHP=par(4);
tcAHP=par(5);
gAMPA=par(6);
tc1AMPA=par(7);
tc2AMPA=par(8);
gfAMPA=gsyn0(1:2);
gNMDA=par(9);
tc1NMDA=par(10);
tc2NMDA=par(11);
gfNMDA=gsyn0(3:4);
gGABA=par(12);
tc1GABA=par(13);
tc2GABA=par(14);
gfGABA=gsyn0(5:6);
Vth=par(15);
ACa=par(16);
tc1Ca=par(17);
tc2Ca=par(18);
gfCa=gsyn0(7:8);
gADP=par(19);
tcADP=par(20);
slope_ADP=par(21);
half_ADP=par(22);

J=zeros(length(V),length(V));

dIahpdV1=-gAHP*exp(-(t-spt_last)/tcAHP);
dIampadV1=-(gsyn(gfAMPA,t-t0,[tc1AMPA tc2AMPA])+gsynAvg(1));
if Flags(1)==1
    a1=par(32:34);
    gs=gsyn(gfNMDA,t-t0,[tc1NMDA tc2NMDA])+gsynAvg(2);
    dInmdadV1=-a1(1)/(1+a1(2)*exp(-a1(3)*V(1)))^2*gs*V(1)*a1(2)*a1(3)*exp(-a1(3)*V(1))-a1(1)/(1+a1(2)*exp(-a1(3)*V(1)))*gs;
else
    dInmdadV1=(gsyn(gfNMDA,t-t0,[tc1NMDA tc2NMDA])+gsynAvg(2))*0.29;
end;
dIgabadV1=-(gsyn(gfGABA,t-t0,[tc1GABA tc2GABA])+gsynAvg(3));
dIadpdV1=-gADP*V(2);
J(1,1)=(-gL+dIahpdV1+dIampadV1+dInmdadV1+dIgabadV1+dIadpdV1)/Cm;
dIadpdV2=gADP*(35-V(1));
J(1,2)=dIadpdV2/Cm;
J(2,2)=-1/tcADP;
Cai=gsyn(gfCa,t-t0,[tc1Ca tc2Ca]);
J(3,3)=-1/10000;
J(4,3)=-2/10*(Cai-V(3));
J(4,4)=-1/10;
J(5,3)=-2/10*(Cai-V(3))/(10000^2);
J(5,5)=-1/10;

% compute synaptic conductances (eqn. 4-5)
function g=gsyn(gf,dt,tc)
c=tc(1)*tc(2)/(tc(2)-tc(1));
g=c*(gf(1)*exp(-dt/tc(2))-gf(2)*exp(-dt/tc(1)));

% time-derivatives of syn. conduc. for Jacobian
function dg=dgsyn(gf,dt,tc)
c=tc(1)*tc(2)/(tc(2)-tc(1));
dg=c*(gf(1)*(-1/tc(2))*exp(-dt/tc(2))-gf(2)*(-1/tc(1))*exp(-dt/tc(1)));

% event detection function for spike-time detection
function [sp,term,direc]=evLIF(t,V,t0,spt_last,gsyn0,Iext,par,gsynAvg,Flags)
Vth=par(15);
sp=ones(length(V),1);
sp(1)=V(1)-Vth;
term=ones(length(V),1);
direc=ones(length(V),1);

% synaptic depression (eq. 6)
function [u,R]=syndep(uprev,Rprev,dt,pdep)
util=pdep([1 4 7]);
trec=pdep([2 5 8]);
tfac=pdep([3 6 9]);
qu=uprev.*exp(-(dt*ones(1,3))./tfac);
qR=exp(-(dt*ones(1,3))./trec);
u=qu+util.*(1-qu);
R=Rprev.*(1-u).*qR+1-qR;
