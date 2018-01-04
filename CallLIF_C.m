function [tall,Vall,T,R,Gav,Cai]=CallLIF_C(p0,tiv,Iaff,TList,mgsyn,NoisePar,s,Flags);
% Copyright (C) 2002,2003 Daniel Durstewitz

% inputs & outputs - same as in LIFNetML.m

% Control parameters
N=length(p0(1,:));  %Number of neurons
CtrPar(1)=N;  %N
CtrPar(2)=tiv(1);   %Tstart
CtrPar(3)=tiv(2);   %Tstop
CtrPar(4)=1e-5; %hmin
CtrPar(5)=0.1;  %hmax
CtrPar(6)=1e8;  %WgtSaveStep
CtrPar(7)=0.1;  %ViewSaveStep
CtrPar(8)=2*(tiv(2)-tiv(1))/CtrPar(7);  %kmax
CtrPar(9)=0.001;    %VthTol
CtrPar(10)=1e-6;    %eps
CtrPar(11)=NoisePar(1); %Tnoise
CtrPar(12)=s;   %rnd-seed
CtrPar(13)=0;   %Warnings (0=suppress)

% Neuron parameters
NeuPar(1:3,:)=p0(1:3,:);    %Cm,gL,EL
NeuPar(4,:)=p0(15,:);   %Vth
NeuPar(5:6,:)=p0(4:5,:);    %gAHP,tcAHP
NeuPar(7,:)=p0(3,:);    %EAHP
NeuPar(8:12,:)=p0(16:20,:); %ACa,tcCa_on,tcCa_off,gADP,tcADP
NeuPar(13,:)=35;    %EADP
NeuPar(14:15,:)=p0(21:22,:);    %slopeADP,halfADP

% List of neuron-types
NPList=(1:N)';

% Synapse-type parameters for recurrent connectivity
STypPar(1,:)=[p0(6,1) p0(9,1) p0(12,1)];    %gAMPA,gNMDA,gGABA
STypPar(2,:)=[p0(7,1) p0(10,1) p0(13,1)];    %tc1AMPA,tc1NMDA,tc1GABA
STypPar(3,:)=[p0(8,1) p0(11,1) p0(14,1)];    %tc2AMPA,tc2NMDA,tc2GABA
STypPar(4,:)=[0 0 p0(3,1)-5];    %Erev
STypPar(5:7,:)=zeros(3,3);
if (length(Flags)==3 & Flags(2)==1) | (length(Flags)==2 & Flags(1)==1)
    STypPar(5:7,2)=p0(32:34); %Mg_gate, Mg_fac, Mg_slope
else
    STypPar(5,:)=[0 -0.29 0];    %Mg_gate
    STypPar(4,2)=24/(-0.29);
end;
STypPar(8:11,:)=0;    %Mg_half,gLF,A_LTP,A_LTD
STypPar(12:13,:)=20;    %tc_LTP,tc_LTD
STypPar(14,:)=0;    %wMIN
STypPar(15,:)=1e6;    %wMAX

% Synapse parameters and connectivity matrix
for i=1:3
    k=(1+(i-1)*N):(i*N);
    SynPar(1,k)=i;  %SynTypeNo
    if i<3
        SynPar(2:4,k)=p0(23:25,:);   %utilE,trecE,tfacE
    else
        SynPar(2:4,k)=p0(29:31,:);   %utilI,trecI,tfacI
    end;
    if STypPar(1,i)>0
        SynPar(5,k)=p0(3+i*3,:)./STypPar(1,i);  % adjusting wgts to gmax
    end;
    SPMtx(:,:,i)=diag(k,0);
end;

% only for configurations with input- and/or read-out neurons:
if length(Flags)==3
    gCR=p0(32,:);
    tc1CR=p0(33,:);
    tc2CR=p0(34,:);
    utilCR=p0(35,:);
    trecCR=p0(36,:);
    tfacCR=p0(37,:);
    switch Flags(1)
        case 1
			% configuration for CS + RS
			%
            k=3*N+(1:2);
            SynPar(1,k)=[2 3];
            SynPar(2:4,k)=[utilCR;trecCR;tfacCR];
            SynPar(5,k)=gCR./STypPar(1,2:3);
            SPMtx(2,1,1)=k(1);
            SPMtx(1,2,1)=k(2);
        case 2
			% configuration for CS + InpN
			%
            k=3*N+(1:2);
            SynPar(1,k)=[1 2];
            SynPar(2:4,k)=[utilCR(2) trecCR(2) tfacCR(2)]'*ones(1,length(k));
            SynPar(5,k)=gCR(2);
            SPMtx(1,2,1:2)=k;
        otherwise
			% configuration for CS + RS + InpN
			%
            k=3*N+(1:4);
            SynPar(1,k)=[2 3 1 2];
            utilCR(4)=utilCR(3);
            trecCR(4)=trecCR(3);
            tfacCR(4)=tfacCR(3);
            SynPar(2:4,k)=[utilCR;trecCR;tfacCR];
            SynPar(5,k)=[gCR(1)/STypPar(1,2) gCR(2)/STypPar(1,3) gCR(3) gCR(3)];
            SPMtx(2,1,1)=k(1);
            SPMtx(1,2,1)=k(2);
            SPMtx(1,3,1:2)=k(3:4);
	end;
end;

% Synapse-type parameters for constant conductance inputs
m=length(SynPar(1,:));
for i=1:3
    k=(4+(i-1)*N):(3+i*N);
    STypPar(1,k)=mgsyn(i,:);    %gAMPA,gNMDA,gGABA
    STypPar(2:15,k)=STypPar(2:15,i)*ones(1,length(k));
    STypPar(2:3,k)=0;   %tcs
    mk=(m+1+(i-1)*N):(m+i*N);
    SynPar(1,mk)=k;
    BackgrdInp(:,i)=mk';
end;
k=find(mgsyn>0);
if (isempty(k)) BackgrdInp=[]; end;

% Synapse parameters for synaptic background noise (require only SynPar)
m=length(SynPar(1,:));
r=length(BackgrdInp);
for i=1:3
    k=(m+1+(i-1)*N):(m+i*N);
    SynPar(1,k)=i;  %SynTypeNo
    if i<3
        SynPar(2,k)=NoisePar(2,:);   %RateE
    else
        SynPar(2,k)=NoisePar(3,:);   %RateI
    end;
    SynPar(5,k)=NoisePar(3+i,:);  % wgts
    BackgrdInp(:,r+i)=k';
end;
k=find(NoisePar(2:3,:)>0);
if (isempty(k)) BackgrdInp=BackgrdInp(:,1:r); end;

% Event matrix & event times
EvtMtx=Iaff';
EvtTimes=TList;
ViewList=(1:N)';

% calling the C LIF-simulator
Gav=[];
[ST,tall,Vall]=LIFNetSim1(CtrPar,NeuPar,NPList,STypPar,SynPar,SPMtx,EvtMtx,EvtTimes,ViewList,BackgrdInp);
if nargout>4
    k=find(ST(1,:)>0);
    Tsp=ST(1,k);
    n=length(k);
    gADP=NeuPar(11,1);
    for i=2:n
        k=find(tall>=Tsp(i-1) & tall<=Tsp(i));
        tx=tall(k);
        gx=Vall(2,k);
        dt=((tx(2:end-1)-tx(1:end-2))+(tx(3:end)-tx(2:end-1)))/2;
        Gav(i-1)=gADP*gx(2:end-1)*dt'/sum(dt);
    end;
end;

R=(1000./diff(ST'))';
T=ST(:,1:end-1);
tall=tall';
Vall=Vall';
