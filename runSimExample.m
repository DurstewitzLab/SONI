%
% run sim. with constant syn. conductance input (no noise)
%
clear all;

load LAConfigPar; 
p0=LAPar1;  % param. config.
tiv=[0 50000];  % start & stop time of sim. in msec
Iaff=1;    % stimulus current
TList=[5000 200]';  % time & duration of stim. curr.
GsynAvg=p0(36:38);  % average synaptic background conductances (AMPA,NMDA,GABAA)
NoisePar=[1e8 0 0 0 0 0]';  % noisy syn. inp. - no noise
s=2;    % random seed (irrelevant if no noise)
Flags=[2 1];    % only for ML sim.: Flags(2)= 1 (use ode23s with Jacobian), 2 (use ode23)
rE=5:2.5:100;    % vector of firing rates for nullclines

p0([6 9])=LAPar1([6 9])*1.04;   % adjust recurr. syn. to produce very slowly climbing activity
[tall1,Vall1,T1,R1,Gav1]=CallLIF_C(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);    % call C simulator
%[tall1,Vall1,T1,R1,Gav1]=LIFNetML(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);   % call ML simulator
[nc_rE1,nc_gADP1]=NullClLIF(rE,p0,p0(36:38));    % compute nullclines

p0([6 9])=LAPar1([6 9])*1.15;   % adjust recurr. syn. to produce faster climbing activity
[tall2,Vall2,T2,R2,Gav2]=CallLIF_C(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);    % call C simulator
%[tall2,Vall2,T2,R2,Gav2]=LIFNetML(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);   % call ML simulator
[nc_rE2,nc_gADP2]=NullClLIF(rE,p0,p0(36:38));    % compute nullclines

save SimResExp1 tall1 Vall1 T1 R1 Gav1 tall2 Vall2 T2 R2 Gav2 ...
    rE nc_rE1 nc_rE2 nc_gADP1 nc_gADP2;

% plot results:
figure(1);
subplot(2,2,1), hold off;
plot(rE,nc_rE1,'g',rE,nc_gADP1,'r','LineWidth',1.2);
hold on, plot(R1,Gav1,'k','LineWidth',0.8);
set(gca,'FontSize',11);
legend('FR-nullcline','<g_A_D_P>-nullcline','trajectory',2);
legend('boxoff');
set(gca,'FontSize',13)
axis([0 100 0 0.12]);
xlabel('Firing rate (Hz)');
ylabel('<g_A_D_P> (mS/cm^2)');
subplot(2,2,2), hold off, plot(T1./1000,R1,'k','LineWidth',1.2);
set(gca,'FontSize',13);
axis([0 45 0 100]); box off;
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
subplot(2,2,3), hold off, plot(rE,nc_rE2,'g',rE,nc_gADP2,'r','LineWidth',1.2);
hold on, plot(R2,Gav2,'k','LineWidth',0.8);
set(gca,'FontSize',11);
legend('FR-nullcline','<g_A_D_P>-nullcline','trajectory',2);
legend('boxoff');
set(gca,'FontSize',13);
axis([0 100 0 0.12]);
xlabel('Firing rate (Hz)');
ylabel('<g_A_D_P> (mS/cm^2)');
subplot(2,2,4), hold off, plot(T2./1000,R2,'k','LineWidth',1.2);
set(gca,'FontSize',13);
axis([0 45 0 100]); box off;
xlabel('Time (s)');
ylabel('Firing rate (Hz)');


%
% run same sim. with noisy syn. inputs
%
clear all;

load LAConfigPar; 
p0=LAPar1;  % param. config.
tiv=[0 50000];  % start & stop time of sim. in msec
Iaff=1;    % stimulus current
TList=[5000 200]';  % time & duration of stim. curr.
GsynAvg=[0 0 0]';  % no constant synaptic background conductances (AMPA,NMDA,GABAA)
% convert constant syn. conductances into noisy inputs with same mean:
Re=20000; Ri=20000; % rate of exc. & inh. inputs (Hz)
tc1AMPA=p0(7);
tc2AMPA=p0(8);
tc1NMDA=p0(10);
tc2NMDA=p0(11);
tc1GABA=p0(13);
tc2GABA=p0(14);
NoisePar=zeros(6,1);
NoisePar(1)=0;  % TnoiseStep (has to be >0 for ML simulator!)
% TnoiseStep=0 means continuous updating (works only with C simulator); 
% to speed up sim. set TnoiseStep>0 (e.g. 0.5 msec).
NoisePar(2:3)=[Re Ri]'; % rates of inputs
NoisePar(4)=p0(36)/(Re*tc1AMPA*tc2AMPA/1000);   % amplitude of backgr. gAMPA
NoisePar(5)=p0(37)/(Re*tc1NMDA*tc2NMDA/1000);   % amplitude of backgr. gNMDA
NoisePar(6)=p0(38)/(Ri*tc1GABA*tc2GABA/1000);   % amplitude of backgr. gGABAA
s=2;    % random seed
Flags=[2 1];    % only for ML sim.: Flags(2)= 1 (use ode23s with Jacobian), 2 (use ode23)

p0([6 9])=LAPar1([6 9])*1.04;   % adjust recurr. syn. to produce very slowly climbing activity
[tall1,Vall1,T1,R1,Gav1]=CallLIF_C(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);    % call C simulator
%[tall1,Vall1,T1,R1,Gav1]=LIFNetML(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);   % call ML simulator

p0([6 9])=LAPar1([6 9])*1.15;   % adjust recurr. syn. to produce faster climbing activity
[tall2,Vall2,T2,R2,Gav2]=CallLIF_C(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);    % call C simulator
%[tall2,Vall2,T2,R2,Gav2]=LIFNetML(p0,tiv,Iaff,TList,GsynAvg,NoisePar,s,Flags);   % call ML simulator

save SimResExp2 tall1 Vall1 T1 R1 Gav1 tall2 Vall2 T2 R2 Gav2;

subplot(2,2,2), hold on, plot(T1./1000,R1,'Color',[0.8 0.8 0.8],'LineWidth',0.8);
subplot(2,2,4), hold on, plot(T2./1000,R2,'Color',[0.8 0.8 0.8],'LineWidth',0.8);



% Copyright (C) 2002,2003 Daniel Durstewitz