% load data for 4 procs
s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1])
f = fittype('a*x^k','options',s)
a_nu=[20 40 60 80 100]';
a_n=a_nu'.^3;
a_nwork_p4=a_nu'.^3/4;
% process 1
a_t1_1_p4=VarName3(1:9:end);
a_e1_1_p4=VarName3(2:9:end);
a_t2_1_p4=VarName3(3:9:end);
a_e2_1_p4=VarName3(5:9:end);
a_t3_1_p4=VarName3(6:9:end);
a_e3_1_p4=VarName3(8:9:end);
a_d1_1_p4=VarName3(4:9:end);
a_d2_1_p4=VarName3(7:9:end);
a_T_1_p4=VarName3(9:9:end);
% process 2
a_t1_2_p4=VarName4(1:9:end);
a_e1_2_p4=VarName4(2:9:end);
a_t2_2_p4=VarName4(3:9:end);
a_e2_2_p4=VarName4(5:9:end);
a_t3_2_p4=VarName4(6:9:end);
a_e3_2_p4=VarName4(8:9:end);
a_d1_2_p4=VarName4(4:9:end);
a_d2_2_p4=VarName4(7:9:end);
a_T_2_p4=VarName4(9:9:end);
% process 3
a_t1_3_p4=VarName5(1:9:end);
a_e1_3_p4=VarName5(2:9:end);
a_t2_3_p4=VarName5(3:9:end);
a_e2_3_p4=VarName5(5:9:end);
a_t3_3_p4=VarName5(6:9:end);
a_e3_3_p4=VarName5(8:9:end);
a_d1_3_p4=VarName5(4:9:end);
a_d2_3_p4=VarName5(7:9:end);
a_T_3_p4=VarName5(9:9:end);
% process 4
a_t1_4_p4=VarName6(1:9:end);
a_e1_4_p4=VarName6(2:9:end);
a_t2_4_p4=VarName6(3:9:end);
a_e2_4_p4=VarName6(5:9:end);
a_t3_4_p4=VarName6(6:9:end);
a_e3_4_p4=VarName6(8:9:end);
a_d1_4_p4=VarName6(4:9:end);
a_d2_4_p4=VarName6(7:9:end);
a_T_4_p4=VarName6(9:9:end);
% plot
figure
loglog(a_nwork_p4,a_t1_1_p4,'k')
hold
loglog(a_nwork_p4,a_t1_2_p4,'b')
loglog(a_nwork_p4,a_t1_3_p4,'r')
loglog(a_nwork_p4,a_t1_4_p4,'g')
title('topos 1')
figure
loglog(a_nwork_p4,a_e1_1_p4,'k')
hold
loglog(a_nwork_p4,a_e1_2_p4,'b')
loglog(a_nwork_p4,a_e1_3_p4,'r')
loglog(a_nwork_p4,a_e1_4_p4,'g')
title('extend 1')
figure
loglog(a_nwork_p4,a_e2_1_p4,'k')
hold
loglog(a_nwork_p4,a_e2_2_p4,'b')
loglog(a_nwork_p4,a_e2_3_p4,'r')
loglog(a_nwork_p4,a_e2_4_p4,'g')
title('extend 2')
figure
loglog(a_nwork_p4,a_e3_1_p4,'k')
hold
loglog(a_nwork_p4,a_e3_2_p4,'b')
loglog(a_nwork_p4,a_e3_3_p4,'r')
loglog(a_nwork_p4,a_e3_4_p4,'g')
title('extend 3')
figure
loglog(a_nwork_p4,a_d1_1_p4,'k')
hold
loglog(a_nwork_p4,a_d1_2_p4,'b')
loglog(a_nwork_p4,a_d1_3_p4,'r')
loglog(a_nwork_p4,a_d1_4_p4,'g')
title('DGM 1')
figure
loglog(a_nwork_p4,a_d2_1_p4,'k')
hold
loglog(a_nwork_p4,a_d2_2_p4,'b')
loglog(a_nwork_p4,a_d2_3_p4,'r')
loglog(a_nwork_p4,a_d2_4_p4,'g')
title('DGM 2')
figure
loglog(a_nwork_p4,a_T_1_p4,'k')
hold
loglog(a_nwork_p4,a_T_2_p4,'b')
loglog(a_nwork_p4,a_T_3_p4,'r')
loglog(a_nwork_p4,a_T_4_p4,'g')
title('total time')
% actual times
a_t1_p4=min([a_t1_1_p4 a_t1_2_p4 a_t1_3_p4 a_t1_4_p4],[],2);
a_e1_p4=min([a_e1_1_p4 a_e1_2_p4 a_e1_3_p4 a_e1_4_p4],[],2);
a_e2_p4=min([a_e2_1_p4 a_e2_2_p4 a_e2_3_p4 a_e2_4_p4],[],2);
a_e3_p4=min([a_e3_1_p4 a_e3_2_p4 a_e3_3_p4 a_e3_4_p4],[],2);
a_d1_p4=min([a_d1_1_p4 a_d1_2_p4 a_d1_3_p4 a_d1_4_p4],[],2);
a_d2_p4=min([a_d2_1_p4 a_d2_2_p4 a_d2_3_p4 a_d2_4_p4],[],2);
a_T_p4=max([a_T_1_p4 a_T_2_p4 a_T_3_p4 a_T_4_p4],[],2);
%a_t1_p4=mean([a_t1_1_p4 a_t1_2_p4 a_t1_3_p4 a_t1_4_p4],2);
%a_e1_p4=mean([a_e1_1_p4 a_e1_2_p4 a_e1_3_p4 a_e1_4_p4],2);
%a_e2_p4=mean([a_e2_1_p4 a_e2_2_p4 a_e2_3_p4 a_e2_4_p4],2);
%a_e3_p4=mean([a_e3_1_p4 a_e3_2_p4 a_e3_3_p4 a_e3_4_p4],2);
%a_d1_p4=mean([a_d1_1_p4 a_d1_2_p4 a_d1_3_p4 a_d1_4_p4],2);
%a_d2_p4=mean([a_d2_1_p4 a_d2_2_p4 a_d2_3_p4 a_d2_4_p4],2);
%a_T_p4=mean([a_T_1_p4 a_T_2_p4 a_T_3_p4 a_T_4_p4],2);
% plots
figure
loglog(a_nwork_p4,a_t1_p4)
title('topos 1')
figure
loglog(a_nwork_p4,a_e1_p4)
title('extend 1')
figure
loglog(a_nwork_p4,a_e2_p4)
title('extend 2')
figure
loglog(a_nwork_p4,a_e3_p4)
title('extend 3')
figure
loglog(a_nwork_p4,a_d1_p4)
title('DGM 1')
figure
loglog(a_nwork_p4,a_d2_p4)
title('DGM 2')
figure
plot(a_nwork_p4,a_T(1:5)./a_T_1_p4/4,'k')
hold
plot(a_nwork_p4,a_T(1:5)./a_T_2_p4/4,'b')
plot(a_nwork_p4,a_T(1:5)./a_T_3_p4/4,'r')
plot(a_nwork_p4,a_T(1:5)./a_T_4_p4/4,'g')
title('eff')
%close all
eff_p4=a_T(1:5)./(4*a_T_p4);
eff0_p4=a_T(1:5)./(a_t1_p4+a_e1_p4+a_e2_p4+a_e3_p4+a_d1_p4+a_d2_p4)/4
a_cs_p4=1-(a_t1_p4+a_e1_p4+a_e2_p4+a_e3_p4+a_d1_p4+a_d2_p4)./a_T_p4;
a_p2_res=[a_t1_p4' a_e1_p4' 0 a_e2_p4' 0 a_e3_p4' a_T_p4'];
a_p2_mod=[a_t1_p4'+a_e1_p4' a_e2_p4' a_e3_p4']
[a_t1n_p4,gof]=fit(a_nwork_p4',a_t1_p4,f,s)
[a_e1n_p4,gof]=fit(a_nwork_p4',a_e1_p4,f,s)
[a_e2n_p4,gof]=fit(a_nwork_p4',a_e2_p4,f,s)
[a_e3n_p4,gof]=fit(a_nwork_p4',a_e3_p4,f,s)
a_a_p4=a_t1n_p4.a+a_e1n_p4.a+a_e2n_p4.a+a_e3n_p4.a
To_p4_Ts=1./eff_p4-a_a_p4/a_a_s