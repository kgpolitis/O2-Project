% load data for 4 procs
s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1])
f = fittype('a*x^k','options',s)
b_nu=[20 40 60 80 100]';
b_n=b_nu'.^3;
b_nwork_p4=b_nu'.^2/4;
b_nwork_t2_p4=2*b_nu.*(2*b_nu/4-3);
b_nwork_t3_p4=2*b_nu.*(2*b_nu/4-7);
% process 1
b_t1_1_p4=VarName3(1:9:end);
b_e1_1_p4=VarName3(2:9:end);
b_t2_1_p4=VarName3(3:9:end);
b_e2_1_p4=VarName3(5:9:end);
b_t3_1_p4=VarName3(6:9:end);
b_e3_1_p4=VarName3(8:9:end);
b_d1_1_p4=VarName3(4:9:end);
b_d2_1_p4=VarName3(7:9:end);
b_T_1_p4=VarName3(9:9:end);
% process 2
b_t1_2_p4=VarName4(1:9:end);
b_e1_2_p4=VarName4(2:9:end);
b_t2_2_p4=VarName4(3:9:end);
b_e2_2_p4=VarName4(5:9:end);
b_t3_2_p4=VarName4(6:9:end);
b_e3_2_p4=VarName4(8:9:end);
b_d1_2_p4=VarName4(4:9:end);
b_d2_2_p4=VarName4(7:9:end);
b_T_2_p4=VarName4(9:9:end);
% process 3
b_t1_3_p4=VarName5(1:9:end);
b_e1_3_p4=VarName5(2:9:end);
b_t2_3_p4=VarName5(3:9:end);
b_e2_3_p4=VarName5(5:9:end);
b_t3_3_p4=VarName5(6:9:end);
b_e3_3_p4=VarName5(8:9:end);
b_d1_3_p4=VarName5(4:9:end);
b_d2_3_p4=VarName5(7:9:end);
b_T_3_p4=VarName5(9:9:end);
% process 4
b_t1_4_p4=VarName6(1:9:end);
b_e1_4_p4=VarName6(2:9:end);
b_t2_4_p4=VarName6(3:9:end);
b_e2_4_p4=VarName6(5:9:end);
b_t3_4_p4=VarName6(6:9:end);
b_e3_4_p4=VarName6(8:9:end);
b_d1_4_p4=VarName6(4:9:end);
b_d2_4_p4=VarName6(7:9:end);
b_T_4_p4=VarName6(9:9:end);
% plot
figure
loglog(b_nwork_p4,b_t1_1_p4,'k')
hold
loglog(b_nwork_p4,b_t1_2_p4,'b')
loglog(b_nwork_p4,b_t1_3_p4,'r')
loglog(b_nwork_p4,b_t1_4_p4,'g')
title('topos 1')
figure
loglog(b_nwork_p4,b_e1_1_p4,'k')
hold
loglog(b_nwork_p4,b_e1_2_p4,'b')
loglog(b_nwork_p4,b_e1_3_p4,'r')
loglog(b_nwork_p4,b_e1_4_p4,'g')
title('extend 1')
figure
loglog(b_nwork_p4,b_e2_1_p4,'k')
hold
loglog(b_nwork_p4,b_e2_2_p4,'b')
loglog(b_nwork_p4,b_e2_3_p4,'r')
loglog(b_nwork_p4,b_e2_4_p4,'g')
title('extend 2')
figure
loglog(b_nwork_p4,b_e3_1_p4,'k')
hold
loglog(b_nwork_p4,b_e3_2_p4,'b')
loglog(b_nwork_p4,b_e3_3_p4,'r')
loglog(b_nwork_p4,b_e3_4_p4,'g')
title('extend 3')
figure
loglog(b_nwork_p4,b_d1_1_p4,'k')
hold
loglog(b_nwork_p4,b_d1_2_p4,'b')
loglog(b_nwork_p4,b_d1_3_p4,'r')
loglog(b_nwork_p4,b_d1_4_p4,'g')
title('DGM 1')
figure
loglog(b_nwork_p4,b_d2_1_p4,'k')
hold
loglog(b_nwork_p4,b_d2_2_p4,'b')
loglog(b_nwork_p4,b_d2_3_p4,'r')
loglog(b_nwork_p4,b_d2_4_p4,'g')
title('DGM 2')
figure
loglog(b_nwork_p4,b_T_1_p4,'k')
hold
loglog(b_nwork_p4,b_T_2_p4,'b')
loglog(b_nwork_p4,b_T_3_p4,'r')
loglog(b_nwork_p4,b_T_4_p4,'g')
title('total time')
% actual times
b_t1_p4=min([b_t1_1_p4 b_t1_2_p4 b_t1_3_p4 b_t1_4_p4],[],2);
b_t2_p4=min([b_t2_1_p4 b_t2_2_p4 b_t2_3_p4 b_t2_4_p4],[],2);
b_t3_p4=min([b_t3_1_p4 b_t3_2_p4 b_t3_3_p4 b_t3_4_p4],[],2);
b_e1_p4=min([b_e1_1_p4 b_e1_2_p4 b_e1_3_p4 b_e1_4_p4],[],2);
b_e2_p4=min([b_e2_1_p4 b_e2_2_p4 b_e2_3_p4 b_e2_4_p4],[],2);
b_e3_p4=min([b_e3_1_p4 b_e3_2_p4 b_e3_3_p4 b_e3_4_p4],[],2);
b_d1_p4=min([b_d1_1_p4 b_d1_2_p4 b_d1_3_p4 b_d1_4_p4],[],2);
b_d2_p4=min([b_d2_1_p4 b_d2_2_p4 b_d2_3_p4 b_d2_4_p4],[],2);
b_T_p4=min([b_T_1_p4 b_T_2_p4 b_T_3_p4 b_T_4_p4],[],2);
%a_t1_p4=mean([a_t1_1_p4 a_t1_2_p4 a_t1_3_p4 a_t1_4_p4],2);
%a_e1_p4=mean([a_e1_1_p4 a_e1_2_p4 a_e1_3_p4 a_e1_4_p4],2);
%a_e2_p4=mean([a_e2_1_p4 a_e2_2_p4 a_e2_3_p4 a_e2_4_p4],2);
%a_e3_p4=mean([a_e3_1_p4 a_e3_2_p4 a_e3_3_p4 a_e3_4_p4],2);
%a_d1_p4=mean([a_d1_1_p4 a_d1_2_p4 a_d1_3_p4 a_d1_4_p4],2);
%a_d2_p4=mean([a_d2_1_p4 a_d2_2_p4 a_d2_3_p4 a_d2_4_p4],2);
%a_T_p4=mean([a_T_1_p4 a_T_2_p4 a_T_3_p4 a_T_4_p4],2);
% plots
figure
loglog(b_nwork_p4,b_t1_p4)
title('topos 1')
figure
loglog(b_nwork_p4,b_e1_p4)
title('extend 1')
figure
loglog(b_nwork_p4,b_e2_p4)
title('extend 2')
figure
loglog(b_nwork_p4,b_e3_p4)
title('extend 3')
figure
loglog(b_nwork_p4,b_d1_p4)
title('DGM 1')
figure
loglog(b_nwork_p4,b_d2_p4)
title('DGM 2')
figure
plot(b_nwork_p4,b_T(1:5)./b_T_1_p4/4,'k')
hold
plot(b_nwork_p4,b_T(1:5)./b_T_2_p4/4,'b')
plot(b_nwork_p4,b_T(1:5)./b_T_3_p4/4,'r')
plot(b_nwork_p4,b_T(1:5)./b_T_4_p4/4,'g')
title('eff')
%close all
b_eff_p4=b_T(1:5)./(4*b_T_p4);
b_eff0_p4=b_T(1:5)./(b_t1_p4+b_e1_p4+b_e2_p4+b_e3_p4+b_t2_p4+b_t3_p4+b_d1_p4+b_d2_p4)/4
b_cs_p4=1-(b_t1_p4+b_e1_p4+b_t2_p4+b_t3_p4+b_e2_p4+b_e3_p4+b_d1_p4+b_d2_p4)./b_T_p4;
b_p4_res=[b_t1_p4 b_e1_p4 b_t2_p4 b_d1_p4 b_e2_p4 b_t3_p4 b_d2_p4 b_e3_p4 b_T_p4];
b_p4_mod=[b_t1_p4+b_e1_p4 b_t2_p4+b_e2_p4+b_d1_p4 b_t3_p4+b_e3_p4+b_d2_p4 b_t1_p4+b_e1_p4+b_t2_p4+b_e2_p4+b_d1_p4+b_t3_p4+b_e3_p4+b_d2_p4 b_cs_p4]
[b_t1n_p4,gof]=fit(b_nwork_p4',b_t1_p4,f,s)
[b_t2n_p4,gof]=fit(b_nwork_t2_p4,b_t2_p4,f,s)
[b_t3n_p4,gof]=fit(b_nwork_t3_p4,b_t3_p4,f,s)
[b_e1n_p4,gof]=fit(b_nwork_p4',b_e1_p4,f,s)
[b_e2n_p4,gof]=fit(b_nwork_p4',b_e2_p4,f,s)
[b_e3n_p4,gof]=fit(b_nwork_p4',b_e3_p4,f,s)
b_a_p4=b_t1n_p4.a+b_t2n_p4.a+b_t3n_p4.a+b_e1n_p4.a+b_e2n_p4.a+b_e3n_p4.a
ideal_eff_p4=b_a_s/b_a_p4
%To_p4_Ts=1./b_eff_p4-b_a_p4/b_a_s