% load data for 2 procs for neighs new : with total time added : case 2
s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1])
f = fittype('a*x^k','options',s)
b_nu=[20 40 60 80 100]';
b_n=b_nu'.^3;
b_nwork_p2=b_nu'.^2/2;
b_nwork_t2_p2=2*b_nu.*(2*b_nu/2-3);
b_nwork_t3_p2=2*b_nu.*(2*b_nu/2-7);
% process 1
b_t1_1_p2=VarName1(1:9:end);
b_e1_1_p2=VarName1(2:9:end);
b_t2_1_p2=VarName1(3:9:end);
b_e2_1_p2=VarName1(5:9:end);
b_t3_1_p2=VarName1(6:9:end);
b_e3_1_p2=VarName1(8:9:end);
b_d1_1_p2=VarName1(4:9:end);
b_d2_1_p2=VarName1(7:9:end);
b_T_1_p2=VarName1(8:9:end);
% process 2
b_t1_2_p2=VarName2(1:9:end);
b_e1_2_p2=VarName2(2:9:end);
b_t2_2_p2=VarName2(3:9:end);
b_e2_2_p2=VarName2(5:9:end);
b_t3_2_p2=VarName2(6:9:end);
b_e3_2_p2=VarName2(8:9:end);
b_d1_2_p2=VarName2(4:9:end);
b_d2_2_p2=VarName2(7:9:end);
b_T_2_p2=VarName2(9:9:end);
% plots
figure
loglog(b_nwork_p2,b_t1_1_p2,'k')
hold
loglog(b_nwork_p2,b_t1_2_p2,'b')
title('topos 1')
figure
loglog(b_nwork_p2,b_e1_1_p2,'k')
hold
loglog(b_nwork_p2,b_e1_2_p2,'b')
title('extend 1')
figure
loglog(b_nwork_p2,b_t2_1_p2,'k')
hold
loglog(b_nwork_p2,b_t2_2_p2,'b')
title('topos 2')
figure
loglog(b_nwork_p2,b_e2_1_p2,'k')
hold
loglog(b_nwork_p2,b_e2_2_p2,'b')
title('extend 2')
figure
loglog(b_nwork_p2,b_t3_1_p2,'k')
hold
loglog(b_nwork_p2,b_t3_2_p2,'b')
title('topos 3')
figure
loglog(b_nwork_p2,b_e3_1_p2,'k')
hold
loglog(b_nwork_p2,b_e3_2_p2,'b')
title('extend 3')
figure
loglog(b_nwork_p2,b_d1_1_p2,'k')
hold
loglog(b_nwork_p2,b_d1_2_p2,'b')
title('DGM 1')
figure
loglog(b_nwork_p2,b_d2_1_p2,'k')
hold
loglog(b_nwork_p2,b_d2_2_p2,'b')
title('DGM 2')
figure
plot(b_nwork_p2,b_T(1:5)./b_T_1_p2/2,'k')
hold
plot(b_nwork_p2,b_T(1:5)./b_T_2_p2/2,'b')
title('eff')
% actual times -> min
%a_t1_p2=a_t1_2_p2 ;
%a_e1_p2=a_e1_2_p2 ;
%a_e2_p2=a_e2_2_p2 ;
%a_e3_p2=a_e3_2_p2 ;
%a_d1_p2=a_d1_2_p2 ;
%a_d2_p2=a_d2_2_p2 ;
%a_T_p2=a_T_2_p2 ;
b_t1_p2=min([b_t1_1_p2 b_t1_2_p2],[],2);
b_e1_p2=min([b_e1_1_p2 b_e1_2_p2],[],2);
b_t2_p2=min([b_t2_1_p2 b_t2_2_p2],[],2);
b_e2_p2=min([b_e2_1_p2 b_e2_2_p2],[],2);
b_t3_p2=min([b_t3_1_p2 b_t3_2_p2],[],2);
b_e3_p2=min([b_e3_1_p2 b_e3_2_p2],[],2);
b_d1_p2=min([b_d1_1_p2 b_d1_2_p2],[],2);
b_d2_p2=min([b_d2_1_p2 b_d2_2_p2],[],2);
b_T_p2=max([b_T_1_p2 b_T_2_p2],[],2);
%a_t1_p2=mean([a_t1_1_p2 a_t1_2_p2],2);
%a_e1_p2=mean([a_e1_1_p2 a_e1_2_p2],2);
%a_e2_p2=mean([a_e2_1_p2 a_e2_2_p2],2);
%a_e3_p2=mean([a_e3_1_p2 a_e3_2_p2],2);
%a_d1_p2=mean([a_d1_1_p2 a_d1_2_p2],2);
%a_d2_p2=mean([a_d2_1_p2 a_d2_2_p2],2);
%a_T_p2=mean([a_T_1_p2 a_T_2_p2],2);
% plots
figure
loglog(b_nwork_p2,b_t1_p2)
title('topos 1')
figure
loglog(b_nwork_p2,b_e1_p2)
title('extend 1')
figure
loglog(b_nwork_p2,b_t2_p2)
title('topos 2')
figure
loglog(b_nwork_p2,b_e2_p2)
title('extend 2')
figure
loglog(b_nwork_p2,b_t3_p2)
title('topos 3')
figure
loglog(b_nwork_p2,b_e3_p2)
title('extend 3')
figure
loglog(b_nwork_p2,b_d1_p2)
title('DGM 1')
figure
loglog(b_nwork_p2,b_d2_p2)
title('DGM 2')
figure
loglog(b_nwork_p2,b_T_p2)
title('total time')
%close all
b_eff_p2=b_T(1:5)./(2*b_T_p2);
b_eff0_p2=b_T(1:5)./(b_t1_p2+b_e1_p2+b_t2_p2+b_e2_p2+b_t3_p2+b_e3_p2+b_d1_p2+b_d2_p2)/2
b_cs_p2=1-(b_t1_p2+b_e1_p2+b_e2_p2+b_e3_p2+b_t2_p2+b_t3_p2+b_d1_p2+b_d2_p2)./b_T_p2;
b_p2_res=[b_t1_p2' b_e1_p2' b_t2_p2' b_e2_p2' b_t3_p2' b_e3_p2' b_T_p2'];
b_p2_mod=[b_t1_p2'+b_e1_p2' b_t2_p2'+b_e2_p2' b_t2_p2'+b_e3_p2']
[b_t1n_p2,gof]=fit(b_nwork_p2',b_t1_p2,f,s)
[b_t2n_p2,gof]=fit(b_nwork_t2_p2,b_t2_p2,f,s)
[b_t3n_p2,gof]=fit(b_nwork_t3_p2,b_t3_p2,f,s)
[b_e1n_p2,gof]=fit(b_nwork_p2',b_e1_p2,f,s)
[b_e2n_p2,gof]=fit(b_nwork_p2',b_e2_p2,f,s)
[b_e3n_p2,gof]=fit(b_nwork_p2',b_e3_p2,f,s)
b_a_p2=b_t1n_p2.a+b_t2n_p2.a+b_t3n_p2.a+b_e1n_p2.a+b_e2n_p2.a+b_e3n_p2.a
ideal_eff_p2=b_a_s/b_a_p2
%To_p2_Ts=1./b_eff_p2-b_a_p2/b_a_s
