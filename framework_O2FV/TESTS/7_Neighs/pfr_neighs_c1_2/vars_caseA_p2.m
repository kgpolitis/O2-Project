% load data for 2 procs for neighs new : with total time added
s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1])
f = fittype('a*x^k','options',s)
a_nu=[20 40 60 80 100]';
a_n=a_nu'.^3;
a_nwork_p2=a_nu'.^3/2;
% process 1
a_t1_1_p2=VarName1(1:9:end);
a_e1_1_p2=VarName1(2:9:end);
a_t2_1_p2=VarName1(3:9:end);
a_e2_1_p2=VarName1(5:9:end);
a_t3_1_p2=VarName1(6:9:end);
a_e3_1_p2=VarName1(8:9:end);
a_d1_1_p2=VarName1(4:9:end);
a_d2_1_p2=VarName1(7:9:end);
a_T_1_p2=VarName1(8:9:end);
% process 2
a_t1_2_p2=VarName2(1:9:end);
a_e1_2_p2=VarName2(2:9:end);
a_t2_2_p2=VarName2(3:9:end);
a_e2_2_p2=VarName2(5:9:end);
a_t3_2_p2=VarName2(6:9:end);
a_e3_2_p2=VarName2(8:9:end);
a_d1_2_p2=VarName2(4:9:end);
a_d2_2_p2=VarName2(7:9:end);
a_T_2_p2=VarName2(9:9:end);
% plots
figure
loglog(a_nwork_p2,a_t1_1_p2,'k')
hold
loglog(a_nwork_p2,a_t1_2_p2,'b')
title('topos 1')
figure
loglog(a_nwork_p2,a_e1_1_p2,'k')
hold
loglog(a_nwork_p2,a_e1_2_p2,'b')
title('extend 1')
figure
loglog(a_nwork_p2,a_e2_1_p2,'k')
hold
loglog(a_nwork_p2,a_e2_2_p2,'b')
title('extend 2')
figure
loglog(a_nwork_p2,a_e3_1_p2,'k')
hold
loglog(a_nwork_p2,a_e3_2_p2,'b')
title('extend 3')
figure
loglog(a_nwork_p2,a_d1_1_p2,'k')
hold
loglog(a_nwork_p2,a_d1_2_p2,'b')
title('DGM 1')
figure
loglog(a_nwork_p2,a_d2_1_p2,'k')
hold
loglog(a_nwork_p2,a_d2_2_p2,'b')
title('DGM 2')
figure
plot(a_nwork_p2,a_T(1:5)./a_T_1_p2/2,'k')
hold
plot(a_nwork_p2,a_T(1:5)./a_T_2_p2/2,'b')
title('eff')
% actual times -> min
%a_t1_p2=a_t1_2_p2 ;
%a_e1_p2=a_e1_2_p2 ;
%a_e2_p2=a_e2_2_p2 ;
%a_e3_p2=a_e3_2_p2 ;
%a_d1_p2=a_d1_2_p2 ;
%a_d2_p2=a_d2_2_p2 ;
%a_T_p2=a_T_2_p2 ;
a_t1_p2=min([a_t1_1_p2 a_t1_2_p2],[],2);
a_e1_p2=min([a_e1_1_p2 a_e1_2_p2],[],2);
a_e2_p2=min([a_e2_1_p2 a_e2_2_p2],[],2);
a_e3_p2=min([a_e3_1_p2 a_e3_2_p2],[],2);
a_d1_p2=min([a_d1_1_p2 a_d1_2_p2],[],2);
a_d2_p2=min([a_d2_1_p2 a_d2_2_p2],[],2);
a_T_p2=max([a_T_1_p2 a_T_2_p2],[],2);
%a_t1_p2=mean([a_t1_1_p2 a_t1_2_p2],2);
%a_e1_p2=mean([a_e1_1_p2 a_e1_2_p2],2);
%a_e2_p2=mean([a_e2_1_p2 a_e2_2_p2],2);
%a_e3_p2=mean([a_e3_1_p2 a_e3_2_p2],2);
%a_d1_p2=mean([a_d1_1_p2 a_d1_2_p2],2);
%a_d2_p2=mean([a_d2_1_p2 a_d2_2_p2],2);
%a_T_p2=mean([a_T_1_p2 a_T_2_p2],2);
% plots
figure
loglog(a_nwork_p2,a_t1_p2)
title('topos 1')
figure
loglog(a_nwork_p2,a_e1_p2)
title('extend 1')
figure
loglog(a_nwork_p2,a_e2_p2)
title('extend 2')
figure
loglog(a_nwork_p2,a_e3_p2)
title('extend 3')
figure
loglog(a_nwork_p2,a_d1_p2)
title('DGM 1')
figure
loglog(a_nwork_p2,a_d2_p2)
title('DGM 2')
figure
loglog(a_nwork_p2,a_T_p2)
title('total time')
%close all
eff_p2=a_T(1:5)./(2*a_T_p2);
eff0_p2=a_T(1:5)./(a_t1_p2+a_e1_p2+a_e2_p2+a_e3_p2+a_d1_p2+a_d2_p2)/2
a_cs_p2=1-(a_t1_p2+a_e1_p2+a_e2_p2+a_e3_p2+a_d1_p2+a_d2_p2)./a_T_p2;
a_p2_res=[a_t1_p2 a_e1_p2 zeros(5,1) a_d1_p2 a_e2_p2 zeros(5,1) a_d2_p2 a_e3_p2 a_T_p2];
a_p2_mod=[a_t1_p2+a_e1_p2 a_e2_p2+a_d1_p2 a_e3_p2+a_d2_p2 a_t1_p2+a_e1_p2+a_e2_p2+a_e3_p2+a_d1_p2+a_d2_p2 a_cs_p2]
[a_t1n_p2,gof]=fit(a_nwork_p2',a_t1_p2,f,s)
[a_e1n_p2,gof]=fit(a_nwork_p2',a_e1_p2,f,s)
[a_e2n_p2,gof]=fit(a_nwork_p2',a_e2_p2,f,s)
[a_e3n_p2,gof]=fit(a_nwork_p2',a_e3_p2,f,s)
a_a_p2=a_t1n_p2.a+a_e1n_p2.a+a_e2n_p2.a+a_e3n_p2.a
ideal_eff_p2=a_a_s/a_a_p2
%To_p2_Ts=1./eff_p2-a_a_p2/a_a_s
