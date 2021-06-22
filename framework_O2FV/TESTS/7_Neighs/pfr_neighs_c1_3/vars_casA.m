% load data for 3 procs for neighs new : with total time added
a_nu=[20 40 60 80 100]';
a_n=a_nu'.^3;
a_nwork_p3=a_nu'.^3/3;
% process 1
a_t1_1_p3=VarName7(1:9:end);
a_e1_1_p3=VarName7(2:9:end);
a_t2_1_p3=VarName7(3:9:end);
a_e2_1_p3=VarName7(5:9:end);
a_t3_1_p3=VarName7(6:9:end);
a_e3_1_p3=VarName7(8:9:end);
a_d1_1_p3=VarName7(4:9:end);
a_d2_1_p3=VarName7(7:9:end);
a_T_1_p3=VarName7(8:9:end);
% process 2
a_t1_2_p3=VarName8(1:9:end);
a_e1_2_p3=VarName8(2:9:end);
a_t2_2_p3=VarName8(3:9:end);
a_e2_2_p3=VarName8(5:9:end);
a_t3_2_p3=VarName8(6:9:end);
a_e3_2_p3=VarName8(8:9:end);
a_d1_2_p3=VarName8(4:9:end);
a_d2_2_p3=VarName8(7:9:end);
a_T_2_p3=VarName8(9:9:end);
% process 3
a_t1_3_p3=VarName9(1:9:end);
a_e1_3_p3=VarName9(2:9:end);
a_t2_3_p3=VarName9(3:9:end);
a_e2_3_p3=VarName9(5:9:end);
a_t3_3_p3=VarName9(6:9:end);
a_e3_3_p3=VarName9(8:9:end);
a_d1_3_p3=VarName9(4:9:end);
a_d2_3_p3=VarName9(7:9:end);
a_T_3_p3=VarName9(9:9:end);
% plots
figure
loglog(a_nwork_p3,a_t1_1_p3,'k')
hold
loglog(a_nwork_p3,a_t1_2_p3,'b')
loglog(a_nwork_p3,a_t1_3_p3,'r')
title('topos 1')
figure
loglog(a_nwork_p3,a_e1_1_p3,'k')
hold
loglog(a_nwork_p3,a_e1_2_p3,'b')
loglog(a_nwork_p3,a_e1_3_p3,'r')
title('extend 1')
figure
loglog(a_nwork_p3,a_e2_1_p3,'k')
hold
loglog(a_nwork_p3,a_e2_2_p3,'b')
loglog(a_nwork_p3,a_e2_3_p3,'r')
title('extend 2')
figure
loglog(a_nwork_p3,a_e3_1_p3,'k')
hold
loglog(a_nwork_p3,a_e3_2_p3,'b')
loglog(a_nwork_p3,a_e3_3_p3,'r')
title('extend 3')
figure
loglog(a_nwork_p3,a_d1_1_p3,'k')
hold
loglog(a_nwork_p3,a_d1_2_p3,'b')
loglog(a_nwork_p3,a_d1_3_p3,'r')
title('DGM 1')
figure
loglog(a_nwork_p3,a_d2_1_p3,'k')
hold
loglog(a_nwork_p3,a_d2_2_p3,'b')
loglog(a_nwork_p3,a_d2_3_p3,'r')
title('DGM 2')
figure
loglog(a_nwork_p3,a_T_1_p3,'k')
hold
loglog(a_nwork_p3,a_T_2_p3,'b')
loglog(a_nwork_p3,a_T_3_p3,'r')
title('total time')
% actual times
a_t1_p3=min([a_t1_1_p3 a_t1_2_p3 a_t1_3_p3],[],2);
a_e1_p3=min([a_e1_1_p3 a_e1_2_p3 a_e1_3_p3],[],2);
a_e2_p3=min([a_e2_1_p3 a_e2_2_p3 a_e2_3_p3],[],2);
a_e3_p3=min([a_e3_1_p3 a_e3_2_p3 a_e3_3_p3],[],2);
a_d1_p3=min([a_d1_1_p3 a_d1_2_p3 a_d1_3_p3],[],2);
a_d2_p3=min([a_d2_1_p3 a_d2_2_p3 a_d2_3_p3],[],2);
a_T_p3=max([a_T_1_p3 a_T_2_p3 a_T_3_p3],[],2);
%a_t1_p3=mean([a_t1_1_p3 a_t1_2_p3 a_t1_3_p3],2);
%a_e1_p3=mean([a_e1_1_p3 a_e1_2_p3 a_e1_3_p3],2);
%a_e2_p3=mean([a_e2_1_p3 a_e2_2_p3 a_e2_3_p3],2);
%a_e3_p3=mean([a_e3_1_p3 a_e3_2_p3 a_e3_3_p3],2);
%a_d1_p3=mean([a_d1_1_p3 a_d1_2_p3 a_d1_3_p3],2);
%a_d2_p3=mean([a_d2_1_p3 a_d2_2_p3 a_d2_3_p3],2);
%a_T_p3=mean([a_T_1_p3 a_T_2_p3 a_T_3_p3],2);
% plots
figure
loglog(a_nwork_p3,a_t1_p3)
title('topos 1')
figure
loglog(a_nwork_p3,a_e1_p3)
title('extend 1')
figure
loglog(a_nwork_p3,a_e2_p3)
title('extend 2')
figure
loglog(a_nwork_p3,a_e3_p3)
title('extend 3')
figure
loglog(a_nwork_p3,a_d1_p3)
title('DGM 1')
figure
loglog(a_nwork_p3,a_d2_p3)
title('DGM 2')
figure
loglog(a_nwork_p3,a_T_p3)
title('total time')