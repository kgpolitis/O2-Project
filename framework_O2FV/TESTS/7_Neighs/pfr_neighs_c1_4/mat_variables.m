% load data for 2 procs
a_nu=[20 40 60 80 100]';
a_n=a_nu'.^3;
a_nwork_p4=a_nu'.^3/4;
% process 1
a_t1_1_p4=VarName3(1:8:40);
a_e1_1_p4=VarName3(2:8:40);
a_t2_1_p4=VarName3(3:8:40);
a_e2_1_p4=VarName3(5:8:40);
a_t3_1_p4=VarName3(6:8:40);
a_e3_1_p4=VarName3(8:8:40);
a_d1_1_p4=VarName3(4:8:40);
a_d2_1_p4=VarName3(7:8:40);
% process 2
a_t1_2_p4=VarName4(1:8:40);
a_e1_2_p4=VarName4(2:8:40);
a_t2_2_p4=VarName4(3:8:40);
a_e2_2_p4=VarName4(5:8:40);
a_t3_2_p4=VarName4(6:8:40);
a_e3_2_p4=VarName4(8:8:40);
a_d1_2_p4=VarName4(4:8:40);
a_d2_2_p4=VarName4(7:8:40);
% process 3
a_t1_3_p4=VarName5(1:8:40);
a_e1_3_p4=VarName5(2:8:40);
a_t2_3_p4=VarName5(3:8:40);
a_e2_3_p4=VarName5(5:8:40);
a_t3_3_p4=VarName5(6:8:40);
a_e3_3_p4=VarName5(8:8:40);
a_d1_3_p4=VarName5(4:8:40);
a_d2_3_p4=VarName5(7:8:40);
% process 4
a_t1_4_p4=VarName6(1:8:40);
a_e1_4_p4=VarName6(2:8:40);
a_t2_4_p4=VarName6(3:8:40);
a_e2_4_p4=VarName6(5:8:40);
a_t3_4_p4=VarName6(6:8:40);
a_e3_4_p4=VarName6(8:8:40);
a_d1_4_p4=VarName6(4:8:40);
a_d2_4_p4=VarName6(7:8:40);
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
% actual times
a_t1_p4=min([a_t1_1_p4 a_t1_2_p4 a_t1_3_p4 a_t1_4_p4],[],2);
a_e1_p4=min([a_e1_1_p4 a_e1_2_p4 a_e1_3_p4 a_e1_4_p4],[],2);
a_e2_p4=min([a_e2_1_p4 a_e2_2_p4 a_e2_3_p4 a_e2_4_p4],[],2);
a_e3_p4=min([a_e3_1_p4 a_e3_2_p4 a_e3_3_p4 a_e3_4_p4],[],2);
a_d1_p4=min([a_d1_1_p4 a_d1_2_p4 a_d1_3_p4 a_d1_4_p4],[],2);
a_d2_p4=min([a_d2_1_p4 a_d2_2_p4 a_d2_3_p4 a_d2_4_p4],[],2);
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