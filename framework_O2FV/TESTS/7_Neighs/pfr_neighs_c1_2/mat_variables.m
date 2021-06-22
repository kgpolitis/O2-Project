% load data for 2 procs
a_nu=[20 40 60 80 100]';
a_n=a_nu'.^3;
a_nwork_p2=a_nu'.^3/2;
% process 1
a_t1_1_p2=VarName1(1:8:40);
a_e1_1_p2=VarName1(2:8:40);
a_t2_1_p2=VarName1(3:8:40);
a_e2_1_p2=VarName1(5:8:40);
a_t3_1_p2=VarName1(6:8:40);
a_e3_1_p2=VarName1(8:8:40);
a_d1_1_p2=VarName1(4:8:40);
a_d2_1_p2=VarName1(7:8:40);
% process 2
a_t1_2_p2=VarName2(1:8:40);
a_e1_2_p2=VarName2(2:8:40);
a_t2_2_p2=VarName2(3:8:40);
a_e2_2_p2=VarName2(5:8:40);
a_t3_2_p2=VarName2(6:8:40);
a_e3_2_p2=VarName2(8:8:40);
a_d1_2_p2=VarName2(4:8:40);
a_d2_2_p2=VarName2(7:8:40);
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
% actual times
a_t1_p2=min([a_t1_1_p2 a_t1_2_p2],[],2);
a_e1_p2=min([a_e1_1_p2 a_e1_2_p2],[],2);
a_e2_p2=min([a_e2_1_p2 a_e2_2_p2],[],2);
a_e3_p2=min([a_e3_1_p2 a_e3_2_p2],[],2);
a_d1_p2=min([a_d1_1_p2 a_d1_2_p2],[],2);
a_d2_p2=min([a_d2_1_p2 a_d2_2_p2],[],2);
% plots
figure
loglog(a_nwork_p4,a_t1_p2)
title('topos 1')
figure
loglog(a_nwork_p4,a_e1_p2)
title('extend 1')
figure
loglog(a_nwork_p4,a_e2_p2)
title('extend 2')
figure
loglog(a_nwork_p4,a_e3_p2)
title('extend 3')
figure
loglog(a_nwork_p4,a_d1_p2)
title('DGM 1')
figure
loglog(a_nwork_p4,a_d2_p2)
title('DGM 2')