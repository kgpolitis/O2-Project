% figure t_{DGM}^k - Case A
figure
plot(n(1:5),a_d1_p2,'b')
hold
plot(n(1:5),a_d1_p4,'g')
plot(n(1:5),a_d2_p2,'b--')
plot(n(1:5),a_d2_p4,'g--')
grid on
legend('t_{DGM}^2 - p=2','t_{DGM}^2 - p=4','t_{DGM}^3 - p=2','t_{DGM}^3 - p=4','location','NorthWest')
%title('Case 1')
xlabel('n_C')
ylabel('t (s)','Rotation',0)
% figure t_{DGM}^k - Case B
figure
plot(n(1:5),b_d1_p2,'b')
hold
plot(n(1:5),b_d1_p4,'g')
plot(n(1:5),b_d2_p2,'b--')
plot(n(1:5),b_d2_p4,'g--')
grid on
legend('t_{DGM}^2 - p=2','t_{DGM}^2 - p=4','t_{DGM}^3 - p=2','t_{DGM}^3 - p=4','location','NorthWest')
%title('Case 2')
xlabel('n_C')
ylabel('t (s)','Rotation',0)
% figure extend serial vs extend parallel: case 1
figure
loglog(n(1:5),a_e2(1:5))
hold
plot(a_nwork_p2',a_e2_p2,'g','linewidth',2)
plot(a_nwork_p4',a_e2_p4,'r')
plot(n(1:5),a_e3(1:5),'--')
plot(a_nwork_p2',a_e3_p2,'g--','linewidth',2)
plot(a_nwork_p4',a_e3_p4,'r--')
grid on
legend('t_{extend}^2 - serial','t_{extend}^2 - p=2','t_{extend}^3 - p=4',...
       't_{extend}^3 - serial','t_{extend}^3 - p=2','t_{extend}^3 - p=4',...
'location','NorthWest')
xlabel('n_{C_i}')
ylabel('t (s)','Rotation',0)
% figure extend serial vs extend parallel:case 2
figure
loglog(nu.^2,b_e2(1:5))
hold
plot(b_nwork_p2',b_e2_p2,'g','linewidth',2)
plot(b_nwork_p4',b_e2_p4,'r')
plot(nu.^2,b_e3(1:5),'--')
plot(b_nwork_p2',b_e3_p2,'g--','linewidth',2)
plot(b_nwork_p4',b_e3_p4,'r--')
grid on
legend('t_{extend}^2 - serial','t_{extend}^2 - p=2','t_{extend}^3 - p=4',...
       't_{extend}^3 - serial','t_{extend}^3 - p=2','t_{extend}^3 - p=4',...
'location','NorthWest')
xlabel('n_{C_i}')
ylabel('t (s)','Rotation',0)
% figure topos serial vs extend parallel
%figure
%loglog(n(1:5),a_t1(1:5))
%hold
%plot(a_nwork_p2',a_t1_p2,'g','linewidth',2)
%plot(a_nwork_p4',a_t1_p4,'r')
%grid on
%legend('t_{extend}^2 - serial','t_{extend}^2 - p=2','t_{extend}^3 - p=4',...
%       't_{extend}^3 - serial','t_{extend}^3 - p=2','t_{extend}^3 - p=4',...
%'location','NorthWest')
%xlabel('n_{C_i}')
%ylabel('t (s)','Rotation',0)
%figure
%plot(a_nwork_p2',2*a_e2_p2./a_e2(1:5),'b')
%hold
%plot(a_nwork_p4',4*a_e2_p4./a_e2(1:5),'g')
%plot(a_nwork_p2',2*a_e3_p2./a_e3(1:5),'--')
%plot(a_nwork_p4',4*a_e3_p4./a_e3(1:5),'g--')
%grid on
speed_up_DGM1_p2=a_d1_p2./b_d1_p2;
speed_up_DGM2_p2=a_d2_p2./b_d2_p2;
speed_up_DGM1_p4=a_d1_p4./b_d1_p4;
speed_up_DGM2_p4=a_d2_p4./b_d2_p4;
figure
plot(n(1:5),speed_up_DGM1_p2,'b')
hold
plot(n(1:5),speed_up_DGM1_p4,'g')
plot(n(1:5),speed_up_DGM2_p2,'b--')
plot(n(1:5),speed_up_DGM2_p4,'g--')
grid on
legend('R_{DMG}^2 - p=2','R_{DGM}^2 - p=4','R_{DGM}^3 - p=2','R_{DGM}^3 - p=4','location','NorthWest')
xlabel('n_C')
ylabel('t (s)','Rotation',0)
figure
plot(n(1:5),eff_p2,'b-')
hold
plot(n(1:5),eff_p4,'g-')
plot(n(1:5),b_eff_p2,'b--')
plot(n(1:5),b_eff_p4,'g--')
xlabel('n_C')
ylabel('Efficiency')
legend('Case 1 - p=2','Case 1 - p=4','Case 2 - p=2','Case 2 - p=4','location','SouthEast')
grid on
figure
plot(n(1:5),2*eff_p2,'b-')
hold
plot(n(1:5),4*eff_p4,'g-')
plot(n(1:5),2*b_eff_p2,'b--')
plot(n(1:5),4*b_eff_p4,'g--')
xlabel('n_C')
ylabel('Speedup')
legend('Case 1 - p=2','Case 1 - p=4','Case 2 - p=2','Case 2 - p=4','location','SouthEast')
grid on