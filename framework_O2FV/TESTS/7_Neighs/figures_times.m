% Generates figures for the results obtained by timing the cell 
% neighborhood  -> use with times.mat 
% construction algorithm in SERIAL 
% Actions performed -> Regression Formula
%                      Plot raw data and regression
s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1])
f = fittype('a*x^k','options',s)
ii=1
jj=5
% figure topos times vs n - results and regression
figure
[t_topos_n_s,gof2] = fit(n(ii:jj),a_t1(ii:jj),f,s)
plot(n(ii:jj),a_t1(ii:jj),'o')
hold
plot(t_topos_n_s,'b')
grid on
legend('t_{topos} - Results','t_{topos}(n_C)=3.002e-6 \cdot n_C^{1.019}','Location','Northwest')
xlabel('n_C','Interpreter','tex')
ylabel('t(s)','Interpreter','tex','rotation',0)
%
% figure extend times vs n - results and regression
figure
[t_extend1_n_s,gof2] = fit(n(ii:jj),a_e1(ii:jj),f,s)
[t_extend2_n_s,gof2] = fit(n(ii:jj),a_e2(ii:jj),f,s)
[t_extend3_n_s,gof2] = fit(n(ii:jj),a_e3(ii:jj),f,s)
a_a_s=t_topos_n_s.a+t_extend1_n_s.a+t_extend2_n_s.a+t_extend3_n_s.a
loglog(n(ii:jj),a_e1(ii:jj),'ob')
hold
plot(t_extend1_n_s,'b')
plot(n(ii:jj),a_e2(ii:jj),'or')
plot(t_extend2_n_s,'r')
plot(n(ii:jj),a_e3(ii:jj),'og')
plot(t_extend3_n_s,'g')
grid on
legend('t_{extend}^1 - Results','t_{extend}^1(n_C)=3.786e-7 \cdot n_C^{0.995}', ...
       't_{extend}^2 - Results','t_{extend}^2(n_C)=2.375e-5 \cdot n_C^{1.018}', ...
       't_{extend}^3 - Results','t_{extend}^2(n_C)=1.369e-4 \cdot n_C^{1.032}', ...
'Location','Northwest')
xlabel('n_C','Interpreter','tex')
ylabel('t(s)','Interpreter','tex','rotation',0)
%
% figure t_topos^2/t_topos^1 t_topos^3/t_topos^2
figure
plot(n(ii:jj),b_t2(ii:jj)./b_t1(ii:jj),'o')
hold
plot(n(ii:jj),b_t3(ii:jj)./b_t2(ii:jj),'go')
legend('t_{topos}^2/t_{topos}^1','t_{topos}^3/t_{topos}^2')
grid on
xlabel('n_C','Interpreter','tex')
% figure R
syms R x
p=1
%R=3/3.8*x.^(2/3)./(9*x.^(1/3)-20*p)
R=x.^(2/3)./(9*x.^(1/3)-20*p)
figure
plot(n(ii:jj),a_t1(ii:jj)./(b_t1(ii:jj)+b_t2(ii:jj)+b_t3(ii:jj)+b_cs(ii:jj).*b_T(ii:jj)),'o')
%plot(n(ii:jj),3.8/3*a_t1(ii:jj)./(b_t1(ii:jj)+b_t2(ii:jj)+b_t3(ii:jj)),'o')
hold
ezplot(R,[n(ii),n(jj)])
legend('Speedup Ratio R - calculated','Speedup Ratio R - theoretical','Location','NorthWest')
grid on
xlabel('n_C','Interpreter','tex')
% figure topos times case 1 vs case 2x
% figure extend times vs n - results and regression
figure
[b_topos_all,gof2] = fit(b_nwork_t(ii:jj),b_t(ii:jj),f,s)
loglog(n(ii:jj),a_t1(ii:jj),'o')
hold
plot(t_topos_n_s,'b')
loglog(b_nwork_t(ii:jj),b_t(ii:jj),'or')
plot(b_topos_all,'r')
grid on
xlabel('n_C - total workload','Interpreter','tex')
ylabel('t(s)','Interpreter','tex','rotation',0)
legend('t_{topos}^1 - Case 1','t_{topos}^1(n_C)=3.002e-6 \cdot n_C^{1.019}', ...
       'sum_k t_{topos}^k - Case 2','t_{topos}(n_C)=3.801e-6 \cdot n_C^{1.011}', ...
'Location','Northwest')
% figure total speedup ratio for level 1, level 2 and level 3 
figure
loglog(n(1:5).^(2/3),(a_t1(1:5)+b_e1(1:5))./(b_t1(1:5)+b_e1(1:5)))
hold
plot(n(1:5).^(2/3),(a_t1(1:5)+b_e1(1:5)+b_e2(1:5))./(b_t1(1:5)+b_t2(1:5)+b_e1(1:5)+b_e2(1:5)),'r')
plot(n(1:5).^(2/3),(a_t1(1:5)+b_e1(1:5)+b_e2(1:5)+b_e3(1:5))./(b_t(1:5)+b_e1(1:5)+b_e2(1:5)+b_e3(1:5)),'g')
grid on
xlabel('n_C - number of cells in input block','Interpreter','tex')
legend('R^1','R^2','R^3', ...
'Location','Northwest')
ii=1:5
[b_t1n_s,gof2] = fit(nu(ii:jj).^2,b_t1(ii:jj),f,s)
[b_t2n_s,gof2] = fit(2*nu(ii:jj).*(2*nu(ii:jj)-3),b_t2(ii:jj),f,s)
[b_t3n_s,gof2] = fit(2*nu(ii:jj).*(2*nu(ii:jj)-7),b_t3(ii:jj),f,s)
[b_e1n_s,gof2] = fit(nu(ii:jj).^2,b_e1(ii:jj),f,s)
[b_e2n_s,gof2] = fit(nu(ii:jj).^2,b_e2(ii:jj),f,s)
[b_e3n_s,gof2] = fit(nu(ii:jj).^2,b_e3(ii:jj),f,s)
b_a_s=b_t1n_s.a+b_t2n_s.a+b_t3n_s.a+b_e1n_s.a+b_e2n_s.a+b_e3n_s.a
