isodata10=importdata('results_vfi_n10.dat')
isodata20=importdata('results_vfi_n20.dat')
isodata25=importdata('results_vfi_n25.dat')
isodata30=importdata('results_vfi_n30.dat')
isodata40=importdata('results_vfi_n40.dat')
isodata50=importdata('results_vfi_n50.dat')
isodata55=importdata('results_vfi_n55.dat')
isodata80=importdata('results_vfi_n80.dat')
isodata90=importdata('results_vfi_n90.dat')
isodata100=importdata('results_vfi_n100.dat')
iplot_1=0
iplot_2=1
iplot_3=1
if (iplot_1==1)
%
loglog(2:4,isodata10(1,:),'o','MarkerFaceColor','b','MarkerSize',7)
hold
loglog(2:4,isodata10(2,:),'go','MarkerFaceColor','g','MarkerSize',7)
loglog(2:4,isodata10(3,:),'ro','MarkerFaceColor','r','MarkerSize',7)
set(gca,'Xtick',[2 3 4])
%
loglog(2:4,isodata20(1,:),'>','MarkerFaceColor','b','MarkerSize',7)
loglog(2:4,isodata20(2,:),'>g','MarkerFaceColor','g','MarkerSize',7)
loglog(2:4,isodata20(3,:),'>r','MarkerFaceColor','r','MarkerSize',7)
%
%loglog(2:4,isodata25(1,:),'<','MarkerFaceColor','b','MarkerSize',7)
%loglog(2:4,isodata25(2,:),'<g','MarkerFaceColor','g','MarkerSize',7)
%loglog(2:4,isodata25(3,:),'<r','MarkerFaceColor','r','MarkerSize',7)
%
%loglog(2:4,isodata30(1,:),'s','MarkerFaceColor','b','MarkerSize',7)
%loglog(2:4,isodata30(2,:),'sg','MarkerFaceColor','g','MarkerSize',7)
%loglog(2:4,isodata30(3,:),'sr','MarkerFaceColor','r','MarkerSize',7)
%
loglog(2:4,isodata40(1,:),'d','MarkerFaceColor','b','MarkerSize',7)
loglog(2:4,isodata40(2,:),'dg','MarkerFaceColor','g','MarkerSize',7)
loglog(2:4,isodata40(3,:),'dr','MarkerFaceColor','r','MarkerSize',7)
%
%loglog(2:4,isodata50(1,:),'*','MarkerFaceColor','b','MarkerSize',7)
%loglog(2:4,isodata50(2,:),'*g','MarkerFaceColor','g','MarkerSize',7)
%loglog(2:4,isodata50(3,:),'*r','MarkerFaceColor','r','MarkerSize',7)
%
%
loglog(2:4,isodata80(1,:),'v','MarkerFaceColor','b','MarkerSize',7)
loglog(2:4,isodata80(2,:),'vg','MarkerFaceColor','g','MarkerSize',7)
loglog(2:4,isodata80(3,:),'vr','MarkerFaceColor','r','MarkerSize',7)
%
loglog(2:4,isodata100(1,:),'^','MarkerFaceColor','b','MarkerSize',7)
loglog(2:4,isodata100(2,:),'^g','MarkerFaceColor','g','MarkerSize',7)
loglog(2:4,isodata100(3,:),'^r','MarkerFaceColor','r','MarkerSize',7)
end
%
fsz_label=16
fsz_leg=16
fsz_ax=14
N=[10 20 30 40 80 90 100]
dx=2./[10 20 30 40 80 90 100]
y=dx.^2
if (iplot_2==1)
lvl=1
figure
loglog(2./dx,[isodata10(lvl,1) isodata20(lvl,1) isodata30(lvl,1) isodata40(lvl,1) isodata80(lvl,1) isodata90(lvl,1) isodata100(lvl,1)])
hold
loglog(2./dx,[isodata10(lvl,2) isodata20(lvl,2) isodata30(lvl,2) isodata40(lvl,2) isodata80(lvl,2) isodata90(lvl,2) isodata100(lvl,2)],'g')
loglog(2./dx,[isodata10(lvl,3) isodata20(lvl,3) isodata30(lvl,3) isodata40(lvl,3) isodata80(lvl,3) isodata90(lvl,3) isodata100(lvl,3)],'r')
loglog(2./dx,y,'k')
grid on
xlabel('$$ \frac{1}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa )  $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx,'fontsize',fsz_ax)
lg=legend('order=2','order=3','order=4')
set(lg,'fontsize',fsz_leg)
title('Level 1 Neighborhoods')
lvl=2
figure
loglog(2./dx,[isodata10(lvl,1) isodata20(lvl,1) isodata30(lvl,1) isodata40(lvl,1) isodata80(lvl,1) isodata90(lvl,1) isodata100(lvl,1)])
hold
loglog(2./dx,[isodata10(lvl,2) isodata20(lvl,2) isodata30(lvl,2) isodata40(lvl,2) isodata80(lvl,2) isodata90(lvl,2) isodata100(lvl,2)],'g')
loglog(2./dx,[isodata10(lvl,3) isodata20(lvl,3) isodata30(lvl,3) isodata40(lvl,3) isodata80(lvl,3) isodata90(lvl,3) isodata100(lvl,3)],'r')
loglog(2./dx,y,'k')
grid on
xlabel('$$ \frac{1}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa )  $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx,'fontsize',fsz_ax)
lg=legend('order=2','order=3','order=4')
set(lg,'fontsize',fsz_leg)
title('Level 2 Neighborhoods')
lvl=3
figure
loglog(2./dx,[isodata10(lvl,1) isodata20(lvl,1) isodata30(lvl,1) isodata40(lvl,1) isodata80(lvl,1) isodata90(lvl,1) isodata100(lvl,1)])
hold
loglog(2./dx,[isodata10(lvl,2) isodata20(lvl,2) isodata30(lvl,2) isodata40(lvl,2) isodata80(lvl,2) isodata90(lvl,2) isodata100(lvl,2)],'g')
loglog(2./dx,[isodata10(lvl,3) isodata20(lvl,3) isodata30(lvl,3) isodata40(lvl,3) isodata80(lvl,3) isodata90(lvl,3) isodata100(lvl,3)],'r')
loglog(2./dx,y,'k')
grid on
xlabel('$$ \frac{2R}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa ) ) $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx,'fontsize',fsz_ax)
lg=legend('order=2','order=3','order=4')
title('Level 3 Neighborhoods')
set(lg,'fontsize',fsz_leg)
lvl=4
figure
loglog(2./dx,[isodata10(lvl,1) isodata20(lvl,1) isodata30(lvl,1) isodata40(lvl,1) isodata80(lvl,1) isodata90(lvl,1) isodata100(lvl,1)])
hold
loglog(2./dx,[isodata10(lvl,2) isodata20(lvl,2) isodata30(lvl,2) isodata40(lvl,2) isodata80(lvl,2) isodata90(lvl,2) isodata100(lvl,2)],'g')
loglog(2./dx,[isodata10(lvl,3) isodata20(lvl,3) isodata30(lvl,3) isodata40(lvl,3) isodata80(lvl,3) isodata90(lvl,3) isodata100(lvl,3)],'r')
loglog(2./dx,y,'k')
grid on
xlabel('$$ \frac{1}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa )  $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx,'fontsize',fsz_ax)
lg=legend('order=2','order=3','order=4')
set(lg,'fontsize',fsz_leg)
title('Level 4 Neighborhoods')
lvl=5
figure
loglog(2./dx,[isodata10(lvl,1) isodata20(lvl,1) isodata30(lvl,1) isodata40(lvl,1) isodata80(lvl,1) isodata90(lvl,1) isodata100(lvl,1)])
hold
loglog(2./dx,[isodata10(lvl,2) isodata20(lvl,2) isodata30(lvl,2) isodata40(lvl,2) isodata80(lvl,2) isodata90(lvl,2) isodata100(lvl,2)],'g')
loglog(2./dx,[isodata10(lvl,3) isodata20(lvl,3) isodata30(lvl,3) isodata40(lvl,3) isodata80(lvl,3) isodata90(lvl,3) isodata100(lvl,3)],'r')
loglog(2./dx,y,'k')
grid on
xlabel('$$ \frac{1}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa )  $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx,'fontsize',fsz_ax)
lg=legend('order=2','order=3','order=4')
set(lg,'fontsize',fsz_leg)
title('Level 5 Neighborhoods')
end
if (iplot_3==1)
figure
ord=1
lvl=2
hold
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'r')
lvl=3
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'g')
lvl=4
loglog(2./dx(3:end),[isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'b')
lvl=5
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'c')
grid on
xlabel('$$ \frac{1}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa )  $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx(3:end),'fontsize',fsz_ax)
lg=legend('level=2','level=3','level=4','level=5')
set(lg,'fontsize',fsz_leg)
title('Order 2')
figure
ord=2
lvl=2
hold
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'r')
lvl=3
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'g')
lvl=4
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'b')
lvl=5
loglog(2./dx(3:end),[ isodata30(lvl,ord) isodata40(lvl,ord) isodata80(lvl,ord) isodata90(lvl,ord) isodata100(lvl,ord)],'c')
grid on
xlabel('$$ \frac{1}{Dx} $$','interpreter','latex','fontsize',fsz_label)
ylabel('$$ L_{max} (\kappa )  $$','interpreter','latex','fontsize',fsz_label)
set(gca,'XTick',2./dx(3:end),'fontsize',fsz_ax)
lg=legend('level=2','level=3','level=4','level=5')
set(lg,'fontsize',fsz_leg)
title('Order 3')
end