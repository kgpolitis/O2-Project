% Standard Figure for wave curvature
x=linspace(-1,1,1000)
z=1/2*sin(2*pi*x)
k=-2*pi^2*sin(2*pi*x)./(sqrt(pi^2*cos(2*pi*x).^2+1).^3)
[haxes,hline1,hline2]=plotyy(x,z,x,k)
set(gca,'fontsize',14);
set(gcf,'Color','w')
set(gcf,'Units','Points');
Pos=get(gcf,'Position');
set(gcf,'Position',[Pos(1) Pos(2) 460 380]);
set(hline1,'color','b','linestyle','--')
axes(haxes(1))
ylim([-1 1])
set(gca,'YColor','k')
xlabel('x','fontsize',16)
ylabel('z','rotation',0,'fontsize',16,'color','k')
axes(haxes(2))
set(gca,'fontsize',14);
ylabel('\kappa','rotation',0,'fontsize',16,'color','k')
set(gca,'YColor','k')
grid
legend('\kappa(x,y)','z(x,y)','fontsize',16,'location','Northwest')
set(hline2,'color','k')
