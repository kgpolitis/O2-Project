% Error Plots Script
% This Script Plots error by a given file
% 
% Input: File and data associations
% The first line of the file must be the name of the variable
% Note that:
%     1.you must specify the variable that will be plotted and 
%       the dependent variable
%     2.you must store the data as columns and each column must begin
%       with the name of the variable 
%     3.the dependant and independant variables are specified by 
%       the indices i_iv and i_dv respectively 
%
clear
% The name of the file you store the data:
% Normal/ Curvature for exact isosurface
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geo_resnormal.txt'
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geo_rescurv.txt'
% Normal/ Curvature for Ci reconstructed isosurface
% filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geoCi_resnormal.txt'
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geoCi_rescurv.txt'
% Normal/ Curvature for Ci-K1 reconstructed isosurface
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geosmCi_resnormal.txt'
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geosmCi_rescurv.txt'
% Normal/ Curvature for Ci-K2 reconstructed isosurface
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geosm2Ci_resnormal.txt'
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geosm2Ci_rescurv.txt'
% Comparison weights/no weights
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geo_wk_comp.txt'
% Compact Stencils calculations
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_normal_fromISO.txt'
% compact stencil curvature comparisons
% filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geo_isop_tperson.txt'
% Normal/ Curvature smooth/no smooth comparisons LSqR
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_sph_geo_no_and_smCi_rescurv_l3o2.txt'
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/sm_errs_wav_geo_no_and_smCi_rescurv_l2o3_nobnd.txt'
% Normal Curvatur smooth/no smooth comparisons volume methods
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/res_sph_comp_smvsnosm_volmeth.tex'
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/res_sph_comp_smvsnosm_volmeth_normal.tex'
% Normal Curvatur smooth/no smooth comparisons all methods
%filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/res_sph_comp_all_vol_surf.txt'
filename='~/O2-project/framework_O2FV/TESTS/15_NC_test/res_sph_comp_all_vol_surf_normal.txt'
% Get data
CD=importdata(filename);
n=size(CD.data,1);
i_dv=[];
% Which column stores the dependant and the independant variable
% Note that idy is the second independant variable if you work want an
% error surface (if not it should be zero). i_skip is the number of
% valued you want to skip (default=1), i_start is the starting index
% (default-1) and i_end is the ending index (default=n)
i_iv=1
i_iv2=0
%----
% Lmax
mi_dv(1)=4
mi_dv(2)=6
mi_dv(3)=8
mi_dv(4)=10
mi_dv(5)=12
mi_dv(6)=14
mi_dv(7)=16
mi_dv(8)=18
mi_dv(9)=20
mi_dv(10)=22
% L1
%mi_dv=mi_dv+1
%----
% all 4 normal
%i_dv=mi_dv
% all 4 curv
%i_dv=mi_dv(1:9)
% isopatch / order 2 or level 1/2
%i_dv(1:4)=mi_dv([1 2 5 8]);
% all expect isopatch
%i_dv(1:9)=mi_dv(2:10)
% order 2/3
%i_dv(1)=mi_dv(2)
%i_dv(2)=mi_dv(3)
%i_dv(3)=mi_dv(5)
%i_dv(4)=mi_dv(6)
%i_dv(5)=mi_dv(8)
%i_dv(6)=mi_dv(9)
% normal level 1
%i_dv=mi_dv(1:4)
% normal level 2 
% i_dv=mi_dv(5:7)
% normal level 3
% i_dv=mi_dv(8:10)
% curv level 1
% i_dv=mi_dv(1:3)
% curv level 2 
% i_dv=mi_dv(4:6)
% curv level 3
% i_dv=mi_dv(7:9)
use_order=0
use_level=0
if (use_order==1) 
% normal/curv order 
is_curv=0 % if curv result this must be 1
order=3 % > for order 2 : order =1 for order 3 order =3 for order 4 order =3
i_dv=mi_dv([2 5 8]-is_curv+order-1)
elseif (use_level==1)
% normal/curv lvl
is_curv=0 % if curv result this must be 1
level=1 % > for order 2 : order =1 for order 3 order =3 for order 4 order =3
i_dv=mi_dv([(level-1)*3+2:level*3+1]-is_curv)
end
i_dv=i_dv+1

% weights/no weights comparison
% level 1
%i_dv(1)=4
%i_dv(2)=6
%i_dv(3)=12
%i_dv(4)=14
% level 2
%i_dv(1)=8
%i_dv(2)=10
%i_dv(3)=16
%i_dv(4)=18

% normals from iso
% here we compare smooth / no smooth approaches
%i_dv(1)=4  % exact
%i_dv(2)=6  % ci
%i_dv(3)=8  % k1(ci)
%i_dv(4)=10 % k2(ci)
%i_dv(3)=12 % L1(ci)
%i_dv(3)=14 % L2(ci)
%i_dv(4)=16 % L4(ci)
%i_dv(5)=12 % L1(ci)
%i_dv(5)=14 % L2(ci)
%i_dv(6)=16 % L4(ci)
% for L1:
%i_dv=i_dv+1

% curvature tpersonal
% exact normals
%i_dv(1)=6
%i_dv(2)=8
%i_dv(3)=10
%i_dv(4)=4
%i_dv=i_dv+1
% approx normals
%i_dv(1)=14
%i_dv(2)=16
%i_dv(3)=18
%i_dv(4)=12
%i_dv=i_dv+1
% comparison all curvature averaging
%i_dv(1)=4
%i_dv(2)=12
%i_dv(3)=20
%i_dv=i_dv+1

% exact isosurface results
%i_dv(1)=4
%i_dv(2)=6
%i_dv(3)=8
%i_dv(4)=10
%i_dv(5)=12
%i_dv(6)=14

% exact isosurface results for LSqR: normals
%i_dv(1)=4
%i_dv(2)=6
%i_dv(3)=8
%i_dv(4)=10
%i_dv(5)=12
%i_dv(6)=14
% for LSqR: curvature 
%i_dv(1)=4
%i_dv(2)=8
%i_dv(3)=10
%i_dv(4)=12
%i_dv(5)=14
%i_dv(6)=16
% for L1
%i_dv=i_dv+1

% comparisons smooth/no smooth for LSqR
% l2o2
%i_dv=[6 8 10 12]
%i_dv=[6 4 12 14]
%i_dv=i_dv+1

% comparisons smooth/no smooth for vol methods
%i_dv=[4 5 6 7 8 9]

% comparisons smooth/no smooth for vol methods
i_dv=[4 6 8]
%i_dv=i_dv+1

i_skip=1
i_start=1
i_end=n
i_use=[i_start:i_skip:i_end];
%i_use=[1 6 16 21 31 36 46] 

% The scale factors for the dependent variables if any ... 
% Note: must be 1 if not used; this is useful when you have only stored
%       the number of discretizations   
s_iv=1/2
s_iv2=1
s_dv(1:length(i_dv))=1

% Options 
% -------
%    iv_minus_1: if 1 redefines the independant variable as:
%                s_iv/column_specified
%                
%    revaxis: if 1 the x-axis orientation is reversed 
%
%    i_tv : column of the ticks variable by default it is equal to the
%           i_iv (the same as the independent variable). If they are 
%           different then this variable of the column i_tv is 
%           plotted in another x axis above
%                                                                                         
%    use_fits: the script finds and plots the fitted curves or surface
%             1-> fit function: a*x^k
%             2-> fit function: polynomial 2nd degree
%             3-> fit function: polynomial 3rd degree
%                                                                           
%    default_lims : change default x axis limit to correspond to the
%                   values defined by x_min_extend and x_max_extend
%    x_min_extend : units to extend the x axis lower limit before min 
%                   value of x (similar for x_max
%   
%    order_line_k : order of line identifying error
%
%------ Defaults: do not change!!
iv_minus_1=0%
revaxis=0%
i_tv=i_iv%
symbol(1)={'ob'}%
symbol(2)={'or'}%
use_fits=1%
default_lims=0%
x_min_ext=3%
x_max_ext=x_min_ext%
order_line_k=2%

%----- You may change these (or just comment them if not required):
iv_minus_1=0
revaxis=0
i_tv=i_iv
symb(1)={'ok'}
symb(2)={'ob'}
symb(3)={'og'}
symb(4)={'or'}
symb(5)={'sqb'}
symb(6)={'sqg'}
symb(7)={'sqr'}
symb(8)={'db'}
symb(9)={'dg'}
symb(10)={'dr'}
% all 4 normal
%symbol=symb
% all 4 curv
%symbol=symb(2:10)
% isopatch order 2
%symbol(1:4)=symb([1 2 5 8])
% all besides isopatch
%symbol(1:9)=symb(2:10)
% order 2/3
%symbol(1)=symb(2)
%symbol(2)=symb(3)
%symbol(3)=symb(5)
%symbol(4)=symb(6)
%symbol(5)=symb(8)
%symbol(6)=symb(9)
% normal level 1
%symbol=symb(2:4)
% normal level 2
%symbol=symb(5:7)
% normal level 3
%symbol=symb(8:10)
% curv level 1
%symbol=symb(1:3)
% curv level 2
%symbol=symb(4:6)
% curv level 3
%symbol=symb(7:8)
if (use_order==1) 
% normal order order
symbol=symb([2 5 8]+order-1)
elseif (use_level==1)
% normal/curv lvl
symbol=symb([(level-1)*3+2:level*3+1])
end

% weights/no weights comparison
% level 1
%symbol(1)={'ob'}
%symbol(2)={'or'}
%symbol(3)={'o-b'}
%symbol(4)={'o-r'}
% level 2
%symbol(1)={'sqb'}
%symbol(2)={'sqr'}
%symbol(3)={'sq-b'}
%symbol(4)={'sq-r'}

% normals from iso
% here we compare smooth / no smooth approaches
%symbol(1)={'-k'}
%symbol(2)={'--k'}
%symbol(3)={'ob'}
%symbol(4)={'sqb'}
%symbol(3)={'om'}
%symbol(4)={'sqm'}
%symbol(5)={'dm'}
%symbol(5)={'om'}
%symbol(6)={'sqm'}
%symbol(7)={'dm'}

% tukoshin
%symbol(1)={'ob'}
%symbol(2)={'sqb'}
%symbol(3)={'db'}
%symbol(4)={'om'}
%symbol(3)={'dk'}
%symbol(3)={'ob'}
%symbol(4)={'sqb'}
% comparison all
%symbol(1)={'ok'}
%symbol(2)={'ob'}
%symbol(3)={'or'}

% exact isosurface results for LSqR: normals
%symbol(1)={'ok'}
%symbol(2)={'ob'}
%symbol(3)={'og'}
%symbol(4)={'sqg'}
%symbol(5)={'sqr'}
%symbol(6)={'dr'}
%exact isosurface results for LSqR: curv
%symbol(1)={'ob'}
%symbol(2)={'og'}
%symbol(3)={'sqg'}
%symbol(4)={'sqr'}
%symbol(5)={'dr'}

% comparisons smooth/no smooth for LSqR
% l2o2
%symbol(1)={'ok'}
%symbol(2)={'ob'}
%symbol(3)={'og'}
%symbol(4)={'or'}
%symbol(1)={'ok'}
%symbol(2)={'ok-'}
%symbol(3)={'or'}
%symbol(4)={'or-'}

% comparisons smooth/no smooth for vol methods
%symbol(1)={'ob'}
%symbol(2)={'sqb'}
%symbol(3)={'og'}
%symbol(4)={'sqg'}
%symbol(5)={'or'}
%symbol(6)={'sqr'}

% comparisons smooth/no smooth for all methods
symbol(1)={'ob'}
symbol(2)={'og'}
symbol(3)={'or'}


use_fits=0
default_lims=0
x_min_ext=3
x_max_ext=x_min_ext
order_line_k=[-2]

%---- You have to change these:
%xname='\Delta\theta'
%labelx=['$$', xname,'=\frac{2\pi}{n_\theta} =2\Delta\phi$$']
labelx=['$$ \frac{2R}{\delta} $$']
yname='\vec{n}'
mlabely(1)={['$$ Isopatch $$']}
mlabely(2)={['$$ Lvl 1 - Ord 2$$']}
mlabely(3)={['$$ Lvl 1 - Ord 3$$']}
mlabely(4)={['$$ Lvl 1 - Ord 4$$']}
mlabely(5)={['$$ Lvl 2 - Ord 2$$']}
mlabely(6)={['$$ Lvl 2 - Ord 3$$']}
mlabely(7)={['$$ Lvl 2 - Ord 4$$']}
mlabely(8)={['$$ Lvl 3 - Ord 2$$']}
mlabely(9)={['$$ Lvl 3 - Ord 3$$']}
mlabely(10)={['$$ Lvl 3 - Ord 4$$']}
% all 4 normal
%labely=mlabely
% all 4 curv
%labely=mlabely(2:10)
% isopatch and order 2
%labely(1:4)=mlabely([1 2 5 8])
% 
% all expect isopatches
%labely(1:9)=mlabely(2:10)
%
% Order 2/3
%labely(1)=mlabely(2)
%labely(2)=mlabely(3)
%labely(3)=mlabely(5)
%labely(4)=mlabely(6)
%labely(5)=mlabely(8)
%labely(6)=mlabely(9)
% normal level 1
%labely=mlabely(2:4)
% normal level 2
%labely=mlabely(5:7)
% normal level 3
%labely=mlabely(8:10)
% curv level 1
%labely=mlabely(1:3)
% curv level 2
%labely=mlabely(4:6)
% curv level 3
%labely=mlabely(7:8)
if (use_order==1)
% normal order order
labely=mlabely([2 5 8]+order-1)
elseif (use_level==1)
% normal/curv lvl
labely=mlabely([(level-1)*3+2:level*3+1])
end

% weights /no weights comparison
% level 1
%labely(1)={['$$ Lvl 1 - Ord 2 - No\,weights$$']}
%labely(2)={['$$ Lvl 1 - Ord 4 - No\,weights$$']}
%labely(3)={['$$ Lvl 1 - Ord 2 - Weights $$']}
%labely(4)={['$$ Lvl 1 - Ord 4 - Weights $$']}
% level 2
%labely(1)={['$$ Lvl 2 - Ord 2 - No\,weights$$']}
%labely(2)={['$$ Lvl 2 - Ord 4 - No\,weights$$']}
%labely(3)={['$$ Lvl 2 - Ord 2 - Weights $$']}
%labely(4)={['$$ Lvl 2 - Ord 4 - Weights $$']}

% normals from iso
% here we compare smooth / no smooth approaches
%labely(1)={['$$ Exact $$']}
%labely(2)={['$$ C_I $$']}
%labely(3)={['$$ K1(C_I) $$']}
%labely(4)={['$$ K2(C_I) $$']}
%labely(3)={['$$ La1(C_I) $$']}
%labely(3)={['$$ La2(C_I) $$']}
%labely(4)={['$$ La4(C_I) $$']}
%labely(5)={['$$ La1(C_I) $$']}
%labely(5)={['$$ La2(C_I) $$']}
%labely(6)={['$$ La4(C_I) $$']}

% Tukoshin
%labely(1)={['$$ \lambda:Bisection $$']}
%labely(2)={['$$ \lambda:Parallel\,to\,LR $$']}
%labely(3)={['$$ \lambda:Circle\,Tangent $$']}
%labely(4)={['$$ Curvature\,Averaging $$']}
%labely(4)={['$$ CA - Approx\, \vec{n} $$']}
%labely(3)={['$$ Simpson Rule $$']}
%labely(1)={['$$ Exact\,Isosurface\,and\,Normals $$']}
%labely(2)={['$$ Exact\,Isosurface $$']}
%labely(3)={['$$ K2(C_I)\, Reconstruction $$']}

% exact isosurface results for LSqR: normals
%labely(1)={['$$ Isopatch $$']}
%labely(2)={['$$ Lvl 1 - Ord 2$$']}
%labely(3)={['$$ Lvl 1 - Ord 4$$']}
%labely(4)={['$$ Lvl 2 - Ord 4$$']}
%labely(5)={['$$ Lvl 2 - Ord 6$$']}
%labely(6)={['$$ Lvl 3 - Ord 6$$']}
% exact isosurface results for LSqR: curv
%labely(1)={['$$ Lvl 1 - Ord 2$$']}
%labely(2)={['$$ Lvl 1 - Ord 4$$']}
%labely(3)={['$$ Lvl 2 - Ord 4$$']}
%labely(4)={['$$ Lvl 2 - Ord 6$$']}
%labely(5)={['$$ Lvl 3 - Ord 6$$']}

% comparisons smooth/no smooth for LSqR
% l2o2
%labely(1)={['$$ Exact $$']}
%labely(2)={['$$ C_I $$']}
%labely(3)={['$$ K1(C_I) $$']}
%labely(4)={['$$ K2(C_I) $$']}
%labely(1)={['$$ Exact $$']}
%labely(2)={['$$ Exact\,:\vec{n}=\vec{k} $$']}
%labely(3)={['$$ K2(C_I) $$']}
%labely(4)={['$$ K2(C_I)\,:\vec{n}=\vec{k} $$']}

% comparisons smooth/no smooth for vol methods
%labely(1)={['$$ L_{\max}:C_I$$']}
%labely(2)={['$$ L_{1}: C_I$$']}
%labely(3)={['$$ L_{\max}:La2(C_I) $$']}
%labely(4)={['$$ L_{1}:La2(C_I) $$']}
%labely(5)={['$$ L_{\max}:K2(C_I) $$']}
%labely(6)={['$$ L_{1}:K2(C_I) $$']}

% comparison smooth/no smooth for all
labely(1)={['$$ CDS-K2(C_I)$$']}
labely(2)={['$$ CA-K2(C_I)$$']}
labely(3)={['$$ LSqR-K2(C_I) $$']}

%----
%labely(1)={['$$ L_{\max} \left(',yname '\right) - Gauss - CDS $$']}
%labely(2)={['$$ L_1 \left(',yname, '\right) - Gauss - CDS $$']}
%labely(3)={['$$ L_{\max} \left(',yname '\right) - LSq3 - Lvl1$$']}
%labely(4)={['$$ L_1 \left(',yname, '\right) - LSq3 - Lvl1 $$']}
%labely(5)={['$$ L_{\max} \left(',yname '\right) - LSq3 - Lvl2 $$']}
%labely(6)={['$$ L_1 \left(',yname, '\right) - LSq3 - Lvl2 $$']}
%labely(1)={['$$ $$']}
%labely(2)={'$$ L_{\max} \left( A_{^s c} \right)$$''}
roty=0
%legloc='NorthWest'
legloc='NorthEast'




%----- Start Work
if (i_iv2==0)  
% variables to be used
if (iv_minus_1==1) 
    iv=s_iv./CD.data(i_use,i_iv);
else
    iv=s_iv*CD.data(i_use,i_iv);
end
m=length(s_dv)
for i=1:m
dv(:,i)=s_dv(i).*CD.data(i_use,i_dv(i));
end

% Plot in loglog
figure
set(gcf,'Color','w')
n_figures=3
at_figure=1
%subplot(n_figures,1,at_figure);
loglog(iv,dv(:,1),char(symbol(1)),'linewidth',2)
mfc=char(symbol(1))
%loglog(iv,dv(:,1),char(symbol(1)),'linewidth',2,'MarkerFaceColor',mfc(end))
hold
grid on
for i=2:length(s_dv)
%   subplot(n_figures,1,at_figure);
   loglog(iv,dv(:,i),char(symbol(i)),'linewidth',2)
mfc=char(symbol(i))
%plot(iv,dv(:,i),char(symbol(i)),'linewidth',2,'MarkerFaceColor',mfc(end))
%   if (mod(i,n_figures)==0)
%       at_figure=at_figure+1
%   end
end

if (default_lims==0)
% limits of x axis
x_min=min(iv)*(1-x_min_ext/10);
x_max=max(iv)*x_max_ext;
set(gca,'xlim',[x_min x_max]);
end

% Change ticks
if (revaxis==1)
    % reverse axis
    set(gca,'xdir','reverse')
end

    
if (~isempty(order_line_k))
% plot parameters for positions and colors for order lines
up=[1 2 0.11];
right=[1 1 1];
line_colors={'b--','g--','r--'};
for (nline=1:length(order_line_k))
if (order_line_k(nline)~=0) 
    labely(m+nline)={['$$ O(x^{',int2str(order_line_k(nline)),'}) $$']}
    %labely(m+1)={['$$ O(',xname,'^',int2str(order_line_k),') $$']}
    % add order line
    [max_dv,i_max_dv]=max(dv);
    yl=max_dv(1)/(iv(i_max_dv(1)).^order_line_k(nline))*iv.^order_line_k(nline);
    loglog(iv/right(nline),yl*up(nline),char(line_colors(nline)),'linewidth',2)
end
end
end

if (use_fits==1) 
    % fits
    s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1]);
    f = fittype('a*x^k','options',s);
    [Err_ivdv,gof2] = fit(iv,dv,f,s)
    plot(Err_ivdv)    
elseif (use_fits==2)
    f = fittype('poly2');
    [Err_ivdv,gof2] = fit(iv,dv,f)
    plot(Err_ivdv)    
elseif (use_fits==3)
    f = fittype('poly4');
    [Err_ivdv,gof2] = fit(iv,dv,f)
    plot(Err_ivdv)    
end


% Add annotations
set(gcf,'Units','Points');
Pos=get(gcf,'Position');
set(gcf,'Position',[Pos(1) Pos(2) 460 380]);
set(gca,'fontsize',14);
xlabel(labelx,'Interpreter','Latex','fontsize',16)
ylabel('$$ L_{\max}(\vec{n}) $$','Interpreter','Latex','fontsize',16,'rotation',roty)
%ylabel('$$ RE(\kappa) $$','Interpreter','Latex','fontsize',16,'rotation',roty)
%ylim([0.1 12])
if (length(order_line_k)~=0 && use_fits==1)
    legend(labely,'2nd Order')
elseif (length(order_line_k)~=0)
    legend('String',labely,'Interpreter','Latex','fontsize',16,'Location',legloc)
end

% add secondary axis
if (i_tv~=i_iv)
    ax1=gca;
    ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'Color','none','YTick',[],'xscale','log');
       linkaxes([ax1 ax2],'x')
       set(ax2,'XTick',sort(unique(iv)))
       set(ax2,'XTicklabel',sort(unique(CD.data(i_use,i_tv))))
end  

% find local n - power at y=a*x^n
ldv=length(dv(:,1))
n_local(1:ldv-1,1:m)=0
for i=1:m
    n_local(:,i)=(log(dv(2:ldv,i))-log(dv(1:ldv-1,i)))./(log(iv(2:ldv))-log(iv(1:ldv-1)));
end
n_local

else

% variables to be used
if (iv_minus_1==1) 
    iv=s_iv./CD.data(i_use,i_iv);
    iv2=s_iv2./CD.data(i_use,i_iv2);
else
    iv=s_iv*CD.data(i_use,i_iv);
    iv2=s_iv2*CD.data(i_use,i_iv2);
end
dv=s_dv*CD.data(i_use,i_dv);
    
% Plot in loglog
figure
set(gcf,'Color','w')
plot3(iv,iv2,dv,symbol)
set(gca,'xscale','log','yscale','log','zscale','log')
hold
grid on

if (use_fits==1) 
    % fits
    s=fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 1 1 1 ]);
    %f = fittype('a*x1^k+b*x2^l+c*x1^m*x2^n','independent',{'x1','x2'},'options',s);
    f = fittype('a*x1^k+b*x2^l','independent',{'x1','x2'},'options',s);
    [Err_ivdv,gof2] = fit([iv iv2],dv,f,s)
    plot(Err_ivdv)    
elseif (use_fits==2)
    f = fittype('poly22');
    [Err_ivdv,gof2] = fit([iv iv2],dv,f,s)
    plot(Err_ivdv)    
elseif (use_fits==3)
    f = fittype('poly33');
    [Err_ivdv,gof2] = fit([iv iv2],dv,f,s)
    plot(Err_ivdv)    
end

end
