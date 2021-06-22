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

% The name of the file you store the data:
filename='~/O2-project/framework_O2FV/TESTS/13_RecEval/rec_errs_sin2.txt'

% Get data
CD=importdata(filename);
n=size(CD.data,1);
i_dv=[];
% Which column stores the dependant and the independant variable
% Note that idy is the second independant variable if you work want an
% error surface (if not it should be zero). i_skip is the number of
% variables you want to skip (default=1), i_start is the starting index
% (default-1) and i_end is the ending index (default=n)
i_iv=1
i_iv2=0
i_dv(1)=5
i_dv(2)=6
i_skip=1
i_start=1
i_end=n
i_use=[i_start:i_skip:i_end];
%i_use=[1 6 16 21 31 36 46] 

% The scale factors for the dependent variables if any ... 
% Note: must be 1 if not used this is useful when you have only stored
%       the number of discretizations   
s_iv=2
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
iv_minus_1=1%
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
iv_minus_1=1
revaxis=0
i_tv=i_iv
symbol(1)={'ob'}
symbol(2)={'sqr'}
use_fits=0
default_lims=1
x_min_ext=3
x_max_ext=x_min_ext
order_line_k=[1 2]

%---- You have to change these:
%xname='\Delta\theta'
%labelx=['$$', xname,'=\frac{2\pi}{n_\theta} =2\Delta\phi$$']
labelx=['$$ \frac{\delta}{L} $$']
yname='d_{^s n}'
labely(1)={['$$ L_{\max} \left(',yname '\right) $$']}
labely(2)={['$$ L_1 \left(',yname, '\right) $$']}
%labely(1)={['$$ $$']}
%labely(2)={'$$ L_{\max} \left( A_{^s c} \right)$$''}
roty=90

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
loglog(iv,dv(:,1),char(symbol(1)),'linewidth',2)
hold
grid on
for i=2:length(s_dv)
    plot(iv,dv(:,i),char(symbol(i)),'linewidth',2)
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
up=[3 0.7];
right=[1.1 1];
line_colors={'b--','r--'};
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
%ylabel(labely,'Interpreter','Latex','fontsize',16,'rotation',roty)

if (length(order_line_k)~=0 && use_fits==1)
    legend(labely,'2nd Order')
elseif (length(order_line_k)~=0)
    legend('String',labely,'Interpreter','Latex','fontsize',16)
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
n_local(1:length(dv)-1,1:m)=0
for i=1:m
    n_local(:,i)=(log(dv(2:length(dv),i))-log(dv(1:length(dv)-1,i)))./(log(iv(2:length(dv)))-log(iv(1:length(dv)-1)));
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
