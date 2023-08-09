clc
clear
close all
load('lgnSphereMu_3.2762_sigma_0.5_n_680_V_1000.mat','R','radii','n','mu','sigma');
%% fitting
edges = 0:0.2:max(radii); 
par_r = lognfit(radii);                            
% r = lognrnd(par_r(1),par_r(2),[1,n]);
% hist_data = histfit(radii,25,'lognormal');
% b = hist_data(2).XData; radii data
% p = hist_data(2).YData; probability data
mu_r = par_r(1); 
sigma_r = par_r(2);

%% nonlinear fit of radii
par = lognfit(radii);
xx = linspace(0,max(radii),numel(radii))';
yy = lognpdf(xx,par(1),par(2));
fun = @(b,radii)lognpdf(radii,b(1),b(2));
b = nlinfit(xx,yy,fun,[par(1),par(2)]);
fb = n*(xx(2)-xx(1))*fun(b,xx);
histogram(radii,xx)
hold on
plot(xx,fb,'r-');

%% Integral
    syms x s t;
phi = @(x) (((x*sigma_r*sqrt(2*pi)).^(-1)).*exp(-(log(x)-mu_r).^2/(2*sigma_r.^2))); 
% pd = phi(r);
    fi = int(phi, x, s, t);         
    s = edges(1:end-1);
    t = edges(2:end);
    f = (double(subs(fi)))';
    f(isnan(f))= 0;
    
    %% Conversion
final = GeneralizedMatrix(length(edges),0.5);
    P = final*f;      
    F = round(1000*P/sum(P));
    F(F<0) = 0; 
bincntrs = (edges(1:end-1) + edges(2:end))/2;
% histogram (R,500)
%     hold on
%     plot(bincntrs, F,'linewidth',2)
    
     %% Parameters
N = sum(F);
mean_F = sum(bincntrs'.*F)/N;
sd_F = sqrt(((N*(sum(F.*bincntrs'.^2)))-(sum(F.*bincntrs'))^2)/(N*(N-1)));
v_F = sd_F^2;
mu_F = log((mean_F^2)/sqrt(v_F+mean_F^2));            % parameter 1
sigma_F = sqrt(log(v_F/(mean_F^2)+1));                % parameter 2
%% Plotting
