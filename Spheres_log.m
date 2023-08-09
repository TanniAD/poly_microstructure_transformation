  %%    Generate Spheres
clc
clear 
close all

n = 510;            % Number of Spheres
v = 1;              % Variance
m = 4;              % Mean

mu = log((m^2)/sqrt(v+m^2));            
sigma = sqrt(log(v/(m^2)+1)); 
R = sort(lognrnd(mu,sigma,[n,1]),'descend');
V = [100, 100, 100];         % Dimension of the Volume
  
M = GenerateRandomSpheres(n, V, R);
limit=V(1);                                         %Limit of the meshgrid
n2 = 500;                                           %Number of meshgrid
coord = linspace(0, limit, n2);
[X, Y, Z] = meshgrid(coord, coord, coord);          %The BOX
    %%   Plot Spheres 
% figure(1);
% axes('NextPlot', 'add', ...
%     'XLim', [0, V(1)], 'YLim', [0, V(2)], 'ZLim', [0, V(3)]);
% view(3); [x, y, z] = sphere(); axis('square');
%     
% for k = 1:n
%     surf(x * R(k) + M(k, 1), y * R(k) + M(k, 2), z * R(k) + M(k, 3));
% end
% xlabel('x-axis')
% ylabel('y-axis')
% zlabel('z-axis')
% shading interp; camlight right; lighting phong; view(60,30);

    %%   Condition to be inside sphere
C = zeros(n2,n2,n2);
for m = 1:n(1)
in_sphere = (X - M(m,1)).^2 + (Y - M(m,2)).^2 + (Z - M(m,3)).^2 <= R(m,1).^2;
C = double(in_sphere)+C;
end

    %%   Slice
xslice = [];   yslice = [];
zslice = 10:2*max(R):V(3); 
radii = slicer(X,Y,Z,C,coord, xslice, yslice, zslice);    % currently slicer is for zslice only
edges = 0:0.2:max(radii);                                     % bin = 0.5 (generalized matrix)
% histogram(R,edges);
% hold on
% histogram(radii,edges);
% legend('Sphere','Circles')
% hold off
    %% Fitting
    figure(1)
    par = lognfit(radii);                       % Parameters of radii lognormal fit
hf = histfit(radii,60,'lognormal');
hold on
xx = linspace(0,max(radii),numel(radii))';
yy = lognpdf(xx,par(1),par(2));
fun = @(b,radii)lognpdf(radii,b(1),b(2));
b = nlinfit(xx,yy,fun,[par(1),par(2)]);
fb = (edges(2)-edges(1))*n*fun(b,xx);
% histogram(radii,edges)
plot(xx,fb,'r-')
legend ('radii','histfit','nlinfit')
    %%   Integral
   s = 0:0.1:xx(end);
e=1;
for ii=1:numel(s)+1
    f(e)=(xx(ii+1)-xx(ii))*fb(ii);
    e=e+1;
end

final = GeneralizedMatrix(numel(xx),0.5);
P = (final)*f;      
F = round(1000*P/sum(P));
F(F<0) = 0; 

bincntrs = (xx(1:end-1) + xx(2:end))/2;    % x coordinate for F plot
figure(2)
plot(bincntrs, F,'LineWidth',2);          
hold on
histogram(R,xx, 'facecolor','r')
legend ('F','R')

N = sum(F);
mean_F = sum(bincntrs.*F)/N;
sd_F = sqrt(((N*(sum(F.*bincntrs.^2)))-(sum(F.*bincntrs))^2)/(N*(N-1)));
