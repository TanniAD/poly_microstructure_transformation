clc, clear
B = h5read('mean 50_350.dream3d','/DataContainers/SyntheticVolumeDataContainer/CellData/FeatureIds');
A = double(squeeze(B(1,:,:,:)));

[m,n]= count_unique(A);
D = nthroot(n,3);
h = 1:99:length(A);
M = cell(length(h),5); 
for i = 1:length(h)
    M{i,1} = slice(A,[],[],h);
    M{i,2} = M{i,1}.CData;
    [M{i,3},M{i,4}] = count_unique(M{1,2}); 
    M{i,5} = sqrt((4/pi)*M{i,4});
end
    
d2d_multislices_fullvolume = cat(1,M{:,5});
save('mean 50_350.mat','A','D','d2d_multislices_fullvolume','-v7.3');
%% Multiple Slices Full Volume
clc,clear
load('Set1_sigma_1_1000.mat','B','Q');
D = Q{2,3};
h = 1:100:length(B);
for i = 1:length(h)
    m{1,i} = slice(B,[],[],h(i));
   m{2,i} = m{1,i}.CData; % Color Data
   m{3,i} = m{2,i}(~isnan(m{2,i}));
    [m{4,i},m{5,i}] = count_unique(m{3,i}); % {i,4} unique values, W{i,5} number of uniques
    x{1,i} = sqrt((4/pi)*m{5,i});   %  EquivDia
    [x{2,i}, x{3,i}, x{4,i}, x{5,i}] = numdist(x{1,i});  % 2D Number dist
[x{6,i}, x{7,i}, x{8,i}, x{9,i}] = areadist(x{1,i});  % 2D Area dist
[x{10,i}, x{11,i}, x{12,i}, x{13,i}] = numdist(D);      % 3D Number dist
[x{14,i}, x{15,i}, x{16,i}, x{17,i}] = grainz(x{1,i}, x{2,i},x{3,i}); % NumDist Transformed
[x{18,i}, x{19,i}, x{20,i}, x{21,i}] = areadist(x{14,i});  % Area dist Transformed
end
% d2d_multislices_fullvolume = cat(1,x{:,6});
% save('Set6_sigma_5_1000.mat','x','-append')
for j = 1:size(x,2)
    ngof(j,1) = x{5,j};
    ngrain(j,1) = numel(x{1,j});
    agof(j,1) = x{9,j};
end
figure(1),scatter(ngrain,ngof,50,'p','MarkerEdgeColor',colors('carmine red'),'LineWidth',1.5), hold on
scatter(ngrain,agof,50,'p','MarkerEdgeColor',colors('blue'),'LineWidth',1.5)

%% 1 slice Full Volume
clc, clear
load('Set6_mean_30_mu_3.1199_sigma_0.75_1000.mat')
h = Dimension(1)/2;
M = cell(length(h),5); 
    M{1,1} = slice(B,[],[],h);
    M{1,2} = M{1,1}.CData;
    [M{1,3},M{1,4}] = count_unique(M{1,2}); % {i,3} unique values, W{i,4} number of uniques
    M{1,5} = sqrt((4/pi)*M{1,4});
    
 d1 = cell2mat(M(1,5));
 save('Set6_mean_30_mu_3.1199_sigma_0.75_700.mat','d1','-append')
 
%% 1 slice of sub volume
clc, clear
load('Set6_mean_30_mu_3.1199_sigma_0.75_700.mat','B','Dimension')
limits = [0.25*Dimension(1)+1,0.75*Dimension(1),0.25*Dimension(1)+1,0.75*Dimension(1),0.25*Dimension(1)+1,0.75*Dimension(1)];
v = subvolume(B,limits);
[m, n] = count_unique(v);
Ds = nthroot((6*n/pi),3);

h = Dimension(1)/2;
w = cell(length(h),5); 
    w{1,1} = slice(v,[],[],h);
    w{1,2} = w{1,1}.CData;
    [w{1,3},w{1,4}] = count_unique(w{1,2}); % {i,3} unique values, W{i,4} number of uniques
    w{1,5} = sqrt((4/pi)*w{1,4});
ds = cell2mat(w(1,5));

save('Set6_mean_30_mu_3.1199_sigma_0.75_700.mat','v','Ds','ds','-append')

%% Multitple slices of sub volume
clc, clear
load('Set6_mean_30_mu_3.1199_sigma_0.75_700.mat','v')

hs = 1:99:length(v);
w = cell(length(hs),5); 
for i = 1:length(hs)
    w{i,1} = slice(v,[],[],hs(i));
    w{i,2} = w{i,1}.CData;
    [w{i,3},w{i,4}] = count_unique(w{i,2}); % {i,3} unique values, W{i,4} number of uniques
    w{i,5} = sqrt((4/pi)*w{i,4});
end
dsm = cat(1,w{:,5});
save('Set6_mean_30_mu_3.1199_sigma_0.75_700.mat','dsm','w','-append')
%% Boundary Reduction
clc, clear
load('Set6_mean_30_mu_2.9012_sigma_1.0_600.mat','B')
[A,D] = BoundaryReduction(B);
save('Set6_mean_30_mu_2.9012_sigma_1.0_600.mat','A','D','-append')
%% Edge Reduction
clc, clear
tic
load('Set1_mean_30_mu_2.9012_sigma_1.0_900.mat','B')
[c,dr,d_astm] = EdgeReduction(B);
toc
save('Set1_mean_30_mu_2.9012_sigma_1.0_900.mat','c','dr','d_astm','-append')

%% Slice with 1% area/volume
clc,clear
tic
load('Set6_mean_30_mu_2.9012_sigma_1.0_900.mat','A')
h = 100:100:length(A);
M = cell(length(h),6); 
for i = 1:length(h)
    M{i,1} = slice(A,[],[],h(i));
    M{i,2} = get(M{i,1},'CData');
    M{i,3} = M{i,2}(~isnan(M{i,2}));
    [M{i,4},M{i,5}] = count_unique(M{i,3}); % {i,3} unique values, W{i,4} number of uniques
    M{i,6} = sqrt((4/pi)*M{i,5});
end
 d100 = cat(1,M{:,6});
 toc
save('Set6_mean_30_mu_2.9012_sigma_1.0_900.mat','d100','M','-append')
%% Individual cropped slices
clc,clear
% load('Set1sigma_25_1000.mat','B')
B = double(h5read('Set6_mean_30_mu_3.2762_sigma_0.5_1000.dream3d','/DataContainers/SyntheticVolumeDataContainer/CellData/FeatureIds'));
B = double(squeeze(B(1,:,:,:)));
m = slice(B,[],[],randi(length(B),1)); c = m.CData;
Q = cell(30,length(B)/50);
for i = 1:length(B)/50
 Q{1,i} = c(1*i:50*i,1*i:50*i);             % croppedArea
[u,v] = count_unique(B);
Q{2,i} = nthroot(6*v/pi,3);                 % 3D Dia
Q{3,i} = numel(u);                          % 3D grain Number
[p,q] = count_unique(Q{1,i});               
Q{4,i} = nthroot(4*q/pi,2);                 % 2D Dia
Q{5,i} = numel(p);                          % 2D Grain Number
[Q{6,i}, Q{7,i}, Q{8,i}] = numdist(Q{2,i});  % 3D Grain NumDist(mu,sigma,L)
[Q{9,i}, Q{10,i},Q{11,i}] = voldist(Q{2,i});  % 3D Grain VolDist
[Q{12,i}, Q{13,i}, Q{14,i}] = numdist(Q{4,i});     % 2D Grain NumDist
[Q{15,i}, Q{16,i}, Q{17,i}, Q{18,i}] = grainz(Q{4,i}, Q{12,i}, Q{13,i}); % NumDist
[Q{19,i}, Q{20,i}, Q{21,i}] = voldist(Q{15,i});     % Transformed Grain VolDist
Q{22,i} = 100*abs(Q{16,i}-Q{6,i})/Q{6,i}; % mu Error
Q{23,i} = 100*abs(Q{17,i}-Q{7,i})/Q{7,i}; % sigma Error
Q{24,i} = 100*abs(Q{18,i}-Q{8,i})/Q{8,i}; % L Error
Q{25,i} = 100*abs(Q{19,i}-Q{9,i})/Q{9,i}; % mu Error
Q{26,i} = 100*abs(Q{20,i}-Q{10,i})/Q{10,i}; % sigma Error
Q{27,i} = 100*abs(Q{21,i}-Q{11,i})/Q{11,i}; % L Error
[Q{28,i}, Q{29,i}, Q{30,i}] = voldist(Q{4,i});     % 2D Grain VolDist
end
save('Set6_sigma_5_1000.mat','B','Q','-append')
%%
clc,clear
load('Set6_sigma_1_1000.mat','B')
limits = [63, 937, 63, 937, 63,937]; b = subvolume(B,limits);
s = slice(b,[],[],randi(length(b),1)); c = s.CData;
k = 35; q = cell(30,length(b)/k);
for i = 1:length(b)/k
 q{1,i} = c(i:k*i, i:k*i);             % croppedArea
[u,v] = count_unique(b);
q{2,i} = nthroot(6*v/pi,3);                 % 3D Dia
q{3,i} = numel(u);                          % 3D grain Number
[m,n] = count_unique(q{1,i});               
q{4,i} = nthroot(4*n/pi,2);                 % 2D Dia
q{5,i} = numel(m);                          % 2D Grain Number
[q{6,i}, q{7,i}, q{8,i}] = numdist(q{2,i});  % 3D Grain NumDist(mu,sigma,L)
[q{9,i}, q{10,i},q{11,i}] = voldist(q{2,i});  % 3D Grain VolDist
[q{12,i}, q{13,i}, q{14,i}] = numdist(q{4,i});     % 2D Grain NumDist
[q{15,i}, q{16,i}, q{17,i}, q{18,i}] = grainz(q{4,i}, q{12,i}, q{13,i}); % NumDist
[q{19,i}, q{20,i}, q{21,i}] = voldist(q{15,i});     % Transformed Grain VolDist
q{22,i} = 100*abs(q{16,i}-q{6,i})/q{6,i}; % mu Error
q{23,i} = 100*abs(q{17,i}-q{7,i})/q{7,i}; % sigma Error
q{24,i} = 100*abs(q{18,i}-q{8,i})/q{8,i}; % L Error
q{25,i} = 100*abs(q{19,i}-q{9,i})/q{9,i}; % mu Error
q{26,i} = 100*abs(q{20,i}-q{10,i})/q{10,i}; % sigma Error
q{27,i} = 100*abs(q{21,i}-q{11,i})/q{11,i}; % L Error
[q{28,i}, q{29,i}, q{30,i}] = voldist(q{4,i});     % 2D Grain VolDist
end
save('Set6_sigma_1_1000.mat','b','q','-append')
%% MSE
clc,clear
load('Set1_sigma_25_1000.mat','d2d_multislices_fullvolume','Q');
k = 100:100:numel(d2d_multislices_fullvolume); D = Q{2,3};
for i = 1:length(k)
y{1,i} = randsample(d2d_multislices_fullvolume,k(i)); % random 2D diameter selection 
[y{2,i}, y{3,i}, y{4,i}, y{5,i}] = numdist(y{1,i});  % 2D Number dist
[y{6,i}, y{7,i}, y{8,i}, y{9,i}] = areadist(y{1,i});  % 2D Area dist
[y{10,i}, y{11,i}, y{12,i}, y{13,i}] = numdist(D);      % 3D Number dist
[y{14,i}, y{15,i}, y{16,i}, y{17,i}] = grainz(y{1,i}, y{2,i}, y{3,i}); % NumDist Transformed
[y{18,i}, y{19,i}, y{20,i}, y{21,i}] = areadist(y{14,i});  % Area dist Transformed
y{22,i} = HellD(y{6,i},y{7,i},y{18,i},y{19,i});             % Distance in area dist
y{23,i} = HellD(y{2,i},y{3,i},y{15,i},y{16,i});             % Distance in number dist
end
save('Set1_sigma_25_1000.mat','y','-append')
%% SSE
clc,clear
load('Set6_sigma_25_1000.mat','d2d_multislices_fullvolume','Q');
k = 100:100:numel(d2d_multislices_fullvolume); D = Q{2,3};
for i = 1:length(k)
z{1,i} = randsample(d2d_multislices_fullvolume,k(i)); % random 2D diameter selection 
[z{2,i}, z{3,i}, z{4,i}, z{5,i}] = numdist(z{1,i});  % 2D Number dist
[z{6,i}, z{7,i}, z{8,i}, z{9,i}] = areadist(z{1,i});  % 2D Area dist
[z{10,i}, z{11,i}, z{12,i}, z{13,i}] = numdist(D);      % 3D Number dist
[z{14,i}, z{15,i}, z{16,i}, z{17,i}] = grainz(z{1,i}, z{2,i}, z{3,i}); % NumDist Transformed
[z{18,i}, z{19,i}, z{20,i}, z{21,i}] = areadist(z{14,i});  % Area dist Transformed
z{22,i} = HellD(z{6,i},z{7,i},z{18,i},z{19,i});             % Distance in area dist
z{23,i} = HellD(z{2,i},z{3,i},z{15,i},z{16,i});             % Distance in number dist
end
save('Set6_sigma_25_1000.mat','z','-append')
%%
clc, clear
load('Set5_sigma_25_1000.mat','B')
surf(squeeze(B(:,:,500)))
view(0,90), shading interp, colormap jet,

C = cell(6,2);
s = [1,1000]; 

i = 5; j = 6;
S1 = slice(B,[],[],s(1));
C{i,1} = S1.CData;
[m,n] = count_unique(C{i,1});
C{i,2}= sqrt((4/pi)*n);

S2 = slice(B,[],[],s(2));
C{j,1} = S2.CData;
[p,q] = count_unique(C{j,1});
C{j,2}= sqrt((4/pi)*q);

data(2,:) = C{2,2};

save('demo_25.mat','C')