clc, clear, close all
% load('Num(MSlice).mat','L_Dn','L_Gn'), load('Vol(MSlice).mat','L_Dv','L_Gv'), load('Area(MSlice).mat','L_Da','L_Ga')
load('mean 50_350.mat','eta_r'), 
% D = Q{2,3}; Dx = datastats(D); v = var(D);
r = 1:500; eta = eta_r(r);
tbl = table(r', eta');
plot(r, eta, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);

%% Exponential
ModelFunction_exp = @(m, x) exp(-2*x./m(1)); 
m0 = mean(r);               % initial assumption
mdl_exp = fitnlm(tbl, ModelFunction_exp, m0);
coeff_exp = mdl_exp.Coefficients;
xFit1 = linspace(0,500,max(r));
yFit1 = ModelFunction_exp(coeff_exp{1,1}, xFit1);
hold on;
plot(xFit1, yFit1, 'r-', 'LineWidth', 2);
legend('Discrete $\eta$','Fitted $e^{-r/L}$','Interpreter','Latex')

