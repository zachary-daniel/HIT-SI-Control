clear;close all;clc;
load("plasma_load_model.mat")
options = optimset('MaxFunEvals', 50000, 'MaxIter', 50000);
Y1 = L2_Current(1:10000);
dT = 1e-7;
RunTime = 10000.*dT;
X = dT:dT:RunTime;

adapter_squared = @(vector) squared_errors(vector(1), vector(2), vector(3),X', Y1);
[opt_vals, fval] = fminsearch(adapter_squared,[15600,15600,1000],options);

plot(X,opt_vals(3).*(log(opt_vals(1).*X).*sin(opt_vals(2).*X) ))
hold on
plot(X,Y1)
%%
X2 = 0:dT:.004;
plot(X2, 20000.*(log(opt_vals(1).*X2).*sin(opt_vals(2).*X2) ) )




%Functions below this line
function [sum_squared_errors] = squared_errors(a, b,c, X, Y)
    f = @(x) c.*(log(a.*X).*sin(b.*X)); 
    sum_squared_errors = sum((Y-f(X)).^2);
