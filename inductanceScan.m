clear;close all;clc;

L1init = (8.0141e-7); %H;
L2init = 2.0462e-6; %H
L1 = L1init;
Lt = (L1init*L2init)/(L1init+L2init);
scalerL2 =@(scalerL1) (scalerL1*Lt*L1init)/ (L2init* (scalerL1*L1init - Lt )  ); % This is how to change L2 as a function of L1 to maintain resonance;

L1scalers = linspace(1,2.5,25);
maxCurrents = zeros(size(L1scalers));
L1vals = zeros(size(maxCurrents));
L2vals = zeros(size(maxCurrents));
for i = 1:length(L1vals)
    L2 = L2init*scalerL2(L1scalers(i));
    L1 = L1init*L1scalers(i);
    L1vals(i) = L1;
    L2vals(i) = L2;
    assignin('base','L1',L1);
    assignin('base', 'L2',L2);
    run('inductanceScanDummy.m');
    maxCurrents(i) = ans.L2CurrentMax.signals.values(end);
end
plot(L1vals, maxCurrents);
[M,I] = max(maxCurrents);
print("The best value of L1 is " + L1vals(I) + " and L2 is " + L2vals(I));