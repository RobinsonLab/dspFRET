% E is un-gamma-corrected efficiency
% TE is from lifetimes
% TE = 1 - (\tau_{DA}/\tau_{DO}

E = 0.7;
tau_DA = 0.11;
tau_DO = 0.20;
TE = 1- tau_DA/tau_DO;
gamma = ((E/TE)-E)/(1-E);
