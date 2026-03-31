function [B_d1] = gateNMDA_b_v(Vmd1)
global conc_exmg k_B
B_d1=1./(1+exp(-k_B*Vmd1)*conc_exmg/3.57);
end

