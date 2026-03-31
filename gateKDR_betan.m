function [betan_d1] = gateKDR_betan(Vmd1)
betan_d1=0.00008*(Vmd1-42)./(exp((Vmd1-42)/10.0)-1.0);
end

