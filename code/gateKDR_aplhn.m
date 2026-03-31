function [alphn_d1] = gateKDR_aplhn(Vmd1)
alphn_d1=0.00049*(Vmd1-32)./(1.0-exp(-(Vmd1-32)/25.0));
end

