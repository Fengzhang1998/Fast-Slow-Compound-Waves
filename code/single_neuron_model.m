function dz =single_neuron_model(z1, z2, z3)

n=z1;

Vmd1=z2(1:n);Vms=z2(n+1:2*n);Vmd2=z2(2*n+1:3*n);
n_d1=z2(3*n+1:4*n);mcaH_d2=z2(4*n+1:5*n);hcaH_d2=z2(5*n+1:6*n);
Ca_in_d=z2(6*n+1:7*n);q_Ca_d=z2(7*n+1:8*n);

Ved1=z3(1:n);Ves=z3(n+1:2*n);Ved2=z3(2*n+1:3*n);
Id1_stim=z3(3*n+1:4*n);Is_stim=z3(4*n+1:5*n);Id2_stim=z3(5*n+1:6*n);

global Cms Cmd
global gms_Leak gmd2_Leak gmd1_Leak 
global gm_KDR_den gm_Ca_den 
global the_mcaH the_hcaH
global K_mcaH K_hcaH
global tao_mcaH tao_hcaH

global gm_KAHP_den gm_NMDA_equ
global v_Ca_in tao_Ca_in
global a_q_Ca tao_q_Ca

global gc_d1_s gc_s_d1 gc_s_d2 gc_d2_s
global E_Leak_d E_Leak_s E_NMDA E_K E_Ca
global Ecoeff_Imd1 Ecoeff_Imd1s Ecoeff_Imd2s Ecoeff_Imd2
global para_NMDA_1 para_NMDA_2
 

%% basal dendrite
Id1_Leak=-gmd1_Leak*(Vmd1-E_Leak_d);

Is_d1_from_intra=gc_s_d1*(Vms-Vmd1);
Is_d1_from_extra=Ecoeff_Imd1*(Ves-Ved1);
Is_d1=Is_d1_from_intra+Is_d1_from_extra;

B_d1=1./(1+exp(-para_NMDA_1*Vmd1)*para_NMDA_2);

Id1_NMDA=-gm_NMDA_equ*B_d1.*(Vmd1-E_NMDA);

Id1_KDR=-gm_KDR_den*(n_d1.^4).*(Vmd1-E_K);

Id1_tot=Id1_Leak+Is_d1+Id1_NMDA+Id1_KDR+Id1_stim;

%% soma
Is_Leak=-gms_Leak*(Vms-E_Leak_s);

Id1_s_from_intra=gc_d1_s*(Vmd1-Vms);
Id1_s_from_extra=Ecoeff_Imd1s*(Ved1-Ves);
Id1_s=Id1_s_from_intra+Id1_s_from_extra;

Id2_s_from_intra=gc_d2_s*(Vmd2-Vms);
Id2_s_from_extra=Ecoeff_Imd2s*(Ved2-Ves);
Id2_s=Id2_s_from_intra+Id2_s_from_extra;

Is_tot=Is_Leak+Id1_s+Id2_s+Is_stim;

%% apical dendrite
Id2_Leak=-gmd2_Leak*(Vmd2-E_Leak_d);

Is_d2_from_intra=gc_s_d2*(Vms-Vmd2);
Is_d2_from_extra=Ecoeff_Imd2*(Ves-Ved2);
Is_d2=Is_d2_from_intra+Is_d2_from_extra;

Id2_Ca=-gm_Ca_den*mcaH_d2.^2.*hcaH_d2.*(Vmd2-E_Ca);

Id2_KAHP=-gm_KAHP_den*q_Ca_d.*(Vmd2-E_K);

Id2_tot=Id2_Leak+Is_d2+Id2_Ca+Id2_KAHP+Id2_stim;

%% 
alphn_d1=gateKDR_aplhn(Vmd1);
betan_d1=gateKDR_betan(Vmd1);

ninf_d1=alphn_d1./(alphn_d1+betan_d1); 
tao_n_d1=1./(alphn_d1+betan_d1);

mcaH_inf_d2=1./(1+exp(-(Vmd2-the_mcaH)/K_mcaH));
hcaH_inf_d2=1./(1+exp(-(Vmd2-the_hcaH)/K_hcaH));

Cainf_in_d=v_Ca_in*Id2_Ca*tao_Ca_in;            
qinf_Ca_d=1./(1+(a_q_Ca./Ca_in_d).^4);          
qinf_Ca_d1=1./(1+(a_q_Ca./Cainf_in_d).^4);

%% 
dVmd1=Id1_tot/Cmd;
dVms=Is_tot/Cms;
dVmd2=Id2_tot/Cmd;

dn_d1=(ninf_d1-n_d1)./(tao_n_d1);
dmcaH_d2=(mcaH_inf_d2-mcaH_d2)./tao_mcaH;
dhcaH_d2=(hcaH_inf_d2-hcaH_d2)./tao_hcaH;

dCa_in_d=v_Ca_in*Id2_Ca-Ca_in_d./tao_Ca_in; 
dq_Ca_d=(qinf_Ca_d-q_Ca_d)./tao_q_Ca;

dz=[dVmd1;dVms;dVmd2;
    dn_d1;dmcaH_d2;dhcaH_d2;
    dCa_in_d;dq_Ca_d;
    ninf_d1;mcaH_inf_d2;hcaH_inf_d2;
    Cainf_in_d;qinf_Ca_d1;];

end

