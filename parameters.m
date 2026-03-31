%%%%%%%%%%%%three-compartmental HPC paramter sets%%%%%%%%%%%%%%
%d1: basal dendritic compartment including NMDAR channels d1是基树突
%d2: apical dendritic compartment including Ca2+ channels d2是顶树突
%s: soma compartment is the passive one
%-------global parameters
global nN nN_z
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
global conc_exmg k_B
global Ecoeff_Imd1 Ecoeff_Imd1s Ecoeff_Imd2s Ecoeff_Imd2
global ts te dt t nt dt05
global para_NMDA_1 para_NMDA_2
 

ts=0;
te=5000;
dt=0.05;
dt05=dt/2;

t=ts:dt:te; %units: ms
nt=length(t);

nN=200;nN_z=18;
SF=20; 

k1=10;     
k2=10;     
k3=15;    


dis_ext=1e-4*2;%units:cm 
Ds=1e-4*17.8;Ls=Ds;%units:cm 

Cms=1;Cmd=1;%units:uF/cm^2
Rms=680;Rmd2=680;Rmd1=780;%units: oumiga*cm^2
Ri=530;
gms_Leak=1/Rms*1e+3;
gmd2_Leak=1/Rmd2*1e+3;
gmd1_Leak=1/Rmd1*1e+3;

Dd=1e-4*17.8;Ld=Dd;%units:cm
D_d1s=1e-4*5.2;%units:cm
D_d2s=1e-4*5.2;%units:cm

L_d1s=1e-4*250;%units:cm
L_d2s=1e-4*500;%units:cm
Lds=[L_d1s;L_d2s];

As=pi*Ds*Ls;%units:cm^2
Ad1=15*As;%units:cm^2
Ad2=9.5*As;%units:cm^2

Rc_d1s=1e-6*4*Ri*L_d1s/(pi*D_d1s^2);
Rc_d2s=1e-6*4*Ri*L_d2s/(pi*D_d2s^2);
gc_d1s=1e-3/Rc_d1s;
gc_d2s=1e-3/Rc_d2s;
gc_d1_s=gc_d1s/As;
gc_s_d1=gc_d1s/Ad1;
gc_d2_s=gc_d2s/As;
gc_s_d2=gc_d2s/Ad2;

Re=300;

E_Leak_d=-58;
E_Leak_s=-58;
E_NMDA=0;
E_K=-60;
E_Ca=10;
gm_KDR_den=200;
gm_Ca_den=42; 
gm_KAHP_den=0.03; 
para_NMDA_1=0.0744; 
para_NMDA_2=0.2801;  
gm_NMDA_equ=6.5;
the_mcaH=-20;the_hcaH=-40;the_c_Ca=-40;
K_mcaH=16;K_hcaH=-7;K_c_Ca=7;
tao_mcaH=10;
tao_hcaH=70;
tao_c_Ca=20; 
v_Ca_in=0.13;tao_Ca_in=13;
a_Ca_in=6;
a_q_Ca=2;tao_q_Ca=300;
conc_exmg=1.0;
k_B=0.0744;

Ecoeff_Imd2=k1*gc_d2s/Ad2*k2;
Ecoeff_Imd2s=k1*gc_d2s/As;
Ecoeff_Imd1s=k1*gc_d1s/As;
Ecoeff_Imd1=k1*gc_d1s/Ad1*k3;



