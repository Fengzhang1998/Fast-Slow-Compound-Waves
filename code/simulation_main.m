clear;
close all;

run parameters.m 

s_c=1;
cn_matrix=generate_connection_matrix(s_c,nN);

Vmd1_rest=-60;
Vms_rest=-60;
Vmd2_rest=-60;

k1=single_neuron_model(1,[Vmd1_rest;Vms_rest;Vmd2_rest;0;0;0;0;0;0;0],[0;0;0;0;0;0]);
n_d1_rest=k1(9);mcaH_d2_rest=k1(10);hcaH_d2_rest=k1(11);

k1=single_neuron_model(1,[Vmd1_rest;Vms_rest;Vmd2_rest;0;mcaH_d2_rest;hcaH_d2_rest;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0]);
Ca_in_d_rest=k1(12);q_Ca_d_rest=k1(13);

X_rest=[Vmd1_rest;Vms_rest;Vmd2_rest;n_d1_rest;mcaH_d2_rest;hcaH_d2_rest;
        Ca_in_d_rest;q_Ca_d_rest];

k1=single_neuron_model(1,X_rest,[0;0;0;0;0;0]);

Id1_intra_rest=-k1(1)*Cmd;
Is_intra_rest=-k1(2)*Cms;
Id2_intra_rest=-k1(3)*Cmd; 

%% 
cor_cell=generate_cor_xyz(nN,nN_z,Lds,Ds,dis_ext);
cor_soma_x=cor_cell(1:nN_z,:);

cor_soma_y=cor_cell(nN_z+1:2*nN_z,:);
cor_soma_z=cor_cell(2*nN_z+1:3*nN_z,:);
cor_den1_x=cor_cell(3*nN_z+1:4*nN_z,:);
cor_den1_y=cor_cell(4*nN_z+1:5*nN_z,:);
cor_den1_z=cor_cell(5*nN_z+1:6*nN_z,:);
cor_den2_x=cor_cell(6*nN_z+1:7*nN_z,:);
cor_den2_y=cor_cell(7*nN_z+1:8*nN_z,:);
cor_den2_z=cor_cell(8*nN_z+1:9*nN_z,:);

wg_matrix=generate_field_coupling(nN,nN_z,cor_cell,Re);
wg_matrix_s_s=wg_matrix(1:nN,:);
wg_matrix_d2_s=wg_matrix(nN+1:2*nN,:);
wg_matrix_d1_s=wg_matrix(2*nN+1:3*nN,:);

wg_matrix_s_d2=wg_matrix(3*nN+1:4*nN,:);
wg_matrix_d2_d2=wg_matrix(4*nN+1:5*nN,:);
wg_matrix_d1_d2=wg_matrix(5*nN+1:6*nN,:);

wg_matrix_s_d1=wg_matrix(6*nN+1:7*nN,:);
wg_matrix_d2_d1=wg_matrix(7*nN+1:8*nN,:);
wg_matrix_d1_d1=wg_matrix(8*nN+1:9*nN,:); 

%% 
X1=zeros(8*nN,nt);

Vm_rest_var=10;
n_d1_rest_var=0.1;mcaH_d2_rest_var=0.1;hcaH_d2_rest_var=0.1;
Ca_in_d_rest_var=0;q_Ca_d_rest_var=0.1;


X1(1:nN,1)=Vmd1_rest-Vm_rest_var/2+Vm_rest_var.*rand(nN,1);
X1(nN+1:2*nN,1)=Vms_rest-Vm_rest_var/2+Vm_rest_var.*rand(nN,1);
X1(2*nN+1:3*nN,1)=Vmd2_rest-Vm_rest_var/2+Vm_rest_var.*rand(nN,1);
X1(3*nN+1:4*nN,1)=n_d1_rest-n_d1_rest_var/2+n_d1_rest_var.*rand(nN,1);
X1(4*nN+1:5*nN,1)=mcaH_d2_rest-mcaH_d2_rest_var/2+mcaH_d2_rest_var.*rand(nN,1);
X1(5*nN+1:6*nN,1)=hcaH_d2_rest-hcaH_d2_rest_var/2+hcaH_d2_rest_var.*rand(nN,1);
X1(6*nN+1:7*nN,1)=Ca_in_d_rest-Ca_in_d_rest_var/2+Ca_in_d_rest_var.*rand(nN,1);
X1(7*nN+1:8*nN,1)=q_Ca_d_rest-q_Ca_d_rest_var/2+q_Ca_d_rest_var.*rand(nN,1);

Stim=zeros(6*nN,nt);

Stim(3*nN+1:4*nN,:)=Id1_intra_rest-0.0;
Stim(4*nN+1:5*nN,:)=Is_intra_rest-0.0;
Stim(5*nN+1:6*nN,:)=Id2_intra_rest-0.0;

%% Id2_in
Ttrig=1000;%units:ms
Ts=1000;%units:ms                     
fs=1/Ts;%units:1/ms
period_pulse=30;%units: ms                      
duty=100*period_pulse/Ts;          
Id2_in_pul_max=15;%units: uA/cm^2                
Id2_in_pul_min=0;%units: uA/cm^2
Id2_in=(Id2_in_pul_max-Id2_in_pul_min)*(1+square(2*pi*fs*(t-Ttrig),duty))/2;%units: mV/mm 
Id2_in=Id2_in_pul_min+Id2_in.*(sign(t-Ttrig)+1)/2;  

Stim(5*nN+20:5*nN+25,:)=Stim(5*nN+20:5*nN+25,:)+ones(6,1)*Id2_in;

%% 
V_ext_d1=zeros(nN,nt);
V_ext_s=zeros(nN,nt);
V_ext_d2=zeros(nN,nt);


%% simulation start
for i=1:nt

    tp=(i-1)*dt
    Vmd1=X1(1:nN,i);Vms=X1(nN+1:2*nN,i);Vmd2=X1(2*nN+1:3*nN,i);   

    Id1_s_sum=gc_s_d1*cn_matrix.*wg_matrix_d1_s*(Vms-Vmd1);  
    Is_s_sum=cn_matrix.*wg_matrix_s_s*(gc_d1_s*(Vmd1-Vms)+gc_d2_s*(Vmd2-Vms));
    Id2_s_sum=gc_s_d2*cn_matrix.*wg_matrix_d2_s*(Vms-Vmd2);
    Id1_d1_sum=gc_s_d1*cn_matrix.*wg_matrix_d1_d1*(Vms-Vmd1);
    Is_d1_sum=cn_matrix.*wg_matrix_s_d1*(gc_d1_s*(Vmd1-Vms)+gc_d2_s*(Vmd2-Vms));
    Id2_d1_sum=gc_s_d2*cn_matrix.*wg_matrix_d2_d1*(Vms-Vmd2);
    Id1_d2_sum=gc_s_d1*cn_matrix.*wg_matrix_d1_d2*(Vms-Vmd1);
    Is_d2_sum=cn_matrix.*wg_matrix_s_d2*(gc_d1_s*(Vmd1-Vms)+gc_d2_s*(Vmd2-Vms));
    Id2_d2_sum=gc_s_d2*cn_matrix.*wg_matrix_d2_d2*(Vms-Vmd2);
    
    Ved1=1*SF*1e-3/(4*pi)*(Is_d1_sum*As+Id2_d1_sum*Ad2+Id1_d1_sum*Ad1);
    Ved2=1*SF*1e-3/(4*pi)*(Is_d2_sum*As+Id2_d2_sum*Ad2+Id1_d2_sum*Ad1);
    Ves=1*SF*1e-3/(4*pi)*(Is_s_sum*As+Id2_s_sum*Ad2+Id1_s_sum*Ad1);
    Vstim_d1=Stim(1:nN,i)+Ved1;
    Vstim_s=Stim(nN+1:2*nN,i)+Ves;
    Vstim_d2=Stim(2*nN+1:3*nN,i)+Ved2;
    V_ext_d1(:,i)=Ved1;
    V_ext_s(:,i)=Ves;
    V_ext_d2(:,i)=Ved2;
    Istim_d1=Stim(3*nN+1:4*nN,i);
    Istim_s=Stim(4*nN+1:5*nN,i);
    Istim_d2=Stim(5*nN+1:6*nN,i);
  
    if i<nt
    %% Runge-Kutta
    k1=single_neuron_model(nN,X1(:,i),[Vstim_d1;Vstim_s;Vstim_d2;Istim_d1;Istim_s;Istim_d2]);
    % 
    X1_inc=X1(:,i)+dt05*k1(1:8*nN);
    Vmd1_inc=X1_inc(1:nN);Vms_inc=X1_inc(nN+1:2*nN);Vmd2_inc=X1_inc(2*nN+1:3*nN);
    Id1_s_sum=gc_s_d1*wg_matrix_d1_s*(Vms_inc-Vmd1_inc);
    Is_s_sum=wg_matrix_s_s*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_s_sum=gc_s_d2*wg_matrix_d2_s*(Vms_inc-Vmd2_inc);
    Id1_d1_sum=gc_s_d1*wg_matrix_d1_d1*(Vms_inc-Vmd1_inc);
    Is_d1_sum=wg_matrix_s_d1*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_d1_sum=gc_s_d2*wg_matrix_d2_d1*(Vms_inc-Vmd2_inc);
    Id1_d2_sum=gc_s_d1*wg_matrix_d1_d2*(Vms_inc-Vmd1_inc);
    Is_d2_sum=wg_matrix_s_d2*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_d2_sum=gc_s_d2*wg_matrix_d2_d2*(Vms_inc-Vmd2_inc);
    Ved1=1*SF*1e-3/(4*pi)*(Is_d1_sum*As+Id2_d1_sum*Ad2+Id1_d1_sum*Ad1);%units: mV
    Ved2=1*SF*1e-3/(4*pi)*(Is_d2_sum*As+Id2_d2_sum*Ad2+Id1_d2_sum*Ad1);%units: mV
    Ves=1*SF*1e-3/(4*pi)*(Is_s_sum*As+Id2_s_sum*Ad2+Id1_s_sum*Ad1);%units: mV
    Vstim_d1=Stim(1:nN,i)+Ved1;
    Vstim_s=Stim(nN+1:2*nN,i)+Ves;
    Vstim_d2=Stim(2*nN+1:3*nN,i)+Ved2;
    k2=single_neuron_model(nN,X1_inc,[Vstim_d1;Vstim_s;Vstim_d2;Istim_d1;Istim_s;Istim_d2]);
    % 
    X2_inc=X1(:,i)+dt05*k2(1:8*nN);
    Vmd1_inc=X2_inc(1:nN);Vms_inc=X2_inc(nN+1:2*nN);Vmd2_inc=X2_inc(2*nN+1:3*nN);
    Id1_s_sum=gc_s_d1*wg_matrix_d1_s*(Vms_inc-Vmd1_inc);
    Is_s_sum=wg_matrix_s_s*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_s_sum=gc_s_d2*wg_matrix_d2_s*(Vms_inc-Vmd2_inc);
    Id1_d1_sum=gc_s_d1*wg_matrix_d1_d1*(Vms_inc-Vmd1_inc);
    Is_d1_sum=wg_matrix_s_d1*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_d1_sum=gc_s_d2*wg_matrix_d2_d1*(Vms_inc-Vmd2_inc);
    Id1_d2_sum=gc_s_d1*wg_matrix_d1_d2*(Vms_inc-Vmd1_inc);
    Is_d2_sum=wg_matrix_s_d2*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_d2_sum=gc_s_d2*wg_matrix_d2_d2*(Vms_inc-Vmd2_inc);
    Ved1=1*SF*1e-3/(4*pi)*(Is_d1_sum*As+Id2_d1_sum*Ad2+Id1_d1_sum*Ad1);%units: mV
    Ved2=1*SF*1e-3/(4*pi)*(Is_d2_sum*As+Id2_d2_sum*Ad2+Id1_d2_sum*Ad1);%units: mV
    Ves=1*SF*1e-3/(4*pi)*(Is_s_sum*As+Id2_s_sum*Ad2+Id1_s_sum*Ad1);%units: mV
    Vstim_d1=Stim(1:nN,i)+Ved1;
    Vstim_s=Stim(nN+1:2*nN,i)+Ves;
    Vstim_d2=Stim(2*nN+1:3*nN,i)+Ved2;
    k3=single_neuron_model(nN,X2_inc,[Vstim_d1;Vstim_s;Vstim_d2;Istim_d1;Istim_s;Istim_d2]);
    % 
    X3_inc=X1(:,i)+dt*k3(1:8*nN);
    Vmd1_inc=X3_inc(1:nN);Vms_inc=X3_inc(nN+1:2*nN);Vmd2_inc=X3_inc(2*nN+1:3*nN);
    Id1_s_sum=gc_s_d1*wg_matrix_d1_s*(Vms_inc-Vmd1_inc);
    Is_s_sum=wg_matrix_s_s*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_s_sum=gc_s_d2*wg_matrix_d2_s*(Vms_inc-Vmd2_inc);
    Id1_d1_sum=gc_s_d1*wg_matrix_d1_d1*(Vms_inc-Vmd1_inc);
    Is_d1_sum=wg_matrix_s_d1*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_d1_sum=gc_s_d2*wg_matrix_d2_d1*(Vms_inc-Vmd2_inc);
    Id1_d2_sum=gc_s_d1*wg_matrix_d1_d2*(Vms_inc-Vmd1_inc);
    Is_d2_sum=wg_matrix_s_d2*(gc_d1_s*(Vmd1_inc-Vms_inc)+gc_d2_s*(Vmd2_inc-Vms_inc));
    Id2_d2_sum=gc_s_d2*wg_matrix_d2_d2*(Vms_inc-Vmd2_inc);
    Ved1=1*SF*1e-3/(4*pi)*(Is_d1_sum*As+Id2_d1_sum*Ad2+Id1_d1_sum*Ad1);%units: mV
    Ved2=1*SF*1e-3/(4*pi)*(Is_d2_sum*As+Id2_d2_sum*Ad2+Id1_d2_sum*Ad1);%units: mV
    Ves=1*SF*1e-3/(4*pi)*(Is_s_sum*As+Id2_s_sum*Ad2+Id1_s_sum*Ad1);%units: mV
    Vstim_d1=Stim(1:nN,i)+Ved1;
    Vstim_s=Stim(nN+1:2*nN,i)+Ves;
    Vstim_d2=Stim(2*nN+1:3*nN,i)+Ved2;
    k4=single_neuron_model(nN,X3_inc,[Vstim_d1;Vstim_s;Vstim_d2;Istim_d1;Istim_s;Istim_d2]);
    % 
    final_step=dt/6*(k1+2*k2+2*k3+k4);
    X1(:,i+1)=X1(:,i)+final_step(1:8*nN); 
    end
end 

Vmd1=X1(1:nN,:);Vms=X1(nN+1:2*nN,:);Vmd2=X1(2*nN+1:3*nN,:);

E_ext_sd1=(V_ext_s-V_ext_d1)./(L_d1s*10);
E_ext_sd2=(V_ext_s-V_ext_d2)./(L_d2s*10);

run generate_LFP.m
