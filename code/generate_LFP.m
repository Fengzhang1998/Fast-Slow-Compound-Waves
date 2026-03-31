% clear all;
multiplyer = 1.0;

dis_1_nN=cor_soma_x(1,nN)-cor_soma_x(1,1); %units: cm

n_x_mea=3;
n_y_mea=2;

dis_x_mea=dis_1_nN/(n_x_mea-1); %ounits: cm 
dis_y_mea=1e-4*210; %units: cm 
dis_z_mea=-1e-4*60;

cor_mea_x=zeros(n_y_mea,n_x_mea);
cor_mea_y=zeros(n_y_mea,n_x_mea);
cor_mea_z=dis_z_mea+zeros(n_y_mea,n_x_mea);
%
for i=1:n_y_mea
    for j=1:n_x_mea
        if j~=1
           cor_mea_x(i,j)=cor_mea_x(i,j-1)+dis_x_mea;
        end
        if i~=1
           cor_mea_y(i,j)=cor_mea_y(i-1,j)+dis_y_mea; 
        end
    end
end

cor_mea_x(1,1)=cor_soma_x(1,floor(nN/4))-cor_soma_x(1,1);
cor_mea_x(1,2)=cor_soma_x(1,floor(nN/2))-cor_soma_x(1,1);
cor_mea_x(1,3)=cor_soma_x(1,ceil(3*nN/4))-cor_soma_x(1,1);
cor_mea_x(2,1)=cor_soma_x(1,floor(nN/4))-cor_soma_x(1,1);
cor_mea_x(2,2)=cor_soma_x(1,floor(nN/2))-cor_soma_x(1,1);
cor_mea_x(2,3)=cor_soma_x(1,ceil(3*nN/4))-cor_soma_x(1,1); 

wg_matrix_s=zeros(n_y_mea*n_x_mea,nN);  
wg_matrix_d1=zeros(n_y_mea*n_x_mea,nN);
wg_matrix_d2=zeros(n_y_mea*n_x_mea,nN);

for i=1:n_y_mea*n_x_mea
        dis_mea_s=sqrt((cor_mea_x(i)-cor_soma_x).^2+(cor_mea_y(i)-cor_soma_y).^2+(cor_mea_z(i)-cor_soma_z).^2);
        dis_mea_d1=sqrt((cor_mea_x(i)-cor_den1_x).^2+(cor_mea_y(i)-cor_den1_y).^2+(cor_mea_z(i)-cor_den1_z).^2);
        dis_mea_d2=sqrt((cor_mea_x(i)-cor_den2_x).^2+(cor_mea_y(i)-cor_den2_y).^2+(cor_mea_z(i)-cor_den2_z).^2); % cor_den2_z大小为(nN_z,nN)
        wg_matrix_s(i,:)=sum(1./(4.*dis_mea_s));
        wg_matrix_d1(i,:)=sum(1./(4.*dis_mea_d1));
        wg_matrix_d2(i,:)=sum(1./(4.*dis_mea_d2));
end

t0=0000;
Is_d1_sum=gc_s_d1*wg_matrix_d1*(Vms-Vmd1); 
Id1d2_s_sum=wg_matrix_s*(gc_d1_s*(Vmd1-Vms)+gc_d2_s*(Vmd2-Vms));
Is_d2_sum=gc_s_d2*wg_matrix_d2*(Vms-Vmd2);

Ve_mea=SF*1e-3*Re/(4*pi)*(Is_d1_sum*Ad1+Id1d2_s_sum*As+Is_d2_sum*Ad2);%units: mV

%% 
Ve_mea_ns=Ve_mea(1:2:length(Ve_mea(:,1)),:);
Ve_mea_nd=Ve_mea(2:2:length(Ve_mea(:,1)),:);
Ve_mea_sd=Ve_mea_ns-Ve_mea_nd;
Ve_mea_ns_mean=mean(Ve_mea_ns);
Ve_mea_nd_mean=mean(Ve_mea_nd);


