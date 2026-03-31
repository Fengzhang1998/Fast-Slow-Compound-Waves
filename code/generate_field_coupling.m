function y = generate_field_coupling(x1,x2,x3,x4)

nN=x1;nN_z=x2;
cor_soma_x=x3(1:nN_z,:);
cor_soma_y=x3(nN_z+1:2*nN_z,:);
cor_soma_z=x3(2*nN_z+1:3*nN_z,:);
cor_den1_x=x3(3*nN_z+1:4*nN_z,:);
cor_den1_y=x3(4*nN_z+1:5*nN_z,:);
cor_den1_z=x3(5*nN_z+1:6*nN_z,:);
cor_den2_x=x3(6*nN_z+1:7*nN_z,:);
cor_den2_y=x3(7*nN_z+1:8*nN_z,:);
cor_den2_z=x3(8*nN_z+1:9*nN_z,:);
Re=x4;

wg_matrix_s_s=zeros(nN,nN);
wg_matrix_d2_s=zeros(nN,nN);
wg_matrix_d1_s=zeros(nN,nN);

wg_matrix_s_d2=zeros(nN,nN);
wg_matrix_d2_d2=zeros(nN,nN);
wg_matrix_d1_d2=zeros(nN,nN);

wg_matrix_s_d1=zeros(nN,nN);
wg_matrix_d2_d1=zeros(nN,nN);
wg_matrix_d1_d1=zeros(nN,nN);

for i=1:nN
        for j=1:nN
           if j~=i
   
               d_ss=0;
               d_d2s=0;
               d_d1s=0;
               
               d_sd2=0;
               d_d2d2=0;
               d_d1d2=0;
               
               d_sd1=0;
               d_d2d1=0;
               d_d1d1=0;
               
               for k=1:nN_z
                   %
                   dr_ss=sqrt((cor_soma_x(k,j)-cor_soma_x(1,i))^2+(cor_soma_y(k,j)-cor_soma_y(1,i))^2+(cor_soma_z(k,j)-cor_soma_z(1,i))^2);
                   dr_d2s=sqrt((cor_den2_x(k,j)-cor_soma_x(1,i))^2+(cor_den2_y(k,j)-cor_soma_y(1,i))^2+(cor_den2_z(k,j)-cor_soma_z(1,i))^2);
                   dr_d1s=sqrt((cor_den1_x(k,j)-cor_soma_x(1,i))^2+(cor_den1_y(k,j)-cor_soma_y(1,i))^2+(cor_den1_z(k,j)-cor_soma_z(1,i))^2);
                   
                   dr_sd2=sqrt((cor_soma_x(k,j)-cor_den2_x(1,i))^2+(cor_soma_y(k,j)-cor_den2_y(1,i))^2+(cor_soma_z(k,j)-cor_den2_z(1,i))^2);
                   dr_d2d2=sqrt((cor_den2_x(k,j)-cor_den2_x(1,i))^2+(cor_den2_y(k,j)-cor_den2_y(1,i))^2+(cor_den2_z(k,j)-cor_den2_z(1,i))^2);
                   dr_d1d2=sqrt((cor_den1_x(k,j)-cor_den2_x(1,i))^2+(cor_den1_y(k,j)-cor_den2_y(1,i))^2+(cor_den1_z(k,j)-cor_den2_z(1,i))^2);
                   
                   dr_sd1=sqrt((cor_soma_x(k,j)-cor_den1_x(1,i))^2+(cor_soma_y(k,j)-cor_den1_y(1,i))^2+(cor_soma_z(k,j)-cor_den1_z(1,i))^2);
                   dr_d2d1=sqrt((cor_den2_x(k,j)-cor_den1_x(1,i))^2+(cor_den2_y(k,j)-cor_den1_y(1,i))^2+(cor_den2_z(k,j)-cor_den1_z(1,i))^2);
                   dr_d1d1=sqrt((cor_den1_x(k,j)-cor_den1_x(1,i))^2+(cor_den1_y(k,j)-cor_den1_y(1,i))^2+(cor_den1_z(k,j)-cor_den1_z(1,i))^2);
                   
                   d_ss=d_ss+Re/(4.*dr_ss);
                   d_d2s=d_d2s+Re/(4.*dr_d2s);
                   d_d1s=d_d1s+Re/(4.*dr_d1s);
                   
                   d_sd2=d_sd2+Re/(4.*dr_sd2);
                   d_d2d2=d_d2d2+Re/(4.*dr_d2d2);
                   d_d1d2=d_d1d2+Re/(4.*dr_d1d2);
                   
                   d_sd1=d_sd1+Re/(4.*dr_sd1);
                   d_d2d1=d_d2d1+Re/(4.*dr_d2d1);
                   d_d1d1=d_d1d1+Re/(4.*dr_d1d1);
                                  
               end
               wg_matrix_s_s(i,j)=d_ss;
               wg_matrix_d2_s(i,j)=d_d2s;
               wg_matrix_d1_s(i,j)=d_d1s;
               
               wg_matrix_s_d2(i,j)=d_sd2;
               wg_matrix_d2_d2(i,j)=d_d2d2;
               wg_matrix_d1_d2(i,j)=d_d1d2;
               
               wg_matrix_s_d1(i,j)=d_sd1;
               wg_matrix_d2_d1(i,j)=d_d2d1;
               wg_matrix_d1_d1(i,j)=d_d1d1;
                             
           end
        end
end
y=[wg_matrix_s_s;wg_matrix_d2_s;wg_matrix_d1_s;
   wg_matrix_s_d2;wg_matrix_d2_d2;wg_matrix_d1_d2;
   wg_matrix_s_d1;wg_matrix_d2_d1;wg_matrix_d1_d1];
end

