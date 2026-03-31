function y = generate_cor_xyz(x1,x2,x3,x4,x5)

nN=x1;
nN_z=x2;
L_d1s=x3(1);
L_d2s=x3(2);
Ds=x4;
dis_ext=x5;

cor_soma_x=zeros(nN_z,nN);
cor_soma_y=zeros(nN_z,nN);
cor_soma_z=zeros(nN_z,nN);
cor_den1_x=zeros(nN_z,nN);
cor_den1_y=zeros(nN_z,nN);
cor_den1_z=zeros(nN_z,nN);
cor_den2_x=zeros(nN_z,nN);
cor_den2_y=zeros(nN_z,nN);
cor_den2_z=zeros(nN_z,nN);

for i=1:nN_z
    for j=1:nN
        cor_soma_y(i,j)=0;
        cor_den2_y(i,j)=L_d2s;
        cor_den1_y(i,j)=-L_d1s;
        if j>=2            
           cor_soma_x(i,j)=cor_soma_x(i,j-1)+Ds+dis_ext;
           cor_den1_x(i,j)=cor_den1_x(i,j-1)+Ds+dis_ext;
           cor_den2_x(i,j)=cor_den2_x(i,j-1)+Ds+dis_ext;           
        end
        if i>=2           
           cor_soma_z(i,j)=cor_soma_z(i-1,j)+Ds+dis_ext;
           cor_den1_z(i,j)=cor_den1_z(i-1,j)+Ds+dis_ext;
           cor_den2_z(i,j)=cor_den2_z(i-1,j)+Ds+dis_ext;           
        end
    end
end
y=[cor_soma_x;cor_soma_y;cor_soma_z;
   cor_den1_x;cor_den1_y;cor_den1_z;
   cor_den2_x;cor_den2_y;cor_den2_z];
end

