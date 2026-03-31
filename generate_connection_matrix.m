function y=generate_connection_matrix(s_c,x)
  n=x;
  if s_c==0
  cn_matrix=zeros(n);
  for i=2:n
      for j=1:n
          if j<i
              cn_matrix(i,j)=1;
          end
      end
  end
  else
  cn_matrix=ones(n);    
  end
  y=cn_matrix;
end

