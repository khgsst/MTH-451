function BS(b,R)
[m,m]=size(R)
for j=m:-1:1
x(:,j)=b(:,j)-x/R(j,j);
end
