A=[[1,0,2,0,0],
   [3,4,0,5,0],
   [0,6,7,0,8],
   [0,0,9,10,0],
   [0,0,0,11,12]];

V=[1,2,3,4,5];

R1=A*V';

DA=[[0,1,2],
    [3,4,5],
    [6,7,8],
    [9,10,-0],
    [11,12,0]];

IOFF=[-1,0,2];

[m,n]=size(DA);
R=zeros(m,1);

for i=1:m
    sum = 0;
    for j=1:n
        index=i+IOFF(j);
        if (index > 0 && index <=m)
            sum = sum + DA(i,j) * V(i+IOFF(j));
        end        
    end
    R(i)=sum;
end