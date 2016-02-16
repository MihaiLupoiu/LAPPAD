AC=[[1,2,0],
    [3,4,5],
    [6,7,8],
    [9,10,0],
    [11,12,0]];

JC=[[1,3,1],
    [1,2,4],
    [2,3,5],
    [3,4,4],
    [4,5,5]];

v=[1,2,3,4,5];

[m,n]=size(AC);
R=zeros(m,1);

for i=1:m
    sum = 0;
    for j=1:n
        sum=sum+AC(i,j)*v(JC(i,j));
    end
    R(i)=sum;
end