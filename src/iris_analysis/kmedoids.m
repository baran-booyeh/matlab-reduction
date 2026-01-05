close all
clear

load IrisDataAnnotated.mat
rng(1)
[n,p]=size(X);

for i =1:p
        for j=i:p
            if (i==j)
               D(i,j)=0;
            else
                  D(i,j)=norm(X(:,i)-X(:,j),2);
            end
        end
end
D=D'
I_m=[2,20,100];
D_m+D(1:end,I_m)

