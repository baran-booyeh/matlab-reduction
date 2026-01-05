
clear
close all


load HandwrittenDigits.mat
I1=find(I==3);
I1l=length(I1);
I2=find(I==6);
I2l=length(I2);
I3=find(I==9);
I3l=length(I3);
J=[I1,I2,I3];
X=X(:,J);
I=[1*ones(1,I1l), 2*ones(1,I2l), 3*ones(1,I3l)];

dim=size(X,1);
[n,m]=size(X);
p=m;
for i=1:p
    for j=i:p
        if (i==j)
            D(i,j)=0;
        else
            D(i,j)=norm(X(:,i)-X(:,j),2);
            %TRY ALSO NORM 1 and NORM 2
            %compare also the iteration number
        end
    end
    %     end
end
D=D+D';

%THIS IS AN EXAMPLE...
I_m=sort(randperm(length(J),3));
% !!!!!pick up 3 random starting medoids indexes !!!!
 

figure()
plot(X(1,:),X(2,:),'bo','MarkerSize',8);
hold on ;
plot(X(1,I_m),X(2,I_m),"xr",'MarkerSize',8);
xlabel("x")
ylabel("y")
legend('Data','Initial medoids')

Err=1;
itmax = 100;
tol=1.0e-10;
iter=0;
while(Err>=tol && iter < itmax) %Stopping criterion

    %% assignment step

    D_m=D(1:end,I_m);%distances w.r.t. medoids submatrix

    [q,I_assign]=min(D_m') %TRANSPOSE or use correctly min()

    %output
    Q=sum(q)

    k=length(I_m); %number of clusters
    %%%%%%%%%%%%%
    oldI_m=I_m;
    clear I_m

    %% updating step
    for ell=1:k

        I_ell = find(I_assign == ell) % Indices to points in the cluster

        D_ell = D(I_ell,I_ell)

        s=sum(D_ell)
        [qq(ell),j] = min(sum(D_ell))

        %from local to global indicees
        I_m(ell) = I_ell(j); %swap (updating) of the medoids indeces


        %final_medoids=I_m
        newI_m=I_m
        oldI_m
    end


    %% recompute the global coherence
    Qnew=sum(qq);


    %% if not converged,  continue from the assignment step
    Err=abs(Q-Qnew)

    Errplot(iter+1)=Err;

    Q=Qnew;

    iter=iter+1

    if (Err < tol)
        flag=0;
    else
        flag=1;
    end
end

disp('flag')
flag

final_medoids=I_m


%% Q plot
figure()
semilogy([1:iter],Errplot,'bo-')
[[1:iter]' Errplot(:)];


%% Final Clusters
figure()
for j=1:k
    X_l{j} = X(:,find(I_assign == j)) ; %k-ith  cluster
end

scatter3(X_l{1}(1,:),X_l{1}(2,:),X_l{1}(3,:),'red' )
hold on
scatter3(X_l{2}(1,:),X_l{2}(2,:),X_l{2}(3,:),'blue' )
hold on
scatter3(X_l{3}(1,:),X_l{3}(2,:),X_l{3}(3,:),'green' )
title('Final clustering')

figure()
scatter(X(1,:),X(2,:));
hold on ;
scatter(X(1,I_m),X(2,I_m),"xr");
xlabel("x")
ylabel("y")
legend('Data','Final medoids')


figure()
cm = confusionchart(I,I_assign)
%The element cm(i,j) is the number of
%times an observation in the ith true class was predicted to be in the
%jth class.


