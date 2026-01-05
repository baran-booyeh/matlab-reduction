% Load the IrisData.mat file
load IrisDataAnnotated.mat

% Compute the mean vector and subtract it from the data
x_c = 1/length(X)*sum(X,2);
X_c = X - x_c*ones(1,length(X));

% Compute the SVD of the centered data
[U,D,V] = svd(X_c);
sigma=diag(D);
figure(1)
plot(sigma)
xlabel('SV index')
ylabel('SV')
set(gca, 'FontSize', 16)
% Project the centered data onto the first two principal components
Z = U(:,1:2)'*X_c;

% Create a scatter plot of the first two principal components
figure(2)
plot(Z(1,:),Z(2,:),'k.','MarkerSize',10);
xlabel('First Principal Component');
ylabel('Second Principal Component');

% Project the centered data onto the first three principal components
Z_1= U(:,1:3)'*X_c;

% Create a scatter plot of the first three principal components
figure(3)
plot3(Z_1(1,:),Z_1(2,:),Z_1(3,:),'k.','MarkerSize',10);
xlabel('First Principal Component');
ylabel('Second Principal Component');
zlabel('Thirs Principal Component');

figure(4)
scatter(Z(1,:),Z(2,:), 25, I, 'filled');
xlabel('First Principal Component');
ylabel('second Principal Component');

I1 = find(I == 1);
I2 = find(I == 2);
I3 = find(I == 3);

%for any subset we have to considered the centered matrices by 
% subtracting the columnwise means of each subset, so we can get
%the within scatter matrices

X1 = X(:,I1);
X2 = X(:,I2);
X3 = X(:,I3);

 

p1 = length(I1);
c1 = 1/p1*sum(X1,2);
X1c = X1 - c1*ones(1,p1);


p2 = length(I2);
c2 = 1/p2*sum(X2,2);
X2c = X2 - c2*ones(1,p2);


p3 = length(I3);
c3 = 1/p3*sum(X3,2);
X3c = X3 - c3*ones(1,p3);

  
%compute the within cluster by summing the scatter matrices of
%the centered subsets

Sw = X1c*X1c'+X2c*X2c'+X3c*X3c';

%Compute the between scatter matrix by using the total mean (c) based on
%the means of each subset

p = p1+p2+p3;
c = 1/p*sum(X,2);
Sb = p1*(c1-c)*(c1-c)'+p2*(c2-c)*(c2-c)'+p3*(c3-c)*(c3-c)';

 

%%

%Check if regularization is required

%visualization of the eignvalues of the within cluster to understand
%the distribution of singualr values and the impact of regularization
 

Xw = [X1c,X2c,X3c];
sXw = svd(Sw);
sXw2 = sXw.^2;
[eig(Sw), flip(sXw2)];

tau = 10^(-10);

%Computes the regularization term 

epsilon = tau*max((eig(Sw)));


figure()
semilogy(sXw2,'*')
hold on
semilogy(sXw2+epsilon,'sm');
legend('eigenvalues in decreasing order','scaled eigenvalues')
title('Sw eigenvalues as Squared singular values of Xw')

 

%%
%Check if  we need regularization by checking 
% whether the determinant of Sw is below 10^-15.
%Regularization is performed by adding epsilon which is a small positive
%value to diagonal elements of Sw. 
%Helps to avoid instability for inversing the Sw.

det(Sw)

if (det(Sw) <=10^(-15))

disp('do regularization')

sSw = eig(Sw);

Sw = Sw+eye(size(Sw))*epsilon;

sSwreg = eig(Sw);

figure()

semilogy(sSw,'r*')

hold on

semilogy(sSwreg,'bd')

legend('eig(Sw','eig(Sw reg)')

disp('det(Sw) after regularization')

det(Sw)

 

end

 

% Cholesky factor of within cluster spread
%To decompose the matrix into the product of an upper triangualr matrix (K)
%and its conjugate transpose

K = chol(Sw);

A = (K'\Sb)/K; %A= inv(k)*Sb*inv(k')

%Solving the generalized eigen value problem to find the
%eigenvectors(W)(Representts the direction that data should be projected to
%maximize seperation between classes)
% and corresponding eigenvalues(D)(The importance of each eigenvector)
% of matrix A

% here I could have also asked just for W = eig(A), which was the first
% element of the list I wanted (otherwise I would have used tilde)


[W,D] = eig(A);

Q = inv(K)*W; 

 

% Q(:,1:2)' : Selects first 2 columns from matrix Q and transpose
% them to have a matrix with eigenvectors as the rows
% (Subset of eigenvectors obtained from step before)


%Project data onto the LDA directions by multiplying Q(:,1:2)' with each
%centered subset to get LDA components

LDA1 = Q(:,1:2)' * X1;
LDA2 = Q(:,1:2)' * X2;
LDA3 = Q(:,1:2)' * X3;

% for PCA recall we had Z = Q(:,1:2)'*X;

%Histogram showing the distribution of LDA components for each subset

figure(5)

histogram(LDA1)

hold on

histogram(LDA2)

histogram(LDA3)

title("Components for LDA")

%Scatterplot of first two LDA components with each 
% subset showing in different colors 

figure(6)

plot(LDA1(1,:),LDA1(2,:),'r.', "MarkerSize",15 )

hold on

plot(LDA2(1,:),LDA2(2,:),'y.', "MarkerSize",15 )

plot(LDA3(1,:),LDA3(2,:),'g.', "MarkerSize",15 )

hold off

set(gca, "FontSize", 20)

title('LDA', "FontSize",20)
