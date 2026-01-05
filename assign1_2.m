load ModelReductionData.mat
% Center the data
x_c= mean(X,2);
X_c = X - x_c*ones(1,length(X));

% Compute the SVD
[U, D, V] = svd(X_c);
%Singualr values
sigma=diag(D);
% Plot singular values
figure()
semilogy(sigma)
set(gca, 'FontSize', 16)

% Plot scatter plots of first few principal components
Z=U(:,1:2)'*X_c;
figure()
plot(Z(1,:),Z(2,:),'k.','MarkerSize',10)
set(gca,'FontSize',20)
