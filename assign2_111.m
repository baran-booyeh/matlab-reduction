close all
clear

load HandwrittenDigits.mat
 
% Extract from X the images that correspond to numbers 0, 1, 3 and 7
for j=[0,4]
    Ij=find(I==j);
    Y= X(:,Ij);

    %Centering the data
    y_c =1/length(Ij)*sum(Y,2);
    Y_c = Y-y_c*ones(1,length(Ij));
    sign

    %Singular value decomposition of Y_c:
    [U,D,V]= svd(Y_c);

    % Select a subset of samples
    Xj= Y_c(:,1:5);
    for l=1:size(Xj,2) %Fix the column of samples
    % Loop over the values of k
    figure;
     k=1;
     for k_values=[5, 10, 15, 20, 25]
       
        % Compute the projection onto the first k principal components
        Z = U(:, 1:k_values)'* Xj(:,l);
        
        % Reconstruct the samples using the projection coefficients and the first k principal components
        ZP  = U(:,1:k_values)*Z;
        
        % Compute the residual for each sample
        residual=Xj(:,l)-ZP;
        
        % Compute the norms of the errors for this value of k
        errors(k_values/5) = norm(residual);
        

           
         
                hold on
                subplot(3, 5, k)
                imagesc(reshape(Xj(:, l), 16, 16)')
                axis('square');
                axis('off');
                title(sprintf('k=%d',k_values));
                subplot(3, 5, k+5)
                imagesc(reshape(ZP , 16, 16)')
                axis('square');
                axis('off');
                subplot(3, 5, k+10)
                imagesc(reshape(residual, 16, 16)')
                colormap(1-gray);
                axis('square');
                axis('off');
                
                k=k+1;
                
            end
    end

    % Plot the norms of the errors as a function of k for this group of samples
    figure();
    plot([5, 10, 15, 20, 25], errors);
    xlabel('k');
    ylabel ('Norm of error');
    % Create a scatter plot of the first two principal components
   
    figure()
plot(Z(1,:),Z(2,:),'k.','MarkerSize',10);
xlabel('First Principal Component');
ylabel('Second Principal Component');

% Project the centered data onto the first three principal components
Z_1= U(:,1:2)'*X_c;

% Create a scatter plot of the first three principal components
figure()
plot3(Z_1(1,:),Z_1(2,:),'k.','MarkerSize',10);
xlabel('First Principal Component');
ylabel('Second Principal Component');


figure()
scatter(Z(1,:),Z(2,:), 25, I, 'filled');
xlabel('First Principal Component');
ylabel('second Principal Component');
    
end 
        
     
 