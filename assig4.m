close all
clear

load Yale_64x64.mat
r1=80;
faceW=64;
faceH=64;
numperLine=5;
showline=1;
Y = zeros(faceH*showline,faceW*numperLine);

% Center the data
x_c =mean(fea,2); 
X_c = fea - x_c*ones(1,length(fea));
[U,D,V] = svds(X_c,r1);
X_r1 = U*D*V';

j=0;
for i=[9,20,33,47,59] %Select 5 random representative faces 
    if j<6
     Y(1:64,j*faceW+1:(j+1)*faceW)= reshape(X_r1(i,:),[faceH,faceW]);
     j=j+1;
    end
  imagesc(Y);
  colormap("gray");
end

%Plot the log of Singular Values of X_c
sigma=diag(D);
figure()
semilogy(sigma, 'k', 'LineWidth', 2)

%Plot the first 5 feature vectors
for a = 1:5
    figure()
    Z = U(:,a)'*X_c;
    imagesc(reshape(Z(1,:),64,64));
    colormap("gray")
end

%approximate the selected images using first k feature vectors
numperLine2=3;
showLine2=2;
Y=zeros(faceH*showLine2,faceW*numperLine2);
for i=[9,20,33,47,59]
  n=0;
  figure()
  for k= [4,8,15]
     k_c =mean(fea,2); 
     K_c = fea - k_c*ones(1,length(fea));
     [U,D,V] = svds(K_c,k);
     X_k = U*D*V';
     Y(1:64,n*faceW+1:(n+1)*faceW)=reshape(X_k(i,:),[faceW,faceH]);
     Y(65:128,n*faceW+1:(n+1)*faceW)=reshape(K_c(i,:)-X_k(i,:),[faceW,faceH]);
     n=n+1;
     imagesc(Y);
     colormap("gray")
  end
end
