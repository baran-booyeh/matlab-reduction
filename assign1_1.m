load ModelReductionData

% Plot scatter plots for all combinations of two dimensions
figure(1)
k = 1;
for i = 1:5
    for j = i+1:6
        subplot(5, 3, k)
        plot(X(i,:), X(j,:), 'k.', 'MarkerSize', 7)
        axis('equal')
        set(gca, 'FontSize', 8)
        xlabel(['X=',num2str(i)])
        ylabel(['Y=',num2str(j)])
        k = k + 1;
    end
end