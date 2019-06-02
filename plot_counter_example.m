% 
% run the experiments
% 
hold on;
[X1 X2] = meshgrid(0.5:0.1:4.5, 4:0.1:6);
Z = (X1-X2).^2/16 + (X2 - 6 - 1/sqrt(5)).^2/4 + sqrt(abs(X1)) + sqrt(abs(X2));
% surf(X1, X2, Z);
contour(X1, X2, Z, 44, 'ShowText','off');

X = [1, 5; 
    1.98143096, 5.19181570; 
    2.96087222, 5.39358042; 
    3.26119870, 5.46791434;
    3.24292043, 5.46411559; 
    3.24286217, 5.46410350];

% text(X(1, 1), X(1, 2), 'x');
% text(X(2, 1), X(2, 2), 'x');
% text(X(3, 1), X(3, 2), 'x');
% text(X(4, 1), X(4, 2), 'x');
% text(X(5, 1), X(5, 2), 'x');
% text(X(6, 1), X(6, 2), '*');
text(X(1, 1)-0.02, X(1, 2)-0.02, 'x');
arrow(X(1,:), X(2,:), 'BaseAngle', 20, 'Length', 7);
arrow(X(2,:), X(3,:), 'BaseAngle', 20, 'Length', 7);
arrow(X(3,:), X(4,:), 'BaseAngle', 20, 'Length', 7);
text(X(4, 1)-0.04, X(4, 2)-0.04, '*');
% arrow(X(4,:), X(5,:));
% arrow(X(5,:), X(6,:));

