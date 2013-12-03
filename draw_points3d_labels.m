function [ ] = draw_points3d_labels( points,labels )
figure;
axis equal;
hold on;
grid on;
[m,n]=size(points);
C=labels; % color matrix
dim=size(points,2);
if dim==2
    scatter(points(:,1),points(:,2),5,C,'filled');
else
    scatter3(points(:,1),points(:,2),points(:,3),5,C,'filled');
end
hold off;
set(gca,'GridLineStyle','-');
set(gcf, 'Color', 'w');
axis equal;
axis([-1 1 -1 1 -1 1]);
axis tight;
end

