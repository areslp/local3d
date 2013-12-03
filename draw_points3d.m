function [ ] = draw_points3d( points,t )
%DRAW_POINTS3D Summary of this function goes here
%   Detailed explanation goes here
figure;
% title(t);
axis equal;
hold on;
daspect([1 1 1])
% view([-70,15]), set(gca,'CameraViewAngle',8);

view([20,15]), set(gca,'CameraViewAngle',8);
grid on;
% scatter3(points(:,1),points(:,2),points(:,3),3,[0,0,0]);
x=points(:,1);
y=points(:,2);
% [X Y]=meshgrid(x,y);
% % Z=points(:,3:size(points,2));
Z=points(:,3);
% mesh(X,Y,Z);
plot3(x,y,Z,'b*');
hold off;
% axis equal;

end

