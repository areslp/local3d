function [  ] = draw_direction3d( X )
%DRAW_DIRECTION3D Summary of this function goes here
%   Detailed explanation goes here
figure;
hold on;
grid on;
daspect([1 1 1]);
view([-70,15]);
dim=size(X,1);
set(gca,'CameraViewAngle',8);
for i=1:size(X,2)
    dir=X(:,i);
    %归一化
    % dir=dir/norm(dir);
    point = zeros(1,dim);
    p_end = point + .5*dir';
    mt = [point;p_end];
%     line(mt(:,1),mt(:,2),'Color',[0,1,0]);
    arrow3(point,p_end,'r',0.9);
end
hold off;
axis equal;

end

