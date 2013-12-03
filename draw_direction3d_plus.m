function [  ] = draw_direction3d_plus( X, X_new, index )
%DRAW_DIRECTION3D Summary of this function goes here
%   Detailed explanation goes here
figure;
hold on;
grid on;
daspect([1 1 1]);
view([-70,15]);
set(gca,'CameraViewAngle',8);
color='r'; % default color
point = [0,0,0]; % original
for i=1:size(X,2)
    dir=X(:,i);
    % normalize
    dir=dir/norm(dir);
    p_end = point + .5*dir';
    if i==index
        color='b';
    end
    arrow3(point,p_end,color,0.9);
    % recovery
    color='r';
end

dir_new=X_new(:,index);
p_end=point+.5*dir_new';
arrow3(point,p_end,'g',0.9);

hold off;
axis equal;

end

