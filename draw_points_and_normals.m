function [ ] = draw_points_and_normals( points, normals,t )

p1=[min(points(:,1)),min(points(:,2)),min(points(:,3))];
p2=[max(points(:,1)),max(points(:,2)),max(points(:,3))];
diag_len=norm(p1-p2);

p_num = size(points,1);
figure_num=figure;
hold on;
scatter3(points(:,1),points(:,2),points(:,3),3,[1,0,0]);
for i=1:p_num
    point = points(i,:);
    p_end = point + diag_len/30*normals(i,:);
    mt = [point;p_end];
    %ARROW(point,p_end);
    line(mt(:,1),mt(:,2),mt(:,3),'Color',[0,1,0]);
end
hold off;
% axis equal;
end

