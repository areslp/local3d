function [ ] = write_xyzn( points_with_normals )
% write xyzn file
num=size(points_with_normals,1);

vertex=points_with_normals(:,1:3)';
normals_new=points_with_normals(:,4:6);
ff=fopen('paper.xyzn','w');
for i=1:num
    fprintf(ff,'%f %f %f %f %f %f\n',vertex(1,i),vertex(2,i),vertex(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
    % fprintf(ff,'%f %f %f %f %f %f\n',vertex1(1,i),vertex1(2,i),vertex1(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
end
fclose(ff);
