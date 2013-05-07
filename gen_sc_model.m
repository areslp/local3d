% fname='test';
fname='fandisk';
mesh=[fname '.off'];
[vertex, face]=read_mesh(mesh);
[normal,normalf] = compute_normal(vertex,face);
ff=fopen([fname '.xyzn'],'w');
vfring = compute_vertex_face_ring(face);
mr=0; % mesh resolution
for i=1:length(vertex)
    v=vertex(:,i);
    fring=vfring{i};
    fidx=fring(1);
    n=normalf(:,fidx);
    fprintf(ff,'%f %f %f %f %f %f\n',v(1),v(2),v(3),n(1),n(2),n(3));
end
fclose(ff);

% compute the mr of the dual graph
for i=1:length(face)
    v0=vertex(:,face(1,i));
    v1=vertex(:,face(2,i));
    v2=vertex(:,face(3,i));
    mr=mr+norm(v0-v1);
    mr=mr+norm(v1-v2);
    mr=mr+norm(v2-v0);
end
mr=mr/(3*length(face));
ff=fopen(['mr.txt'],'w');
fprintf(ff,'%f\n',mr);
fclose(ff);

