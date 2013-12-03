function [] = extract_sharp_points(in,out)
[pn1]=load(in);
[pn2]=load(out);

num=size(pn1,1);

flag=zeros(num,1);
index=1;
tau=30/180*pi;
for i=1:num
    a=pn1(i,4:6);
    b=pn2(i,4:6);
    costheta = dot(a,b)/(norm(a)*norm(b));
    theta = acos(costheta);
    if theta>tau
        flag(index)=i;
        index=index+1;
    end
end
flag(find(flag==0))=[];

% output
ff=fopen('sharp_points.txt','w');
for i=1:length(flag)
    fprintf(ff,'%d %d\n',flag(i)-1,1); % meshlab
end
fclose(ff);
