function [ ] = viewas1(id)
% set views of all figures the same as the first one
set(0, 'currentfigure', id);
fh=findall(0,'type','figure');
num=length(fh);
[az el]=view;
[axis_]=axis;
for i=1:num
    set(0, 'currentfigure', fh(i));
    view(az,el);
    axis(axis_);
end

