
%This script takes the outline from a combined KPP + PP profiles
%figure and swaps in new data from other figures.

f_from = gcf;
h_from = zeros(6,2);
a_from = get(f_from,'children');
cnt = 1;
for a=1:length(a_from)
    pos = get(a_from(a),'Position');
    if (pos(1)<0.5)
        h_from(cnt,1) = a_from(a);
        h_from(cnt,2) = pos(2);
        cnt = cnt+1;
    end
end
[Y,I] = sort(h_from(:,2),1);
h_from = h_from(I,1);

f_to   = gcf;
h_to = zeros(6,1);
a_to = get(f_to,'children');
cnt = 1;
for a=1:length(a_to)
    pos = get(a_to(a),'Position');
    if (pos(1)>0.35 & pos(1) < 0.6)
% $$$     if (pos(1)<0.35)
        h_to(cnt) = a_to(a);
        h_to(cnt,2) = pos(2);
        cnt = cnt+1;
    end
end
[Y,I] = sort(h_to(:,2),1);
h_to = h_to(I,1);

for a=1:length(h_to)
delete(findobj(h_to(a),'type','line'));
delete(findobj(h_to(a),'type','surface'));
delete(findobj(h_to(a),'type','hggroup'));
copyobj(findobj(h_from(a),'type','surface'),h_to(a));
copyobj(findobj(h_from(a),'type','hggroup'),h_to(a));
copyobj(findobj(h_from(a),'type','line'),h_to(a));
uistack(findobj(h_to(a),'type','text'),'top');
end
