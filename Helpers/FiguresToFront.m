%% Bring all figures to front. 

f = findobj(0,'type','figure');
for i=[f.Number]
    figure(i);
end