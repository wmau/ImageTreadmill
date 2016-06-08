[~,a] = find(graphData.A);
a = unique(a);
ratio = [];
TMoffsets = [];
for i=a'
    [r,t] = TMvsCellOnset(MD(260),graphData,i);
    ratio = [ratio r]; TMoffsets = [TMoffsets t];
end
scatter(TMoffsets,ratio,'.');