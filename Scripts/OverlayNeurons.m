base = MD(215);
rotateme = MD(217);
load(fullfile(rotateme.Location,'ProcOut.mat'),'Xdim','Ydim');

%In 11_30_2015 folder.
load(fullfile(MD(215).Location,...
    'RegistrationInfo-GCamp6f_45_treadmill-12_01_2015-session4.mat'));

% %In 12_01_2015 folder.
% load(fullfile(MD(217).Location,...
%     'RegistrationInfo-GCamp6f_45_treadmill-11_30_2015-session4.mat'));

x = y2;
y = x2; 
%R = RegistrationInfoX.tform.T; 
R = RegistrationInfoX.tform.T(1:2,1:2);
aligned = cell(length(x),1); 

for n=1:length(x)
%     center = repmat([Xdim/2; Ydim/2; 0],1,length(x{n}));
%     v = [x{n}'; y{n}'; ones(1,length(x{n}))];

    center = repmat([Ydim/2; Xdim/2],1,length(x{n}));
    v = [x{n}'; y{n}'];
    t = repmat(RegistrationInfoX.tform.T(3,1:2)',1,length(x{n}));
    
    aligned{n} = R*(v-center) + center + t;
end

figure;
PlotNeurons(base,1:1202,[1 0 0 .5],1);
hold on;
for n=1:length(x)
    plot(aligned{n}(1,:),aligned{n}(2,:),'color',[0 0 1 .5],'linewidth',1);
end
title('x = y2');