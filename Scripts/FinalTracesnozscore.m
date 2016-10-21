NumFrames = 30849;
NumNeurons = 1202;
load('FinalOutput.mat','NeuronPixels');
moviefile = 'SLPDF.h5';
info = h5info(moviefile,'/Object');
trace = zeros(1,NumFrames);
ROI = NeuronPixels{1};

%parpool('local');
resol = 2;                                  % Percent resolution for progress bar, in this case 5%
update_inc = round(NumFrames/(100/resol));  % Get increments for updating ProgressBar
p = ProgressBar(100/resol);
parfor i=1:NumFrames
    tempFrame = loadframe(moviefile,i,info);
    tempFrame = tempFrame(:);
    
    trace(i) = mean(tempFrame(ROI));
    
    if round(i/update_inc) == (i/update_inc)
        p.progress;
    end
end
p.stop; % Terminate progress bar
delete(gcp);
