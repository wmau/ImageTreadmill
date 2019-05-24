%Load the sessions.
loadMD;

%% Figure 1. 
%a. Schematic of maze. Not coded.

%b,c. Maximum projection and traces. 
GenerateTraces;

%d. Time cell rasters.
plotTimeCells(MD(292),'singletraces','LPtrdmll','timecells',906);
plotTimeCells(MD(292),'singletraces','LPtrdmll','timecells',2);
plotTimeCells(MD(292),'singletraces','LPtrdmll','timecells',113);
plotTimeCells(MD(292),'singletraces','LPtrdmll','timecells',431);

%e. Run-averaged time lapse images of example time cell.
TCStills(MD(292),'movie','BPDFF.h5','neurons',431)

%f. Time cell ensemble.
ensemble = PastalkovaPlot(MD(292)); 
PlopLogLineonPastalkovaPlot;

%% Figure 2. 
%a. Posterior probability plots of individual trials. Difficult to
%reproduce exactly because of the randomized nature of trial selection for
%training and test. Other figures like this can be made by: 
[Mdl,~,testX] = TimeDecoder(MD(293));
[~,postProbs] = PredictTime(Mdl,testX,'plotit',false);
TimeDecoderPosteriorProbabilityPlot(postProbs(:,:,N)); %Where N is a trial.

%b. Average posterior probabilities of decoder predictions. Also difficult
%to reproduce exactly but one similar to the one in the manuscript can be
%made by:
PredictTime(Mdl,testX);

%c,d. Decoder error as a function of elapsed time and overall error
%compared to shuffle.
FullDataDecodeError;

%% Figure 3.
%a. Within-session exiting time cell.
plotTimeCells(MD(292),'timecells',585);
TrialBreakdownFluor_time(MD(292),585,6);

%b. Within-session entering time cell.
plotTimeCells(MD(292),'timecells',379);
TrialBreakdownFluor_time(MD(292),379,6);

%c,d. Correlation matrix and correlation as a function of lag. 
%Change codingCells variable to 'timecells'.
TimeCorrMatrix;

%e. Sorted trial bias profiles.
trialDensityMap(MD(292));

%f. Trial decoder performance. 
FullDataTrialDecodeError;

%% Figure 4. 
%a. Fields of view. 
FOVAlignmentPlot;

%b. Time cells sized by field location across days. 
PlotMultipleTimeCellGradients;

%c. ROIs over days. 
plotROIAvgOverDays(MD(292),MD(293:295),[294 240]);

%d. Ensemble plots for one animal's recording sessions.
PastalkovaPlot(MD(292));
PastalkovaPlot(MD(293));
PastalkovaPlot(MD(294));
PastalkovaPlot(MD(295));

%e. Time cell ensemble overlaps.
ComputeTCOverlaps;

%f. Elapsed day time decoding for one session pair. 
[decodedTime,~,Mdl] = ElapsedDayTimeDecoder(MD(292),MD(293),'plotit',true);

%g. Error as a function of elapsed days. 
FullDataDecodeError_DayElapse;

%% Figure 5.
%a. Example time cells over days. 
%From left to right are cells #240, 170, 556 for G45 on 11/30/2016.
%ROIs...
OverlayNeurons;         %Edit neurons vector to include the above cells.
%Tuning curves...
msPlotTimeCells(MD(292:295),'timecells',240);
msPlotTimeCells(MD(292:295),'timecells',107);
msPlotTimeCells(MD(292:295),'timecells',556);

%b. Time cell ensemble over days. 
msPastalkovaPlot(MD(292),[MD(293:295)],10*ones(1,5),1);

%c,d. Correlation matrix and correlation as a function of lag. 
%Change codingCells variable to 'timecells'.
TimeCorrMatrix; 

%e. Proportion of cells exhibiting each across-days behavior. 
PropStabilityAll;       %Change cellType variable to 'time'.

%f. Day deocder error.
FullDataDayDecodeError;

%% Figure S1.
%a,b. Minimum projection and histology. Images saved to disk. 

%c. Time cells.
plotTimeCells(MD(292),'singletraces','LPtrdmll','timecells',240);
plotTimeCells(MD(296),'singletraces','LPtrdmll','timecells',133);
plotTimeCells(MD(300),'singletraces','LPtrdmll','timecells',68);
plotTimeCells(MD(304),'singletraces','LPtrdmll','timecells',173);

%d. Place cells.
LinearizedPF_raster(MD(292),'neurons',84,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(296),'neurons',8,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(301),'neurons',201,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(304),'neurons',8,'plotTrials',true,'noTreadmill',true);

%% Figure S2. 
%Distribution of within-session trial bias scores and overlaid control.
%Change cellType variable to 'placecells' and rasterType variable to
%'place'.
TrialSkew;

%% Figure S3. 
%a,b. Overlaid two registered sessions' neurons.
OverlayNeurons2; 

%c. Pairwise correlations between all cell ROIs between two sessions.
AllPairwiseMaskSpatialCorr(MD(292),MD(293));

%d. ROI centroid and orientation drifts. 
getRegistrationDrifts;

%e. Compare ROI displacements of stable, exiting, entering cells. Change
%cellType variable to 'time' for time cells or 'place' for place cells. 
StabilityReg2; %Something broke here.