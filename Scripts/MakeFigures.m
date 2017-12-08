%% First things first. Load all the session entries. 
loadMD;

%% Figure 1. 
%a. Schematic of maze. Not coded.

%b,c. Maximum projection and traces. 
GenerateTraces;

%d. Time cell rasters.
plotTimeCells(MD(292),'singletraces','DFDTTrdmll','timecells',906);
plotTimeCells(MD(292),'singletraces','DFDTTrdmll','timecells',431);

%e. Time cell ensemble.
ensemble = PastalkovaPlot(MD(292)); 
PlopLogLineonPastalkovaPlot;

%f. Place cell rasters. 
LinearizedPF_raster(MD(292),'neurons',152,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(292),'neurons',190,'plotTrials',true,'noTreadmill',true);

%g. Place cell ensemble.
LinearizedPFs_treadmill(MD(292));

%% Figure 2. 
%a. Time and place raster of a joint-coding cell.
TempSpatRasters(MD(293),664);

%b. Correlate spatial and temporal info for joint-coding cells. 
corrInfo(MD(292:309),'dual');

%c. Make model spatial curve on treadmill.
DecodeTimeWithSpatialLookup(MD(293),'neurons',664,...
    'plotit',true,'neuralActivity','transient');

%d. Make model temporal curve on track. 
DecodePlaceWithTemporalLookup(MD(293),'neurons',664,...
    'plotit',true,'neuralActivity','transient');

%e. Correlate DS. 
DS_Corr; 

%% Figure 3. 
%a. Within-session exiting time cell.
plotTimeCells(MD(292),'timecells',585);
TrialBreakdownFluor_time(MD(292),585,6);

%b. Within-session entering time cell.
plotTimeCells(MD(292),'timecells',379);
TrialBreakdownFluor_time(MD(292),379,6);

%c. Distribution of within-session trial bias scores and overlaid control.
%Change cellType variable to 'timecells' and rasterType varaible to 'time'.
TrialSkew;

%d,e. Correlation matrix and correlation as a function of lag. 
%Change codingCells variable to 'timecells'.
TimeCorrMatrix; 

%% Figure 4. 
%a. Within-session exiting place cell.
LinearizedPF_raster(MD(292),'neurons',279,'plotTrials',false,'noTreadmill',true);
TrialBreakdownFluor_place(MD(292),279,6);

%b. Within-session entering place cell.
LinearizedPF_raster(MD(292),'neurons',890,'plotTrials',false,'noTreadmill',true);
TrialBreakdownFluor_place(MD(292),890,6);

%c. Distribution of within-session trial bias scores and overlaid control.
%Change cellType variable to 'placecells' and rasterType variable to
%'place'.
TrialSkew;

%d,e. Correlation matrix and correlation as a function of lag. 
%Change codingCells variable to 'placecells'.
TimeCorrMatrix; 

%% Figure 5.
%a. Joint-coding cell with different temporal and spatial onsets.
TempSpatRasters(MD(292),895);

%b. Correlate bias scores for track and treadmill. 
PvT_TrialSkew;

%% Figure 6. 
%a. Example time cells over days. 
%From left to right are cells #240, 170, 556 for G45 on 11/30/2016.
%ROIs...
OverlayNeurons;         %Edit neurons vector to include the above cells.
%Tuning curves...
msPlotTimeCells(MD(292:295),'timecells',240);
msPlotTimeCells(MD(292:295),'timecells',107);
msPlotTimeCells(MD(292:295),'timecells',556);

%b. Proportion of cells exhibiting each across-days behavior. 
PropStabilityAll;       %Change cellType variable to 'time'.

%c. Time cell ensemble over days. 
msPastalkovaPlot(MD(292),[MD(293:295)],10*ones(1,5),1);

%d. Correlation matrix and correlation as a function of lag. 
%Change codingCells variable to 'timecells'.
TimeCorrMatrix; 

%% Figure 7. 
%a. Example time cells over days. 
%From left to right are cells #152, 357, 58 for G45 on 11/30/2016.
%ROIs...
OverlayNeurons;         %Edit neurons vector to include the above cells.
%Tuning curves...
msLinearizedPF_raster(MD(292:295),'neurons',152);
msLinearizedPF_raster(MD(292:295),'neurons',357);
msLinearizedPF_raster(MD(292:295),'neurons',58);

%b. Proportion of cells exhibiting each across-days behavior. 
PropStabilityAll;       %Change cellType variable to 'place'.

%c. Place cell ensemble over days.
msLinearizedPFs(MD(292),MD(293:295));

%d. Correlation matrix and correlation as a function of lag. 
%Change codingCells variable to 'placecells'.
TimeCorrMatrix; 

%% Figure 8. 
%a. Joint-coding cell with different temporal and spatial stabilities. 
OverlaNeurons;      %Edit neurons vector to be 145.
msPlotTimeCells(MD(292:295),'timecells',145)
msLinearizedPF_raster(MD(292:295),'neurons',145);

%b,d. Information as a function of stability. 
InfoStabilityANOVA;

%c,e. SVM. 
SVMperformance; 

%% Supplementary Figure 1. 
%a,b. Minimum projection and histology. Images saved to disk. 

%c. Example time cells. 
%From left to right are cells from: G45, G48, Bellatrix, Polaris. 
plotTimeCells(MD(292),'singletraces','DFDTTrdmll','timecells',240);
plotTimeCells(MD(296),'singletraces','DFDTTrdmll','timecells',133);
plotTimeCells(MD(300),'singletraces','DFDTTrdmll','timecells',68);
plotTimeCells(MD(305),'singletraces','DFDTTrdmll','timecells',66);

%d. Example place cells. 
LinearizedPF_raster(MD(292),'neurons',84,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(296),'neurons',8,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(301),'neurons',201,'plotTrials',true,'noTreadmill',true);
LinearizedPF_raster(MD(305),'neurons',372,'plotTrials',true,'noTreadmill',true);

%% Suplementary Figure 2. 
%a. Scalar fit for time cells.
scale_inv_distr_time_cells;

%b. Scalar fit for place cells.
scale_inv_distr_place_cells;

%% Supplementary Figure 3. 
%a. Time cell gradient. 
plotTimeCellGradient(MD(292));

%b. Place cell gradient. 
plotPlaceCellGradient(MD(292)); 

%c. Time cell topology.
TimeCellTopology(MD(292));

%d. Place cell topology.
PlaceCellTopology(MD(292));

%e. Place vs time cell topology vs random
PvT_Topology2;

%% Supplementary Figure 4.
%Check for bleaching.
bleachANOVA;

%% Supplementary Figure 5.
%a,b. Overlaid two registered sessions' neurons.
OverlayNeurons2; 

%c. Closeups of overlaid neurons over 4 days. I've lost what exact cells
%were in the figure, but any random 3 should be fine. 
OverlayNeurons; 

%d. ROI centroid and orientation drifts. 
getRegistrationDrifts;

%e. Compare ROI displacements of stable, exiting, entering cells. Change
%cellType variable to 'time' for time cells or 'place' for place cells. 
CellStabilityStatus; 