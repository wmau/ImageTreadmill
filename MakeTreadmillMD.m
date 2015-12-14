%% Script for creating directory of where treadmill log lives. 

TreadmillMD = struct;
MasterDirectory = 'C:\MasterData';
CurrDir = pwd; 
cd(MasterDirectory); 

i=1; 

TreadmillMD(i).Animal = 'GCamp6f_45_treadmill';
TreadmillMD(i).Path = 'E:\Imaging Data\Endoscope\GCaMP6f_45\Treadmill';
TreadmillMD(i).File = fullfile('E:\Imaging Data\Endoscope\GCaMP6f_45\Treadmill','g45_treadmill.mat'); 

i=i+1;
TreadmillMD(i).Animal = 'GCamp6f_48_treadmill';
TreadmillMD(i).Path = 'E:\Imaging Data\Endoscope\GCaMP6f_45\Treadmill';
TreadmillMD(i).File = fullfile('E:\Imaging Data\Endoscope\GCaMP6f_48\Treadmill','g48_treadmill.mat'); 

save('TreadmillMD.mat','TreadmillMD'); 

cd(CurrDir); 