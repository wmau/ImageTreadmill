%% Script for creating directory of where treadmill log lives. 

TreadmillMD = struct;
MasterDirectory = 'C:\MasterData';
CurrDir = pwd; 
cd(MasterDirectory); 

i=1; 

TreadmillMD(i).Animal = 'GCamp6f_45_treadmil';
TreadmillMD(i).Path = 'E:\Imaging Data\Endoscope\GCaMP6f_45\Treadmill';
TreadmillMD(i).File = fullfile(Treadmill(i).Path,'g45_treadmill.mat'); 

TreadmillMD(i).Animal = 'GCamp6f_48_treadmil';
TreadmillMD(i).Path = 'E:\Imaging Data\Endoscope\GCaMP6f_45\Treadmill';
TreadmillMD(i).File = fullfile(Treadmill(i).Path,'g48_treadmill.mat'); 

save('TreadmillMD.mat','TreadmillMD'); 

cd(CurrDir); 