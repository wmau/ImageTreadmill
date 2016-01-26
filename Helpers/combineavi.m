function combineavi
% COMBINEAVI lets you select two avi movies of the same length and height
% and makes a new movie showing both movies side by side. 
%
% This is just a short snipplet I use to show video recordings made from
% 2 different perspectives simultaniously. Enhance to your liking.
%
%  
% 2 comments for Matlab versions below 7 and older Linux mplayer versions:
%  
% COMMENT 1:
%   
% comment for older Linux (older than 2005; i.e. SuSE8.2): mplayer has a bug
% that will not let you play matlab.avi file. You can either use xanim/xine
% to play the movie or trick mplayer to play the movie by using stdin:
%
%    cat movie.avi | mplayer -
%
% This tip was proposed by Herbert Ramoser in Message-ID:
% <slrnbsk2tg.1m9.jlb17@chaos.egr.duke.edu>
%  
%
% COMMENT 2:
%  
% Matlab version below 7 might need this:
%  
% this is needed due to a bug in the older versions of the `aviinfo'
% function:
%  
%  warning off MATLAB:mir_warning_variable_used_as_function
%
%
% END COMMENTS
%
%  
% Patrick Drechsler <patrick.drechsler at gmx dot net>
% 2005-06-03  
%
  
  % select two files:
  [filename1,pathname1] = uigetfile('.avi','pick first AVI file');
  [filename2,pathname2] = uigetfile('.avi','pick second AVI file');
  file1 = fullfile(pathname1,filename1);
  file2 = fullfile(pathname2,filename2);  
  pdMovie1 = aviread(file1);
  pdMovie2 = aviread(file2);
  fileinfo1 = aviinfo(file1);
  fileinfo2 = aviinfo(file2);
  
  % check if AVI files have the same length and height:
  if fileinfo1.NumFrames~=fileinfo2.NumFrames || ...
        fileinfo1.Height~=fileinfo2.Height
    errordlg('files are not compatible!')
  else
    % inspired by Herbert Ramoser in Message-ID:
    % <art0c0$l9fip$1@ID-148798.news.dfncis.de>
    for i=1:size(pdMovie1,2)
      output(i).cdata = [pdMovie1(i).cdata, pdMovie2(i).cdata];
      output(i).colormap = pdMovie1(i).colormap;
    end;
    
    % name of the new avi file:
    [pathstr,name,ext,versn] = fileparts(filename1);
    newmoviename = [pathname1,name,'_combined', ...
                    num2str(fileinfo1.FramesPerSecond;),ext];

    % create the avi file:
    movie2avi(output, newmoviename, ...
              'fps', fileinfo1.FramesPerSecond;, ...
              'compression', 'none');
    close
  end
