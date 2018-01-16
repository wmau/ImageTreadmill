loadMD;

day1 = imread(fullfile(MD(292).Location,'ICmovie_min_proj.tif'));
day2 = imread(fullfile(MD(293).Location,'ICmovie_min_proj.tif'));
regInfo = image_registerX(MD(292).Animal,MD(292).Date,MD(292).Session,MD(293).Date,...
    MD(293).Session,'suppress_output',true);
day2 = imwarp_quick(day2,regInfo);
figure; imshow(day1,[]);
figure; imshow(day2,[]);
figure; imshow(day1-day2,[]);
