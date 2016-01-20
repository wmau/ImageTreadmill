for i=1:length(session)
    cd(fullfile(session(i).Location,'ICmovie_smoothed-Objects'));
    movie = dir('*.h5');
    infile = fullfile(pwd,movie.name);     
    cd ..
    Tenaspis(infile,'animal_id',session(i).Animal,'sess_date',session(i).Date,'sess_num',session(i).Session);
end