loadMD;
session = MD(292);
cd(session.Location);
load('Pos_align.mat', 'RawTrace');
TCs = getTimeCells(session);
%PSAbool = RawTrace(TCs,:);
PSAbool = RawTrace;
PSAbool(PSAbool<0) = 0;

splitN = floor(size(PSAbool,2)*.75);
trainNEURAL = PSAbool(:,1:splitN);
testNEURAL = PSAbool(:,(splitN+1):end);

X = trainNEURAL;
K = 5;
L = 200; 
lambda = 0.005;

[W, H, cost, loadings, power] = seqNMF(X, 'K', K, 'L', L, ...
    'maxiter',50,'lambda',lambda);

p = 0.05;

[p_vals, is_significant] = test_significance(testNEURAL,W,p);