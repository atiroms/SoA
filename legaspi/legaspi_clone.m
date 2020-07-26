% Published: August 14, 2019
�� Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Create the noisy sensory input signals arriving at various times taoA and tao0
% Cleaners
clear % Clears all variables from the workspace
cle % clears the command window
close all % closes all figures (such as plots)
%{
Only if random generation should be repeatable
rng(1,'twister');
%}
% Simulation conditions
taoInstances = 35000; % Number of taoA and taoO instances to be generated
ExpR = 1; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Voluntary, Involuntary, Sham)
% Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Intermediate, High)
tAp=0; dist_tAt0=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
for CondBO = 1:numCond
% Initialize baseline parameters with reported empirical data
[muA, sigmaA, mu0, sigma0] = soa_IBexperiment(ExpR, CondBO);
% Generate from Gaussian distributions
Vec_taoA = normrnd(tAp + muA, sigmaA, [1 taoInstances]); % Alterntaive: sigmaA .* randn(1,numInstances) + (tAp + muA);
Vec_taoO = normrnd(tOp + mu0, sigma0, [1 taoInstances]); % sSigmaO .* randn(1,numInstances) + (tOp + muO);
%{
% Check the generated data
figure; normplot(Vec_taoA);
figure; normplot(Vec_tao0);
figure; histfit(Vec_taoA);
figure; histfit(Vec_tao0);
%}
% Generate statistics
sizeVec_taoA = numel(Vec_taoA);
sizeVec_taoO = numel(Vec_tao0);
taoA_min = min(Vec_taoA); taoA_max = max(Vec_taoA);
tao0O_min = min(Vec_tao0); taoO_max = max(Vec_tao0);
uVectaoA = mean(Vec_taoA); stdVectaoA = std(Vec_taoA);
uVectaoO = mean(Vec_tao0); stdVectaoO = std(Vec_tao0);
% Store generated simulation data
fprintf ('\n============== tao DataSet Exp %d Cond %d ================\n' ExpR, CondBO);
fprintf('taoA [%0.2f, %?@.2f] taoO [%0.2f, %@.2f]\n, taoA_min, taoA_max, tao0_min, taoO_max);
fprintf('tao statistics: %O@.2f (%@.2f)\t %@.2f (%@.2f)\n, uVectaoA, stdVectaoA, uVectao0, stdVectao0);
fprintf('taoA elements: %d tao0 elements: %d\n, sizeVec_taoA, sizeVec_tao0);
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
fnametaoO = sprintf (�eExp%dCond%d_Vec_tao0.csv',ExpR, CondBO) ;
dlmwrite(fnametaoA, Vec_taoA,;delimiter',',');
dlmwrite(fnametao0O, Vec_tao0,delimiter',',');
end





% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Find the optimal value for the free parameter muAO using Haggard et al.'s data,
% with sigmaAO equal to 1@ms.
% Cleaners
clear % Clears all variables from the workspace
clc % Clears the command window
close all % closes all figures (such as plots)
% Simulation Conditions
taoInstances = 35000; % Number of taoA and taoO instances to be generated
ExpR = 1; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
tAp=0; dist_tAt0O=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
SigmaAO = 10; % To obtain discernible perceptual shifts, sigmaAOQ should be small
for muAO = [190 200 210 220 230 240 250]
fprintf('Action and outcome perceptual shifts per condition given muA0=%d\n' muA0);
sumError = Q;
for CondBO = 1:numCond
% Read from files taoA and taoO values derived from a Gaussian distribution
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR,CondBO) ;
fnametaoO = sprintf ('Exp%dCond%d_Vec_tao0.csv',ExpR,CondBO);
Vec_taoA = dlmread(fnametaoA) ;
Vec_taoO = dlmread(fnametao0);
% Get the reported empricial baseline parameters
[muA, sigmaA, mu0, sigmaO] = soa_IBexperiment(ExpR, CondBO);
% Compute for sigma_Tot
SigmaTot2 = sigmaA*2 + sigma0*2 + sigmaA0%2;
% Compute the action and oputcome perceptual shifts
Vec_PrcShftA = (sigmaA*2 / sigmaTot2) *�� (Vec_tao0O - Vec_taoA - muAQ);
Vec_PrcShftO = ?- (sigma0*2 / sigmaTot2) * (Vec_tao0 -? Vec_taoA ?- muAQ);
uVec_PrcShftA = mean(Vec_PrcShftA); sdVec_PrcShftA = std(Vec_PrcShftA);
uVec_PrcShftO = mean(Vec_PrcShft0); sdVec_PrcShft0O = std(Vec_PrcShft0O) ;
ruVec_PrcShftA = round(uVec_PrcShftA); rsdVec_PrcShftA = round(sdVec_PrcShftA);
ruVec_PrcShftO = round(uVec_PrcShft0); rsdVec_PrcShftO = round(sdVec_PrcShft0);
% Compute for model estimation errors given the reported empirical results
[targPrcShftA, targPrcShft0] = soa_IBTargets(ExpR,CondB0) ;
errPrcShft = abs(ruVec_PrcShftA - targPrcShftA);
%6{
NOTE: Even when we consider the average of action and outcome estimation errors,
the optimal result is the same.
errPrcShft = (abs(ruVec_PrcShftA - targPrcShftA) + abs(ruVec_PrcShftO - targPrcShft0O))/2;
%}
SumError = sumError + errPrcShft;
fprintf('Condition %d:\t %0.1f(%@.2f)\t %@.1f(%@.2f)\n,; CondBO, ruVec_PrcShftA, rsdVec_PrcShftA, ruVec_PrcShft0, rsdVec_PrcShft0);
end
modelEE = sumError/numCond;
fprintf('model estimation error:\t%?.2f:\n\n,modelEE) ;
if muAO==198 || (muA0~=190 && modelEE < min_modelEE)
    
    

    
    
    
min_modelEE = modelEE;
opt_muAO = muAO;
end
end
fprintf('Optimal muAO is %d ms.\n', opt_muA0);
%{
Notes to METHODS:
? Estimates of the perceptual shift in action timing alone was sufficient to indicate
the optimal muAO.
? The optimal muAO for Experiment 1 (Haggard et al.) is 230 ms.
?- Retain this same muAO value for Experiment 2 (Wolpe et al.).
%s}



% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Fit Haggard et al's data to find the optimal value for the free parameter P(Xi=1),
% muAOQ is 230ms and sigmaAO 1s 10ms.
% Plot the action and outcome perceptual shifts given P(Xi=1) in the range
95 [@.0,1.0@] with @.1 increments
% Cleaners
clear �� Clears all variables from the workspace
cle % Clears the command window
close all % closes all figures (such as plots)
% Graph display fonts
fontsize = 20;
sizeBin = 200;
% Simulation Conditions
taoInstances = 35000; % Number of taoA and taoO instances to be generated
ExpR = 1; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
% Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=0; dist_tAt0=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
% Optimal condition-independent parameters
muAO = 230;
SigmaAO = 10;
% Interval length in consideration
T = 250; % Large enough but finite constant
% Data Matrices
LB = 0.0; INC = 0.1; UB = 1.0;
arrPXil = LB: INC:UB;
size_pXil = numel(arrPXil);
arrPrcShftA = zeros(numCond,size_pXil);
arrPrcShftO = zeros(numCond,size_pXil);
arrAOBinding = zeros(numCond,size_pXil);
for CondBO = 1:numCond
% Read from files taoA and taoO values derived from a Gaussian distribution
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
fnametaoO = sprintf ('Exp%dCond%d_Vec_tao0.csv',ExpR, CondBO);
Vec_taoA = dlmread(fnametaoA) ;
Vec_tao0O = dlmread(fnametao0) ;
indxPXil = size_pXil + 1;
for PXi_1 = UB:-INC:LB
PXi_@ = 1 - PXi_1;
% Matrices to track optimal action and outcome estimates
Vec_PrcShftA = soa_InitMatrix(1,taoInstances);
Vec_PrcShftO = soa_InitMatrix(1,taoInstances);
Vec_AOBinding = soa_InitMatrix(1,taoInstances);
for indx_tao = 1:taoInstances
% Do for each pair of taoA and tao0
taoA = Vec_taoA(indx_tao);
tao0O = Vec_tao0(indx_tao);





% Get the reported empricial baseline parameters
[muA, sigmaA, mu0, sigmaO] = soa_IBexperiment(ExpR, CondBO);
�� Compute for the posterior-ratio
Z1 = sqrt(2x*pi)*xsigmaAOxT;
Z? = T%2;
Theta = log( (PXi_1*Z0)/(PXi_@*Z1));
SigmaTot2 = sigmaA*2 + sigma0*2 + sigmaA0%2;
r = exp(Theta - ((tao0-taoA-muA0)*2/(2*sigmaTot2) ));
�� Compute for strength of temporal binding
if r>1% Causal
tAhat = taoA + (sigmaA*2/sigmaTot2)x(tao0-taoA-muA0) ;
tOhat = taoO - (sigma0*2/sigmaTot2)x*(tao0-taoA-muA0) ;
Xihat =1;
else % Acausal
tAhat = taoA;
tOhat = tao0;
Xihat=0;
end
Vec_PrcShftA(1, indx_tao) = tAhat ?- taoA;
Vec_PrcShft0(1, indx_tao) = tOhat - tao0;
Vec_AOBinding(1, indx_tao) = 25@ + (tOhat-tao0) - (tAhat-taoA);
end
uPrcShftA = mean(Vec_PrcShftA(:)); sdPrcShftA = std(Vec_PrcShftA(:));
uPrcShftO = mean(Vec_PrcShft0O(:)); sdPrcShftO = std(Vec_PrcShft0(:));
uAOBinding = mean(Vec_AOBinding(:)); sdAOBinding = std(Vec_AOBinding(:));
% Compute for model estimation errors given the reported empirical results
[targPrcShftA, targPrcShft0] = soa_IBTargets(ExpR,CondB0O) ;
ruVec_PrcShftA = round(uPrcShftA) ;
ruVec_PrcShftO = round(uPrcShft0O);
errPrcShftA = abs(ruVec_PrcShftA -? targPrcShfta);
errPrcShftO = abs(ruVec_PrcShft0O - targPrcShft0);
fprintf('Condition %d \t P(Xi=1): %@.2f\n, CondBO, PXi_1);
fprintf('uPercShfts \t %O.2f(%O.2f)\t %@.1f(%@.1f)\n, uPrcShftA, sdPrcShftA, uPrcShft0,sdPrcShft0) ;
fprintf('Error in action perceptual shift: %@.2f\n,; errPrcShftA);
forintf('Error in outcome perceptual shift: %@.2f\n\n, errPrcShft0O);
indxPXi1 = indxPXil - 1;
arrPrcShftA(CondBO, indxPXi1l) = uPrcShftA;
arrPrcShft0(CondBO, indxPXil) = uPrcShft0;
arrAOBinding(CondBO, indxPXi1) = uAOBinding;
end
fprintf('\n');
end
% Plot and store the perceptual shifts and action-outcome binding
soa_plotPrcShts(ExpR, arrPrcShftA, arrPrcShft0O, arrPXil, fontsize);
fnamePrcShft = sprintf ('Exp%d_PXisPrcShfts.png',ExpR);
saveas(gcf, fnamePrcShft) ;
soa_plotBehaviors(ExpR, arrAOBinding, arrPXil, fontsize, 1);
fnamePrcShft = sprintf ('Exp%d_PXisPrcShfts.png',ExpR);
saveas(gcf, fnamePrcShft);
*s Store the perceptual shifts




fnamePrcShftA = sprintf ('Exp%d_arrPrcShftA.csv',ExpR);
fnamePrcShftO = sprintf ('Exp%d_arrPrcShft0.csv',ExpR) ;
dlmwrite(fnamePrcShftA, arrPrcShftA, 'delimiter',',');
dlmwrite(fnamePrcShft0O, arrPrcShft0, 'delimiter',',');
%{
Notes to METHODS:
- Estimates of the perceptual shift in action timing alone was sufficient to indicate
the optimal P(Xi=1) value. However, note the following.
? Although the optimal P(Xi=1) value for the voluntary and involuntary conditions is 1.@,
the result is saturated. Hence, report P(Xi=1)=@.9 for both conditions.
? Report P(Xi=1)=0.1 for the sham condition, assuming causality is less frequently detected.
- Although the ptimal P(Xi=1) value for the intermediate tone uncertainty condition is @.5,
the outcome binding behavior is not consistent with the reported outcome binding.
Report P(Xi=1)=0.6 for the intermediate tone uncertainty condition since it also best minimized
the estimation error while reproducing the reported action and outcome bindings.
%}




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Compute the trial-to-trial temporal binding and repulsion effects,
% as well as the baseline and operant temporal bindings,
% as functions of temporal disparity, i.e., tao0-taoA
% Cleaners
clear % Clears all variables from the Workspace
clc % clears the Command Window
close all
�� Graph display fonts
fontsize = 20;
SsizeBin = 200;
% Simulation Conditions
taoInstances = 35000; % Number of taoA and taoQ instances to be generated
ExpR = 2; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
% Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=@; dist_tAt0O=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
% Optimal condition-independent parameters
muAO = 230;
SigmaAO = 10;
% Interval length in consideration
T = 250; % Large enough but finite constant
% Data Matrices
Vec_PrcShftA = zeros(numCond,taoInstances);
Vec_PrcShftO = zeros(numCond,taoInstances);
Vec_taol = zeros(numCond,taoInstances);
Vec_OpPrcShfts = zeros(numCond, taoInstances) ;
Vec_BsPrcShfts = zeros(numCond, taoInstances) ;
for CondBO = 1:numCond
% Read from files taoA and taoO values derived from a Gaussian distribution
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
fnametaoO = sprintf ('Exp%dCond%d_Vec_tao0.csv',ExpR, CondBO) ;
Vec_taoA = dlmread(fnametaoA) ;
Vec_tao0O = dlmread(fnametao0);
% Simulated using the fitted P(Xi=1) optimal values
if ExpR==1
if CondBO == 1
PXi_1 = Q.9;
elseif CondBO == 2
PXi_1 = 0.9;
else
PXi_1 = @.1;
end
elseif ExpR==2
if CondBO == 1
PXi_1 = 0.9;
elseif CondBO == 2
PXi_1 = 0.6;
else PXi_1 = Q.5;
end
end




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Compute the trial-to-trial temporal binding and repulsion effects,
% as well as the baseline and operant temporal bindings,
% as functions of temporal disparity, i.e., tao0-taoA
% Cleaners
clear % Clears all variables from the Workspace
clc % clears the Command Window
close all
�� Graph display fonts
fontsize = 20;
SsizeBin = 200;
% Simulation Conditions
taoInstances = 35000; % Number of taoA and taoQ instances to be generated
ExpR = 2; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
% Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=@; dist_tAt0O=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
% Optimal condition-independent parameters
muAO = 230;
SigmaAO = 10;
% Interval length in consideration
T = 250; % Large enough but finite constant
% Data Matrices
Vec_PrcShftA = zeros(numCond,taoInstances);
Vec_PrcShftO = zeros(numCond,taoInstances);
Vec_taol = zeros(numCond,taoInstances);
Vec_OpPrcShfts = zeros(numCond, taoInstances) ;
Vec_BsPrcShfts = zeros(numCond, taoInstances) ;
for CondBO = 1:numCond
% Read from files taoA and taoO values derived from a Gaussian distribution
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
fnametaoO = sprintf ('Exp%dCond%d_Vec_tao0.csv',ExpR, CondBO) ;
Vec_taoA = dlmread(fnametaoA) ;
Vec_tao0O = dlmread(fnametao0);
% Simulated using the fitted P(Xi=1) optimal values
if ExpR==1
if CondBO == 1
PXi_1 = Q.9;
elseif CondBO == 2
PXi_1 = 0.9;
else
PXi_1 = @.1;
end
elseif ExpR==2
if CondBO == 1
PXi_1 = 0.9;
elseif CondBO == 2
PXi_1 = 0.6;
else PXi_1 = Q.5;
end
end




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Compute the Confidence on Causal Estimate (CCE) along P(Xi=1)
% values in the range [@.0,1.@] with increments of @.1
% Cleaners
clear % clears all variables from the Workspace
cle % clears the Command Window
close all
�� Graph display fonts
fontsize = 20;
sizeBin = 200;
% Simulation Conditions
taoInstances = 35000; % Number of taoA and taoO instances to be generated
ExpR = 2; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
% Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=0; dist_tAt0O=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
% Optimal condition-independent parameters
muAQ = 230;
SigmaAO = 10;
% Interval length in consideration
T = 250; % Large enough but finite constant
% Data Matrices
LB = 0.0; INC = 0.1; UB = 1.0;
arrPXil = LB: INC:UB;
size_pXil = numel(arrPXil);
arrCCE = zeros(numCond,size_pXil);
for CondBO = 1:numCond
% Read from files taoA and taoO values derived from a Gaussian distribution
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
fnametaoO = sprintf ('Exp%dCond%d_Vec_tao0.csv',ExpR, CondBO) ;
Vec_taoA = dlmread(fnametaoA) ;
Vec_tao0 = dlmread(fnametao0);
indxPXil = size_pXil + 1;
for PXi_1 = UB:-INC:LB
PXi_@ = 1 - PXi_1;
Vec_CCE = soa_InitMatrix(1,taoInstances) ;
for indx_tao = 1:taoInstances
% Do for each pair of taoA and tao0
taoA = Vec_taoA(indx_tao);
tao0O = Vec_tao0(indx_tao);
% Get the reported empricial baseline parameters
[muA, sigmaA, mu0, sigma0] = soa_IBexperiment(ExpR, CondBO);
% Compute CCE
Z1 = sqrt(2x*pi)*sigmaAOxT;




Z? = T%2;
Theta = log((PXi_1*Z?)/(PXi_@*Z1));
SigmaTot2 = sigmaA*2 + sigmaQ*2 + sigmaAQ%2;
X = Theta - ((tao0-taoA-muA0)*2/(2*sigmaTot2)) + log(sigmaA0/sqrt(sigmaTot2) );
Vec_CCE(1,indx_tao) = (sqrt(sigmaTot2)/(2*pixsigmaA*sigma0*sigmaAQO)) x..
( 1 / (1 + exp(-X)));
end
uCCE = mean(Vec_CCE(1,:)); sdCCE = std(Vec_CCE(1,:));
fprintf('Condition %d\t P(Xi=1): %@.2f\n} CondBO, PXi_1);
fprintf('CCE �e\t %@.2e(%@.2e)\n,  uCCE, sdCCE);
indxPXil = indxPXil - 1;
arrCCE(CondBO, indxPXi1l) = uCCE;
end
end
�f*5 Plot and store CCE as function of causal prior strength
soa_plotBehaviors(ExpR, arrCCE, arrPXil, fontsize, 1);
fnameCCEPXi = sprintf ('Exp%d_CCEPXi.png',ExpR);
saveas(gcf, fnameCCEPXi) ;




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Compute the trial-to-trial Confidence on Causal Estimate (CCE)
�� Cleaners
clear �� Clears all variables from the Workspace
clc % clears the Command Window
close all
% Graph display fonts
fontsize = 20;
sizeBin = 200;
% Simulation Conditions
taoInstances = 3500Q; % Number of taoA and tao0 instances to be generated
ExpR = 1; numCond = 3; % Experimental set-up
% Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
% Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=@; dist_tAt0=250; tOp=tAp+dist_tAto; % Actual physical stimulus timings
% Optimal condition-independent parameters
muAO = 230;
SigmaAO = 10;
% Interval length in consideration
T = 250; % Large enough but finite constant
% Data Matrices
Vec_CCE = soa_InitMatrix(numCond, taoInstances);
Vec_taoI = soa_InitMatrix(numCond, taoInstances) ;
Vec_Pc = soa_InitMatrix(numCond,taoInstances) ;
for CondBO = 1:numCond
% Simulated using the fitted P(Xi=1) optimal values
if ExpR==1
if CondBO == 1
PXi_1 = Q@.9;
elseif CondBO == 2
PXi_1 = Q.9;
else
PXi_1 = Q@.1;
end
elseif ExpR==2
if CondBO == 1
PXi_1 = 0.9;
elseif CondBO == 2
PXi_1 = 0.6;
else PXi_1 = 0.5;
end
end
PXi_@ = 1 - PXi_1;
% Read from files taoA and tao0O values derived from a Gaussian distribution
fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
fnametaoO = sprintf (' Exp%dCond%d_Vec_tao0.csv',ExpR, CondBO) ;
Vec_taoA = dlmread(fnametaoA) ;
Vec_tao0 = dlmread(fnametao0);




for indx_tao = 1:taoInstances
�� Do for each pair of taoA and tao0
taoA = Vec_taoA(indx_tao);
tao0 = Vec_tao0(indx_tao);
* Get the reported empricial baseline parameters
[muA, sigmaA, mu0, sigmaO] = soa_IBexperiment(ExpR, CondBO);
��. Compute CCE
Z1 = sqrt(2x*pi)*xsigmaAOxT;
Z? = T%2;
Theta = log((PXi_1*Z@)/(PXi_0*Z1) );
SigmaTot2 = sigmaA*2 + sigma0*2 + sigmaAQ0%2;
X = Theta - ((tao0-taoA-muA0)*2/(2*sigmaTot2)) + log(sigmaA0O/saqrt(sigmaTot2) );
Vec_CCE(CondBO, indx_tao) = (sqrt(sigmaTot2)/(2*xpixsigmaAxsigma0*sigmaAO)) x..
(1 / (1 + exp(-X)));
Vec_taoI(CondBO,indx_ tao) = tao0-taoA;
end
end
*5 Plot and store trial-to-trial CCE as function of temporl disparity
sortedtaolI = Vec_taol;
[sortedtaol(1,:), sortIndx1] = sort(Vec_taoI(1,:));
[sortedtaoI(2,:), sortIndx2] = sort(Vec_taoI(2,:));
[sortedtaol(3,:), sortIndx3] = sort(Vec_taoI(3,:));
sortedCCE = soa_sortMatrices(Vec_CCE, sortIndx1, sortIndx2, sortIndx3);
soa_plotErrorBars(ExpR, sortedtaol, sortedCCE, fontsize, 1, sizeBin);
fnameCCE = sprintf ('Exp%d_perTrialCCE.png',ExpR);
saveas(gcf, fnameCCE) ;




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Return the means and standard deviations of the reported
% baseline action and outcome timing Judgment errors
function [mu_A, sigma_A, mu_0, sigma_0O] = soa_IBexperiment(experiment_case, condition)
if experiment_case == 1
% Haggard et al., 2002 (Nat Neurosci): Seminal intentional binding experiment
% Different keypress (i.e., the action) conditions
if condition == 1
mu_A = 6; sigma_A = 66;
elseif condition == 2
mu_A = 83; sSigma_A = 83;
elseif condition == 3
mu_A = 32; sigma_A = 78;
end
mu_O = 15; sigma_O = 72;
elseif experiment_case == 2
% Wolpe et al. 2013 (Exp Brain Res): Uncertainty is with the outcome
% Different tone (i.e., the outcome) conditions
mu_A = -8; sSigma_A = 75;
if condition == 1
mu_O = 35; sigma_O = 61;
elseif condition == 2
mu_O = 46; sigma_O = 66;
elseif condition == 3
mu_O = 95; sigma_O = 9@;
end
end




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Return the reported action and outcome perceptual shifts
% for the operant conditions
function [trgtPrcShftA, trgtPrcShft0] = soa_IBTargets(experiment, condition)
if experiment == 1
% Haggard et al., 2002 (Nat Neurosci): Seminal intentional binding experiment
% Different keypress (i.e., the action) conditions
if condition == 1
trgtPrcShftA = 15;
trgtPrcShftO = ?46;
elseif condition == 2
trgtPrcShftA = -27;
trgtPrcShftO = 31;
elseif condition == 3
trgtPrcShftA = -7;
trgtPrcShftO = -8;
end
elseif experiment == 2
% Wolpe et al. 2013 (Exp Brain Res): Uncertainty is with the outcome
95 Different tone (i.e., the outcome) conditions
if condition == 1
trgtPrcShftA = 39;
trgtPrcShftO = ?-51;
elseif condition == 2
trgtPrcShftA = 31;
trgtPrcShftO = -65;
elseif condition == 3
trgtPrcShftA = 32;
trgtPrcShftO = -105;
end
end




% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Instantiate a rows-by-cols matrix with zero values
function [matrix_] = soa_InitMatrix(rows_, cols_)
matrix _ = zeros(rows_, cols_);




% Published: August 14, 2019
* Copyright
% Lab for Neural Computation and Adaptation
%5 RIKEN Center for Brain Science
% Objective: Sort the cotents of the matrices
function [sortedMatrix] = soa_sortMatrices(vectMatrix, sortIndx1, sortIndx2, sortIndx3)
sortedMatrix = vectMatrix;
sortedMatrix(1,:) = vectMatrix(1,sortIndx1);
sortedMatrix(2,:) = vectMatrix(2,sortIndx2);
sortedMatrix(3,:) = vectMatrix(3,sortIndx3);




% Published: August 14, 2019

% Copyright

% Lab for Neural Computation and Adaptation

% RIKEN Center for Brain Science

% Objective: Graph the action and outcome perceptual shifts as functions of

% the strenth of the causal prior

function F = soa_plotPrcShts(experiment, arrPrcShftA, arrPrcShft0, arrPXil, fontsize)

F = figure;

linewidth = 2;

if experiment == 1
% Haggard et al., 2002 (Nat Neurosci): Seminal intentional binding experiment
% Different keypress (i.e., the action) conditions
plot(arrPXil,arrPrcShftA(1,:),;b-', arrPXil,arrPrcShft0(1,:),'b--',...

arrPXil,arrPrcShftA(2,:),;r-', arrPXil,arrPrcShft0(2,:),'r?-',��e.
arrPXil,arrPrcShftA(3,:),;k-', arrPXil,arrPrcShft0(3,:),k--', 'Linewidth', Linewidth) ;

%legend('Voluntary Action�f, ' Tone', �eInvoluntary MEP�f, ' Tone',' Sham TMS',' Tone', �eLocation', 'northwest');
lgnd = legend('Voluntary action',' and tone', �eInvoluntary action',' and tone','Sham',' and tone', 'Location', 'northwest', 'Orientation','vertical');
lgnd.FontSize = 18;
set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1; 1; 1; @.8]));

elseif experiment == 2
% Wolpe et al. 2013 (Exp Brain Res): Uncertainty is with the outcome
95 Different tone (i.e., the outcome) conditions
hold on;
plot(arrPXil,arrPrcShftA(1,:),;Color', [@ @ 250/255], 'LineStyle','-', 'Linewidth', linewidth);
plot(arrPXil,arrPrcShft0(1,:);Color', [@ @ 250/255], 'LineStyle', '--', 'Linewidth', linewidth) ;
plot(arrPXil,arrPrcShftA(2,:),;Color', [@ 140/255 255/255], 'LineStyle','-', 'Linewidth', Linewidth);
plot(arrPXil,arrPrcShft0(2,:),;Color', [@ 140/255 255/255], 'LineStyle', '--', 'Linewidth', Linewidth);
plot(arrPXil,arrPrcShftA(3,:),;Color', [@ 240/255 255/255], 'LineStyle','-', 'Linewidth', linewidth);
plot(arrPXil,arrPrcShft0(3,:);Color', [@ 240/255 255/255], 'LineStyle', '--', 'Linewidth', linewidth) ;
lgnd = legend('Action', ' and low uncertainty tone�f, 'Action', ' and intermediate uncertainty tone, 'Action', ' and high uncertainty tone�f, 'Location', 'northwest');
lgnd.FontSize = 18;
set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1; 1; 1; @.8]));

end

set(gca,'FontSize', fontsize);

set(gca, 'Box', '�eon');

set(lgnd, 'Color', 'none');

set(gca,'color', 'white');

set(gca,'FontSize', fontsize);




% Published: August 14, 2019
* Copyright
% Lab for Neural Computation and Adaptation
%5 RIKEN Center for Brain Science
% Objective: Plot the optimal baheviors used in the figures of the paper with ERROR BARS displayed
function F = soa_plotErrorBars(experiment, arrAxes, arrBehavior, fontsize, flag, sizeBin)
F = figure;
markerSize = 27;
LineWidth = @.1;
if experiment == 1
% Haggard et al., 2002 (Nat Neurosci): Seminal intentional binding experiment
% Different keypress (i.e., the action) conditions
hold all
if flag == 1
errorbar(mean( reshape(arrAxes(1,:),sizeBin, []),1), mean(reshape(arrBehavior(1,:),sizeBin,[]),1), std(reshape(arrBehavior(1,:),sizeBin,[]),2b,', 'LineWidth', LineWidth,'MarkerSize', \
markerSize);
errorbar(mean( reshape(arrAxes(2,:),SizeBin, []),1), mean(reshape(arrBehavior(2,:),sizeBin,[]),1), std(reshape(arrBehavior(2,:),sizeBin,[]),1y,', 'LineWidth', LineWidth,'MarkerSize', &
markerSize);
errorbar(mean(reshape(arrAxes(3,:),SsizeBin, []),1), mean(reshape(arrBehavior(3,:),SizeBin,[]),1), std(reshape(arrBehavior(3,:),sizeBin, []),2k,', 'LineWidth', LineWidth,'MarkerSize', 2
markerSize);
lgnd = legend('Voluntary condition','Involuntary condition','Sham condition', 'Location', 'northwest');
elseif flag == 2
errorbar(mean(reshape(arrAxes(2,:),SsizeBin, []),1), mean(reshape(arrBehavior(2,:),sizeBin,[]),1), std(reshape(arrBehavior(2,:),sizeBin,[]),1I},', '�eLineWidth', LineWidth,'MarkerSize', /
markerSize) ;
errorbar(mean( reshape(arrAxes(1,:),sizeBin, []),1), mean(reshape(arrBehavior(1,:),sizeBin,[]),1), std(reshape(arrBehavior(1,:),sizeBin,[]),2b,', 'LineWidth', LineWidth,'MarkerSize', 4
markerSize);
lgnd = legend('Voluntary condition','Involuntary condition', 'Location', 'northwest');
end
set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1; 1; 1; @.8]));
set(gca, 'Box', �eon');
hold off
elseif experiment == 2
% Wolpe et al. 2013 (Exp Brain Res): Uncertainty is with the outcome
% Different tone (i.e., the outcome) conditions
hold all
if flag == 1
errorbar(mean( reshape(arrAxes(1,:),sizeBin, []),1), mean(reshape(arrBehavior(1,:),sizeBin,[]),1), std(reshape(arrBehavior(1,:),sizeBin,[]),1Xplor', [@ @ 250/255] uv
'LineStyle','none', 'Marker', '.', 'LineWidth', LineWidth,'MarkerSize', markerSize);
errorbar(mean( reshape(arrAxes(2,:),SizeBin, []),1), mean(reshape(arrBehavior(2,:),sizeBin,[]),1), std(reshape(arrBehavior(2,:),sizeBin,[]),1fplor', [@ 140/255 255/255] �h
'LineStyle','none', 'Marker', '.', 'LineWidth', LineWidth, 'MarkerSize',markerSize) ;
errorbar(mean( reshape(arrAxes(3,:),sizeBin, []),1), mean(reshape(arrBehavior(3,:),sizeBin,[]),1), std(reshape(arrBehavior(3,:),sizeBin,[]),2X��plor', [0 240/255 255/255] �h@
'LineStyle','none', 'Marker', '.', 'LineWidth', LineWidth,'MarkerSize',markerSize) ;
end
lgnd = legend('Low uncertainty', �eIntermediate uncertainty�f, 'High uncertainty', 'Location', 'northwest');
set(lgnd.BoxFace, 'ColorType', �etruecoloralpha', 'ColorData', uint8(255*[1; 1; 1; 0.8]));
set(gca, 'Box', 'on');
hold off




end
set(gca,'FontSize', fontsize);





% Function to graph the SoA related measures
% Added 09/06/2017
function F = soa_plotBehaviors(experiment, arrBehavior, arrPXil, fontsize, behavior)
F = figure;
linewidth = 2;
if experiment == 1
if behavior == 1
plot( arrPXil,arrBehavior(1,:)};b', arrPXil,arrBehavior(2,:),'r', arrPXil,arrBehavior(3,:),'k', �eLinewidth', linewidth);
%{
hold all
plot(arrPXil,arrBehavior(1,:),'b', 'Linewidth',3);
plot(arrPXil, arrBehavior(2,:),'r', 'Linewidth',3);
plot(arrPXil,arrBehavior(3,:),'Color', [0 @ @]+0.05*13, 'Linewidth',3);
hold off
%}
elseif behavior == 3
plot( arrPXil,arrBehavior(2,:),r', arrPXil,arrBehavior(1,:),'b', arrPXil,arrBehavior(3,:),'k', �eLinewidth', linewidth);
elseif behavior == 0
plot( arrPXil,arrBehavior(1,:)}b', arrPXil,arrBehavior(2,:),'r', 'Linewidth', linewidth);
elseif behavior == 2
ylim([@.@ 1.0]);
plot( arrPXil,arrBehavior(1,:);b', arrPXil,arrBehavior(2,:),'r', 'Linewidth', linewidth) ;
end
elseif experiment == 2
hold on;
plot(arrPXil,arrBehavior(1,:);Color', [@ @ 250/255], 'LineStyle','-', 'Linewidth', Linewidth) ;
plot(arrPXil,arrBehavior(2,:),Color', [@ 140/255 255/255], 'LineStyle','-', 'Linewidth', Linewidth);
plot(arrPXil,arrBehavior(3,:),Color', [@ 240/255 255/255], 'LineStyle','-', 'Linewidth', linewidth);
hold off;
elseif experiment == 3
hold on;
plot(arrPXil,arrBehavior(1,:),;Color', [@ @ 250/255], 'LineStyle','-', 'Linewidth', linewidth) ;
plot(arrPXil,arrBehavior(2,:);Color', [@ 140/255 255/255], 'LineStyle','-', 'Linewidth', linewidth) ;
hold off;
end
%{
Xlabel('P(\xi=1) of Prior');
if behavior == 1
ylabel('Feeling of Agency');
elseif behavior == 2
ylabel('Judgment of Agency');
elseif behavior == 3
ylabel({'Bias in Action Estimates');
elseif behavior == 4
ylabel('Bias in Outcome Estimates');
end
%}
if experiment == 1
if behavior == 1
lgnd = legend(�eVoluntary condition','Involuntary condition','Sham condition', 'Location', 'northwest');
elseif behavior ==3
lgnd = legend('Voluntary condition','Involuntary condition','Sham condition', 'Location', 'northwest');
else





lgnd = legend(�eVoluntary condition','Involuntary condition','Sham condition', 'Location', 'northwest');
end
elseif experiment == 2
lgnd = legend('Low uncertainty condition, �eIntermediate uncertainty condition, �eHigh uncertainty condition, 'Location', 'northeast');
elseif experiment == 3
lgnd = legend('Active, Instructed', 'Passive, Instructed', 'Location', 'northeast');
end
set(gca,'FontSize', fontsize);
set(gca, 'Box', 'on');
%set(lgnd, 'Color', 'none');
set(gca,'color', 'white');
set(lgnd.BoxFace, 'ColorType', �etruecoloralpha', 'ColorData', uint8(255*(1; 1; 1; @.8]));
% If you want to bold
%ylabel('Feeling of Agency', 'FontWeight', �ebold');