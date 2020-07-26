% Published: August 14, 2019
% Copyright
%   Lab for Neural Computation and Adaptation
%   RIKEN Center for Brain Science
%
% Objective: Fit Haggard et al's data to find the optimal value for the free parameter P(Xi=1),
%               muAO is 230ms and sigmaAO 1s 10ms.
%            Plot the action and outcome perceptual shifts given P(Xi=1) in the range
%               [@.0,1.0@] with @.1 increments

% Cleaners
clear           % clears all variables from the workspace
cle             % clears the command window
close all       % closes all figures (such as plots)

% Graph display fonts
fontsize = 20;
sizeBin = 200;

% Simulation Conditions
taoInstances = 35000;                       % Number of taoA and taoO instances to be generated
ExpR = 1; numCond = 3;                      % Experimental set-up
                                            %   Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
                                            %   Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=0; dist_tAt0=250; tOp=tAp+dist_tAto;    % Actual physical stimulus timings

% Optimal condition-independent parameters
muAO = 230;
SigmaAO = 10;

% Interval length in consideration
T = 250;                                    % Large enough but finite constant

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
    fnametaoO = sprintf ('Exp%dCond%d_Vec_taoO.csv',ExpR, CondBO);
    Vec_taoA = dlmread(fnametaoA) ;
    Vec_taoO = dlmread(fnametaoO) ;
    
    indxPXi1 = size_pXi1 + 1;
    for PXi_1 = UB:-INC:LB
        PXi_0 = 1 - PXi_1;
        
        % Matrices to track optimal action and outcome estimates
        Vec_PrcShftA = soa_InitMatrix(1,taoInstances);
        Vec_PrcShftO = soa_InitMatrix(1,taoInstances);
        Vec_AOBinding = soa_InitMatrix(1,taoInstances);
        
        for indx_tao = 1:taoInstances
            % Do for each pair of taoA and taoO
            taoA = Vec_taoA(indx_tao);
            taoO = Vec_taoO(indx_tao);

            % Get the reported empricial baseline parameters
            [muA, sigmaA, mu0, sigmaO] = soa_IBexperiment(ExpR, CondBO);
            
            % Compute for the posterior-ratio
            Z1 = sqrt(2*pi)*sigmaAO*T;
            Z0 = T^2;
            Theta = log( (PXi_1*Z0)/(PXi_0*Z1));
            SigmaTot2 = sigmaA^2 + sigmaO^2 + sigmaAO^2;
            r = exp(Theta - ((taoO-taoA-muA0)*2/(2*sigmaTot2) ));
            
            % Compute for strength of temporal binding
            if r>1      % Causal
                tAhat = taoA + (sigmaA^2/sigmaTot2)*(taoO-taoA-muAO) ;
                tOhat = taoO - (sigma0^2/sigmaTot2)*(taoO-taoA-muAO) ;
                Xihat =1;
            else        % Acausal
                tAhat = taoA;
                tOhat = taoO;
                Xihat=0;
            end
            
            Vec_PrcShftA(1, indx_tao) = tAhat - taoA;
            Vec_PrcShftO(1, indx_tao) = tOhat - taoO;
            Vec_AOBinding(1, indx_tao) = 250 + (tOhat-taoO) - (tAhat-taoA);
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
        
        fprintf('Condition %d \t P(Xi=1): %@.2f\n', CondBO, PXi_1);
        fprintf('uPercShfts \t %O.2f(%O.2f)\t %0.1f(%0.1f)\n', uPrcShftA, sdPrcShftA, uPrcShft0,sdPrcShft0) ;
        fprintf('Error in action perceptual shift: %0.2f\n', errPrcShftA);
        forintf('Error in outcome perceptual shift: %0.2f\n\n', errPrcShftO);

        indxPXi1 = indxPXi1 - 1;

        arrPrcShftA(CondBO, indxPXi1) = uPrcShftA;
        arrPrcShftO(CondBO, indxPXi1) = uPrcShft0;
        arrAOBinding(CondBO, indxPXi1) = uAOBinding;
    end
fprintf('\n');
end

% Plot and store the perceptual shifts and action-outcome binding

soa_plotPrcShts(ExpR, arrPrcShftA, arrPrcShftO, arrPXi1, fontsize);
fnamePrcShft = sprintf ('Exp%d_PXisPrcShfts.png',ExpR);
saveas(gcf, fnamePrcShft) ;
soa_plotBehaviors(ExpR, arrAOBinding, arrPXil, fontsize, 1);
fnamePrcShft = sprintf ('Exp%d_PXisPrcShfts.png',ExpR);
saveas(gcf, fnamePrcShft);

% Store the perceptual shifts

fnamePrcShftA = sprintf ('Exp%d_arrPrcShftA.csv',ExpR);
fnamePrcShftO = sprintf ('Exp%d_arrPrcShft0.csv',ExpR) ;
dlmwrite(fnamePrcShftA, arrPrcShftA, 'delimiter',',');
dlmwrite(fnamePrcShft0O, arrPrcShftO, 'delimiter',',');

%{
    Notes to METHODS:
    - Estimates of the perceptual shift in action timing alone was sufficient to indicate
        the optimal P(Xi=1) value. However, note the following.
    - Although the optimal P(Xi=1) value for the voluntary and involuntary conditions is 1.@,
        the result is saturated. Hence, report P(Xi=1)=@.9 for both conditions.
    - Report P(Xi=1)=0.1 for the sham condition, assuming causality is less frequently detected.
    - Although the ptimal P(Xi=1) value for the intermediate tone uncertainty condition is @.5,
        the outcome binding behavior is not consistent with the reported outcome binding.
        Report P(Xi=1)=0.6 for the intermediate tone uncertainty condition since it also best minimized
        the estimation error while reproducing the reported action and outcome bindings.
%}