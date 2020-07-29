% Published: August 14, 2019
% Copyright
%   Lab for Neural Computation and Adaptation
%   RIKEN Center for Brain Science
%
% Objective: Fit Haggard et al's data to find the optimal value for the free parameter P(Xi=1),
%               muAO is 230ms and sigmaAO 1s 10ms.
%            Plot the action and outcome perceptual shifts given P(Xi=1) in the range
%               [0.0,1.0] with 0.1 increments

% Cleaners
clear           % clears all variables from the workspace
clc             % clears the command window
close all       % closes all figures (such as plots)

% Graph display fonts
fontsize = 20;
sizeBin = 200;

% Simulation Conditions
taoInstances = 35000;                       % Number of taoA and taoO instances to be generated
ExpR = 1; numCond = 3;                      % Experimental set-up
                                            %   Haggard et al. (2002): ExpR = 1; NumCond = 3; (Vol, Invol, Sham)
                                            %   Wolpe et al. (2013) : ExpR = 2; NumCond = 3; (Low, Int, High)
tAp=0; dist_tAtO=250; tOp=tAp+dist_tAtO;    % Actual physical stimulus timings

% Optimal condition-independent parameters
muAO = 230;
sigmaAO = 10;

% Interval length in consideration
T = 250;                                    % Large enough but finite constant

% Data Matrices
LB = 0.0; INC = 0.1; UB = 1.0;
% LB = 0.2; INC = 0.1; UB = 0.6;
arrPXi1 = LB: INC:UB;
size_pXi1 = numel(arrPXi1);

% matrices for output
arrPrcShftA = zeros(numCond,size_pXi1);
arrPrcShftO = zeros(numCond,size_pXi1);
arrAOBinding = zeros(numCond,size_pXi1);

% Visualization params added by SM
hm_inc = 10;


for CondBO = 1:numCond
    % Read from files taoA and taoO values derived from a Gaussian distribution
    % These are tau A/O (noisy sensory inputs) in the article
    % Used in Monte Carlo simulation to update prior distribution of
    % actual action and outcome timings modeled as vetroloquism effect
    % CondBO = 1 (voluntary) is action-only baseline condition
    % outcome prior distribution uses outcome-only baselinecondition
    fnametaoA = sprintf ('Exp%dCond%d_Vec_taoA.csv',ExpR, CondBO) ;
    fnametaoO = sprintf ('Exp%dCond%d_Vec_taoO.csv',ExpR, CondBO);
    Vec_taoA = dlmread(fnametaoA) ;
    Vec_taoO = dlmread(fnametaoO) ;
    
    % added by SM
    % Monte Carlo results of posterior distribution
    Mat_PrcShftA = soa_InitMatrix(size_pXi1,taoInstances);
    Mat_PrcShftO = soa_InitMatrix(size_pXi1,taoInstances);
    Mat_AOBinding = soa_InitMatrix(size_pXi1,taoInstances);
    
    indxPXi1 = size_pXi1;
    
    % Loop over prior distribution of Xi (prior for causality)
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
            [muA, sigmaA, muO, sigmaO] = soa_IBexperiment(ExpR, CondBO);
            
            % Compute for the posterior-ratio
            % these calculation appear in Methods
            Z1 = sqrt(2*pi)*sigmaAO*T;
            Z0 = T^2;
            Theta = log((PXi_1*Z0)/(PXi_0*Z1));
            sigmaTot2 = sigmaA^2 + sigmaO^2 + sigmaAO^2;
            r = exp(Theta - ((taoO-taoA-muAO)^2/(2*sigmaTot2)));
            
            % Compute for strength of temporal binding
            % Here, MAP estimates (tAhat, tOhat, Xihat) are calculated
            % in classic Bayesian inference, Xi posterior is expressed as
            % probability distribution as in prior, but binary MAP estimate
            % is used instead
            if r>1      % Causal
                tAhat = taoA + (sigmaA^2/sigmaTot2)*(taoO-taoA-muAO) ;
                tOhat = taoO - (sigmaO^2/sigmaTot2)*(taoO-taoA-muAO) ;
                Xihat =1;
            else        % Acausal
                tAhat = taoA;
                tOhat = taoO;
                Xihat=0;
            end
            
            %{
            Vec_PrcShftA(1, indx_tao) = tAhat - taoA;
            Vec_PrcShftO(1, indx_tao) = tOhat - taoO;
            Vec_AOBinding(1, indx_tao) = 250 + (tOhat-taoO) - (tAhat-taoA);
            %}
            
            % Added by SM
            % Record Mote Carlo results of posterior in perceived action
            % and outcome timing results
            Mat_PrcShftA(indxPXi1,indx_tao) = tAhat - taoA;
            Mat_PrcShftO(indxPXi1,indx_tao) = tOhat - taoO;
            Mat_AOBinding(indxPXi1,indx_tao) = 250 + (tOhat-taoO) - (tAhat-taoA);
        end
        
        %{
        % Added by SM
        if CondBO==1
            figure(indxPXi1)
            % histogram(Vec_AOBinding,30)
            histogram(Vec_PrcShftA,30)
        end
        %}
        
        %{
        uPrcShftA = mean(Vec_PrcShftA(:)); sdPrcShftA = std(Vec_PrcShftA(:));
        uPrcShftO = mean(Vec_PrcShftO(:)); sdPrcShftO = std(Vec_PrcShftO(:));
        uAOBinding = mean(Vec_AOBinding(:)); sdAOBinding = std(Vec_AOBinding(:));
        %}
        
        % Added by SM
        % Mean and SD of perceived shifts in action and outcome timings (means are used for article figures)
        uPrcShftA = mean(Mat_PrcShftA(indxPXi1,:)); sdPrcShftA = std(Mat_PrcShftA(indxPXi1,:));
        uPrcShftO = mean(Mat_PrcShftO(indxPXi1,:)); sdPrcShftO = std(Mat_PrcShftO(indxPXi1,:));
        uAOBinding = mean(Mat_AOBinding(indxPXi1,:)); sdAOBinding = std(Mat_AOBinding(indxPXi1,:));
        
        % Update matrices for graphical output
        arrPrcShftA(CondBO, indxPXi1) = uPrcShftA;
        arrPrcShftO(CondBO, indxPXi1) = uPrcShftO;
        arrAOBinding(CondBO, indxPXi1) = uAOBinding;
        
        % Compute for model estimation errors given the reported empirical results
        [targPrcShftA, targPrcShftO] = soa_IBTargets(ExpR,CondBO) ;

        % Calculation for command window output
        ruVec_PrcShftA = round(uPrcShftA) ;
        ruVec_PrcShftO = round(uPrcShftO);
        errPrcShftA = abs(ruVec_PrcShftA - targPrcShftA);
        errPrcShftO = abs(ruVec_PrcShftO - targPrcShftO);
        
        fprintf('Condition %d \t P(Xi=1): %0.2f\n', CondBO, PXi_1);
        fprintf('uPercShfts \t %O.2f(%O.2f)\t %0.1f(%0.1f)\n', uPrcShftA, sdPrcShftA, uPrcShftO,sdPrcShftO) ;
        fprintf('Error in action perceptual shift: %0.2f\n', errPrcShftA);
        fprintf('Error in outcome perceptual shift: %0.2f\n\n', errPrcShftO);

        % indxPXi1 = indxPXi1 - 1;


        
        indxPXi1 = indxPXi1 - 1;
    end % end of looping over P(Xi)
    
    % Added by SM
    % Bin posterior distributions of perceived action/outcome timing shifts
    % and A-O binding
    
    
fprintf('\n');
end % end of looping over experiment conditions (voluntary, etc...)

% Plot and store the perceptual shifts and action-outcome binding

soa_plotPrcShts(ExpR, arrPrcShftA, arrPrcShftO, arrPXi1, fontsize);
fnamePrcShft = sprintf ('Exp%d_PXisPrcShfts1.png',ExpR);
saveas(gcf, fnamePrcShft) ;
soa_plotBehaviors(ExpR, arrAOBinding, arrPXi1, fontsize, 1);
fnamePrcShft = sprintf ('Exp%d_PXisPrcShfts2.png',ExpR);
saveas(gcf, fnamePrcShft);

% Store the perceptual shifts

fnamePrcShftA = sprintf ('Exp%d_arrPrcShftA.csv',ExpR);
fnamePrcShftO = sprintf ('Exp%d_arrPrcShftO.csv',ExpR) ;
dlmwrite(fnamePrcShftA, arrPrcShftA, 'delimiter',',');
dlmwrite(fnamePrcShftO, arrPrcShftO, 'delimiter',',');

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