function [] = LearningSys4_13_2018Script(varargin)
%% Script used to generate thesis data. Works in conjunction with LearningSysPlotScript.m to generate figures.

%% Initialization
addpath('F:\\Dropbox\\Science\\2017 paper\\')
close all
set(0,'DefaultAxesFontSize',24); 
saveBGdata = 1;
loadBGdata = 0;
rng(137)

%% Set up options

% First argument determines which experiments to run
if (nargin() >= 1)
    plotToggle = varargin{1};
else
      plotToggle =  [0 0 0 0 0  ...
                 0 0 0 0 0 ...
                 0 0 0 0 0 ...
                 0 0 0 1];
end

if (nargin() >= 2)
    fpath = varargin{2};
else
    fpath = 'F:\\Dropbox\\Science\\Figures\\';
end

if (nargin() >= 3)
    numSubjects = varargin{3};
else
    numSubjects=1;
end


%% Set trial-invariant parameters 
fixSecondary = 0; %Increase secondary weights to begin under secondary control
VIschedule = 0; %Is this experiment subject to Woods & Bouton VI schedule?
padamount = 19; % How long are windows used to determine per-subject averages?

blockLengthA= 500;         %Trial block lengths
blockLengthB= 50;
blockLengthC= 100;
blockLengthD= 50;
blockLengthE= 20;
blockLengthF= 40;
blockLengthG= 20;
blockLengthH= 1000;
blockLengthI = 2000;
blockLengthLong = 1500;

% Parameters used for Woods & Bouton experiments
WBsegment = 45;
WBtrainingBlock = 1000;
WBextBlock = WBsegment*8;
WBreacBlock = WBsegment*2;
prfChance = 0.035; %How likely is a random reinforcement under partial schedule?
VIscale = 0.5; % How much does VI schedule reduce learning?
VIchance = 1; % How rare is reinforcement per trial on VI schedule (unused)

% Parameters used for savings effect experiments
SVtrainingBlock = 25*4;
SVextBlock = 25*4;
SVreacBlock = 25*4;



stimPeriod = 4; %Time steps per trial before secondary stimuli are processed
stimSize = 2; %Number of unique stimuli
respSize = 3; %Number of unique responses
contextCells = 2*stimSize; %How many stimuli are accessible to secondary only
aInit = 0.25; %Initial activation (all neurons)
DtInit = 0.025; %Timestep
S = 200; %Steps per trial before a decision must be made
activityInitialNoise = 0.1; %Std dev of noise at beginning of each trial
activityNoise = 0.1; %Std dev. of noise at each time step
sVal = 1; %Strength of stimulus input to CTX 
uVal = 0.0; %Strength of broad urgency signal in CTX
predictionLength = 50;
history = zeros(1,predictionLength);
historyLong = zeros(1,3*predictionLength);
interrupt = 0;

% Weight adjustment parameters
cocaineBias = 1; %Multiplies positive and divides negative rewards 
wOffset = 0.2; %Baseline weight value in primary
wDblOffset = 0.2; %As above but for the secondary weights
wPriSpread = 0.02; %Difference between primary GO and NOGO weight init
wSecSpread = 0.1; %As above, but for secondary
wMaxHi = 0.9; %Maximum weight value (High compartment)
wMaxLo = 0.5; %Max weight value (Lo compartment)
wCtxBaseline = 1; %" weights onto cortex
weightsInitialNoiseDev = 0.01; %Stdev of initial noise in weights
weightsInitialNoise = weightsInitialNoiseDev*randn(respSize,contextCells); %initial noise in all weights
wDecayRate = 0.000; %Rate that weights return to baseline
rDeltaHi = 0.0045; %Learning rate in the secondary cpt
rDeltaLo = 0.0036; %Learning rate in the primary cpt
ForwardCouplingTerm = 1; %How much reinforcement passes from secondary to primary
BackwardCouplingTerm = 0.0; %As above, but from primary to secondary
NGSprimary = 5.0; %Learning rate scaling factor for the primary cpt NOGO pathway
NGSsecondary = 5.0; %Learning rate scaling factor for the 2ndary cpt NOGO pathway
NGPunishmentAdded = 0; %Amount added to negative reinforcement for the NOGO pathway
punishmentVal = 0.8; %Is negative reward valued differently from positive?
SecondaryPunishmentScale = 1.0; % Multiply punishmentVal by this, but only for the Secondary compartment. Unused.
nonresponseR = 0; %Reinforcement applied if the model doesn't produce a response due to timeout

% Parameters for drug effect simulation variation
%Is there an input bias affecting cortical input? 
secondaryStimBias = 1;
primaryStimBias = 1;

%Beta parameter affecting DA efficacy
DAbetaSecondary = 0.0;
DAbetaPrimary = 0.0;
%Auto-inhibition (decay rate) parameters
secondaryGABAbias = 0.05;
primaryGABAbias = 0.05;
%Steepness of RPE parametric curve
RPEPrimaryExponent = 1/2;
RPESecondaryExponent = 3/2;


% snapshotTrials = [1 88 280 741 959 1248 1373 1766 1876]; %Which trials to save ALL data for
snapshotTrials =[];
SATparadigm = 0; %Toggle used in SAT simulations


%% Save the workspace so parameters are kept
paramPath = sprintf('%s%sWorkspace.mat',fpath,date);
save(paramPath);


%% Iterate the entire script for each figure to be generated
numFigs = sum(plotToggle);

for figInd = 1:numFigs  % For each experiment

clearvars -except figInd subjInd numFigs paramPath plotToggle 
load(paramPath,'-regexp','^(?!plotToggle)\w')
rng(1337*figInd); %Experiment-dependent RNG seed
plotDraw = [0 0 0 0 0 ...
            0 0 0 0 0 ...
            0 0 0 0 0 ...
            0 0 0];


%% Learning parameters
% Defaults used for key parameters. May be varied per-experiment.

L = wMaxLo;
H = wMaxHi;
HR = rDeltaHi;
LR = rDeltaLo;
dTh = 0.75;

% Params: [thGamma Secondary, thGamma Primary, rDelta Primary, rDelta Secondary, 
%           decision Threshold]
Params = [H L LR HR dTh 0.75 1;... %Expected best performance
    ];

%% Reversal figure
if (plotToggle(1) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [blockLengthA,1;...
                 blockLengthA,2;...
                 blockLengthA,1;...
                  blockLengthA,2;...                 
                  blockLengthA,1;...          
                  blockLengthA,2;...        
                ];      
            
            chanceMatrix = ruleMatrix;
            
plotToggle(1) = 0;
plotDraw (1) = 1;
end

%% Matching / probabilistic 
if (plotToggle(2) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 
ruleMatrix =    [blockLengthA,10;...
                 blockLengthA,10;...
                 blockLengthA,10;...
                 blockLengthA,10;...
                ];           

chanceMatrix = [0.5 0.75;...
    0.7 0.75;...
    0.9 0.75;...
    1 0.75;...
    ];
plotToggle(2) = 0;
plotDraw (2) = 1;

end

%% Speed-accuracy tradeoff probe
if (plotToggle(3) ~= 0 && sum(plotDraw) == 0)
    
ruleMatrix = [...
                 2000,12;...
%                  500,12;...
                 ];
             
chanceMatrix = [0.5 0.5; 
                0.5 0.5];           
Params = [ H L LR HR dTh 0.5];
    
plotToggle(3) = 0;
plotDraw (3) = 1;
SATparadigm = 1;
end

%% Plus-maze imitation
if (plotToggle(4) ~= 0 && sum(plotDraw) == 0)  

% Duration, Contingency 
ruleMatrix =    [blockLengthB,3;...
                 blockLengthC,4;...
                 blockLengthD,5;...
                 blockLengthB,3;...
                 blockLengthH,4;...
                 blockLengthD,5;...               
                 blockLengthB,3;...
                 blockLengthI,4;...
                 blockLengthD,5;...   
                ];
chanceMatrix = ruleMatrix;
plotToggle(4) = 0;
plotDraw (4) = 1;
end


%% Extended acquisition
if (plotToggle(5) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [blockLengthLong,14;...
                 blockLengthLong,13;...
                ];      
chanceMatrix = [1 0.5;...
0.75 0.5;...
];          
plotToggle(5) = 0;
plotDraw (5) = 1;
end


%% Extinction
if (plotToggle(6) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [blockLengthE,1;...
                 blockLengthF,0;...
                 blockLengthG,1;...        
                ];      
chanceMatrix = [1 0.5;...
0.75 0.5;...
0.5 0.5;...
];          
plotToggle(6) = 0;
plotDraw (6) = 1;
end

%% Cocaine effect figure
if (plotToggle(7) ~= 0 && sum(plotDraw) == 0)  

Params = [H*1.2 L*1.5 LR*1.5 HR*1.2 dTh 0];     

ruleMatrix =    [blockLengthA,1;... 
                 blockLengthA,2;...
                 blockLengthA,1;...
                 blockLengthA,2;... 
                 blockLengthA,1;...
                 blockLengthA,2;...
                ];
            
chanceMatrix = ruleMatrix;

plotToggle(7) = 0;
plotDraw(7) = 1;
end


%% Ethanol effect figure
if (plotToggle(8) ~= 0 && sum(plotDraw) == 0)  

Params = [...
          H L LR HR dTh 0;...
         ];
     
secondaryStimBias = 1;
primaryStimBias = 1;

secondaryGABAbias = 0.05;
primaryGABAbias = 0.05;

ruleMatrix =    [500,1;... 
                 500,2;...
                 500,1;...
                 500,2;... 
                 500,1;...
                 500,2;...
                ];
            
chanceMatrix = ruleMatrix;

plotToggle(8) = 0;
plotDraw(8) = 1;
end

%% Amphetamine effect figure
fignum = 9;
if (plotToggle(fignum) ~= 0 && sum(plotDraw) == 0)  

Params = [...
          H*1.75 L*1.1 LR*1.75 HR*1.5 dTh 0;...
         ];
     
secondaryStimBias = 1;
primaryStimBias = 1;

secondaryGABAbias = 0.05;
primaryGABAbias = 0.05;


ruleMatrix =    [500,1;... 
                 500,2;...
                 500,1;...
                 500,2;... 
                 500,1;...
                 500,2;...
                ];
chanceMatrix = [500,1;... 
                 500,2;...
                 500,1;...
                 500,2;... 
                 500,1;...
                 500,2;...
                ];

plotToggle(fignum) = 0;
plotDraw(fignum) = 1;
end


%% Nicotine effect figure
fignum = 10;
if (plotToggle(fignum) ~= 0 && sum(plotDraw) == 0)  

% Elevated phasic relase preferentially in ventral    
Params = [...
          H*1.0 L*1.0 LR*1.1 HR*1.2 dTh 0;...
         ];
     

%Reduced tonic DA preferentially in dorsal
DAbetaSecondary = -0.05;
DAbetaPrimary = -0.1;


predictability = 1;
chanceMatrix = [predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                ];

ruleMatrix =    [blockLengthA,1;...
                 blockLengthA,2;...
                 blockLengthA,1;...
                  blockLengthA,2;...                 
                  blockLengthA,1;...          
                   blockLengthA,2;...  
                ];
plotToggle(fignum) = 0;
plotDraw(fignum) = 1;
end

%% Nicotene withdrawal effect figure
fignum = 11;
if (plotToggle(fignum) ~= 0 && sum(plotDraw) == 0)  

loadBGdata = 0;
saveBGdata = 1;

Params = [...
          H*1.0 L*1.0 LR*1 HR*1 dTh 0;...
         ];
     
%Reduced tonic DA preferentially in NAcc
DAbetaSecondary = -0.05;
DAbetaPrimary = -0.1;


predictability = 1;
chanceMatrix = [predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                ];

ruleMatrix =    [blockLengthA,1;...
                 blockLengthA,2;...
                 blockLengthA,1;...
                  blockLengthA,2;...                 
                  blockLengthA,1;...          
                   blockLengthA,2;...  
                ];
            
plotToggle(fignum) = 0;
plotDraw(fignum) = 1;
end

%% Bupropion + nicotene withdrawal effect figure
fignum = 12;
if (plotToggle(fignum) ~= 0 && sum(plotDraw) == 0)  

loadBGdata = 0;
saveBGdata = 1;

%Elevated learning rate, both cpts, preferentially ventral
Params = [...
          H*1.0 L*1.0 LR*1.0 HR*1.0 dTh 0;...
         ];
     
%Reduced tonic DA 
DAbetaSecondary = +0.1;
DAbetaPrimary = +0.05;

%Preferential increase in D2 (NOGO)
punishmentVal = 1.2;

% Increase feedforward
ForwardCouplingTerm = 1;

predictability = 1;
chanceMatrix = [predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                ];

ruleMatrix =    [blockLengthA,1;...
                 blockLengthA,2;...
                 blockLengthA,1;...
                  blockLengthA,2;...                 
                  blockLengthA,1;...          
                   blockLengthA,2;...  
                ];
            
plotToggle(fignum) = 0;
plotDraw(fignum) = 1;
end

%P.E. of dopamine figure - vary only predictability
fignum = 13;
if (plotToggle(fignum) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [blockLengthA,15;...
                 blockLengthA,16;...
                 blockLengthA,15;...
                  blockLengthA,16;...                 
                  blockLengthA,15;...          
                   blockLengthA,16;...   
                ];      
     
predictability = 0.75;
chanceMatrix = [predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                ];
plotToggle(fignum) = 0;
plotDraw (fignum) = 1;

DAbetaSecondary = 0.0;
DAbetaPrimary = 0.0;
end

%Performance effect of dopamine figure - vary predictability and beta
fignum = 14;
if (plotToggle(fignum) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [blockLengthA,15;...
                 blockLengthA,16;...
                 blockLengthA,15;...
                  blockLengthA,16;...                 
                  blockLengthA,15;...          
                   blockLengthA,16;...   
                ];      
            
predictability = 1;
chanceMatrix = [predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                predictability predictability; ...
                ];
plotToggle(fignum) = 0;
plotDraw (fignum) = 1;

DAbetaSecondary = 0.0;
DAbetaPrimary = 0.0;
end


%% Woods Bouton Figure - Prf2
if (plotToggle(15) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [WBtrainingBlock,19;...
                 WBextBlock,15;...
                 WBreacBlock,18;...     
                ];      
chanceMatrix = [1 1;...
prfChance 1;...
VIchance 1;...
];          
plotToggle(15) = 0;
plotDraw (15) = 1;
woodsBoutonFig = 1;
end

%% Woods Bouton Figure - Prf8
if (plotToggle(16) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 
     
VIschedule = 8;
ruleMatrix =    [WBtrainingBlock,19;...
                 WBextBlock,15;...
                 WBreacBlock,18;...     
                ];    
            
chanceMatrix = [1 1;...
prfChance 1;...
VIchance 1;...
];          
plotToggle(16) = 0;
plotDraw (16) = 1;
woodsBoutonFig = 1;
end

%% Woods Bouton Figure - Ext2
if (plotToggle(17) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [WBtrainingBlock,19;...
                 WBextBlock,15;...
                 WBreacBlock,18;...     
                ];     
chanceMatrix = [...
1 1;...
0 1;...
VIchance 1;...
];          
plotToggle(17) = 0;
plotDraw (17) = 1;

woodsBoutonFig = 1;
end


%% Woods Bouton Figure - Ext8
if (plotToggle(18) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

VIschedule = 8;
ruleMatrix =    [WBtrainingBlock,19;...
                 WBextBlock,15;...
                 WBreacBlock,18;...     
                ];       
            
chanceMatrix = [1 1;...
0 1;...
VIchance 1;...
];          

plotToggle(18) = 0;
plotDraw (18) = 1;

woodsBoutonFig = 1;
end

%% Savings figure
if (plotToggle(19) ~= 0 && sum(plotDraw) == 0)
% Duration, Contingency 

ruleMatrix =    [SVtrainingBlock,18;...
                 SVextBlock,17;...
                 SVreacBlock,18;...     
                ];       
            
chanceMatrix = [1 1;...
1 1;...
1 1;...
];          

plotToggle(19) = 0;
plotDraw (19) = 1;
savingsFig = 1;
end

% Iterate over the number of subjects for this experiment
for subjInd = 1:numSubjects
    
%Iterate over parameter sets per subject
[Psets, ~] = size(Params);
for pSweep = 1:Psets    
    

%% Generate stimuli
T = sum(ruleMatrix(:,1)); %Number of trials
stimVector = zeros(stimSize,T);
contextVector = zeros(stimSize,T);
[Len, ~] = size(ruleMatrix);
for block = 1:Len
    [sV, cV] = stimulusByRule( ruleMatrix(block,2), ruleMatrix(block,1), chanceMatrix(block,2));
    
    if (block == 1)
        startInd = 1;
    else
        startInd = 1+sum(ruleMatrix(1:block-1,1));
    end
    endInd = startInd + ruleMatrix(block,1)-1;
    
    stimVector(:,startInd : endInd)...
        = sV;
    contextVector(:,startInd : endInd)...
        = cV;
end
clear block

% Initialize compartment structures
if(loadBGdata == 0)
SINGLE_BG = struct(...
    'GO',aInit*ones(respSize,contextCells),...
    'NOGO',aInit*ones(respSize,contextCells),...
    'MC',aInit*ones(respSize,contextCells),...
    'S_IN',aInit*ones(respSize,contextCells),...
    'W_GO',wOffset + wPriSpread + weightsInitialNoise,...
    'W_NOGO',wOffset -wPriSpread + weightsInitialNoise,...
    'W_CTX',wCtxBaseline + weightsInitialNoise,...
    'thGamma',0,...
    'rDelta',0,...
    'dt',DtInit,...
    'STOP',aInit,...
    'noise',activityNoise,...
    'stimulusStrength',sVal,...
    'urgencyStrength',uVal,...
    'NOGOscale',NGSprimary,...
    'NOGOAddedPunishment',NGPunishmentAdded,...
    'weightDecay',wDecayRate,...
    'weightInit',wOffset,...
    'decayParameter',primaryGABAbias,...
    'DA_beta',DAbetaPrimary...
);

DOUBLE_BG = struct(...
    'GO',aInit*ones(respSize,contextCells),...
    'NOGO',aInit*ones(respSize,contextCells),...
    'MC',aInit*ones(respSize,contextCells),...
    'S_IN',aInit*ones(respSize,contextCells),...
    'W_GO',wDblOffset +wSecSpread + weightsInitialNoise,...
    'W_NOGO',wDblOffset -wSecSpread +weightsInitialNoise,...
    'W_CTX',wCtxBaseline + weightsInitialNoise,...
    'thGamma',0,...
    'rDelta',0,...
    'dt',DtInit,...
    'STOP',aInit,...
    'noise',activityNoise,...
    'stimulusStrength',sVal,... 
    'urgencyStrength',uVal,...
    'NOGOscale',NGSsecondary,...
    'NOGOAddedPunishment',NGPunishmentAdded,...
    'weightDecay',wDecayRate,...
    'weightInit',wOffset,...
    'decayParameter',secondaryGABAbias,...
    'DA_beta',DAbetaSecondary...
);
else
    load BGdata.mat
end

%Storage variables
rVector = zeros(T,5);
wVector = zeros(T,6,respSize,contextCells,length(Params));
tVector = zeros(T,1);
MCvector = zeros(T,S,2,respSize,contextCells);
GOVector = zeros(T,S,2,respSize,contextCells);
NOGOVector = zeros(T,S,2,respSize,contextCells);

% Reset system to begin each experiment
ruleCount = 0;
counter = 1;
SINGLE_BG = resetWeights(SINGLE_BG,weightsInitialNoiseDev,wOffset,wPriSpread);
DOUBLE_BG = resetWeights(DOUBLE_BG,weightsInitialNoiseDev,wDblOffset,wSecSpread);

% Fix secondary compartment to begin in control, if desired
if(fixSecondary == 1) 
DOUBLE_BG.W_GO = 1.5*DOUBLE_BG.W_GO;
DOUBLE_BG.W_NOGO = 1.5*DOUBLE_BG.W_NOGO; 
else
    if(fixSecondary == 2)
    SINGLE_BG.W_GO = 0.5 * (SINGLE_BG.W_GO);
    SINGLE_BG.W_NOGO = 1.5 * (SINGLE_BG.W_NOGO);
    DOUBLE_BG.W_GO = 1.5.*(DOUBLE_BG.W_GO);
    DOUBLE_BG.W_NOGO = 0.5.*DOUBLE_BG.W_NOGO; 
    end
end

%% Trials
    for t = 1:T % Iterate over each trial
        

     %Reset values for new trial
     SINGLE_BG = resetActivity(SINGLE_BG,aInit,activityInitialNoise);
     DOUBLE_BG = resetActivity(DOUBLE_BG,aInit,activityInitialNoise);     
     SINGLE_BG.thGamma = Params(pSweep,1);
     DOUBLE_BG.thGamma = Params(pSweep,2);
     SINGLE_BG.rDelta = Params(pSweep,3);
     DOUBLE_BG.rDelta = Params(pSweep,4);
     decisionTh = Params(pSweep,5);

     stimulus = [stimVector(:,t) stimVector(:,t) stimVector(:,t)]';
     trackWinner = 0.5;
     trackDecision = 0;
      
     % Use the trial/rule matrix to choose the correct one for this trial
     if(ruleCount >= ruleMatrix(counter,1))
         ruleCount = 0;
         counter = counter + 1;
         rule = ruleMatrix(counter,2);
     else
         rule = ruleMatrix(counter,2);
     end
     ruleCount = ruleCount + 1;
        
     % If VI schedule is in effect, adjust learning rate
     if(t > sum(ruleMatrix(1:2,1)) && VIschedule == 8)
         SINGLE_BG.rDelta = VIscale*Params(pSweep,3);
         DOUBLE_BG.rDelta = VIscale*Params(pSweep,4);
     else
         SINGLE_BG.rDelta = Params(pSweep,3);
         DOUBLE_BG.rDelta = Params(pSweep,4);
     end

    % Determine contingencies to be rewarded
    rAct = contingencyByRule(rule,stimVector(:,t),contextVector(:,t),chanceMatrix(counter,1),punishmentVal);
          
%% Stimulus loop        
    for step = 1:S % A trial has begun

        %Adjust reward contingencies with timeouts
%         if (plotDraw(3) == 1) %During SAT figure generation
%             if (t >= windowOn)                    
%                 win = shutoff; %Time after which reward cannot be obtained
%             else
%                 win = S;
%             end
%             if(step > win)
%             rAct = [-1 -1 -1];
%            end
%         end
                    
        SINGLE_BG.S_IN = zeros(respSize,contextCells);
        SINGLE_BG.S_IN(:,1:2) = primaryStimBias*stimulus(:,1:2);
        DOUBLE_BG.S_IN = zeros(respSize,contextCells);
            
        if (step >= stimPeriod)            
            DOUBLE_BG.S_IN = secondaryStimBias*[stimulus(:,1:stimSize), contextVector(1,t)*ones(respSize,1),contextVector(2,t)*ones(respSize,1)];
        else    
            DOUBLE_BG.S_IN = zeros(respSize,contextCells);
        end
        
            SINGLE_BG = ModelBG(SINGLE_BG);
            DOUBLE_BG = ModelBG(DOUBLE_BG);
            [single_vals , ~]= max(SINGLE_BG.MC);
            [double_vals , ~] = max(DOUBLE_BG.MC);
            single_winner = [0 0 0];
            double_winner = [0 0 0];
            single_val = 0;
            double_val = 0;


if (SATparadigm == 1)
    finalStep = 40+ randi(20,1);
    if (step >= finalStep)
        interrupt = 1;
    end
end
        %End loop if either:
        % a) Time is up (S)
        % b1) One compartment has an activity value greater than the other
        % compartment and the decision threshold, and
        % b2) The winning compartment has settled on one winner (the
        % winning value is unique)
        if (step == S ...
                || (max(single_vals) > decisionTh  && ...
                   max(single_vals) > max(double_vals) && ...
                    ~isequal(single_vals(1),single_vals(2),single_vals(3)))...
                || (max(double_vals) > decisionTh && ...
                   max(double_vals) > max(single_vals) && ...
                   ~isequal(double_vals(1),double_vals(2),double_vals(3)))...
                || (interrupt == 1)...
                )
            
            %Last step - pick a winner    
            [single_vals, single_winner]= max(SINGLE_BG.MC);
            [double_vals, double_winner] = max(DOUBLE_BG.MC);
            [single_val, single_ind] = max(single_vals);
            [double_val, double_ind] = max(double_vals);
            break
            
        end


    end

tVector(t) = step; %Indicates how long it took to make a decision

%% Check result, apply learning and update weights
R1 = zeros(respSize,contextCells);
R2 = zeros(respSize,contextCells);

%% Calculate reward
        if ((~isequal(single_vals(1),single_vals(2),single_vals(3)) ||  ~isequal(double_vals(1),double_vals(2),double_vals(3))) && step < S)
            % The primary compartment won the competition
            if (single_val > double_val &&  ~isequal(single_vals(1),single_vals(2),single_vals(3)))
                
                R1 = zeros(respSize,contextCells);
                R2 = zeros(respSize,contextCells);

                if (single_val >= 0)
                R_raw = rewardCalc(single_winner(single_ind),rAct);
                R = R_raw*computeRPE(R_raw,history)^(RPESecondaryExponent);
                R_primary = R_raw*computeRPE(R_raw,history)^(RPEPrimaryExponent);
                

                if (R<0)
                    R1(single_winner(single_ind),single_ind) = R_primary/cocaineBias;
                    R2(single_winner(single_ind),single_ind) = BackwardCouplingTerm*R/cocaineBias; 
                else
                    R1(single_winner(single_ind),single_ind) = R_primary*cocaineBias;
                    R2(single_winner(single_ind),single_ind) = BackwardCouplingTerm*R*cocaineBias; 
                end
                trackWinner = 1;
                trackDecision = single_winner(single_ind);
                else
                R_raw = 0; 
                R1(single_winner(single_ind),single_ind) = 0;
                R2(single_winner(single_ind),single_ind) = 0;
                trackWinner = -1;
                trackDecision = single_winner(single_ind);
                end


            else 
                %The secondary compartment won the competition
                if (double_val >= single_val &&  ~isequal(double_vals(1),double_vals(2),double_vals(3)))

                    R1 = zeros(respSize,contextCells);
                    R2 = zeros(respSize,contextCells); 

                    if (double_val >= 0)
                    R_raw = rewardCalc(double_winner(double_ind),rAct); 
                    R = R_raw*computeRPE(R_raw,history)^(RPESecondaryExponent);
                    R_primary = R_raw*computeRPE(R_raw,history)^(RPEPrimaryExponent);
                    
                    if (R < 0)
                        R1(double_winner(double_ind),:) = ForwardCouplingTerm*R_primary/cocaineBias; %Asymmetrical "spiral" learning
                        R2(double_winner(double_ind),double_ind) = SecondaryPunishmentScale*R/cocaineBias;
                    else
                        R1(double_winner(double_ind),:) = ForwardCouplingTerm*R_primary*cocaineBias; %Asymmetrical "spiral" learning
                        R2(double_winner(double_ind),double_ind) = R*cocaineBias;
                    end
                    trackDecision = double_winner(double_ind);
                    trackWinner = 0;
                    else
                    R_raw = 0; 
                    R1(double_winner(double_ind),double_ind) = 0;
                    R2(double_winner(double_ind),double_ind) = 0;
                    trackWinner = -1;
                    trackDecision = double_winner(double_ind);
                    end
                end
            end  
        else
            %Timeout - neither won
            R_raw = 0; 
            R1 = nonresponseR*ones(respSize,contextCells);
            R2 = nonresponseR*ones(respSize,contextCells);
            trackWinner = -1;
            trackDecision = 0;
        end
        
        %Update RPE history
        history = [R_raw history(1:predictionLength-1)];
        historyLong = [R_raw historyLong(1:predictionLength-1)];
        
%% Update weights
        SINGLE_BG = LearnBG(SINGLE_BG,R1);
        DOUBLE_BG = LearnBG(DOUBLE_BG,R2);
        
%% Record for plotting        
      %  Save weights,stimuli, and outcomes
        for WrInd = 1:respSize    
            for WsInd = 1:contextCells
            wVector(t,1,WrInd,WsInd,pSweep) = SINGLE_BG.W_GO(WrInd,WsInd);
            wVector(t,2,WrInd,WsInd,pSweep) = SINGLE_BG.W_NOGO(WrInd,WsInd);
            wVector(t,3,WrInd,WsInd,pSweep) = DOUBLE_BG.W_GO(WrInd,WsInd);
            wVector(t,4,WrInd,WsInd,pSweep) = DOUBLE_BG.W_NOGO(WrInd,WsInd);   
            end
        end
            
        if (trackDecision ~= 0)
         rVector(t,1) = rAct(trackDecision);
        else
            rVector(t,1) = 0;
        end
        rVector(t,2) = single_winner(1);
        rVector(t,3) = double_winner(1);
        rVector(t,4) = trackWinner;
        rVector(t,5) = trackDecision;
        
      % Take a snapshot on certain trials and save it 
      if (ismember(t,snapshotTrials))
          fileString = sprintf('Snapshot_trial%d',t);
          save(fileString);
      end
    end
    
%% Generate percentage data

pctPrimaryDecision = zeros(1,T);
pctSuccesses = zeros(1,T);
avgTimes = zeros(1,T);


pad = ones(padamount,1);


    
pastResponses = zeros(3,T);
    
%Iterate over entire experiment's set of trial blocks
    for k = 1:length(ruleMatrix)  


        if k > 1
            interval = sum(ruleMatrix(1:k-1,1))+1:sum(ruleMatrix(1:k,1));
            paddedVector = [paddedVector(end-padamount+1:end); (rVector(interval, 4))]; 
            % RVector is formatted as 2 for secondary, 1 for primary, -1 for abstention
            % Want format of -1 for abstention, 1 for primary, 0 for secondary
            paddedVector(paddedVector == 2) = 0;
            paddedSuccessVector = [paddedSuccessVector(end-padamount+1:end); (rVector(interval, 1) > 0)];
            paddedTimeVector = [paddedTimeVector(end-padamount+1:end); (tVector(interval))];
            paddedActionVector = [paddedActionVector(end-padamount+1:end); rVector(interval,5)];
          else
            interval = 1:ruleMatrix(1,1);
            paddedVector = [0*pad; (rVector(interval, 4))]; 
            paddedVector(paddedVector == 2) = 0;
            paddedSuccessVector = [0*pad; (rVector(interval, 1) > 0)];
            paddedTimeVector = [tVector(1)*pad; (tVector(interval))];
            paddedActionVector = [0*pad; rVector(interval,5)];
        end       

        %Per-block statistics
        for j = interval 
            
            pastKdecisions =  paddedVector(j-interval(1)+1:j+padamount-interval(1));
            pctPrimaryDecision(j) = mean(pastKdecisions); 
            pastKsuccesses = paddedSuccessVector(j-interval(1)+1:j+padamount-interval(1));
            pastKresponses = paddedActionVector(j-interval(1)+1:j+padamount-interval(1));
            pctSuccesses(j) = mean(pastKsuccesses);
            avgTimes(j) = mean(paddedTimeVector(j-interval(1)+1:j+padamount-interval(1)));
            
            for respInd = 1:3
             pastResponses(respInd,j) = sum((pastKresponses == respInd))/(padamount);
            end
        end
    end
    
    %Used only for probabilistic reversal figure - shows "ideal" (P=1)
    %performance
    if(length(ruleMatrix) == 6)
        
        leftWindow = 0;
        rightWindow = 500;
        if (ruleMatrix(1,1) > leftWindow)
        indicesBlock1 = ruleMatrix(1,1)-leftWindow:ruleMatrix(1,1)+rightWindow;
        else
        indicesBlock1 = 1:ruleMatrix(1,1)+rightWindow+1;
        end
        indicesBlock2 = sum(ruleMatrix(1:3))-leftWindow:sum(ruleMatrix(1:3))+rightWindow;
        indicesBlock3 = sum(ruleMatrix(1:5))-leftWindow:sum(ruleMatrix(1:5))+rightWindow;
    %     indicesBlock4 = sum(ruleMatrix(1:7))-leftWindow:sum(ruleMatrix(1:7))+rightWindow;

        desiredResponses = [pastResponses(1,1:ruleMatrix(1,1))...
            pastResponses(2,ruleMatrix(1,1)+1:sum(ruleMatrix(1:2))) ...
             pastResponses(1,sum(ruleMatrix(1:2))+1:sum(ruleMatrix(1:3))) ...
          pastResponses(2,sum(ruleMatrix(1:3))+1:sum(ruleMatrix(1:4))) ...   
          pastResponses(1,sum(ruleMatrix(1:4))+1:sum(ruleMatrix(1:5))) ...
          pastResponses(2,sum(ruleMatrix(1:5))+1:sum(ruleMatrix(1:6)))...
         ];
    end
    
%% Find crossing points between compartments
    crossing = zeros(1,6);
    for blockInd = 1:length(ruleMatrix)
    
        if (blockInd == 1)
            relevantIndices = 1:ruleMatrix(1);
        else
            relevantIndices = sum(ruleMatrix(1:blockInd-1)):sum(ruleMatrix(1:blockInd));
        end
        
        criterion = 0.5;
        [above, aboveInds] = find(pctPrimaryDecision(relevantIndices) >= criterion + 0.05);    
        [below, belowInds] = find(pctPrimaryDecision(relevantIndices) <= criterion - 0.05);

        if (~isempty(belowInds))
            [nextAbove, nextInd] = find(aboveInds > belowInds(end),1,'first');
        else
            crossing(blockInd) = NaN;
            continue
        end
        
        if (~isempty(nextInd))
            crossing(blockInd) = aboveInds(nextInd) + relevantIndices(1);
        else
            crossing(blockInd) = NaN;
        end

       
    end
    

    % This is the file used later by LearningSysPlotScript
    if (saveBGdata == 1)
        exptext = sprintf('%s%sExperiment%dSubject%d.mat',fpath,date,figInd,subjInd);
        save(exptext);
    end

    
 
end %End of param loop
end %End of subject loop
end %End figure generating loop

end %End script

%% Model differential equations
function [state] = ModelBG(state)
A = state.decayParameter; %Decay rate
B = 1; %Upper activity bound
Sval = state.stimulusStrength; %Impact of stimulus on activity
noiseInput = state.noise*randn(size(state.MC));
urgencySignal = state.urgencyStrength; %Cortical input bias, if any

state.GO = state.GO + state.dt.*(-A.*state.GO + (B - state.GO).*state.W_GO.*state.MC*(1+state.DA_beta)); % GO pathway
state.NOGO = state.NOGO + state.dt.*(-A.*state.NOGO + (B - state.NOGO).*state.W_NOGO.*state.MC*(1-state.DA_beta)); % NOGO pathway
state.MC = state.MC + state.dt.*(noiseInput -A.*state.MC + ...
                                (B - state.MC).*(urgencySignal + Sval.*state.S_IN + state.W_CTX.*state.GO) ...
                                 - state.MC.*state.W_CTX.*state.NOGO); % "Motor cortex"

end



%% Adjust weights based on reward
function [state] = LearnBG(state,R)
    
    %Linear learning rule      
    state.W_GO = state.W_GO + (state.rDelta.*R);
    state.W_NOGO = state.W_NOGO + (-state.NOGOscale*state.rDelta.*R.*( ones(size(R)) + (R<0).*(state.NOGOAddedPunishment)));

    
     % Weight decay to equilibrium
     state.W_GO = state.W_GO - state.weightDecay*(state.W_GO-state.weightInit);
     state.W_NOGO = state.W_NOGO - state.weightDecay*(state.W_NOGO-state.weightInit);

    % Thresholding
    state.W_GO = state.thGamma .* (state.W_GO >= state.thGamma) + state.W_GO .* (state.W_GO < state.thGamma);
    state.W_NOGO = state.thGamma .* (state.W_NOGO >= state.thGamma) + state.W_NOGO .* (state.W_NOGO < state.thGamma);
    state.W_GO = 0 .* (state.W_GO <= 0) + state.W_GO .* (state.W_GO > 0);
    state.W_NOGO = 0 .* (state.W_NOGO <= 0) + state.W_NOGO .* (state.W_NOGO>0);

end

%% Determine whether the contingency was followed
function[R_out] = rewardCalc(winner,rAct)
%Based on contingency between stimulus and rewarded outcome, determine
%reward
R_out = 0;
    if (winner(1) ~= 0)
        R_out = rAct(winner(1));
    end
end

%% Reset activities for a new trial
function [COMP] = resetActivity(COMP,aInit,aNoise)
COMP.MC = aInit*ones(size(COMP.MC))+aNoise*randn(size(COMP.MC));
COMP.GO = aInit*ones(size(COMP.GO));
COMP.NOGO = aInit*ones(size(COMP.NOGO));
COMP.STOP = aInit;
end

%% Reset weights between parameters
function [COMP] = resetWeights(COMP,weightsInitialNoise,offset,spread)

COMP.W_GO = offset + spread + weightsInitialNoise*randn(size(COMP.W_GO));
COMP.W_NOGO = offset -spread + weightsInitialNoise*randn(size(COMP.W_NOGO));
end

function [stimVector,contextVector] = stimulusByRule(rule,trials,param)

stimVector = zeros(2,trials);
contextVector = zeros(2,trials);
    switch rule
        case {0 1 2 7 14 15 16 17 18 19 20} %Present S1 only
            stimVector = [ones(1,trials); zeros(1,trials)];
            contextVector = [zeros(1,trials); zeros(1,trials)];        
        case 3 %Present mixture of S3 and S4
            randVec = rand(1,trials);
            stimVector = [zeros(1,trials); zeros(1,trials)];
            contextVector = [randVec>0.5; randVec<=0.5];
        case 4 %Present S1+S3
            stimVector = [ones(1,trials); zeros(1,trials)];
            contextVector = [ones(1,trials); zeros(1,trials)];
        case 5 %Present S1+S4
            stimVector = [ones(1,trials); zeros(1,trials)];
            contextVector = [zeros(1,trials); ones(1,trials)];
        case 6 %Present either S1 OR S2, plus random S3.
            stimAlt = randi(2,1,trials)-1;
             stimVector = [stimAlt; (stimAlt == 0)];
            contextVector = [randi(2,1,trials)-1; zeros(1,trials)];
        case 8 %Present S1, S1+S3, or S3
            randVec = rand(1,trials);
            stimVector = [randVec<0.33; zeros(1,trials)];
            contextVector = [randVec<0.66; zeros(1,trials)];
        case 9 %Present S1 with param chance of co-presenting S2
            randVec = rand(1,trials);
            stimVector = [ones(1,trials); randVec<param];
            contextVector = [zeros(1,trials); zeros(1,trials)];

        case 10 %Present S1 with param chance of co-presenting S3
            randVec = rand(1,trials);
            stimVector = [ones(1,trials); zeros(1,trials)];
            contextVector = [randVec<param; zeros(1,trials)];

        case 11 %Present S3 or S4 and S2
            randVec = rand(1,trials);
            randVec2 = rand(1,trials);
            stimVector = [zeros(1,trials); ones(1,trials)];
            contextVector = [randVec2 >= 0.5; randVec2 < 0.5];

        case 12 
            randVec = rand(1,trials);
            randVec2 = rand(1,trials);
            stimVector = [randVec >= 0.7; randVec < 0.3];
            contextVector = [randVec2 >= 0.5; randVec2 < 0.5];

        case 13 %S3 only
            stimVector = [zeros(1,trials); zeros(1,trials)];
            contextVector = [ones(1,trials); zeros(1,trials)];

    end

end

% Compute a reward-prediction error signal given reward history
function [RPE] = computeRPE(reward,rewardHistory)

predictionLength = length(rewardHistory);
matches = sum(rewardHistory == reward);
mismatches = sum(rewardHistory ~= reward);
RPE = 0.25*((mismatches+predictionLength)/(matches+predictionLength))^2;

end   
 
%%Determine stimulus-response pairings given integer rule and stimulus
function [rAct] = contingencyByRule(rule,stimVector,contextVector,chance,omission)
    rVal = omission;
    rAct = [-rVal -rVal -rVal];
    RandomVar = rand(); %Random number between 0 and 1 for probabilistic reward
            switch rule
                case 0 %Extinction
                    rAct = [-1 -1 -1];

                case 1 %Response 1

                        switch stimVector(1)
                        case 1
                            switch stimVector(2)
                                case 1
                                    rAct(1) = 1;
                                case 0
                                    rAct(1) = 1;
                            end
                        case 0
                            switch stimVector(2)
                                case 1
                                case 0
                                    rAct = [0 0 0];
                            end
                        end
                        
                case 2 %Response 2 only

                        switch stimVector(1)
                        case 1
                            switch stimVector(2)
                                case 1
                                    rAct(2) = 1;
                                case 0
                                    rAct(2) = 1;
                            end
                        case 0
                            switch stimVector(2)
                                case 1
                                    rAct(2) = 0;
                                case 0
                                    rAct = [0 0 0];
                            end
                        end
                        
 
                case 3 %Plus maze initialization: S3->R1, S4->R2

                    switch contextVector(1)
                        case 1
                            rAct(1)=1;
                        case 0
                    end      
                    switch contextVector(2)
                        case 1
                            rAct(2)=1;
                        case 0
                    end                            
                         
                 case 4 %Plus maze training block: S1+S3->R1
                    switch stimVector(1)
                        case 1
                           switch contextVector(1)
                                case 1
                                    rAct(1) = 1;
                            end           
                    end
             case 5 %Plus maze testing block: S1+S4->R2
                    switch stimVector(1)
                        case 1
                            switch contextVector(2)
                                case 1
                                    rAct(2) = 1;
                            end
                                              
                    end
                       

                case 6
                    switch stimVector(1)
                        case 1
                            switch stimVector(2)
                                case 1          
                                case 0
                                    
                                    switch contextVector(1)
                                        case 1
                                            if (RandomVar<chance)
                                            else
                                            end
                                                  rAct(3) = 1;
                                        case 0
                                            
                                            if (RandomVar<chance)
                                            else
                                            end
                                           rAct(1) = 1;
                                    end                         
                            end
                        case 0
                            switch stimVector(2)
                                case 1
                                    switch contextVector(1)
                                        case 1
                                            if (RandomVar<chance)
                                            else
                                            end
                                                  rAct(2) = 1;
                                        case 0
                                            
                                            if (RandomVar<chance)
                                            else
                                            end
                                           rAct(2) = 1;
                                    end      
                                    
                                    switch contextVector(2)
                                        case 1
                                            if(RandomVar<3)                                                 
                                               
                                            else
                                             
                                            end
                                        case 0
                                            if (RandomVar<chance)
                                            else
                                              
                                            end
                                    end
                          
                                case 0
                                    rAct = [0 0 0];
                            end
                    end
                    
               case 7 %Response 2 only

                        switch stimVector(1)
                        case 1
                            switch stimVector(2)
                                case 1
                                    rAct(3) = 1;
                                case 0
                                    rAct(3) = 1;
                            end
                        case 0
                            switch stimVector(2)
                                case 1
                                    rAct(3) = 1;
                                case 0
                                    rAct = [0 0 0];
                            end
                        end
                        
                case 8 %Match, random if conjunction
                    
                    switch stimVector(1)
                        case 1
                            switch contextVector(1)
                                case 1
                                    if (rand() < chance)
                                        rAct(3) = 1;
                                    else
                                        rAct(1) = 1;
                                    end
                                case 0
                                    rAct(1)=1;
                            end                         
                        case 0
                            switch contextVector(1)
                                case 1
                                    rAct(3) = 1;
                            end
                           
                    end
                
             case 9 %S1->R1, chance of R2 if S2
                   
                    switch stimVector(1)
                        case 0
                        case 1
                            switch stimVector(2)    
                                case 1
                                    if (rand()<chance)
                                        rAct(2) = 1;
                                    else
                                        rAct(1) = 1;
                                    end
                                case 0
                                    rAct(1)=1;
                            end
                    end
               case 10 %S1->R1 or chance of R3 if S3
                   
                    switch stimVector(1)
                        case 0 
                        case 1
                               switch contextVector(1)
                                    case 1
                                        if (rand()<chance)
                                            rAct(3) = 1;
                                        else
                                            rAct(1) = 1;
                                        end
                                    case 0
                                        rAct(1) = 1;
                               end    
                    end                  
                     
               case 11 
                   

                    if(contextVector(1) == 1)
                         rAct(1) = 1;
                    else
                        rAct(2) = 1;
                    end

               
                case 12 
                   
                    
                    if(contextVector(1) == 1)
                        rAct(3) = 1;
                    else
                        rAct(1) = 1;
                    end
                    
                    if(stimVector(1) == 1)
                        if (rand() < 0.9)
                            
                            
                            
                        else
                            rAct(2) = 1;
                        end
                    else
                        if (rand() < 0.9)
                            
                        else
                            rAct(2) = 1;
                        end
                    end
                    
                case 13
                   rAct(1)=1;
                   
                case 14 
                   rAct(2)=1;
                   
                case 15 %Extinction in WB
                   
                    rAct = [-1 0 -1]; 
                    
                    if rand() < chance
                            rAct(1) = 1;
                    end


                        
                case 16
                    if (stimVector(1) == 1 && rand() < chance)
                        rAct(2)=1;
                    else
                        rAct(2)=0;
                    end
                    
                case 17 %Resp 1 extinction in savings paradigm
                    rAct = -1*[1 1 0]; 
                    
                    if (rand() < chance)
                    rAct(1) = -1;
                    end
                    
                case 18 %Stim 1 reacquisition
                    rAct = -1*[1 1 1]; 
                    
                    if(rand() < chance)
                    rAct(1) = 1;
                    end                   
                    
                case 19 %Stim 1 initial acquisition
                    rAct = [1 0.0 0] + 0.0*[0 rand() rand()];
%                     rAct(1) = 1;
                    
                case 20 %Stim 1 training
                    rAct = [1 0.1 0.1];
            end
            
            
            
end   
    