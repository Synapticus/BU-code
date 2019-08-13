function [] = LearningSysPlotScript(varargin)
%Script used to generate figures for thesis. Works using the data generated
%by LearningSysScript.

addpath('F:\\Dropbox\\Science\\2017 paper\\')
num = 1;
subjects = 20;


%Set file path
if (nargin > 0)
    fpath = varargin{1};
else
    fpath = 'F:\\Dropbox\\Science\\Figures\\';
end

%Set number of experiments
if (nargin > 1)
    num = varargin{2};
end

%Set number of subjects
if (nargin > 2)
    subjects = varargin{3};
end

%Check for figure-specific toggle flags
if(exist('savingsFig','var') == 0)
    savingsFig = 0;
end

if(exist('reversalPlot','var') == 0)
    reversalPlot = 0;
end

if(exist('woodsBoutonFig','var') == 0)
    woodsBoutonFig = 0;
end

plotWeights = 1;
plotPercentages = 0;
plotTimes = 0;

% Load data
destpath = sprintf('%s%sExperiment%dSubject%d.mat',fpath,date,num,subjects);
load(destpath,'-regexp','^(?!fpath$|varargin$|num$|subjects$).')

colorSpan = 0:0.1:1;
blacks = [zeros(size(colorSpan)); zeros(size(colorSpan)); colorSpan]';
colors = [colorSpan; zeros(size(colorSpan)); 1-colorSpan]';
cmap = [blacks; colors];

%Generate sample statistics
if(subjects > 1)
    
 if(exist('crossing','var') == 0)
    crossing = zeros(6,length(rVector(:,4)));
 end
         
    for expInd = 1:num

        
        sampleAvgPath = sprintf('%s%ssampleAvgExperiment%d.mat',fpath,date,expInd);
        sampleAvgResponses = zeros(size(pastResponses));
         sampleAvgWeightsPGO = zeros(size(wVector(:,1,1,1,pSweep)));
         sampleAvgWeightsPNOGO = zeros(size(wVector(:,1,1,1,pSweep)));
         sampleAvgWeightsSGO = zeros(size(wVector(:,1,1,1,pSweep)));
         sampleAvgWeightsSNOGO = zeros(size(wVector(:,1,1,1,pSweep)));
         
         sampleAvgWeightsPGOAlt = zeros(size(wVector(:,1,2,1,pSweep)));
         sampleAvgWeightsPNOGOAlt = zeros(size(wVector(:,1,2,1,pSweep)));
         sampleAvgWeightsSGOAlt = zeros(size(wVector(:,1,2,1,pSweep)));
         sampleAvgWeightsSNOGOAlt = zeros(size(wVector(:,1,2,1,pSweep)));
         

         sampleAvgCrossing = zeros(1,6);
         sampleAvgCpt = zeros(size(rVector(:,4)));
         
         eligibleSubjs = zeros(1,6);
            for subjInd = 1:subjects
                
                subjpath = sprintf('%s%sExperiment%dSubject%d.mat',fpath,date,expInd,subjInd);
                load(subjpath,'-regexp','^(?!fpath$|varargin$|num$|subjects$).')


                
                cptTemp = rVector(:,4);
                cptTemp(cptTemp == 2) = 0;
                
                sampleAvgResponses = sampleAvgResponses + pastResponses;
                sampleAvgWeightsPGO = sampleAvgWeightsPGO + wVector(:,1,1,1,pSweep); 
                sampleAvgWeightsPNOGO = sampleAvgWeightsPNOGO + wVector(:,2,1,1,pSweep); 
                sampleAvgWeightsSGO = sampleAvgWeightsSGO + wVector(:,3,1,1,pSweep); 
                sampleAvgWeightsSNOGO = sampleAvgWeightsSNOGO + wVector(:,4,1,1,pSweep); 
                
                sampleAvgWeightsPGOAlt = sampleAvgWeightsPGOAlt + wVector(:,1,2,1,pSweep); 
                sampleAvgWeightsPNOGOAlt = sampleAvgWeightsPNOGOAlt + wVector(:,2,2,1,pSweep); 
                sampleAvgWeightsSGOAlt = sampleAvgWeightsSGOAlt + wVector(:,3,2,1,pSweep); 
                sampleAvgWeightsSNOGOAlt = sampleAvgWeightsSNOGOAlt + wVector(:,4,2,1,pSweep); 
                
                for blockInd = 1:6
                    if (~isnan(crossing(blockInd)))
                     sampleAvgCrossing(blockInd) = sampleAvgCrossing(blockInd) + crossing(blockInd);
                     eligibleSubjs(blockInd) = eligibleSubjs(blockInd)+1;
                    end
                end
                sampleAvgCpt = sampleAvgCpt + cptTemp;
            end
        sampleAvgResponses = sampleAvgResponses ./ subjects;
        sampleAvgWeightsPGO = sampleAvgWeightsPGO ./ subjects;
        sampleAvgWeightsPNOGO = sampleAvgWeightsPNOGO ./ subjects;
        sampleAvgWeightsSGO = sampleAvgWeightsSGO ./ subjects;
        sampleAvgWeightsSNOGO = sampleAvgWeightsSNOGO ./ subjects;
        
        sampleAvgWeightsPGOAlt = sampleAvgWeightsPGOAlt ./ subjects;
        sampleAvgWeightsPNOGOAlt = sampleAvgWeightsPNOGOAlt ./ subjects;
        sampleAvgWeightsSGOAlt = sampleAvgWeightsSGOAlt ./ subjects;
        sampleAvgWeightsSNOGOAlt = sampleAvgWeightsSNOGOAlt ./ subjects;
        
        sampleAvgCrossing = floor(sampleAvgCrossing ./ eligibleSubjs);
        sampleAvgCpt = sampleAvgCpt ./ subjects;
        save(sampleAvgPath)
    end

end

M = padamount+1;
successCriterion = 0.9;
transferCriterion = 0.9;


%% Plot weight trajectories
if (plotWeights ~= 0)
                     %Draw markers on the weight that was in control of
                     %beheavior
                     primaryControlledCorrectTrials = find((rVector(:,4) == 1) ...
                         & rVector(:,2) == rVector(:,5) & (rVector(:,1) > 0));
                     primaryControlledIncorrectTrials = find((rVector(:,4) == 1) ...
                         & rVector(:,2) == rVector(:,5) & (rVector(:,1) <= 0));
                     secondaryControlledCorrectTrials = find((rVector(:,4) == 2)...
                         & (rVector(:,3) == rVector(:,5)) & (rVector(:,1) > 0));
                     secondaryControlledIncorrectTrials = find((rVector(:,4) == 2) ...
                         & (rVector(:,3) == rVector(:,5)) & (rVector(:,1) <= 0));
                     
                    pCCT = EveryNthElement(primaryControlledCorrectTrials,M,2);
                    pCIT = EveryNthElement(primaryControlledIncorrectTrials,M,2);
                    sCCT = EveryNthElement(secondaryControlledCorrectTrials,M,2);
                    sCIT = EveryNthElement(secondaryControlledIncorrectTrials,M,2);
                    
                    mSize = 8;
%                     minDeflection = 0.0;                        
                    h = zeros(1,4);      
    respInds = [1 2];                
    stimInds = [1 2 3 4];
    for respInd = respInds
        for stimInd = stimInds
            figure();
            hold on
            imageH = image(ones(T,T,3));
            yMaxVal = wMaxHi+0.2;
            set(imageH,'xdata',[1 T],'ydata',[-0.2 yMaxVal]);
            set(gca,'xlim',[1 T],'ylim',[-0.2 yMaxVal]);
            setBackground(imageH,ruleMatrix);

            for weightLineInd = 1:4
                      activePCCT = pCCT(rVector(pCCT,2) == respInd);
                      activePCIT = pCIT((rVector(pCIT,2) == respInd));
                      activeSCCT = sCCT((rVector(sCCT,3) == respInd));
                      activeSCIT = sCIT((rVector(sCIT,3) == respInd));
                      mString = strcat('k',markerStringByResponse(respInd));
                      switch weightLineInd
                        case 1 %Primary compartment active GO weights
                            
                            hold on
                            h(weightLineInd) = plot(1:T,wVector(1:T,weightLineInd,respInd,stimInd,pSweep),'k-','LineWidth',3);%Weight trajectories
                            set(h(weightLineInd),'Color','r');
                            set(h(weightLineInd),'LineStyle','-');
                            
                            plot(activePCCT,wVector(activePCCT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',mSize);
                            plot(activePCIT,wVector(activePCIT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',mSize);
                        case 2 %Primary compartment active NOGO weights
                            
                            hold on
                            h(weightLineInd) = plot(1:T,wVector(1:T,weightLineInd,respInd,stimInd,pSweep),'k-','LineWidth',3);%Weight trajectories
                            set(h(weightLineInd),'Color','r');
                            set(h(weightLineInd),'LineStyle','-.');       
                            plot(activePCCT,wVector(activePCCT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',mSize);
                            plot(activePCIT,wVector(activePCIT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',mSize);
                        case 3 %Secondary compartment active GO weights   
                            
                            hold on
                       
                            h(weightLineInd) = plot(1:T,wVector(1:T,weightLineInd,respInd,stimInd,pSweep),'k-','LineWidth',3);%Weight trajectories
                            plot(activeSCCT,wVector(activeSCCT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',mSize);
                            plot(activeSCIT,wVector(activeSCIT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',mSize);     
                            set(h(weightLineInd),'Color','b');
                            set(h(weightLineInd),'LineStyle','-');
                        case 4 %Secondary compartment active NOGO weights
                            
                            hold on
                       
                            h(weightLineInd) = plot(1:T,wVector(1:T,weightLineInd,respInd,stimInd,pSweep),'k-','LineWidth',3);%Weight trajectories
                            plot(activeSCCT,wVector(activeSCCT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',mSize);
                            plot(activeSCIT,wVector(activeSCIT,weightLineInd,respInd,stimInd,pSweep),...
                            mString,'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',mSize);          
                            set(h(weightLineInd),'Color','b');
                            set(h(weightLineInd),'LineStyle','-.'); 
                      end
                     
            end
        
            titleTxt = sprintf('Weight trajectories: Response %d Stimulus %d',respInd,stimInd);
            title(titleTxt)
            ylabel('Weight value');
            xlabel('Trials');
            
            % Make marker legend entries
            m(1) = plot(-10,-10);
            set(m(1),'Marker','o'); 
            set(m(1),'MarkerFaceColor','green');
            set(m(1),'MarkerEdgeColor','black');
            m(2) = plot(-10,-10);
            set(m(2),'Marker','o'); 
            set(m(2),'MarkerEdgeColor','black');
            set(m(2),'MarkerFaceColor','red');           

            legend([h m],'Go, Primary','NOGO, Primary','GO, Secondary','NOGO, Secondary','Controlled rewarded decision','Controlled un-rewarded decision','Location','SouthOutside','FontSize',16);
            
            [rmL , ~] = size(ruleMatrix);
            p = zeros(rmL,1);
            t = 1;
            for j = 1:rmL
                p(j) = ruleMatrix(j,1);
            tAdj = t+0.1*p(j);
            descTxt = descTextByRule(ruleMatrix(j,2),chanceMatrix(j,2),chanceMatrix(j,1));
            text('units','data','position',[tAdj yMaxVal],'FontSize',16,'FontWeight','bold','VerticalAlignment','top','string',descTxt);

            t=t+p(j);
            end
                ax= gca;
            ax.YTick = [0, wOffset, wMaxLo, wMaxHi];
            ax.YTickLabel = {'0',sprintf('Initial weight: %0.2f',wDblOffset),sprintf('Secondary weight limit: %0.1f',wMaxLo),sprintf('Primary weight limit: %0.1f',wMaxHi)};
            % Print to file
            ptext = sprintf('%s%sWeightsResponseStim%dResp%dparam%d.png',fpath,date,stimInd,respInd,pSweep);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            print(gcf,'-dpng',ptext,'-r0');
            hold off
            close
        end
    end


end


 if(plotPercentages == 1)   
    colorSpan = 0:0.1:1;
    blacks = [zeros(size(colorSpan)); zeros(size(colorSpan)); colorSpan]';
    colors = [colorSpan; zeros(size(colorSpan)); 1-colorSpan]';
    cmap = [blacks; colors];
%% Percentage of responses produced figure
    figure()
    hold on
    imageH = image(ones(T,T,3));
    set(imageH,'xdata',[1 T],'ydata',[-10 110]);
    set(gca,'xlim',[1 T],'ylim',[-10 110]);
    setBackground(imageH,ruleMatrix);

    h= plot(1:T,100*pastResponses(1,:),1:T,100*pastResponses(2,:),1:T,100*pastResponses(3,:));
    legend(h,'Response 1','Response 2','Response 3','Location','SouthOutside','FontSize',20)

    
    xlabel('Trial')
    ylabel('Percentage of responses')
    title(sprintf('Responses produced (%d-trial moving window)',padamount+1))
    ax= gca;
    ax.YTick = [0 25 50 75 100];
    set(h,'LineWidth',3);
    % Print to file
    ptext = sprintf('%s%sPctResponsesFig_ParamSet%d.png',fpath,date,pSweep);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    close
    
    %% Color coded compartment figure
    figure()
    hold on
    imageH = image(ones(T,T,3));
    set(imageH,'xdata',[1 T],'ydata',[-10 110]);
    set(gca,'xlim',[1 T],'ylim',[-10 110]);
    setBackground(imageH,ruleMatrix);
    surface([1:length(100*pctPrimaryDecision);1:length(100*pctPrimaryDecision)],[100*pctPrimaryDecision;100*pctPrimaryDecision],[zeros(size(100*pctPrimaryDecision));zeros(size(100*pctPrimaryDecision))],[pctPrimaryDecision;pctPrimaryDecision],...
    'facecol','no',...
    'edgecol','interp',...
    'LineWidth',4);
	colormap(colors)
    caxis([0 1]);
    colorbar('Ticks',[0.1, 0.9],...
        'TickLabels',{'Secondary','Primary'});
    title(sprintf('Transition of control to primary compartment (%d-trial moving window)',padamount+1))
    ylabel('% of trials controlled by primary compartment')
    
      ptext = sprintf('%s%sCptControlColor_ParamSet%d.png',fpath,date,pSweep);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    close
%% Percentage correct decisions figure
    figure()
    hold on
    imageH = image(ones(T,T,3));
    YMaxVal = 125;
    set(imageH,'xdata',[1 T],'ydata',[-20 YMaxVal]);
    set(gca,'xlim',[1 T],'ylim',[-20 YMaxVal]);
    setBackground(imageH,ruleMatrix);
    % Plot line
%     plot(100*pctSuccesses,'LineWidth',5); 
    surface([1:length(100*pctSuccesses);1:length(100*pctSuccesses)],[100*pctSuccesses;100*pctSuccesses],[zeros(size(100*pctSuccesses));zeros(size(100*pctSuccesses))],[pctPrimaryDecision;pctPrimaryDecision],...
        'facecol','no',...
        'edgecol','interp',...
        'LineWidth',4);
    colormap(cmap);
    caxis([-1 1]);
    h= colorbar('Ticks',[-0.9, 0, 0.5, 0.9],...
        'TickLabels',{'Indecisive','Secondary-biased','Unbiased','Primary-biased'});
%     set(h,'YDir','reverse');
    
    [rmL , ~] = size(ruleMatrix);
    p = zeros(rmL,1);
    t = 1;
    for j = 1:rmL
       p(j) = ruleMatrix(j,1);
       tAdj = t+0.1*p(j);
        belowCriterionIndex = find(pctSuccesses(t:t+p(j)-1) < (1-successCriterion), 1, 'last');

        [minVal, minIndex] = min(pctSuccesses(t:t+p(j)-1));
            pastCriterionIndex = find(pctSuccesses(t+minIndex:t+p(j)-1) > successCriterion,1)+minIndex;
               critTxt = sprintf('%0.1f%% accuracy \n %d trials after rule change',100*successCriterion, pastCriterionIndex);
               
               text('units','data','position',[tAdj -10],'FontSize',14,'FontWeight','bold','string',critTxt);


%         end
        
    txt = descTextByRule(ruleMatrix(j,2),chanceMatrix(j,2),chanceMatrix(j,1));      
    text('units','data','position',[tAdj YMaxVal],'FontSize',12,'FontWeight','bold','VerticalAlignment','top','string',txt); 
        t = t+p(j);
    end
    title('Model Performance')
    ylabel(sprintf('Percentage of past %d trials ending in reward',padamount));
    xlabel('Trial number')
    ax= gca;
    ax.YTick = [0, 25, 50, 75, 100];
    % Print to file
    ptext = sprintf('%s%sPerformanceFig_ParamSet%d.png',fpath,date,pSweep);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    close
    
%% Figure showing compartment % of trials in control
    figure()
    set(gcf,'PaperPositionMode','auto');    
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    hold on
    clear h
    axis([1 T -1 1]);
    imageH = image(ones(T,T,3));
    set(imageH,'xdata',[1 T],'ydata',[-1.25 1.75]);
    set(gca,'xlim',[1 T],'ylim',[-1.25 1.75]);
    setBackground(imageH,ruleMatrix);    
    
    [rmL, ~] = size(ruleMatrix);
    p = zeros(rmL,1);
    t=1;
    for j = 1:rmL
    p(j) = ruleMatrix(j,1);
    tAdj = t+0.1*p(j);
    
    [~, aboveCriterionIndex] = find(pctPrimaryDecision(t:(t+p(j)-1)) > transferCriterion , 1,'last');
    [~, belowCriterionIndex] = find( pctPrimaryDecision(t:t+p(j)-1) < (1-transferCriterion), 1, 'last');
    [~, zeroIndex] = find ( abs( pctPrimaryDecision(t:t+p(j)-1)-0.5) < 0.05 , 1, 'last');
        
    txt = descTextByRule(ruleMatrix(j,2),chanceMatrix(j,2),chanceMatrix(j,1));    
    text('units','data','position',[tAdj 1.5],'FontSize',12,'FontWeight','bold','VerticalAlignment','top','string',txt);

    
     t = t + p(j);
    end
    
    plot(pctPrimaryDecision,'LineWidth',5)
    titleTxt = sprintf('Sliding window average (N=%d) of controller by trial',padamount+1);
    title(titleTxt)
    ax= gca;
    ax.YTick = [-1, 0, 0.5, 1];
    ax.YTickLabel = {'Indecisive','Secondary-biased','Unbiased','Primary-biased'};
    ylabel('Control of behavior by compartment')
    xlabel('Trials')
    hold off
    
    % Print to file
    ptext = sprintf('%s%sCompartmentRatioparam%d.png',fpath,date,pSweep);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    
    plotResponseAvg = 0;
    
    if (plotResponseAvg == 1)
        
            figure()
            hold on
            for trialBlock = 1:rmL
                if (trialBlock == 1)
                    responses = rVector(1:ruleMatrix(trialBlock,1),5);
                else
                    responses = rVector(sum(ruleMatrix(1:trialBlock-1,1))+1:sum(ruleMatrix(1:trialBlock,1)),5);
                end
                subplot(1,rmL,trialBlock)
                h = histogram(responses,'Normalization','probability');
                h.NumBins = 3;
                ax = gca;
                ax.XTick = 1:3;
                ax.YLim = [0 1];
                xlabel('Response number')
                titleText = sprintf('Trial block %d',trialBlock);
         
                title(titleText)
            end
            
            % Print to file
            ptext = sprintf('%s%sResponseHistogramparam%d.png',fpath,date,pSweep);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            print(gcf,'-dpng',ptext,'-r0');
            hold off
            close
    end
end

if (plotDraw(6) == 1)
    figure()
    acquisition = 100*pctSuccesses(1:ruleMatrix(1,1));
    reacquisition = 100*pctSuccesses((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
    h = plot((1:ruleMatrix(1,1))-1,acquisition,'--',((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)))-sum(ruleMatrix(1:2,1))-1,reacquisition,'-');
    set(h,'LineWidth',3);
    legend('Initial acquisition','Reacquisition','Location','SouthOutside')
    xlabel('Trials after rule change')
    ylabel('% conditioned response')
    titleTxt = sprintf('Reacquisition following extinction');
    title(titleTxt); 
    ptext = sprintf('%s%sSavingsComparisonFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    % Weights
    figure()
    acquisitionGO = wVector(1:ruleMatrix(1,1),1,1,1,pSweep);
    acquisitionNOGO = wVector(1:ruleMatrix(1,1),2,1,1,pSweep);
    reacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
    reacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);
    subplot(211)
    h = plot((1:ruleMatrix(1,1))-1,acquisitionGO,'--',((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)))-sum(ruleMatrix(1:2,1))-1,reacquisitionGO,'-');
    ylabel('GO Weight')
    legend('Initial acquisition','Reacquisition','Location','EastOutside')
    titleTxt = sprintf('Reacquisition following extinction');
    title(titleTxt); 
    set(h,'LineWidth',3);
    set(gca,'ylim',[0 0.4]);
    subplot(212)
    h2 = plot((1:ruleMatrix(1,1))-1,acquisitionNOGO,'--',((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)))-sum(ruleMatrix(1:2,1))-1,reacquisitionNOGO,'-');
    set(h2,'LineWidth',3);
    legend('Initial acquisition','Reacquisition','Location','EastOutside')
    xlabel('Trials after rule change')
    ylabel('NOGO Weight')
    set(gca,'ylim',[0 0.4]);
    
    ptext = sprintf('%s%sSavingsComparisonWeightsFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
end

    leftWindow = 0;
    rightWindow = 500;
    indicesBlock1 = ruleMatrix(1,1)-leftWindow:ruleMatrix(1,1)+rightWindow;
    indicesBlock2 = sum(ruleMatrix(1:3))-leftWindow:sum(ruleMatrix(1:3))+rightWindow;
    indicesBlock3 = sum(ruleMatrix(1:5))-leftWindow:sum(ruleMatrix(1:5))+rightWindow;
    
    
if (reversalPlot ==1)
    figure()

    

    subplot(311)
    hold on
    line([0 0],[0 100],'LineStyle',':','LineWidth',3,'Color','black');
    surface([-leftWindow:rightWindow;-leftWindow:rightWindow],...
        [100*pctSuccesses(indicesBlock1);100*pctSuccesses(indicesBlock1)],...
        [zeros(size(pctSuccesses(indicesBlock1)));zeros(size(pctSuccesses(indicesBlock1)))],...
        [pctPrimaryDecision(indicesBlock1);pctPrimaryDecision(indicesBlock1)],...
    'facecol','no',...
    'edgecol','interp',...
    'LineWidth',4);
    colormap(colors)
    caxis([0 1]);
    ylabel(sprintf('%% of past %d trials ending in reward',M))
    title('Reinstatement after repeated reversals')
    
    subplot(312)
    
    hold on
    line([0 0],[0 100],'LineStyle',':','LineWidth',3,'Color','black');
            surface([-leftWindow:rightWindow;-leftWindow:rightWindow],...
            [100*pctSuccesses(indicesBlock2);100*pctSuccesses(indicesBlock2)],...
            [zeros(size(pctSuccesses(indicesBlock2)));zeros(size(pctSuccesses(indicesBlock2)))],...
            [pctPrimaryDecision(indicesBlock2);pctPrimaryDecision(indicesBlock2)],...
        'facecol','no',...
        'edgecol','interp',...
        'LineWidth',4);
        colormap(colors)
    caxis([0 1]);
    ylabel(sprintf('%% of past %d trials ending in reward',M))
    
    subplot(313)
    
    hold on
    line([0 0],[0 100],'LineStyle',':','LineWidth',3,'Color','black');
            surface([-leftWindow:rightWindow;-leftWindow:rightWindow],...
            [100*pctSuccesses(indicesBlock3);100*pctSuccesses(indicesBlock3)],...
            [zeros(size(pctSuccesses(indicesBlock3)));zeros(size(pctSuccesses(indicesBlock3)))],...
            [pctPrimaryDecision(indicesBlock3);pctPrimaryDecision(indicesBlock3)],...
        'facecol','no',...
        'edgecol','interp',...
        'LineWidth',4);
        colormap(colors)
    caxis([0 1]);
    xlabel('Trials')
    ylabel(sprintf('%% of past %d trials ending in reward',M))


    ptext = sprintf('%s%sReversalFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    

end

if (reversalPlot >= 2)
    figure()



    subplot(311)

    
    if (reversalPlot == 3)
        if (subjects > 1)
            subjpathTemp = sprintf('%s%ssampleAvgExperiment1.mat',fpath,date);
            clear sampleAvg* crossing
            load(subjpathTemp,'sampleAvg*');
            crossing = sampleAvgCrossing;
            pctSuccesses = sampleAvgResponses(2,:);
            pctPrimaryDecision = sampleAvgCpt;
            
        else
            destpathTemp = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
            clear pct* crossing
            load(destpathTemp,'pct*','crossing');
        end

            hold on
    % Plot baseline accuracy
         plot(-leftWindow:rightWindow,100*pctSuccesses(indicesBlock1),'k:','LineWidth',2) 

    % Plot baseline crossings
        line([crossing(1)-ruleMatrix(1,1) crossing(1)-ruleMatrix(1,1)],[0 100],'LineStyle','-.','LineWidth',2,'Color',[0.75 0.75 0.75]);
        line([crossing(2)-ruleMatrix(1,1) crossing(2)-ruleMatrix(1,1)],[0 100],'LineStyle','-.','LineWidth',2,'Color',[0.75 0.75 0.75]);

        plot(crossing(1)-ruleMatrix(1,1), 105,'kv','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
        plot(crossing(2)-ruleMatrix(1,1), 105,'kv','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
        

    else

    end

    if (subjects > 1)
    subjpathTemp = sprintf('%s%ssampleAvgExperiment2.mat',fpath,date);
    load(subjpathTemp,'sampleAvg*');
    crossing = sampleAvgCrossing;
    pctSuccesses = sampleAvgResponses(2,:);
    pctPrimaryDecision = sampleAvgCpt;

    else
    destpathTemp = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
    clear pct* crossing
    load(destpathTemp,'pct*','crossing');
    end
        
    h1 = surface([-leftWindow:rightWindow;-leftWindow:rightWindow],...
        [100*pctSuccesses(indicesBlock1);100*pctSuccesses(indicesBlock1)],...
        [zeros(size(pctSuccesses(indicesBlock1)));zeros(size(pctSuccesses(indicesBlock1)))],...
        [pctPrimaryDecision(indicesBlock1)';pctPrimaryDecision(indicesBlock1)'],...
    'facecol','no',...
    'edgecol','interp',...
    'LineWidth',4);
    colormap(colors)
    caxis([0 1]);
    
    h3 = plot(crossing(1)-ruleMatrix(1,1), 105,'kv','MarkerFaceColor','k','MarkerSize',10);
    h4 = plot(crossing(2)-ruleMatrix(1,1), 105,'kv','MarkerFaceColor','k','MarkerSize',10);
 
    title('Reinstatement after repeated reversals')
    axis([-leftWindow rightWindow 0 105]);
    
    line([crossing(1)-ruleMatrix(1,1) crossing(1)-ruleMatrix(1,1)],[0 100],'LineStyle',':','LineWidth',2,'Color','black');
    line([crossing(2)-ruleMatrix(1,1) crossing(2)-ruleMatrix(1,1)],[0 100],'LineStyle',':','LineWidth',2,'Color','black');
    line([crossing(3)-sum(ruleMatrix(1:3)) crossing(3)-sum(ruleMatrix(1:3))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(4)-sum(ruleMatrix(1:3)) crossing(4)-sum(ruleMatrix(1:3))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(5)-sum(ruleMatrix(1:5)) crossing(5)-sum(ruleMatrix(1:5))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(6)-sum(ruleMatrix(1:5)) crossing(6)-sum(ruleMatrix(1:5))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    

    
    
    subplot(312)
    
    hold on
    
        if (reversalPlot == 3)
            
        if (subjects > 1)
            subjpathTemp = sprintf('%s%ssampleAvgExperiment1.mat',fpath,date);
            clear sampleAvg* crossing
            load(subjpathTemp,'sampleAvg*');
            crossing = sampleAvgCrossing;
            pctSuccesses = sampleAvgResponses(2,:);
            pctPrimaryDecision = sampleAvgCpt;
        else
            destpathTemp = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
            clear pct* crossing
            load(destpathTemp,'pct*','crossing');
        end    
            
            

        hold on
        plot(-leftWindow:rightWindow,100*pctSuccesses(indicesBlock2),'k:','LineWidth',2) 

        line([crossing(3)-sum(ruleMatrix(1:3)) crossing(3)-sum(ruleMatrix(1:3))],[0 100],'LineStyle','-.','LineWidth',2,'Color',[0.75 0.75 0.75]);
        line([crossing(4)-sum(ruleMatrix(1:3)) crossing(4)-sum(ruleMatrix(1:3))],[0 100],'LineStyle','-.','LineWidth',2,'Color',[0.75 0.75 0.75]);

        plot(crossing(3)-sum(ruleMatrix(1:3)), 105,'kv','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
        plot(crossing(4)-sum(ruleMatrix(1:3)), 105,'kv','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);

        else
        h2 = plot(-leftWindow:rightWindow,100*desiredResponses(indicesBlock1),'k--');
        set(h2,'LineWidth',2);
        end
    

        if (subjects > 1)
            subjpathTemp = sprintf('%s%ssampleAvgExperiment2.mat',fpath,date);
            clear sampleAvg* crossing
            load(subjpathTemp,'sampleAvg*');
            crossing = sampleAvgCrossing;
            pctSuccesses = sampleAvgResponses(2,:);
            pctPrimaryDecision = sampleAvgCpt;
        else
            destpathTemp = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
            clear pct* crossing
            load(destpathTemp,'pct*','crossing');
        end    
        
            surface([-leftWindow:rightWindow;-leftWindow:rightWindow],...
            [100*pctSuccesses(indicesBlock2);100*pctSuccesses(indicesBlock2)],...
            [zeros(size(pctSuccesses(indicesBlock2)));zeros(size(pctSuccesses(indicesBlock2)))],...
            [pctPrimaryDecision(indicesBlock2)';pctPrimaryDecision(indicesBlock2)'],...
        'facecol','no',...
        'edgecol','interp',...
        'LineWidth',4);
        colormap(colors)
    caxis([0 1]);

    plot(crossing(3)-sum(ruleMatrix(1:3)), 105,'kv','MarkerFaceColor','k','MarkerSize',10);
    plot(crossing(4)-sum(ruleMatrix(1:3)), 105,'kv','MarkerFaceColor','k','MarkerSize',10);

    line([crossing(1)-ruleMatrix(1,1) crossing(1)-ruleMatrix(1,1)],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(2)-ruleMatrix(1,1) crossing(2)-ruleMatrix(1,1)],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(3)-sum(ruleMatrix(1:3)) crossing(3)-sum(ruleMatrix(1:3))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0 0 0]);
    line([crossing(4)-sum(ruleMatrix(1:3)) crossing(4)-sum(ruleMatrix(1:3))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0 0 0]);
    line([crossing(5)-sum(ruleMatrix(1:5)) crossing(5)-sum(ruleMatrix(1:5))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(6)-sum(ruleMatrix(1:5)) crossing(6)-sum(ruleMatrix(1:5))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    
    ylabel(sprintf('%% of past %d trials ending in reward',M))
    axis([-leftWindow rightWindow 0 105]);
    
    
    subplot(313)
    
    hold on
    if (reversalPlot == 3)
            
        if (subjects > 1)
            subjpathTemp = sprintf('%s%ssampleAvgExperiment1.mat',fpath,date);
            load(subjpathTemp,'sampleAvg*');
            crossing = sampleAvgCrossing;
            pctSuccesses = sampleAvgResponses(2,:);
            pctPrimaryDecision = sampleAvgCpt;
        else
            destpathTemp = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
            clear pct* crossing
            load(destpathTemp,'pct*','crossing');
        end    
        
        hold on
         plot(-leftWindow:rightWindow,100*pctSuccesses(indicesBlock3),'k:','LineWidth',2) 

        line([crossing(5)-sum(ruleMatrix(1:5)) crossing(5)-sum(ruleMatrix(1:5))],[0 100],'LineStyle','-.','LineWidth',2,'Color',[0.75 0.75 0.75]);
        line([crossing(6)-sum(ruleMatrix(1:5)) crossing(6)-sum(ruleMatrix(1:5))],[0 100],'LineStyle','-.','LineWidth',2,'Color',[0.75 0.75 0.75]);

        plot(crossing(5)-sum(ruleMatrix(1:5)), 105,'kv','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
        plot(crossing(6)-sum(ruleMatrix(1:5)), 105,'kv','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
        
    end

        if (subjects > 1)
            subjpathTemp = sprintf('%s%ssampleAvgExperiment2.mat',fpath,date);
            load(subjpathTemp,'sampleAvg*');
            crossing = sampleAvgCrossing;
            pctSuccesses = sampleAvgResponses(2,:);
            pctPrimaryDecision = sampleAvgCpt;
        else
            destpathTemp = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
            clear pct* crossing
            load(destpathTemp,'pct*','crossing');
        end  
        
    surface([-leftWindow:rightWindow;-leftWindow:rightWindow],...
            [100*pctSuccesses(indicesBlock3);100*pctSuccesses(indicesBlock3)],...
            [zeros(size(pctSuccesses(indicesBlock3)));zeros(size(pctSuccesses(indicesBlock3)))],...
            [pctPrimaryDecision(indicesBlock3)';pctPrimaryDecision(indicesBlock3)'],...
        'facecol','no',...
        'edgecol','interp',...
        'LineWidth',4);
        colormap(colors)
    caxis([0 1]);
    xlabel('Trials')
    
    h3 = plot(crossing(5)-sum(ruleMatrix(1:5)), 105,'kv','MarkerFaceColor','k','MarkerSize',10);
    h4 = plot(crossing(6)-sum(ruleMatrix(1:5)), 105,'kv','MarkerFaceColor','k','MarkerSize',10);

    line([crossing(1)-ruleMatrix(1,1) crossing(1)-ruleMatrix(1,1)],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(2)-ruleMatrix(1,1) crossing(2)-ruleMatrix(1,1)],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(3)-sum(ruleMatrix(1:3)) crossing(3)-sum(ruleMatrix(1:3))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(4)-sum(ruleMatrix(1:3)) crossing(4)-sum(ruleMatrix(1:3))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    line([crossing(5)-sum(ruleMatrix(1:5)) crossing(5)-sum(ruleMatrix(1:5))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0 0 0]);
    line([crossing(6)-sum(ruleMatrix(1:5)) crossing(6)-sum(ruleMatrix(1:5))],[0 100],'LineStyle',':','LineWidth',2,'Color',[0 0 0]);
    
    
    axis([-leftWindow rightWindow 0 105]);

    ptext = sprintf('%s%sReversalFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    
     
  
end


  figure()
        
    clear pct* crossing
    subjpathTemp = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
    load(subjpathTemp,'crossing');
    crBaseline = crossing([2 4 6])-[500 1500 2500];
    
    clear pct* crossing
    subjpathTemp = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
    load(subjpathTemp,'crossing');
    crComparison = crossing([2 4 6])-[500 1500 2500];
        
    b=bar([0 1 2],[crBaseline; crComparison]');
    b(1).FaceColor = 'black';
    b(2).FaceColor = 'white';
    xlabel('Initial acquisition and two reversals');
    ylabel('Trials to habit transition');
    
    title('Cocaine accelerates transitions to habitual control')
    legend('Baseline performance','Asymmetric learning rate enhancement (cocaine)','Location','northeast')
    ptext = sprintf('%s%sReversalBarFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close


if (woodsBoutonFig == 1)
        %Savings figure formatted to match Woods & Bouton 2007

        if (subjects > 1)
            
            clear pctSuccesses ruleMatrix pastResponses destpath   
            sampleAvgPath = sprintf('%s%ssampleAvgExperiment1.mat',fpath,date);
            load(sampleAvgPath)
            
            %Response 1
            PartialExtinction = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            PartialReacquisition = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
            
            %Response 2
            PartialExtinctionAlt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            PartialReacquisitionAlt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
            
            clear pctSuccesses ruleMatrix pastResponses destpath   
            sampleAvgPath = sprintf('%s%ssampleAvgExperiment2.mat',fpath,date);
            load(sampleAvgPath)            
            PartialExtinction8 = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            PartialReacquisition8 = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
            
            PartialExtinction8Alt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            PartialReacquisition8Alt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
    
            clear pctSuccesses ruleMatrix pastResponses destpath   
            sampleAvgPath = sprintf('%s%ssampleAvgExperiment3.mat',fpath,date);
            load(sampleAvgPath)
            FullExtinction = 100*sampleAvgResponses(1,ruleMatrix(1,1)+1:sum(ruleMatrix(1:2,1)) );
            FullReacquisition = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
            
            FullExtinctionAlt = 100*sampleAvgResponses(2,ruleMatrix(1,1)+1:sum(ruleMatrix(1:2,1)) );
            FullReacquisitionAlt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
    
            clear pctSuccesses ruleMatrix pastResponses destpath   
            sampleAvgPath = sprintf('%s%ssampleAvgExperiment4.mat',fpath,date);
            load(sampleAvgPath)
            FullExtinction8 = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            FullReacquisition8 = 100*sampleAvgResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
            
            FullExtinction8Alt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            FullReacquisition8Alt = 100*sampleAvgResponses(2,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
        else
             clear pctSuccesses ruleMatrix pastResponses destpath    
            destpath = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
            load(destpath)

            PartialExtinction = 100*pastResponses(1,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            PartialReacquisition = 100*pastResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

            clear pctSuccesses ruleMatrix pastResponses destpath    
            destpath = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
            load(destpath)

            PartialExtinction8 = 100*pastResponses(1,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            PartialReacquisition8 = 100*pastResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

            clear pctSuccesses ruleMatrix pastResponses destpath
            destpath = sprintf('%s%sExperiment3Subject1.mat',fpath,date);
            load(destpath)

            FullExtinction = 100*pastResponses(1,ruleMatrix(1,1)+1:sum(ruleMatrix(1:2,1)) );
            FullReacquisition = 100*pastResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

            clear pctSuccesses ruleMatrix pastResponses destpath    
            destpath = sprintf('%s%sExperiment4Subject1.mat',fpath,date);
            load(destpath)

            FullExtinction8 = 100*pastResponses(1,(sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
            FullReacquisition8 = 100*pastResponses(1,(sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
        end

    % First block is initial acquisition

    leftBlock = WBextBlock;
    rightBlock = WBreacBlock;
    
    leftBin = floor(leftBlock/8);
    rightBin = floor(rightBlock/2);
    
    figure()
    hold on
    group = 0;
    binnedAvg = zeros(4,11);
    for vec = [FullExtinction8; PartialExtinction8; FullExtinction; PartialExtinction]'
        group = group + 1;
        for bin = 1:8
            binnedAvg(group,bin) = mean(vec((bin-1)*leftBin+1:bin*leftBin));

        end
        binnedAvg(group,11) = vec(1); %"Bin 0" showing end of acquisition
    end
    
    
    
    group = 0;
    for vec = [FullReacquisition8; PartialReacquisition8; FullReacquisition; PartialReacquisition]'
        group = group+1;
        for bin = 1:2
            binnedAvg(group,bin+8) = mean(vec((bin-1)*rightBin+1:bin*rightBin));
        end
    end
    
    L1 = plot(1:10,binnedAvg(1,1:10),'ks-','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    L2 = plot(1:10,binnedAvg(2,1:10),'ks:','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    L3 = plot(1:10,binnedAvg(3,1:10),'ko-','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    L4 = plot(1:10,binnedAvg(4,1:10),'ko:','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    plot(0:1,[binnedAvg(1,11) binnedAvg(1,1)],'ks-','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    plot(0:1,[binnedAvg(2,11) binnedAvg(2,1)],'ks:','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    plot(0:1,[binnedAvg(3,11) binnedAvg(3,1)],'ko-','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);
    plot(0:1,[binnedAvg(4,11) binnedAvg(4,1)],'ko:','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',3);

    
    
    line([8 8],[0 100],'Color','black','LineStyle','--','LineWidth',3);
    
    
    axis([0 10 0 100])
    
    
    legend([L1 L2 L3 L4],'Ext8','Prf8','Ext2','Prf2','Location','North')
    xlabel('Trial blocks')
    ylabel('% conditioned response')
    titleTxt = sprintf('Reacquisition following extinction');
    title(titleTxt); 
    ptext = sprintf('%s%sWoodsBoutonBinnedComparisonFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    
    
    figure()
    hold on
    plot(1:leftBlock,FullExtinction(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'k-','LineWidth',3);
    plot(leftBlock+1:leftBlock+rightBlock,FullReacquisition(1:rightBlock),'k-','LineWidth',3);        
    plot(1:leftBlock,PartialExtinction(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'k:','LineWidth',3);
    plot(leftBlock+1:leftBlock+rightBlock,PartialReacquisition(1:rightBlock),'k:','LineWidth',3);
    
    plot(1:leftBlock,FullExtinction8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'k-','LineWidth',3);
    plot(leftBlock+1:leftBlock+rightBlock,FullReacquisition8(1:rightBlock),'k-','LineWidth',3);        
    plot(1:leftBlock,PartialExtinction8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'k:','LineWidth',3);
    plot(leftBlock+1:leftBlock+rightBlock,PartialReacquisition8(1:rightBlock),'k:','LineWidth',3);

    L3 = plot(1:10:leftBlock,FullExtinction(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,FullReacquisition(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    L4 =plot(1:10:leftBlock,PartialExtinction(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,PartialReacquisition(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    L1 = plot(1:10:leftBlock,FullExtinction8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,FullReacquisition8(1:10:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    L2 = plot(1:10:leftBlock,PartialExtinction8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,PartialReacquisition8(1:10:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    line([leftBlock leftBlock],[0 100],'Color','black','LineStyle','--','LineWidth',3);
    
    axis([0 leftBlock+rightBlock 0 100])
    legend([L1 L2 L3 L4],'Ext8','Prf8','Ext2','Prf2','Location','North')
    xlabel('Trials')
    ylabel('% conditioned response')
    titleTxt = sprintf('Reacquisition following extinction');
    title(titleTxt); 
    ptext = sprintf('%s%sWoodsBoutonComparisonFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
%
%   
  
    %Weights comparison
    figure()
    % First block is initial acquisition

    if (subjects > 1)
        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment1.mat',fpath,date);
        load(destpath)

        partialExtinctionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        partialExtinctionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        partialReacquisitionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        partialReacquisitionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryPartialExtinctionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryPartialExtinctionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryPartialReacquisitionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryPartialReacquisitionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));


        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment2.mat',fpath,date);
        load(destpath)

        partialExtinctionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        partialExtinctionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        partialReacquisitionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        partialReacquisitionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryPartialExtinctionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryPartialExtinctionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryPartialReacquisitionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryPartialReacquisitionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment3.mat',fpath,date);
        load(destpath)

        fullExtinctionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        fullExtinctionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        fullReacquisitionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        fullReacquisitionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryFullExtinctionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryFullExtinctionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryFullReacquisitionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryFullReacquisitionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment4.mat',fpath,date);
        load(destpath)

        fullExtinctionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        fullExtinctionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        fullReacquisitionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        fullReacquisitionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryFullExtinctionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryFullExtinctionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryFullReacquisitionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryFullReacquisitionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
    
    else
        
       clear pctSuccesses pastResponses ruleMatrix destpath wVector
        destpath = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
        load(destpath)

        partialExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        partialExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        partialReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        partialReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);

        secondaryPartialExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryPartialExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryPartialReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryPartialReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);



        clear pctSuccesses pastResponses ruleMatrix destpath wVector  
        destpath = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
        load(destpath)

        partialExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        partialExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        partialReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        partialReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);


        secondaryPartialExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryPartialExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryPartialReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryPartialReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);


        clear pctSuccesses pastResponses ruleMatrix destpath wVector   
        destpath = sprintf('%s%sExperiment3Subject1.mat',fpath,date);
        load(destpath)

        fullExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        fullExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        fullReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        fullReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);

        secondaryFullExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryFullExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryFullReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryFullReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);

        destpath = sprintf('%s%sExperiment4Subject1.mat',fpath,date);
        load(destpath)

        fullExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        fullExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        fullReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        fullReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);

        secondaryFullExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryFullExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryFullReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryFullReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);
    end
    
    hold on
          
    line([leftBlock leftBlock],[0 1],'Color','black','LineStyle',':','LineWidth',3);
    
    H1 = plot(1:leftBlock,fullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'g-','LineWidth',3)
    H2 = plot(1:leftBlock,fullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r-','LineWidth',3)
    H3 = plot(1:leftBlock,partialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'g--','LineWidth',3)
    H4 = plot(1:leftBlock,partialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r--','LineWidth',3)

    plot(1:leftBlock,secondaryFullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b-','LineWidth',3)
    plot(1:leftBlock,secondaryFullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'y-','LineWidth',3)
    plot(1:leftBlock,secondaryPartialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b--','LineWidth',3)
    plot(1:leftBlock,secondaryPartialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'y--','LineWidth',3)

    plot(1:leftBlock,secondaryFullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b-','LineWidth',3)
    plot(1:leftBlock,secondaryFullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'y-','LineWidth',3)
    plot(1:leftBlock,secondaryPartialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b--','LineWidth',3)
    plot(1:leftBlock,secondaryPartialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'y--','LineWidth',3)

    
    plot(1:leftBlock,fullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'g-','LineWidth',3)
    plot(1:leftBlock,fullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r-','LineWidth',3)
    plot(1:leftBlock,partialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'g--','LineWidth',3)
    plot(1:leftBlock,partialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r--','LineWidth',3)
    
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO(1:rightBlock),'g-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionNOGO(1:rightBlock),'r-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO(1:rightBlock),'g--','LineWidth',3) 
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionNOGO(1:rightBlock),'r--','LineWidth',3)
    
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO(1:rightBlock),'b-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionNOGO(1:rightBlock),'y-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO(1:rightBlock),'b--','LineWidth',3) 
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO(1:rightBlock),'y--','LineWidth',3)
    
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO8(1:rightBlock),'b-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionNOGO8(1:rightBlock),'y-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO8(1:rightBlock),'b--','LineWidth',3) 
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO8(1:rightBlock),'y--','LineWidth',3)
    
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO8(1:rightBlock),'g-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionNOGO8(1:rightBlock),'r-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO8(1:rightBlock),'g--','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionNOGO8(1:rightBlock),'r--','LineWidth',3)
        
    plot(1:10:leftBlock,fullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,fullReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,partialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,partialReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    L1 = plot(1:10:leftBlock,fullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,fullReacquisitionNOGO8(1:10:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    L2 = plot(1:10:leftBlock,partialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,partialReacquisitionNOGO8(1:10:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    L3 = plot(1:10:leftBlock,fullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,fullReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    L4 = plot(1:10:leftBlock,partialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,partialReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    plot(1:10:leftBlock,fullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,fullReacquisitionGO8(1:10:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,partialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,partialReacquisitionGO8(1:10:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    plot(1:10:leftBlock,fullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,fullReacquisitionGO(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,partialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,partialReacquisitionGO(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    % Secondary
    plot(1:10:leftBlock,secondaryFullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryFullReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,secondaryPartialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    plot(1:10:leftBlock,secondaryFullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryFullReacquisitionNOGO8(1:10:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,secondaryPartialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO8(1:10:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    plot(1:10:leftBlock,secondaryFullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryFullReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,secondaryPartialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    plot(1:10:leftBlock,secondaryFullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryFullReacquisitionGO8(1:10:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,secondaryPartialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryPartialReacquisitionGO8(1:10:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    plot(1:10:leftBlock,secondaryFullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryFullReacquisitionGO(1:10:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
    plot(1:10:leftBlock,secondaryPartialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:10:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    plot(leftBlock+1:10:leftBlock+rightBlock,secondaryPartialReacquisitionGO(1:10:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    
    legend([L1 L2 L3 L4],'Ext8','Prf8','Ext2','Prf2','Location','North')
    axis([0 leftBlock+rightBlock 0 1])
    set(gca,'YMinorGrid','on')
    xlabel('Trials')
    ylabel('Synaptic weight')
    titleTxt = sprintf('GO and NOGO weight trajectories');
    title(titleTxt); 
    ptext = sprintf('%s%sWoodsBoutonWeightsComparisonFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    
        %Weights comparison
    figure()
    % First block is initial acquisition
    
    if (subjects > 1)
        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment1.mat',fpath,date);
        load(destpath)

        partialExtinctionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        partialExtinctionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        partialReacquisitionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        partialReacquisitionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryPartialExtinctionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryPartialExtinctionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryPartialReacquisitionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryPartialReacquisitionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
        
        partialExtinctionGOAlt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        partialExtinctionNOGOAlt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        partialReacquisitionGOAlt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        partialReacquisitionNOGOAlt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryPartialExtinctionGOAlt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryPartialExtinctionNOGOAlt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryPartialReacquisitionGOAlt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryPartialReacquisitionNOGOAlt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));


        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment2.mat',fpath,date);
        load(destpath)

        partialExtinctionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        partialExtinctionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        partialReacquisitionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        partialReacquisitionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryPartialExtinctionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryPartialExtinctionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryPartialReacquisitionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryPartialReacquisitionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
        
        partialExtinctionGO8Alt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        partialExtinctionNOGO8Alt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        partialReacquisitionGO8Alt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        partialReacquisitionNOGO8Alt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryPartialExtinctionGO8Alt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryPartialExtinctionNOGO8Alt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryPartialReacquisitionGO8Alt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryPartialReacquisitionNOGO8Alt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment3.mat',fpath,date);
        load(destpath)

        fullExtinctionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        fullExtinctionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        fullReacquisitionGO = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        fullReacquisitionNOGO = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryFullExtinctionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryFullExtinctionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryFullReacquisitionGO = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryFullReacquisitionNOGO = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
        
        
        fullExtinctionGOAlt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        fullExtinctionNOGOAlt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        fullReacquisitionGOAlt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        fullReacquisitionNOGOAlt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryFullExtinctionGOAlt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryFullExtinctionNOGOAlt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryFullReacquisitionGOAlt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryFullReacquisitionNOGOAlt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        clear pctSuccesses pastResponses ruleMatrix destpath sampleAvgWeightsPGO sampleAvgWeightsPNOGO sampleAvgWeightsSGO sampleAvgWeightsSNOGO   
        destpath = sprintf('%s%ssampleAvgExperiment4.mat',fpath,date);
        load(destpath)

        fullExtinctionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        fullExtinctionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        fullReacquisitionGO8 = sampleAvgWeightsPGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        fullReacquisitionNOGO8 = sampleAvgWeightsPNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryFullExtinctionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryFullExtinctionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryFullReacquisitionGO8 = sampleAvgWeightsSGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryFullReacquisitionNOGO8 = sampleAvgWeightsSNOGO((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
        
        fullExtinctionGO8Alt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        fullExtinctionNOGO8Alt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        fullReacquisitionGO8Alt = sampleAvgWeightsPGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        fullReacquisitionNOGO8Alt = sampleAvgWeightsPNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));

        secondaryFullExtinctionGO8Alt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));   
        secondaryFullExtinctionNOGO8Alt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)));
        secondaryFullReacquisitionGO8Alt = sampleAvgWeightsSGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));   
        secondaryFullReacquisitionNOGO8Alt = sampleAvgWeightsSNOGOAlt((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
    
    else
        
       clear pctSuccesses pastResponses ruleMatrix destpath wVector
        destpath = sprintf('%s%sExperiment1Subject1.mat',fpath,date);
        load(destpath)

        partialExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        partialExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        partialReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        partialReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);

        secondaryPartialExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryPartialExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryPartialReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryPartialReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);



        clear pctSuccesses pastResponses ruleMatrix destpath wVector  
        destpath = sprintf('%s%sExperiment2Subject1.mat',fpath,date);
        load(destpath)

        partialExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        partialExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        partialReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        partialReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);


        secondaryPartialExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryPartialExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryPartialReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryPartialReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);


        clear pctSuccesses pastResponses ruleMatrix destpath wVector   
        destpath = sprintf('%s%sExperiment3Subject1.mat',fpath,date);
        load(destpath)

        fullExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        fullExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        fullReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        fullReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);

        secondaryFullExtinctionGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryFullExtinctionNOGO = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryFullReacquisitionGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryFullReacquisitionNOGO = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);

        destpath = sprintf('%s%sExperiment4Subject1.mat',fpath,date);
        load(destpath)

        fullExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),1,1,1,pSweep);   
        fullExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),2,1,1,pSweep);
        fullReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),1,1,1,pSweep);   
        fullReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),2,1,1,pSweep);

        secondaryFullExtinctionGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),3,1,1,pSweep);   
        secondaryFullExtinctionNOGO8 = wVector((sum(ruleMatrix(1:1,1))+1):sum(ruleMatrix(1:2,1)),4,1,1,pSweep);
        secondaryFullReacquisitionGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),3,1,1,pSweep);   
        secondaryFullReacquisitionNOGO8 = wVector((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)),4,1,1,pSweep);
    end
    
    MFreq = 20;
    titleTxt = sprintf('GO - NOGO weight trajectories');
    title(titleTxt); 
    
    %Ext8
    subplot(221)
    title('Ext8')
    hold on
    line([leftBlock leftBlock],[-1 1],'Color','black','LineStyle',':','LineWidth',3);

    
    % Ext8 Extinction
    plot(1:leftBlock,fullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -fullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r-','LineWidth',3)
    plot(1:leftBlock,fullExtinctionGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -fullExtinctionNOGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r--','LineWidth',3)
    % Ext8 Reacquisition
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO8(1:rightBlock)...
        -fullReacquisitionNOGO8(1:rightBlock),'r-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO8Alt(1:rightBlock)...
        -fullReacquisitionNOGO8Alt(1:rightBlock),'r--','LineWidth',3)
        
%     Ext8 Secondary
     plot(1:leftBlock,secondaryFullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryFullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b-','LineWidth',3)
         plot(1:leftBlock,secondaryFullExtinctionGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryFullExtinctionNOGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b--','LineWidth',3)
    
     plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO8(1:rightBlock)...
        -secondaryFullReacquisitionNOGO8(1:rightBlock),'b-','LineWidth',3)
     plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO8Alt(1:rightBlock)...
        -secondaryFullReacquisitionNOGO8Alt(1:rightBlock),'b--','LineWidth',3)
    
    
    %Ext8 Markers
%     plot(1:MFreq:leftBlock,fullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -fullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionGO8(1:MFreq:rightBlock)...
%         -fullReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1); 
%   
%     plot(1:MFreq:leftBlock,fullExtinctionGO8Alt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -fullExtinctionNOGO8Alt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'s','MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionGO8Alt(1:MFreq:rightBlock)...
%         -fullReacquisitionNOGO8Alt(1:MFreq:rightBlock),'s','MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
%  
%     %Secondary Ext8 Markers
%     plot(1:MFreq:leftBlock,secondaryFullExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -secondaryFullExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryFullReacquisitionGO8(1:MFreq:rightBlock)...
%         -secondaryFullReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        

    axis([0 leftBlock+rightBlock -1 1])
    set(gca,'YMinorGrid','on','YTick',-1:0.2:1,'MinorGridAlpha',1,'XMinorGrid','on')
    xlabel('Trials')
    ylabel('GO-NOGO difference')
    
    %Prf8
    subplot(222)
    hold on
    line([leftBlock leftBlock],[-1 1],'Color','black','LineStyle',':','LineWidth',3);
    title('Prf8')
        

    % Prf8 Extinction
    plot(1:leftBlock,partialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
    -partialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r-','LineWidth',3)
    plot(1:leftBlock,partialExtinctionGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
    -partialExtinctionNOGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r--','LineWidth',3)

    
    %Prf8 Reacquisition 
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO8(1:rightBlock)...
        -partialReacquisitionNOGO8(1:rightBlock),'r-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO8Alt(1:rightBlock)...
        -partialReacquisitionNOGO8Alt(1:rightBlock),'r--','LineWidth',3)
    
%     Prf8 Secondary
    plot(1:leftBlock,secondaryPartialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryPartialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b-','LineWidth',3)
    plot(1:leftBlock,secondaryPartialExtinctionGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryPartialExtinctionNOGO8Alt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b--','LineWidth',3)
    
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO8(1:rightBlock)...
        -secondaryPartialReacquisitionNOGO8(1:rightBlock),'b-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO8Alt(1:rightBlock)...
        -secondaryPartialReacquisitionNOGO8Alt(1:rightBlock),'b--','LineWidth',3)
    
        %Prf8 markers
%     plot(1:MFreq:leftBlock,partialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -partialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionGO8(1:MFreq:rightBlock)...
%         -partialReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
% 
%     plot(1:MFreq:leftBlock,partialExtinctionGO8Alt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -partialExtinctionNOGO8Alt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'s','MarkerFaceColor',[0 0.3 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionGO8Alt(1:MFreq:rightBlock)...
%         -partialReacquisitionNOGO8Alt(1:MFreq:rightBlock),'s','MarkerFaceColor',[0 0.3 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
% 
%     L2 = plot(1:MFreq:leftBlock,secondaryPartialExtinctionGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -secondaryPartialExtinctionNOGO8(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryPartialReacquisitionGO8(1:MFreq:rightBlock)...
%         -secondaryPartialReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
    axis([0 leftBlock+rightBlock -1 1])
    set(gca,'YMinorGrid','on','YTick',-1:0.2:1,'MinorGridAlpha',1,'XMinorGrid','on')
    xlabel('Trials')
    ylabel('GO-NOGO difference')
    
    %Ext2
    subplot(223)
    hold on
    line([leftBlock leftBlock],[-1 1],'Color','black','LineStyle',':','LineWidth',3);
    title('Ext2')
    
    plot(1:leftBlock,fullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -fullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r-','LineWidth',3)
    plot(1:leftBlock,fullExtinctionGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -fullExtinctionNOGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r--','LineWidth',3)
    
    % Ext2 reacquisition
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO(1:rightBlock)...
        -fullReacquisitionNOGO(1:rightBlock),'r-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGOAlt(1:rightBlock)...
        -fullReacquisitionNOGOAlt(1:rightBlock),'r--','LineWidth',3)
        
    %Ext2 Secondary
    plot(1:leftBlock,secondaryFullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryFullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b-','LineWidth',3)
    plot(1:leftBlock,secondaryFullExtinctionGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryFullExtinctionNOGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b--','LineWidth',3)
    
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO(1:rightBlock)...
        -secondaryFullReacquisitionNOGO(1:rightBlock),'b-','LineWidth',3)
    
        plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGOAlt(1:rightBlock)...
        -secondaryFullReacquisitionNOGOAlt(1:rightBlock),'b--','LineWidth',3)
    
        %Ext2 markers
%     plot(1:MFreq:leftBlock,fullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -fullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionGO(1:MFreq:rightBlock)...
%         -fullReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);  
%     
%     plot(1:MFreq:leftBlock,fullExtinctionGOAlt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -fullExtinctionNOGOAlt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'o','MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionGOAlt(1:MFreq:rightBlock)...
%         -fullReacquisitionNOGOAlt(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor',[0 0.8 0],'MarkerSize',15,'LineWidth',1); 
%     % Secondary Ext2
%     L3 = plot(1:MFreq:leftBlock,secondaryFullExtinctionGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -secondaryFullExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryFullReacquisitionGO(1:MFreq:rightBlock)...
%         -secondaryFullReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        

        axis([0 leftBlock+rightBlock -1 1])
    set(gca,'YMinorGrid','on','YTick',-1:0.2:1,'MinorGridAlpha',1,'XMinorGrid','on')
    xlabel('Trials')
    ylabel('GO-NOGO difference')
    
    %Prf2
    subplot(224)
    hold on
    line([leftBlock leftBlock],[-1 1],'Color','black','LineStyle',':','LineWidth',3);      
    title('Prf2')
    
    plot(1:leftBlock,partialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -partialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r-','LineWidth',3)
        plot(1:leftBlock,partialExtinctionGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -partialExtinctionNOGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'r--','LineWidth',3)
    
    plot(1:leftBlock,secondaryPartialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryPartialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b-','LineWidth',3)
    plot(1:leftBlock,secondaryPartialExtinctionGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1))...
        -secondaryPartialExtinctionNOGOAlt(ruleMatrix(2,1)-leftBlock+1:ruleMatrix(2,1)),'b--','LineWidth',3)
    
    % Prf2 reacquisition
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO(1:rightBlock)...
        -partialReacquisitionNOGO(1:rightBlock),'r-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGOAlt(1:rightBlock)...
        -partialReacquisitionNOGOAlt(1:rightBlock),'r--','LineWidth',3)

    %Prf2 Secondary
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO(1:rightBlock)...
        -secondaryPartialReacquisitionNOGO(1:rightBlock),'b-','LineWidth',3)
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGOAlt(1:rightBlock)...
        -secondaryPartialReacquisitionNOGOAlt(1:rightBlock),'b--','LineWidth',3)
    
    %Prf2 markers
%     plot(1:MFreq:leftBlock,partialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -partialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);    
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionGO(1:MFreq:rightBlock)...
%         -partialReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%   
%     plot(1:MFreq:leftBlock,partialExtinctionGOAlt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -partialExtinctionNOGOAlt(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'o','MarkerFaceColor',[0 0.2 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionGOAlt(1:MFreq:rightBlock)...
%         -partialReacquisitionNOGOAlt(1:MFreq:rightBlock),'o','MarkerFaceColor',[0 0.2 0],'MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);        
%  
%      %Secondary Prf2
%     L4 = plot(1:MFreq:leftBlock,secondaryPartialExtinctionGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1))...
%         -secondaryPartialExtinctionNOGO(ruleMatrix(2,1)-leftBlock+1:MFreq:ruleMatrix(2,1)),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
%     plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryPartialReacquisitionGO(1:MFreq:rightBlock)...
%         -secondaryPartialReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    axis([0 leftBlock+rightBlock -1 1])
    set(gca,'YMinorGrid','on','YTick',-1:0.2:1,'MinorGridAlpha',1,'XMinorGrid','on')
    xlabel('Trials')
    ylabel('GO-NOGO difference')

    ptext = sprintf('%s%sWoodsBoutonWeightsSpreadComparisonFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
   
%% Reacquisition weight trajectories for Woods/Bouton sims
    figure()

    subplot(221) %Primary GO
    hold on
    xlabel('Trials')
    ylabel('Weight')
    title('Primary GO weights')
    axis([leftBlock+1,leftBlock+rightBlock,0.65 0.85])  
    set(gca,'YMinorGrid','on','YTick',0:0.1:1,'MinorGridAlpha',1,'XMinorGrid','on')
    
    %Prf2
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionGO(1:MFreq:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Prf8
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionGO8(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Ext2
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionGO(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    
    %Ext8
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionGO8(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    
    subplot(222) %Primary NOGO
    hold on
    xlabel('Trials')
    ylabel('Weight')
    title('Primary NOGO weights')
    axis([leftBlock+1,leftBlock+rightBlock,0.3 0.9])  
    set(gca,'YMinorGrid','on','YTick',0:0.1:1,'MinorGridAlpha',1,'XMinorGrid','on')
    
    %Prf2
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionNOGO(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Prf8
    plot(leftBlock+1:leftBlock+rightBlock,partialReacquisitionNOGO8(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,partialReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Ext2
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionNOGO(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Ext8
    plot(leftBlock+1:leftBlock+rightBlock,fullReacquisitionNOGO8(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,fullReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    
    subplot(223) %Secondary GO
   hold on
    xlabel('Trials')
    ylabel('Weight')
    title('Secondary GO weights')
    axis([leftBlock+1,leftBlock+rightBlock,0.4 0.5])  
    set(gca,'YMinorGrid','on','YTick',0:0.1:1,'MinorGridAlpha',1,'XMinorGrid','on')

      %Prf2
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO(1:rightBlock),'k-','LineWidth',3)
     plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryPartialReacquisitionGO(1:MFreq:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
   
    %Prf8
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionGO8(1:rightBlock),'k-','LineWidth',3)       
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryPartialReacquisitionGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Ext2
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryFullReacquisitionGO(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
   
    %Ext8
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionGO8(1:rightBlock),'k-','LineWidth',3)
    plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryFullReacquisitionGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

      
    subplot(224) %Secondary NOGO
    hold on
    xlabel('Trials')
    ylabel('Weight')
    title('Secondary NOGO weights')
    axis([leftBlock+1,leftBlock+rightBlock,0 0.35])  
    set(gca,'YMinorGrid','on','YTick',0:0.1:1,'MinorGridAlpha',1,'XMinorGrid','on')
      
    %Prf2
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO(1:rightBlock),'k-','LineWidth',3)
    L4=plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Prf8
    plot(leftBlock+1:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO8(1:rightBlock),'k-','LineWidth',3)      
    L2=plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryPartialReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    %Ext2
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionNOGO(1:rightBlock),'k-','LineWidth',3)
    L3=  plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryFullReacquisitionNOGO(1:MFreq:rightBlock),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);
  
    %Ext8
    plot(leftBlock+1:leftBlock+rightBlock,secondaryFullReacquisitionNOGO8(1:rightBlock),'k-','LineWidth',3)
    L1 = plot(leftBlock+1:MFreq:leftBlock+rightBlock,secondaryFullReacquisitionNOGO8(1:MFreq:rightBlock),'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',15,'LineWidth',1);

    
    ax=axes('Units','Normal','Position',[.025 .025 .95 .95],'Visible','off');
    set(get(ax,'Title'),'Visible','on','Units','normalized','Position',[0.5 0.985])
    title('Reacquisition weight trajectories')
    lh = legend([L1 L2 L3 L4],'Ext8','Prf8','Ext2','Prf2');
    set(lh,'Position',[0 0.5 0.1 0.1])
    ptext = sprintf('%s%sWoodsBoutonReacqWeightsComparisonFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
end

if (savingsFig == 1)
    
    acqIndices = 1:SVtrainingBlock;
    extIndices = SVtrainingBlock+1:SVtrainingBlock+SVextBlock;
    reacqIndices = SVtrainingBlock + SVextBlock+1:SVtrainingBlock + SVextBlock+SVreacBlock;
    
    if (subjects == 1)
        acquisition = 100*pastResponses(1,acqIndices);
        extinction = 100*pastResponses(1,extIndices);
        reacquisition = 100*pastResponses(1,reacqIndices);
        
        acqPGOweights = wVector(acqIndices,1,1,1,pSweep);   
        acqPNOGOweights = wVector(acqIndices,2,1,1,pSweep);
        reacqPGOweights = wVector(reacqIndices,1,1,1,pSweep);
        reacqPNOGOweights =  wVector(reacqIndices,2,1,1,pSweep);
        acqSGOweights = wVector(acqIndices,3,1,1,pSweep);   
        acqSNOGOweights = wVector(acqIndices,4,1,1,pSweep);
        reacqSGOweights = wVector(reacqIndices,3,1,1,pSweep);
        reacqSNOGOweights =  wVector(reacqIndices,4,1,1,pSweep);
    else
        acquisition = 100*sampleAvgResponses(1,acqIndices);
        extinction = 100*sampleAvgResponses(1,extIndices);
        reacquisition = 100*sampleAvgResponses(1,reacqIndices);
        
        acqPGOweights = sampleAvgWeightsPGO(acqIndices);
        acqPNOGOweights = sampleAvgWeightsPNOGO(acqIndices);
        reacqPGOweights = sampleAvgWeightsPGO(reacqIndices);
        reacqPNOGOweights = sampleAvgWeightsPNOGO(reacqIndices);
        acqSGOweights = sampleAvgWeightsSGO(acqIndices);
        acqSNOGOweights =  sampleAvgWeightsSNOGO(acqIndices);
        reacqSGOweights =  sampleAvgWeightsSGO(reacqIndices);
        reacqSNOGOweights =   sampleAvgWeightsSNOGO(reacqIndices);
    end
    
    figure()
    hold on
    h1 = plot(acqIndices,acquisition,'k','LineWidth',3);
    h2 = plot(acqIndices,reacquisition,'r','LineWidth',3);
    
    legend([h1 h2],'Acquisition','Reacquisition','Location','SouthEast')
    axis([1 SVtrainingBlock 0 100])
    set(gca,'YMinorGrid','on','LineWidth',2,'MinorGridLineStyle','-','MinorGridColor',[.6 .6 .6])
    xlabel('Trials')
    ylabel('% CRs')
    titleTxt = sprintf('Exploitation of reinstated contingency is faster following extinction');
    title(titleTxt); 
    ptext = sprintf('%s%sSavingsFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    
    figure()
    hold on
    h1 = plot(acqIndices,acquisition,'k','LineWidth',3);
    h2 = plot(extIndices,extinction,'k:','LineWidth',3);
    h3 = plot(reacqIndices,reacquisition,'r','LineWidth',3);
    
    legend([h1 h2 h3],'Acquisition','Extinction','Reacquisition','Location','SouthEast')
    axis([1 SVtrainingBlock+SVextBlock+SVreacBlock 0 100])
    set(gca,'YMinorGrid','on','LineWidth',2,'MinorGridLineStyle','-','MinorGridColor',[.6 .6 .6])
    xlabel('Trials')
    ylabel('% CRs')
    titleTxt = sprintf('Response production tracks reward contingency');
    title(titleTxt); 
    ptext = sprintf('%s%sSavingsFullFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
    
    f1 = figure()
    
    axis1 = axes();
    
    hold on
    h1 = plot(acqIndices,acqPGOweights,'g','LineWidth',3);
    h2 = plot(acqIndices,acqPNOGOweights,'r','LineWidth',3);
    h3 = plot(acqIndices,reacqPGOweights,'g:','LineWidth',3);
    h4 = plot(acqIndices,reacqPNOGOweights,'r:','LineWidth',3);


    axis([1 SVtrainingBlock 0 1])

%     box on
    axis2 = copyobj(axis1,f1);
    
    ylabel('Synaptic weights')
    xlabel('Trials')
    titleTxt = sprintf('Savings in GO/NOGO weight difference accelerates reacquisition');
    title(titleTxt); 
    set(axis1,'YMinorGrid','on','LineWidth',2,'MinorGridLineStyle','-','MinorGridColor',[.6 .6 .6]) 
    
    legendY = 0.775;
    legendX = 0.2;
    legendXwidth = 0.3;
    legendYwidth = 0.1;
    
    l1 = legend([h1 h2 h3 h4],{'Acquisition Primary GO','Acquisition Primary NOGO',...
         'Reacquisition Primary GO','Reacquisition Primary NOGO'},...
         'Position',[legendX legendY legendXwidth legendYwidth],'Color','w','EdgeColor','none');
    
   
    delete(get(axis2,'Children'));
    set(axis2,'Color', 'none');
    set(axis2,'YAxisLocation','Right');

    
    h5 = plot(acqIndices,acqSGOweights,'b','LineWidth',3);
    h6 = plot(acqIndices,acqSNOGOweights,'k','LineWidth',3);
    h7 = plot(acqIndices,reacqSGOweights,'b:','LineWidth',3);
    h8 = plot(acqIndices,reacqSNOGOweights,'k:','LineWidth',3);
    l2 = legend(axis2,[h5 h6 h7 h8],{'Acquisition Secondary GO','Acquisition Secondary NOGO',...
         'Reacquisition Secondary GO','Reacquisition Secondary NOGO'},...
         'Position',[legendX+legendXwidth legendY legendXwidth legendYwidth],'Color','w','EdgeColor','none');
     
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
   
    ptext = sprintf('%s%sSavingsWeightFig.png',fpath,date);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
end

if (plotTimes == 1)
    
    colorSpan = 0:0.1:1;
    blacks = [zeros(size(colorSpan)); zeros(size(colorSpan)); colorSpan]';
    colors = [colorSpan; zeros(size(colorSpan)); 1-colorSpan]';
    cmap = [blacks; colors];
           
    
    normalizedTimes = avgTimes / max(avgTimes);
    
    figure()
    hold on
    imageH = image(ones(T,T,3));
    set(imageH,'xdata',[1 T],'ydata',[min(normalizedTimes) max(normalizedTimes)]);
    set(gca,'xlim',[1 T],'ylim',[min(normalizedTimes) max(normalizedTimes)]);
    setBackground(imageH,ruleMatrix);
    hold on

    surface([1:length(normalizedTimes);1:length(normalizedTimes)],...
        [normalizedTimes;normalizedTimes],...
        [zeros(size(normalizedTimes));zeros(size(normalizedTimes))],...
        [pctPrimaryDecision;pctPrimaryDecision],...
        'facecol','no',...
        'edgecol','interp',...
        'LineWidth',4);
    colormap(colors);
    caxis([0 1]);
    h= colorbar('Ticks',[0.1, 0.9],...
        'TickLabels',{'Secondary','Primary'});
    set(h,'YDir','reverse');
    title('Per-trial time elapsed to decision')
    xlabel('Trial number')
    ylabel('Time elapsed')
    [rmL , ~] = size(ruleMatrix);
    p = zeros(rmL,1);
    t = 1;
    for j = 1:rmL
        p(j) = ruleMatrix(j,1);
        tAdj = t+0.1*p(j);
        txt = descTextByRule(ruleMatrix(j,2),chanceMatrix(j,2),chanceMatrix(j,1));
        text('units','data','position',[tAdj 1+max(tVector)],'FontSize',14,'FontWeight','bold','string',txt);
        t=t+p(j);
    end

    % Print to file
    ptext = sprintf('%s%sTimeElapsedFig_param%d.png',fpath,date,pSweep);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    print(gcf,'-dpng',ptext,'-r0');
    hold off
    close
end


%% SAT figure
if(plotDraw(3) == 1)
primaryTrials = (rVector(:,4) == 1);
secondaryTrials = (rVector(:,4) == 2);

primaryTimes = tVector(primaryTrials);
secondaryTimes = tVector(secondaryTrials);

primarySpeeds = 1./primaryTimes;
secondarySpeeds = 1./secondaryTimes;

primarySuccesses = rVector(primaryTrials,1).*(rVector(primaryTrials,1) > 0);
secondarySuccesses = rVector(secondaryTrials,1).*(rVector(secondaryTrials,1) > 0);

window = 49;

primaryAccuracy = zeros(length(primarySuccesses),1);
for i = 1:length(primarySuccesses)
    if (i<=window)
        primaryAccuracy(i) = mean(primarySuccesses(1:i));
    else
        primaryAccuracy(i) = mean(primarySuccesses(i-window:i));
    end
end

secondaryAccuracy = zeros(length(secondarySuccesses),1);
for i = 1:length(secondarySuccesses)
     if (i<=window)
        secondaryAccuracy(i) = mean(secondarySuccesses(1:i));
    else
        secondaryAccuracy(i) = mean(secondarySuccesses(i-window:i));
    end
end

[sortedPrimarySpeeds sortIndices] = sort(primarySpeeds);
[uniquePrimarySpeeds IA IC] = unique(sortedPrimarySpeeds); 
primaryAccuraciesSortedBySpeed = primaryAccuracy(sortIndices);
primaryAccuracySamples = zeros(1,length(uniquePrimarySpeeds));
meanPrimaryAccuracyBySpeed = zeros(1,length(uniquePrimarySpeeds));
for i = 1:length(uniquePrimarySpeeds)
    if (i < length(uniquePrimarySpeeds))
    meanPrimaryAccuracyBySpeed(i) = mean(primaryAccuraciesSortedBySpeed(IA(i):IA(i+1)-1));
    stdPrimaryAccuracyBySpeed(i) = mean(primaryAccuraciesSortedBySpeed(IA(i):IA(i+1)-1));
    primaryAccuracySamples(i) = length(IA(i):IA(i+1)-1);
    else
    meanPrimaryAccuracyBySpeed(i) = mean(primaryAccuraciesSortedBySpeed(IA(i):end));
    stdPrimaryAccuracyBySpeed(i) = mean(primaryAccuraciesSortedBySpeed(IA(i):end));
  
    primaryAccuracySamples(i) = length(IA(i):IA(end));
    end
end

[sortedSecondarySpeeds sortSecondaryIndices] = sort(secondarySpeeds);
[uniqueSecondarySpeeds IA IC] = unique(sortedSecondarySpeeds); 
secondaryAccuraciesSortedBySpeed = secondaryAccuracy(sortSecondaryIndices);
secondaryAccuracySamples = zeros(1,length(uniqueSecondarySpeeds));
meanSecondaryAccuracyBySpeed = zeros(1,length(uniqueSecondarySpeeds));
for i = 1:length(uniqueSecondarySpeeds)
    if(i < length(uniqueSecondarySpeeds))
        meanSecondaryAccuracyBySpeed(i) = mean(secondaryAccuraciesSortedBySpeed(IA(i):IA(i+1)-1));
        stdSecondaryAccuracyBySpeed(i) = std(secondaryAccuraciesSortedBySpeed(IA(i):IA(i+1)-1));
        secondaryAccuracySamples(i) = length(IA(i):IA(i+1)-1);
    else
        meanSecondaryAccuracyBySpeed(i) = mean(secondaryAccuraciesSortedBySpeed(IA(i):end));
        stdSecondaryAccuracyBySpeed(i) = std(secondaryAccuraciesSortedBySpeed(IA(i):end));
        secondaryAccuracySamples(i) = length(IA(i):IA(end));
    end
end



figure()
h = plot(1./uniquePrimarySpeeds,meanPrimaryAccuracyBySpeed,'r',1./uniqueSecondarySpeeds,meanSecondaryAccuracyBySpeed,'b');
xlabel('Decision Time')
ylabel('Accuracy')
title('Speed-Accuracy Tradeoff')
legend('Primary','Secondary','Location','SouthOutside');
% axis([15 50 0 1]);
set(h,'LineWidth',3)
%  axis([1 50 0 1])


figure()
plot(primarySpeeds,primaryAccuracy,'rx',secondarySpeeds,secondaryAccuracy,'bx')

aggregateTimes = [primaryTimes; secondaryTimes];
aggregateAccuracy = [primaryAccuracy; secondaryAccuracy];

fastTimes = aggregateTimes(aggregateTimes < median(aggregateTimes));
slowTimes = aggregateTimes(aggregateTimes >= median(aggregateTimes));

fastAccuracy = aggregateAccuracy(aggregateTimes < median(aggregateTimes));
slowAccuracy = aggregateAccuracy(aggregateTimes >= median(aggregateTimes));

figure();
bar([mean(slowTimes) mean(fastTimes)],[mean(slowAccuracy) mean(fastAccuracy)]);


[sortedSlowTimes sortSlowTimeIndices] = sort(slowTimes);
[uniqueSlowTimes IA IC] = unique(sortedSlowTimes); 
slowAccuraciesSortedByTime = slowAccuracy(sortSlowTimeIndices);
for i = 1:length(uniqueSlowTimes)
    if(i < length(uniqueSlowTimes))
        meanSlowAccuracyByTime(i) = mean(slowAccuraciesSortedByTime(IA(i):IA(i+1)-1));
        steSlowAccuracyByTime(i) = std(slowAccuraciesSortedByTime(IA(i):IA(i+1)-1))/sqrt(length(IA(i):IA(i+1)-1));
    else
        meanSlowAccuracyByTime(i) = mean(slowAccuraciesSortedByTime(IA(i):end));
        steSlowAccuracyByTime(i) = std(slowAccuraciesSortedByTime(IA(i):IA(end)))/sqrt(length(IA(i):length(sortedSlowTimes)));
    end
end


[sortedFastTimes sortFastTimeIndices] = sort(fastTimes);
[uniqueFastTimes IA IC] = unique(sortedFastTimes); 
fastAccuraciesSortedByTime = slowAccuracy(sortFastTimeIndices);
for i = 1:length(uniqueFastTimes)
    if(i < length(uniqueFastTimes))
        meanFastAccuracyByTime(i) = mean(fastAccuraciesSortedByTime(IA(i):IA(i+1)-1));
        steFastAccuracyByTime(i) = std(fastAccuraciesSortedByTime(IA(i):IA(i+1)-1))/sqrt(length(IA(i):IA(i+1)-1));
    else
        meanFastAccuracyByTime(i) = mean(fastAccuraciesSortedByTime(IA(i):end));
        steFastAccuracyByTime(i) = std(fastAccuraciesSortedByTime(IA(i):end))/sqrt(length(IA(i):length(sortedFastTimes)));
    end
end

figure()
bar([uniqueSlowTimes; uniqueFastTimes],[meanSlowAccuracyByTime meanFastAccuracyByTime]);

[sortedUniqueTimes I] = sort([uniqueSlowTimes; uniqueFastTimes]);
meanAccuraciesByTimes = [meanSlowAccuracyByTime meanFastAccuracyByTime];
steAccuraciesByTimes = [steSlowAccuracyByTime steFastAccuracyByTime];

figure()
hold on
h1= plot(sortedUniqueTimes, meanAccuraciesByTimes(I),'k',...
         sortedUniqueTimes,meanAccuraciesByTimes(I) + 1.96*steAccuraciesByTimes(I),'k:',...
         sortedUniqueTimes,meanAccuraciesByTimes(I) - 1.96*steAccuraciesByTimes(I),'k:',...
         'LineWidth',3);
h2 = plot(1./uniquePrimarySpeeds,meanPrimaryAccuracyBySpeed,'r--',1./uniqueSecondarySpeeds,meanSecondaryAccuracyBySpeed,'b--','LineWidth',3);
title('Speed-Accuracy Tradeoff')
legend([h2(2), h1(1), h2(1)],'Secondary-controlled','Overall performance','Primary-controlled','Location','SouthOutside');
% plot(1./primarySpeeds,primaryAccuracy,'rx',1./secondarySpeeds,secondaryAccuracy,'bx')
xlabel('Decision time')
ylabel('Accuracy')


ptext = sprintf('%s%sSATCurveFig.png',fpath,date);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
print(gcf,'-dpng',ptext,'-r0');
hold off
close




end 

if(plotDraw(4) == 1)

testBlockA = ((sum(ruleMatrix(1:2,1))+1):sum(ruleMatrix(1:3,1)));
testBlockB = ((sum(ruleMatrix(1:5,1))+1):sum(ruleMatrix(1:6,1)));
testBlockC = ((sum(ruleMatrix(1:8,1))+1):sum(ruleMatrix(1:9,1)));

hold on
h1 = plot(0:length(testBlockA)-1,100*pastResponses(1,testBlockA),'b-');
h2 = plot(0:length(testBlockB)-1,100*pastResponses(1,testBlockB),'--');
set(h2,'Color',[0.5 0 0.5])
h3 = plot(0:length(testBlockC)-1,100*pastResponses(1,testBlockC),'r:');
set(h1,'LineWidth',3);
set(h2,'LineWidth',3);
set(h3,'LineWidth',3);
xlabel('Trials after rule change')
ylabel('% incorrect response')
legend([h1 h2 h3],'Brief training','Extended training','Overtraining','Location','SouthOutside');
   
% Print to file
ptext = sprintf('%s%sPlusMazeFigparam%d.png',fpath,date,pSweep);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
print(gcf,'-dpng',ptext,'-r0');
hold off
close
end 

end

function [clippedVector] = EveryNthElement(fullVector,N,Dim)
     %Filter down to every N'th element
    clippedVector = [];
    for mIndex = 1:N:length(fullVector)
         clippedVector = cat(Dim,clippedVector,fullVector(mIndex));
    end
end

function [markerString] = markerStringByResponse(response)
    switch response
        case 1
            markerString = 'o';
        case 2
            markerString = 's';
        case 3
            markerString = 'd';
        otherwise
            markerString = '';
    end

end

function [imageH] = setBackground(imageH,ruleMatrix)
backGcolors = [1 1 1;0.7 0.7 0.7;...
1 1 1;0.7 0.7 0.7;...
1 1 1;0.7 0.7 0.7;...
1 1 1;0.7 0.7 0.7;...
1 1 1;0.7 0.7 0.7;...
1 1 1;0.7 0.7 0.7;...
               ];
        
        backg = get(imageH,'cdata');
        [rmL, ~] = size(ruleMatrix);
        p = zeros(rmL,1);
        t = 1;
        for j = 1:rmL
           p(j) = ruleMatrix(j,1);
           backg(:,t:t+p(j),1) = backGcolors(j,1);
           backg(:,t:t+p(j),2) = backGcolors(j,2);
           backg(:,t:t+p(j),3) = backGcolors(j,3);
           t = t + p(j);
        end
        set(imageH,'cdata',backg);
end

%% Give text contingency descriptor as a function of rule - unused
function [descText] = descTextByRule(rule,param1,param2)
     descText = 'Function deprecated';
end