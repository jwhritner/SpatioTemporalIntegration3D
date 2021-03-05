%% Analysis code for Spatiotemporal Integration 3D
% jake.whritner@utexas.edu 

% This code will load in cell arrays with all of the individual data, as
% well as aggregate data structures that include fitted values for the
% aggregate. It will then plot the results figures for the temporal and
% spatial experiments presented in the paper.

clear; close all;

%%%% Params %%%% 
runBootstrapAnalysis = 0; % by default, we'll just load the previous dataStruct 
% Can flip it on to run new simulations or just poke around the code. Keep
% in mind that it takes a few minutes due to the large number of
% repetitions/large matrices 

%%%% Loading the data %%%%
%% Load it
load('data/allSubjsT.mat') % load temporal data, 8 cells: 4 subj x 2 conditions (cd odds, iovd evens) 
% temporal info:
% 3200 trials per condition per subject
% manipulated Duration (allSubjT{i}.Duration) - 16 durs, note they differ
% for cd and iovd 
load('data/aggDataT.mat'); % this struct has some info from the aggregate data, 
% such as percent correct for the concatenated subject data, fitted values 
load('data/allSubjsS.mat') % load spatial data, 6 cells: 3 subj (1-3 from temporal) x 2 conditions (cd odds, iovd evens) 
% spatial info: 
% 1600 trials per condition per subject
% manipulated Alpha (allSubjS{i}.Alpha) - 4 alphas (22.5, 45, 90, 180) 
% might not need these. just use the means of the bootstrapped fits? 
load('data/aggDataS.mat'); % this struct has some info from the aggregate data, 
% such as percent correct for the concatenated subject data, fitted values


%%%% Getting responses, computing percent correct %%%%
%% Temporal 
% each subject ran 200 trials per duration, per condition
% sort data into 200 x 16 x 4 (nTrials x nDurs x nSubjs) for each iovd and cd
% bootstrap will do 10 k reps of 200 trials x 16 durs x 4 subjs
%allSubjsT{[1 3 5 7]}  %(1 = jawCD, 3 = achCD, 5 = chgCD, 7 = lkcCD)
%allSubjsT{[2 4 6 8]}  %(2 = jawIOVD, 4 = achIOVD, 6 = chgIOVD, 8 = lkcIOVD)
for ii = 1:size(allSubjsT,2)
    % getting performance for individual subjects
    Dur(:,ii) = unique(allSubjsT{ii}.Duration);
    trialsPerDur(:,ii) = hist(allSubjsT{ii}.Duration,unique(Dur(:,ii)))';    
    for kk = 1:length(Dur)
        currentDur = round(Dur(kk,ii),3);
        whichTrials = find(round(allSubjsT{ii}.Duration,3) == currentDur);
        responsesT(:,kk,ii) = allSubjsT{ii}.Response(whichTrials) ==  allSubjsT{ii}.Direction(whichTrials);       
    end
end

% compute pct corrects per duration here (will use for plotting) 
pctCrctPerDur = squeeze(sum(responsesT,1)/trialsPerDur(1)'.*100);
pctCrctPerDur_cd = pctCrctPerDur(:,[1 3 5 7]); % just nice to seperate for ease of use later 
pctCrctPerDur_iovd = pctCrctPerDur(:,[2 4 6 8]);

%% Spatial

for ii = 1:size(allSubjsS,2)
    % getting performance for individual subjects
    Alpha = unique(allSubjsS{ii}.Alpha);
    trialsPerAlpha(:,ii) = hist(allSubjsS{ii}.Alpha,unique(Alpha(:,1)))';
    for kk = 1:length(Alpha)
        currentAlpha = round(Alpha(kk,1),2);
        whichTrials = find(round(allSubjsS{ii}.Alpha,3) == currentAlpha);
        responsesS(:,kk,ii) = allSubjsS{ii}.Response(whichTrials) ==  allSubjsS{ii}.Direction(whichTrials);
    end
end
% compute pct corrects per alpha here (will use for plotting) 
pctCrctPerA = squeeze(sum(responsesS,1)/trialsPerAlpha(1)'.*100);
pctCrctPerA_cd = pctCrctPerA(:,[1 3 5]); % just nice to seperate for ease of use later 
pctCrctPerA_iovd = pctCrctPerA(:,[2 4 6]);
% responses -- 400 (trials per alpha) x 4 alphas x 6 subjs (3 cd (odds) and 3
% iovd (evens))


%% Bootsrapping 
% here, either load in a struct that has saved BS info or run the function
% to get that info

% now I have 200 (nTrialsPerDur) x 16 (nDurs) x 8 (nSubjs--odds are cd, evens iovd) respones
% keep in mind that odd and even columns will have different Durs (match
% up with Dur later)
if runBootstrapAnalysis == 1
    bsDataStruct = analysis_bootstrapCIs(responsesT,responsesS,Dur,Alpha)
else 
    load('data/bsDataStruct.mat')
end 

%% Fitting saturating exponential functions to the individual data -- put this in a seprate function 
indDataStacked = reshape(responsesT,3200,8); % want responses in one long vector for feeding into the fitting algorithm (I think)
nTrialsPerDur = (size(indDataStacked,1)/length(Dur));
repDurs_cd = sort(repmat(Dur(:,1),nTrialsPerDur,1));
repDurs_iovd = sort(repmat(Dur(:,2),nTrialsPerDur,1));
durs_cd = unique(Dur(:,1));
durs_iovd = unique(Dur(:,2));
indDataStacked_cd = indDataStacked(:,[1 3 5 7]);
indDataStacked_iovd = indDataStacked(:,[2 4 6 8]);

% fitting for individual data, disparity-based first 
f_cd = @(b,durs_cd) b(1).*exp(-1*durs_cd./b(2))+b(3);
options = optimset('MaxFunEvals',1e5);
% fit our function 
for ii = 1:size(indDataStacked_cd,2)
    fitParams_cd(:,ii) = fminsearch(@(b) norm(indDataStacked_cd(:,ii) - f_cd(b,repDurs_cd)), [0; 0; 0], options);
    indFitVals_cd(:,ii) = (f_cd(fitParams_cd(:,ii),durs_cd))*100;
    corrCoeff_cd(ii) = corr2(indFitVals_cd(:,ii), pctCrctPerDur_cd(:,ii));
    
end
tau_cd = fitParams_cd(2,:);
expRsquare_cd = power(corrCoeff_cd,2);

% now velocity-based
f_iovd = @(b,durs_iovd) b(1).*exp(-1*durs_iovd./b(2))+b(3);
options = optimset('MaxFunEvals',1e5);
% fit our function 
for ii = 1:size(indDataStacked_iovd,2)
    fitParams_iovd(:,ii) = fminsearch(@(b) norm(indDataStacked_iovd(:,ii) - f_cd(b,repDurs_iovd)), [0; 0; 0], options);
    indFitVals_iovd(:,ii) = (f_iovd(fitParams_iovd(:,ii),durs_iovd))*100;
    corrCoeff_iovd(ii) = corr2(indFitVals_iovd(:,ii), pctCrctPerDur_iovd(:,ii));
    
end
tau_iovd = fitParams_iovd(2,:);
expRsquare_iovd = power(corrCoeff_iovd,2);

%% Plotting individual subjects w/ bootstrapped error bars

% grab the SEMS from our bsDataStruct for our CI error bars  
bsSEMt_cd = bsDataStruct.semT (:,[1 3 5 7]); 
bsSEMt_iovd = bsDataStruct.semT (:,[2 4 6 8]); 

figure()
PaperPositionMode = 'manual';
orient('landscape')
for ii = 1:size(pctCrctPerDur,2)/2
    subplot(size(pctCrctPerDur,2)/4,size(pctCrctPerDur,2)/4,ii)
%     errorbar(durs_cd.*1000, pctCrctPerDur_cd(:,ii),bsSEMt_cd(:,ii),'.','MarkerSize',30,'CapSize',0,'Color',[0 0.4470 0.7410]); % 68% conf
    errorbar(durs_cd.*1000, pctCrctPerDur_cd(:,ii),bsSEMt_cd(:,ii)*1.96,'.','MarkerSize',30,'CapSize',0,'Color',[0 0.4470 0.7410]); % *1.96 for 95% conf
    hold on
    plot(durs_cd.*1000, indFitVals_cd(:,ii),'Color',[0 0.4470 0.7410],'LineWidth',2) % adding the fitted line from the bootstrapped
%     errorbar(durs_iovd.*1000, pctCrctPerDur_iovd(:,ii),bsSEMt_iovd(:,ii),'.','MarkerSize',30,'CapSize',0,'Color',[0.9100 0.4100 0.1700]); % just bsSEM if 68% (+-1 std)
    errorbar(durs_iovd.*1000, pctCrctPerDur_iovd(:,ii),bsSEMt_iovd(:,ii)*1.96,'.','MarkerSize',30,'CapSize',0,'Color',[0.9100 0.4100 0.1700]); % just bsSEM if 68% (+-1 std) if 95, multiply by 1.96
    plot(durs_iovd.*1000, indFitVals_iovd(:,ii),'Color',[0.9100 0.4100 0.1700],'LineWidth',2) % adding the fitted line from the bootstrapped
    % fix iovd color here.
    grid off
    xlim([0 1300])
    xticks([0:250:1200])
    ylim([40 110])
    yticks([50:10:100])
    formatFigure('Stim duration (ms)', 'Percent correct', "Subj"+ii, 18, 14);
    set(gcf,'units','points','position',[50,50,800,600])
    box off
end

%% Aggregate plot temporal 

% fits for the aggregate data are in the aggDataT struct
% can see the process of bootstrapping the fits (to get CIs) in the
% analysis_bootstrapCIs function 

% pull bootstrapped SEMs from the bsDataStruct 
bsAggSEMt_cd = bsDataStruct.aggSEMt(:,1);
bsAggSEMt_iovd = bsDataStruct.aggSEMt(:,2);

% plot our psychometric w/ error bars and fits (main temporal results figure #) 
f = figure();
f.PaperPositionMode = 'manual';
orient(f,'landscape')
% c = errorbar(durs_cd.*1000, aggDataT.CD.pctCrctPerDur*100,bsAggSEMt_cd,'.','MarkerSize',30,'Capsize',0); % plot w/ BS error bars (68% conf intervals)
c = errorbar(durs_cd.*1000, aggDataT.CD.pctCrctPerDur*100,bsAggSEMt_cd*1.96,'.','MarkerSize',30,'Capsize',0); % plot w/ BS error bars (95% conf intervals)
c.Color = [0 0.4470 0.7410];
c.CapSize = 0;
hold on
plot(durs_cd.*1000, aggDataT.CD.fittedValues*100,'Color',[0 0.4470 0.7410],'LineWidth',2)
% i = errorbar(durs_iovd.*1000.*1000, aggDataT.IOVD.pctCrctPerDur*100,bsAggSEMt_iovd,'.','MarkerSize',30,'Capsize',0); % plot with errorbars
i = errorbar(durs_iovd.*1000, aggDataT.IOVD.pctCrctPerDur*100,bsAggSEMt_iovd*1.96,'.','MarkerSize',30,'Capsize',0); % plot with errorbars
i.Color = [0.9100 0.4100 0.1700];
i.CapSize = 0;
plot(durs_iovd.*1000, aggDataT.IOVD.fittedValues*100, 'Color',[0.9100 0.4100 0.1700],'LineWidth',2)
grid off
xlim([0 1200])
ylim([45 100])
ax = gca;
ax.FontSize = 12;
legend('Disparity-based',' ','Velocity-based',' ')
formatFigure('Stim duration (ms)', 'Percent correct', 'Aggregate data (n = 4), Trials per condition = 12,800', 18, 14);
set(gcf,'units','points','position',[50,50,800,600])
box off

%% Spatial -  plotting individual subjects w/ bootstrapped error bars
% grab the SEMS from our bsDataStruct for our CI error bars  
bsSEMs_cd = bsDataStruct.semS(:,[1 3 5]); 
bsSEMs_iovd = bsDataStruct.semS(:,[2 4 6]); 

figure()
PaperPositionMode = 'manual';
orient('landscape')
for ii = 1:size(pctCrctPerA,2)/2
    subplot(1,size(pctCrctPerA,2)/2,ii)
    %     errorbar(Alpha(:,1), pctCrctPerA_cd(:,ii),bsSEMs_cd(:,ii),'.','MarkerSize',20,'CapSize',0); % just bsSEM if 68% (+-1 std) if 95, multiply by 1.96 (note, bsSEM in either direction is the same as plotting lci and upperci)
    errorbar(Alpha(:,1), pctCrctPerA_cd(:,ii),bsSEMs_cd(:,ii)*1.96,'.--','MarkerSize',30,'CapSize',0,'LineWidth',2); % *1.96 for 95% conf
    hold on
    %     errorbar(Alpha(:,1).*1000, pctCrctPerA_iovd(:,ii),bsSEMs_iovd(:,ii),'.','MarkerSize',20,'CapSize',0); % just bsSEM if 68% (+-1 std) if 95, multiply by 1.96
    errorbar(Alpha(:,1), pctCrctPerA_iovd(:,ii),bsSEMs_iovd(:,ii)*1.96,'.--','MarkerSize',30,'CapSize',0,'LineWidth',2); % just bsSEM if 68% (+-1 std) if 95, multiply by 1.96
    grid off
    xlim([0 200])
    xticks([22.5 45 90 180])
    ylim([45 100])
    yticks([50:10:100])
    formatFigure('Alpha (degrees)', 'Percent correct', "Subj"+ii, 18, 14);
    %     legend('Disparity-based','Velocity-based','Location','Northeastoutside')
    set(gcf,'units','points','position',[50,50,1400,600])
    box off
    axis square
end

%% plot spatial aggregate -- figure #

bsAggSEM_cd = bsDataStruct.aggSEMs(:,1);
bsAggSEM_iovd = bsDataStruct.aggSEMs(:,2);

h = figure();
h.PaperPositionMode = 'manual';
orient(h,'landscape')
c = errorbar(Alpha(:,1), aggDataS.CD.pctCrctPerAlpha*100,bsAggSEM_cd*1.96,'.--','MarkerSize',30,'Capsize',0,'LineWidth',2); % * 1.96 for 95% conf intervals
c.Color = [0 0.4470 0.7410];
c.CapSize = 0;
hold on
i = errorbar(Alpha(:,1), aggDataS.IOVD.pctCrctPerAlpha*100,bsAggSEM_iovd*1.96,'.--','MarkerSize',30,'Capsize',0,'LineWidth',2); % * 1.96 for 95% conf intervals
i.Color = [0.9100 0.4100 0.1700];
i.CapSize = 0;
grid off
xlim([0 190])
xticks(round(Alpha)) %fix this
ylim([50 100])
yticks([50:10:100])
formatFigure('Alpha (deg)', 'Percent correct', 'Bootstrapped aggregate data - spatial', 18, 14);
legend('Disparity-based','Velocity-based','Location','Northeastoutside')
set(gcf,'units','points','position',[50,50,800,600])
box off
