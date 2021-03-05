function [bsDataStruct] = analysis_bootstrapCIs(responsesT,responsesS,Dur,Alpha)
%analysis_bootstrapCIs Runs bootstrap simulations for individual and
%aggregate data to provide CIs for the data and fits 

nBootReps = 10000;

%% bootstrap analysis for temporal experiment

% bootstrapping CIs for the individual data
bsDataT = zeros(200,16,8,nBootReps);
% 10k reps of our experiment (pulling 200 new trials from the data for each
% dur/subj 10k times)
h = waitbar(0,'Please wait...');
for kk = 1:nBootReps
    for jj = 1:size(responsesT,2)
        for ii = 1:size(responsesT,3)
            bsDataT(:,jj,ii,kk) = randsample(responsesT(:,jj,ii),length(responsesT),true); % sampling w/ replacement from our data
            % bsSampledMeans(kk,1) = mean(bsData); % store our means for each of our 10k replications
            bsDataPctCrctT(jj,ii,kk) = sum(bsDataT(:,jj,ii,kk) == 1)/size(bsDataT,1);     % find pct correct per duration per subject for this rep
        end
    end
    waitbar(kk/nBootReps)
end
close(h)

bsDataPctCrctT = bsDataPctCrctT * 100; % get into percentages
% evaluate/visualize the distribution/pct corrects
bsDataStruct.meanPctCrctT = mean(bsDataPctCrctT,3); % make sure these are basically equal to the pct corrects from our data for each duration/subject. should be very close if not equal
bsStdPctCrctT = std(bsDataPctCrctT,0,3); % very small std expected
% want standard error/68 and 95% conf intervals for the means of the pct
% corrects (10 k reps)
bsDataStruct.semT = bsStdPctCrctT; % std of the sampling distribution is the standard error


%% bootstrapping CIs for aggregate (stacking all subjects on top of each other so we have 800 trials instead of 200 from 4 subjects
bsAggDataT = zeros(800,16,2,nBootReps); % 800 trials x 16 durs x 2 conditions (will split after)
% need to reshape responses to have 800 x 16 x 2 (stacking subjects on top
% can do this in fewer steps but wanted to be explicit
responsesT_cd = responsesT(:,:,[1 3 5 7]);
responsesT_iovd = responsesT(:,:,[2 4 6 8]);
aggResponsesT(:,:,1) = [responsesT_cd(:,:,1);responsesT_cd(:,:,2);responsesT_cd(:,:,3);responsesT_cd(:,:,4)];
aggResponsesT(:,:,2) = [responsesT_iovd(:,:,1);responsesT_iovd(:,:,2);responsesT_iovd(:,:,3);responsesT_iovd(:,:,4)];
% aggPctCrcts = squeeze(sum(aggResponses,1)./length(aggResponses))*100; % just confirming these match our aggDataT.IOVD.pctCrctPerDur*100

% pulling 800 new trials from the data for each dur/condition 10k times
h = waitbar(0,'Please wait...');
for kk = 1:nBootReps
    
    for jj = 1:size(aggResponsesT,2)
        for ii = 1:size(aggResponsesT,3)
            bsAggDataT(:,jj,ii,kk) = randsample(aggResponsesT(:,jj,ii),length(aggResponsesT),true); % sampling w/ replacement from our data
            bsAggDataPctCrctT(jj,ii,kk) = sum(bsAggDataT(:,jj,ii,kk) == 1)/size(bsAggDataT,1);     % find pct correct per duration per subject for this rep
            
        end
    end
    waitbar(kk/nBootReps)
end
close(h)

bsAggDataPctCrctT = bsAggDataPctCrctT * 100;
bsDataStruct.aggMeanPctCrctT = mean(bsAggDataPctCrctT,3); % make sure these are basically equal to the pct corrects from our data for each duration/subject. should be very close if not equal
bsAggStdPctCrctT = std(bsAggDataPctCrctT,[],3); %
bsDataStruct.aggSEMt = bsAggStdPctCrctT; % std of the sampling distribution is the standard error
% 95% ci
uci95 = bsDataStruct.aggMeanPctCrctT + 1.96.*  bsDataStruct.aggSEMt; % approximate upper 95 pct bound
lci95 = bsDataStruct.aggMeanPctCrctT - 1.96.*  bsDataStruct.aggSEMt; % approximate lower bound
% 68% ci
uci68 = bsDataStruct.aggMeanPctCrctT + 1.*  bsDataStruct.aggSEMt; % approximate upper 68% pct bound
lci68 = bsDataStruct.aggMeanPctCrctT - 1.*  bsDataStruct.aggSEMt; % approximate lower bound

% Have a look at the bootstrapped sampling distributions
figure()
for ii = 1:size(bsAggDataPctCrctT,1)
    subplot(4,4,ii)
    histogram(bsAggDataPctCrctT(ii,1,:),30)
    xline(bsDataStruct.aggMeanPctCrctT(ii,1),'k','LineWidth',2); % black xline at mean
    xline(uci95(ii,1),'r','LineWidth',2); % red xlines at lci and uci
    xline(lci95(ii,1),'r','LineWidth',2); % red xlines at lci and uci
    hold on
    histogram(bsAggDataPctCrctT(ii,2,:),30)
    xline(bsDataStruct.aggMeanPctCrctT(ii,2),'k','LineWidth',2); % black xline at mean
    xline(uci95(ii,2),'r','LineWidth',2); % red xlines at lci and uci
    xline(lci95(ii,2),'r','LineWidth',2); % red xlines at lci and uci
    title("Dur"+ii)
end
sgtitle("Boostrapped distributions of pct correct 16 durations (10 k reps)")

%% Bootstrapping tau fits/cis for aggregate

bsAggDataTStacked = reshape(bsAggDataT,12800,2,10000); % want responses in one long vector for feeding into the fitting algorithm
nTrialsPerDur_agg = (size(bsAggDataTStacked,1)/length(Dur));
repDurs_cd_agg = sort(repmat(Dur(:,1),nTrialsPerDur_agg,1));
repDurs_iovd_agg = sort(repmat(Dur(:,2),nTrialsPerDur_agg,1));
durs_cd = unique(Dur(:,1));
durs_iovd = unique(Dur(:,2));

bsAggDataT_cd = squeeze(bsAggDataTStacked(:,1,:));
clear indFitVals_cd_agg;
clear fitParams_cd_agg;
indFitVals_cd_agg = zeros(16,10000);
fitParams_cd_agg = zeros(3,10000);

%%fitting for disparity-based agg data
f_cd = @(b,durs_cd) b(1).*exp(-1*durs_cd./b(2))+b(3);
options = optimset('MaxFunEvals',1e5);

h = waitbar(0,'Please wait...'); % this one takes a few minutes

for ii = 1:size(bsAggDataT_cd,2)
    fitParams_cd_agg(:,ii) = fminsearch(@(b) norm(bsAggDataT_cd(:,ii) - f_cd(b,repDurs_cd_agg)), [0; 0; 0], options);
    indFitVals_cd_agg(:,ii) = f_cd(fitParams_cd_agg(:,ii),durs_cd);
    waitbar(ii/size(bsAggDataT_cd,2))
end
close(h)

clear tau_cd_agg;
tau_cd_agg(:,:)  = fitParams_cd_agg(2,:);
tau_cd_agg = round(tau_cd_agg,4);
weirdTauInds_cd = find(tau_cd_agg > 1); % need to reject taus greater than 1
tau_cd_agg(weirdTauInds_cd) = nan; % need to reject taus greater than 1
meanTau_cd_agg = mean(tau_cd_agg,2,'omitnan') %
% std(tau_cd,3,'omitnan') % can i get this with nans?
meanIndFitVals_cd_agg = mean(indFitVals_cd_agg,2)*100;
% conf intervals for fitted values
bsAggStdTau_cd = std(tau_cd_agg,[],2); %
bsAggSEMtau_cd = bsAggStdTau_cd; % std of sampling distribution is sem
% 95% ci
uci95tau_cd = meanTau_cd_agg + 1.96.* bsAggSEMtau_cd; % approximate upper 95 pct bound
lci95tau_cd = meanTau_cd_agg - 1.96.* bsAggSEMtau_cd; % approximate lower bound

% same thing for iovd (could do this in one loop it just takes a whle so
% w.e)
bsAggDataT_iovd = squeeze(bsAggDataTStacked(:,2,:));
clear indFitVals_iovd_agg;
clear fitParams_iovd_agg;
indFitVals_iovd_agg = zeros(16,10000);
fitParams_iovd_agg = zeros(3,10000);
%%fitting for agg data
f_iovd = @(b,durs_iovd) b(1).*exp(-1*durs_iovd./b(2))+b(3);
options = optimset('MaxFunEvals',1e5);

h = waitbar(0,'Please wait...'); % again, takes a while because the bootstrapped data set is large
for ii = 1:size(bsAggDataT_iovd,2)
    fitParams_iovd_agg(:,ii) = fminsearch(@(b) norm(bsAggDataT_iovd(:,ii) - f_cd(b,repDurs_iovd_agg)), [0; 0; 0], options);
    indFitVals_iovd_agg(:,ii) = f_iovd(fitParams_iovd_agg(:,ii),durs_iovd);
    waitbar(ii/size(bsAggDataT_iovd,2))
end
close(h)

clear tau_iovd_agg;
tau_iovd_agg(:,:)  = fitParams_iovd_agg(2,:);
tau_iovd_agg = round(tau_iovd_agg,4);
weirdTauInds_iovd = find(tau_iovd_agg > 1); % need to reject taus greater than 1
tau_iovd_agg(weirdTauInds_iovd) = nan; % need to reject taus greater than 1
meanTau_iovd_agg = mean(tau_iovd_agg,2,'omitnan') %
meanIndFitVals_iovd_agg = mean(indFitVals_iovd_agg,2)*100; % these are basically the same as the actualy pctCrct (confirm from saved aggData fitted values)
%  just need CIs now...
bsAggStdTau_iovd = std(tau_iovd_agg,[],2); %
bsAggSEMtau_iovd = bsAggStdTau_iovd;
% 95% ci
uci95tau_iovd = meanTau_iovd_agg + 1.96.* bsAggSEMtau_iovd; % approximate upper 95 pct bound
lci95tau_iovd = meanTau_iovd_agg - 1.96.* bsAggSEMtau_iovd; % approximate lower bound

% Histogram of taus w/ confidence intervals
figure()
PaperPositionMode = 'manual';
orient('landscape')
histogram(tau_cd_agg,'Normalization','probability','FaceColor',[0 0.4470 0.7410]) % ok for iovd looks good except for subj 3
hold on
xline(meanTau_cd_agg,'k','LineWidth',2);
xline(uci95tau_cd,'r','LineWidth',2);
xline(lci95tau_cd,'r','LineWidth',2);
histogram(tau_iovd_agg,'Normalization','probability','FaceColor',[0.9100 0.4100 0.1700]) % ok for iovd looks good except for subj 3
hold on
xline(meanTau_iovd_agg,'k','LineWidth',2);
xline(uci95tau_iovd,'r','LineWidth',2);
xline(lci95tau_iovd,'r','LineWidth',2);
formatFigure('Tau', 'Proportion', "Dist of taus, temporal agg data (n = 4)", 18, 14);
% throw on 95% conf intervals. will report these
% notice larger spread for the CD distribution

%%  spatial stuff (same aggregate analysis, no tau fits (only 4 points))

% spatial individual data bootstrapping for CIs 
bsDataS = zeros(400,4,6,nBootReps);

% 10k reps of our experiment (pulling 200 new trials from the data for each
% dur/subj 10k times)
for kk = 1:nBootReps
    
    for jj = 1:size(responsesS,2)
        for ii = 1:size(responsesS,3)
            bsDataS(:,jj,ii,kk) = randsample(responsesS(:,jj,ii),length(responsesS),true); % sampling w/ replacement from our data
            bsDataPctCrctS(jj,ii,kk) = sum(bsDataS(:,jj,ii,kk) == 1)/size(bsDataS,1);     % find pct correct per duration per subject for this rep     
        end
    end   
end

bsDataPctCrctS = bsDataPctCrctS * 100; % get into percentages
bsDataStruct.meanPctCrctS = mean(bsDataPctCrctS,3); % make sure these are basically equal to the pct corrects from our data for each duration/subject. should be very close if not equal
bsStdPctCrctS = std(bsDataPctCrctS,0,3); % very small std expected
bsDataStruct.semS = bsStdPctCrctS; 
 
%% bootstrapping error bars for agg spatial data
bsAggDataS = zeros(1200,4,2,nBootReps);  % should be 1200 trials x 4 alphas x 2 conditions (iovd/cd) x nBootReps

% need to reshape responses to have 1200 x 4 x 2 (stacking subjects on top
% can do this in fewer steps but wanted to be explicit for now. fix later
responsesS_cd = responsesS(:,:,[1 3 5]);
responsesS_iovd = responsesS(:,:,[2 4 6]);

aggResponsesS(:,:,1) = [responsesS_cd(:,:,1);responsesS_cd(:,:,2);responsesS_cd(:,:,3)];
aggResponsesS(:,:,2) = [responsesS_iovd(:,:,1);responsesS_iovd(:,:,2);responsesS_iovd(:,:,3)];
aggPctCrctsS = squeeze(sum(aggResponsesS,1)./length(aggResponsesS))*100; % looks good (cd is first col, iovd 2nd)

% ok these match our aggData.IOVD.pctCrctPerDur*100

% pulling 1200 new trials from the data for each dur/condition 10k times
for kk = 1:nBootReps
    
    for jj = 1:size(aggResponsesS,2)
        for ii = 1:size(aggResponsesS,3)
            bsAggDataS(:,jj,ii,kk) = randsample(aggResponsesS(:,jj,ii),length(aggResponsesS),true); % sampling w/ replacement from our data
            bsAggDataPctCrctS(jj,ii,kk) = sum(bsAggDataS(:,jj,ii,kk) == 1)/size(bsAggDataS,1);     % find pct correct per duration per subject for this rep
            
        end
    end
    
end
bsAggDataPctCrctS = bsAggDataPctCrctS * 100;
bsDataStruct.aggMeanPctCrctS = mean(bsAggDataPctCrctS,3); % make sure these are basically equal to the pct corrects from our data for each duration/subject. should be very close if not equal
bsAggStdPctCrctS = std(bsAggDataPctCrctS,[],3);
bsDataStruct.aggSEMs = bsAggStdPctCrctS; 
% 95% ci
uci95 = bsDataStruct.aggMeanPctCrctS + 1.96.* bsDataStruct.aggSEMs; % approximate upper 95 pct bound
lci95 = bsDataStruct.aggMeanPctCrctS - 1.96.* bsDataStruct.aggSEMs; % approximate lower bound
% 68% ci
uci68 = bsDataStruct.aggMeanPctCrctS + 1.* bsDataStruct.aggSEMs; % approximate upper 68% pct bound
lci68 = bsDataStruct.aggMeanPctCrctS - 1.* bsDataStruct.aggSEMs; % approximate lower bound
% split them up for ease of id/use later
%cd
bsAggMeanPctCrct_cd = bsDataStruct.aggMeanPctCrctS(:,1);
% uci_cd = uci68(:,1);
% lci_cd = lci68(:,1);
uci_cd = uci95(:,1);
lci_cd = lci95(:,1);

%iovd
bsAggMeanPctCrct_iovd = bsDataStruct.aggMeanPctCrctS(:,2);
% uci_iovd = uci68(:,2);
% lci_iovd = lci68(:,2);
uci_iovd = uci95(:,2);
lci_iovd = lci95(:,2);

% let's look at the bootstrap data and see if they look ok
figure()
for ii = 1:size(bsAggDataPctCrctS,1)
    subplot(1,4,ii)
    histogram(bsAggDataPctCrctS(ii,1,:),30)
    xline(bsAggMeanPctCrct_cd(ii),'k','LineWidth',2); % black xline at mean
    xline(uci_cd(ii),'r','LineWidth',2); % red xlines at lci and uci
    xline(lci_cd(ii),'r','LineWidth',2); % red xlines at lci and uci
    hold on
    histogram(bsAggDataPctCrctS(ii,2,:),30)
    xline(bsAggMeanPctCrct_iovd(ii),'k','LineWidth',2); % black xline at mean
    xline(uci_iovd(ii),'r','LineWidth',2); % red xlines at lci and uci
    xline(lci_iovd(ii),'r','LineWidth',2); % red xlines at lci and uci
    title("Alpha"+ii)
end
sgtitle("Boostrapped distributions of pct correct 4 alphas (10 k reps)")


%% save bsDataStruct
save('bsDataStruct.mat','bsDataStruct')

end
