function p=suicide_eyetracking(fname,figstoplot)  
% reads edf data from Eyelink eye-tracker, cleans/deblinks data, and calculates
% many eye-tracking indices 
% usage: p=suicide_eyetracking(fname,figstoplot)  
% where fname is the name of the .edf file
% figstoplot is optional - all figures are plotted by default
% p is a structure containing the following fields:
%    .... ** SAY WHAT THE FIELDS ARE AND WHAT THEY MEAN HERE ***
if nargin<1, fname='245fv.edf'; end   %% fname is the filename with a postfix of ".edf"
if nargin<2, figstoplot=1:100; end

fprintf('Working on file %s\n',fname);
% read the raw edf eye-tracking data into matlab by a toolbox (https://github.com/uzh/edf-converter)
pr=Edf2Mat(fname); 
%%%%%%%%%%%%%%%%%%%%
% pr.Samples has the X, Y, and pupil size data
if ismember(1,figstoplot)
  figure(1); clf;  %% clf: clean current figure window
  plot(pr.Samples.posX,pr.Samples.posY)  %% plot scanpath
  xlabel('X');
  ylabel('Y');
end

if ismember(2,figstoplot)
  figure(2); clf;
  plot(pr.Samples.time,pr.Samples.pupilSize)  %% plot pupile size
  xlabel('samples');
  ylabel('Pupil size');
end

%%%%%%%%%%%%%
% clean the pupil data and identify where blinks occurred.
% store the cleaned data in a structure called "p"
p.filename=pr.filename;
% function gotham: filtering blink and smooth it with output rescaleData,BlinkTimes,Noblinks,NoblinksUnsmoothed input should mm-unit
% This function could be seen in the repository
p=gotham_stublinks(p,ismember(3,figstoplot),0); 
% add the X and Y data to p
p.X=pr.Samples.posX;
p.Y=pr.Samples.posY;
% add saccade values that were previously unmeasured
p.Xrs=p.X;
p.Xrs(find(isnan(p.Xrs)))=-9999;
p.Xrs=interpmissing(p.Xrs,-9999);

p.Yrs=p.Y;
p.Yrs(find(isnan(p.Yrs)))=-9999;
p.Yrs=interpmissing(p.Yrs,-9999);
 % Add regions of interests in spatial coordinates
p.rois(1).Xlb=0; p.rois(1).Xub=860; 
p.rois(1).Ylb=-10; p.rois(1).Yub=515; 
p.rois(2).Xlb=1075; p.rois(2).Xub=1935;
p.rois(2).Ylb=-10; p.rois(2).Yub=515; 
p.rois(3).Xlb=1075; p.rois(3).Xub=1935;
p.rois(3).Ylb=570; p.rois(3).Yub=1100; 
p.rois(4).Xlb=0; p.rois(4).Xub=860;
p.rois(4).Ylb=575; p.rois(4).Yub=1100; 
for roi=1:4  
  p.rois(roi).InRegionSamples=find((p.X>p.rois(roi).Xlb) & (p.X<p.rois(roi).Xub) & (p.Y>p.rois(roi).Ylb) & (p.Y<p.rois(roi).Yub)); 
end  
if ismember(4,figstoplot)
  figure(4); clf;
  for ct=1:4
    subplot(2,2,ct);
    plot(p.X,p.Y);
    xlabel('X');
    ylabel('Y');
    hold on;
    plot(p.X(p.rois(ct).InRegionSamples),p.Y(p.rois(ct).InRegionSamples),'r.');
    title(ct);
  end   % visulization of ROI 
end

% Add label messages from orignal edf-files
p.time=pr.Samples.time;
p.RescaleFactor=1000;  
p.EventLabels={'Fixation gap','fixation_timer','Stimulations','stimulation-time'};
goodmessages=[];
p.EventCodes=[];
p.EventTimes=[];
for msg=1:length(pr.Events.Messages.info) 
  for lab=1:length(p.EventLabels) 
    if ~isempty(strfind(pr.Events.Messages.info{msg},p.EventLabels{lab})) 
      goodmessages(end+1)=msg;
      p.EventCodes(end+1)=lab;
      p.EventTimes(end+1)=pr.Events.Messages.time(msg);
    end
  end
end
for ct=1:length(p.EventTimes)
  p.EventTicks(ct)=position(p.EventTimes(ct),p.time);
end
% for each trial, start at 1/10th second prior to Stimulation
p.TrialStarts=p.EventTicks(find(p.EventCodes==3))-100; 
p.TrialEnds=p.TrialStarts+25100; % 25100=trial length(25000msec)+prior period (100msec) 
if p.TrialEnds(end)>length(p.X) 
    p.TrialStarts=p.TrialStarts(1:end-1);
    p.TrialEnds=p.TrialEnds(1:end-1);
end
p.TrialLengths=p.TrialEnds-p.TrialStarts;
p.StimLatencies=100.*ones(size(p.TrialStarts)); 
p.NumTrials=length(p.TrialStarts); 
p.TrialTypesNoDrops=ones(size(p.TrialStarts));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get dwell time in each ROI for each trial
p.inroi=zeros(size(p.XTrials,1),4,size(p.XTrials,2));
for ct=1:size(p.XTrials,1)
  for roi=1:4
    indsinroi=find((p.XTrials(ct,:)>p.rois(roi).Xlb) & (p.XTrials(ct,:)<p.rois(roi).Xub) & (p.YTrials(ct,:)>p.rois(roi).Ylb) & (p.YTrials(ct,:)<p.rois(roi).Yub));
    p.inroi(ct,roi,indsinroi)=1;
    p.msinroi(ct,roi)=length(indsinroi);
  end
end
        
    p.FixatTrials = zeros(size(p.YTrials));
    for trl=1:size(p.XTrials,1)
            % get rate of change
            xchg=[0 (p.XTrials(trl,2:end)-p.XTrials(trl,1:end-1))];
            ychg=[0 (p.YTrials(trl,2:end)-p.YTrials(trl,1:end-1))];
            xyvel=sqrt(xchg.^2+ychg.^2);

       for ct=1:size(p.XTrials,2)-(minfixtime+1)
         xint=(p.XTrials(trl,ct):p.XTrials(trl,ct+minfixtime));
         xrange=(max(xint)-min(xint));
         yint=(p.YTrials(trl,ct):p.YTrials(trl,ct+minfixtime));
         yrange=(max(yint)-min(yint));
         maxvel=max(xyvel(ct:ct+minfixtime));
         if ((xrange<minfixpoints) & (yrange<minfixpoints) & (maxvel<maxfixvel))
           p.FixatTrials(trl,ct:ct+minfixtime)=1;
         end
       end   
    end
   p.inroi=zeros(size(p.XTrials,1),4,size(p.XTrials,2));
   for ct=1:size(p.XTrials,1)
     for roi=1:4
       indsinroi=find((p.FixatTrials(ct,:)==1) & (p.XTrials(ct,:)>p.rois(roi).Xlb) & (p.XTrials(ct,:)<p.rois(roi).Xub) & (p.YTrials(ct,:)>p.rois(roi).Ylb) & (p.YTrials(ct,:)<p.rois(roi).Yub));
       p.fix_inroi(ct,roi,indsinroi)=1;
       p.fix_msinroi(ct,roi)=length(indsinroi); 
       p.fix_meanpupinroi(ct,roi)=mean(rescaleoutliers(p.NormedPupTrials(ct,indsinroi)')); 
       if length(indsinroi)>0
         p.latency_first_fixation(ct,roi)=indsinroi(1)./p.RescaleFactor;
         p.gaze_duration_first2500ms(ct,roi)=length(find(indsinroi<2500))./p.RescaleFactor; % define early-stage as the initial 2500 msec of each trial.
         notinroi=find(1-p.fix_inroi(ct,roi,indsinroi(1):end));
         p.latency_saccade_away_from_first_fixation(ct,roi)=(notinroi(1))./p.RescaleFactor;
       else
         p.latency_first_fixation(ct,roi)=999; % lab missing value
         p.gaze_duration_first2500ms(ct,roi)=0; % if there is no fixation in this roi label it as value "0"
         p.latency_saccade_away_from_first_fixation(ct,roi)=-999;
       end
       p.num_fixations_in_roi(ct,roi)=length(find(diff(p.fix_inroi(ct,roi,:))==1)); 
     end
     for roi=1:4     % Percentage of the first fixation location
        p.initial_fix_is_this_roi(ct,roi)=(p.latency_first_fixation(ct,roi)==min(p.latency_first_fixation(ct,:))); 
        p.dwell_time_other(ct,roi)=length(find(p.fix_inroi(ct,setdiff(1:4,roi),:)))./length(find(p.fix_inroi(ct,:,:)));  
     end
   end
% regression calculate (strict definition: count only across different ROIs)
if ismember(9,figstoplot)
 for trl=1:20  
    whichroi(trl,:)=[squeeze(p.fix_inroi(trl,1,:))+2.*squeeze(p.fix_inroi(trl,2,:))+3.*squeeze(p.fix_inroi(trl,3,:))+4.*squeeze(p.fix_inroi(trl,4,:))]'; 
    lastroi(trl,:)=whichroi(trl,:);
    lastroinonzero=find(whichroi(trl,:)>0);
    lastroizeros=find(lastroi(trl,:)==0);
    for ct=1:length(lastroizeros)
        inroipoints=lastroinonzero(find(lastroinonzero<lastroizeros(ct))); 
        if length(inroipoints)>0 
          lastinroi=inroipoints(end);    
          lastroi(trl,lastroizeros(ct))=whichroi(trl,lastinroi); 
        else
          lastroi(trl,lastroizeros(ct))=0;
        end
    end
    whichroidiff(trl,:)=diff(whichroi(trl,:)); 
    for roinum=1:4
        fixintoroi(trl,roinum)=length(find((whichroi(trl,2:end)==roinum) & (lastroi(trl,1:end-1)~=roinum) & (whichroidiff(trl,:)~=0)));
    end
 end
end
%write output
if ismember(10,figstoplot)
  dat=[[1:length(p.msinroi)]' p.msinroi p.meanpupinroi p.fix_msinroi p.latency_first_fixation p.gaze_duration_first2500ms p.dwell_time_other p.initial_fix_is_this_roi p.num_fixations_in_roi p.latency_saccade_away_from_first_fixation  fixintoroi]; 
  outfname=sprintf('%s_scanpath_results.csv',p.filename(1:end-4));  % '%s_dwell_and_pupil.csv'--> _scanpath_results.csv
  fprintf('Saving file %s\n',outfname);
  gdlmwrite(outfname,dat,',',0,0,5,['Trial,Dwell_ROI1,Dwell_ROI2,Dwell_ROI3,Dwell_ROI4,Fix_Dwell_ROI1,Fix_Dwell_ROI2,Fix_Dwell_ROI3,Fix_Dwell_ROI4,Fix_Pup_ROI1,Fix_Pup_ROI2,Fix_Pup_ROI3,Fix_Pup_ROI4,' ...
      'latency_FF_R1,latency_FF_R2,latency_FF_R3,latency_FF_R4,gaze2500_R1,gaze2500_R2,gaze2500_R3,gaze2500_R4,otherDW_R1,otherDW_R2,otherDW_R3,otherDW_R4,FF_distri_R1,FF_distri_R2,FF_distri_R3,FF_distri_R4,'...
      'Count_RI,Count_R2,Count_R3,Count_R4,latency_avoid_R1,latency_avoid_R2,latency_avoid_R3,latency_avoid_R4,count_rg_r1,count_rg_r2,count_rg_r3,count_rg_r4']); 
end
