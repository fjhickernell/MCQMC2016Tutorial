function ConstructMCQMCTutorialFigures(reDoData)
% This function recomputes the examples and recreates the plots for the
% MCQMC 2016 talk.
%
% reDoData = true deletes all the data

if nargin < 1
   reDoData = false;
end
if reDoData
   disp('Are you sure that you want to <strong>delete</strong> all old data?')
   confirm = input('Type ''Y'' (upper or lower case) to confirm: ','s');
   reDoData = false;
   if numel(confirm)
      reDoData = strcmpi(confirm(1),'y');
   end
end
if reDoData %deleting all data files
   disp('Deleting all data files')
   delete AsianCallExampleAllData.mat
   delete discrepancyData.mat
   delete MVNProbExampleData.mat
end

% Now re-creating all plots
AsianArithmeticMeanOptionExample %Re-compute Asian Call
AsianCallExamplePlots %and plot results
MVNExample %Re-compute multivariate normal probability
MVNExamplePlots %and plot results
PlotDiscrepancy %Plot discrepancy figures