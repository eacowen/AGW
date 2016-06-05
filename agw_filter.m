function [Data,Time]=agw_filter(Data, Time, DataMin, DataMax)
%
% [DataF,TimeF]=agw_filter(Data, Time, DataMax, DataMin)
%
% Function to threshold and adaptive Gaussian filter data that is a
% function of one variable (assumed to be time here but could be one
% spatial coordinate).
%
% DataF - returned matrix of filtered data
% TimeF - returned vecotr of times of filtered data
% Data  - Data matrix to be filtered - assumed to be a function of 1
%         variable and that dependence coincides with the longer dimension of the
%         matrix
% DataMax - Maximum allowable value in the data set
% DataMin - Minimum allowable value in the data set
%
% Note: Cowen encourages setting DataMin/DataMax as wide as possible, it is
% really there to remove non-physical spurious data - the AGW portion will 
% remove the rest.
%
% Written by Edwin A. Cowen III
% Copyright reserved 2006
% Contact Cowen at eac20@cornell.edu before distributing this software.
%
% Note on May 6,2016 Cowen edited input argument order so that DataMin
% comesbefore DataMax meeting user expectations.  Older versions have this
% swapped.

[NI, NJ]=size(Data);
%Test if data has time access running in column direction 
%Note: assume that longer axis is time axis
%Take transpose if test not met
if NI>NJ
    Data=Data';
    [NI, NJ]=size(Data);
end
Nstart=NJ;
%Threshold filter the data to lie between DataMin and DataMax
valid=(Data <= DataMax & Data >= DataMin);
if NI>1
  ind=reshape(valid,NI,NJ);
 %Remove data at a given time that does not meet threshold in any component
  valid=(sum(ind)==NI);
else
  ind=valid;
end
%Return to "Data" filtered data and time of valid data
Data=Data(:,find(valid));
Time=Time(find(valid));
[NI, NJ]=size(Data);
disp(sprintf('Started with %d, filtered to %d, threshold',Nstart,NJ))
% Adaptive Gaussian Filter

coef=[2 1.25 1.15 1.1 1.06 1.03 1.01]; %Likely overkill - can stop at 1.15

% First pass through data based on robust statistics - median and IQR-baed
% filter width
for i=1:length(coef)
  DOF=NJ-1; %degrees of freedom
  p=1/(2*NJ); %two-sided probability of exceedence event given record size
  filtFact=-tinv(p,DOF); %coefficient of dynamic filter based on Student's t
  t_fact=tinv(0.75,DOF); %coefficient for normalzing width of IQR based on student's-t
  DataMedian=median(Data')';
  s=iqr(Data')'/(2*t_fact); %estimate standard deviation based on IQR
  HiLim=DataMedian+coef(i)*filtFact*s; %set high limit for data
  LoLim=DataMedian-coef(i)*filtFact*s; %set low limit for data
  HiLim=repmat(HiLim,1,NJ);
  LoLim=repmat(LoLim,1,NJ);
  valid=(Data <= HiLim & Data >= LoLim); %test for data that is within limits
  
  if NI>1
    ind=reshape(valid,NI,NJ);
    %Remove data at a given time that does not meet threshold in any component
    valid=(sum(ind)==NI);
  else
    ind=valid;
  end
  
  %Return to "Data" filtered data and time of valid data
  Data=Data(:,find(valid));
  Time=Time(find(valid));
  [NI, NJ]=size(Data);
  disp(sprintf('Started with %d, filtered to %d',Nstart,NJ))
end

%continue iteratively filtering based on Gaussian statistics now until
%NJ=NJprev, i.e., until number of data points in time passing filter test
%remains unchanged.

NJprev=NJ+1; %initialize NJprev to be larger than NJ
while NJprev>NJ
  DOF=NJ-1; %degrees of freedom
  p=1/(2*NJ); %two-sided probability
  filtFact=-tinv(p,DOF); %coefficient of dynamic filter
  DataMean=mean(Data')';  %note switched to mean in loop
  s=std(Data')'; %note switched to standard deviation in loop
  HiLim=DataMean+coef(i)*filtFact*s; %set high limit for data
  LoLim=DataMean-coef(i)*filtFact*s; %set low limit for data
  HiLim=repmat(HiLim,1,NJ);
  LoLim=repmat(LoLim,1,NJ);
  valid=(Data <= HiLim & Data >= LoLim); %test for data that is within limits
  
  if NI>1
    ind=reshape(valid,NI,NJ);
    %Remove data at a given time that does not meet threshold in any component
    valid=(sum(ind)==NI);
  else
    ind=valid;
  end
  
  %Return to "Data" filtered data and time of valid data
  Data=Data(:,find(valid));
  Time=Time(find(valid));
  NJprev=NJ;
  [NI, NJ]=size(Data);
  disp(sprintf('Started with %d, filtered to %d, in iterative loop',Nstart,NJ))
end
