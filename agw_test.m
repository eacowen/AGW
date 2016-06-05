% Example code to test agw_filter.
% create two very simple data vectors a and b based on the normal
% distribution with some uniform noise added in.  The noise in a appears in
% the final 200 points.  The noise in b appeas in the first 100 and final
% 100 points.  The normal data is 1000 points long.  Data is assumed
% temporal with time vector t

clear all
close all

%Create a and b data vectors and tiem vector t
a=[randn(1000,1); (rand(200,1)-0.5)*40];
b=[(rand(100,1)-0.5)*40; randn(1000,1); (rand(100,1)-0.5)*40 ]
t=0:(length(a)-1);

%Plot raw data
figure
plot(t,a,'r',t,b,'g')
legend ('a','b','Location','Best')
xlabel('t (time)')
ylabel('a,b (value)')

%Histograms of raw data
figure
subplot(121)
hist(a,201);
xlabel('a')
ylabel('Occurrence count')
subplot(122)
hist(b,201);
xlabel('b')
ylabel('Occurrence count')

%Set just inside the maximum value - this is arbitrary, for PIV I usually
%set at +/- half the subwindow length scale if there is no initial guess.
%It shoudl be sould very wide - just to remove something the user absoluely
%knows to be bad data but corrupts the calcuation of the robust statistics
%(median and IQR)
DataMin=-19;
DataMax=19; 

%Build a single Data array with both data vectors as a monotonic function
%of time
Data=[a, b];
%Call agw_filter to filter the data
[DataF,TimeF]=agw_filter(Data, t, DataMin, DataMax);

%retrieve individual filtered data vectors
aF=DataF(1,:)';
bF=DataF(2,:)';

%Histograms of filtered data
figure
subplot(121)
hist(aF,101);
xlabel('a')
ylabel('Occurrence count')
subplot(122)
hist(bF,101);
xlabel('b')
ylabel('Occurrence count')

%example of interpolating data in

aI = interp1(TimeF,aF,t,'linear');
bI = interp1(TimeF,bF,t,'linear');

%Plot filtered data and interpolated data
figure
plot(t,aI,'ro',TimeF,aF,'r+',t,bI,'go',TimeF,bF,'g+')
legend ('aI','aF','bI','bF','Location','Best')
xlabel('t (time)')
ylabel('a,b (value)')
%Note that all of the noise has been filtered out other than those within 
%about +/- 3 std.  The missing data has been interpolated in linearly and 
%the open circles represent the interpolated data, hence the linear
%segements at the beginning and end of the temporal record.

%Plot filtered data and raw data
figure
plot(TimeF,aF,'r+',TimeF,bF,'g+',t,a,'ro',t,b,'go')
legend ('aF','bF','a','b','Location','Best')
xlabel('t (time)')
ylabel('a,b (value)')
%Note that all of the noise has been filtered out other than those within 
%about +/- 3 std.  The open circles represent the removed 'outlier' data.