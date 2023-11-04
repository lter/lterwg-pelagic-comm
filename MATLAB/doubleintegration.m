function [int1,int2,correlations] = doubleintegration(driver,bio,timespan)

%This code takes as input:
%        a physical driver (2-column matrix, first column is time and second column is magnitude of index)
%        a biological time-series  (2-column matrix, first column is time and second column is properly transformed measured of abundance of biomass)
%        the representative timespan of the organism (units should be same as the units of your time-series (most likely years)

%This code returns as output:
%        the first temporal integration (really a temporal average) of the physical driver over the timespan of the organism
%        the second temporal integration (really a temporal average) of the physical driver over the timespan of the organism
%        correlations between the biological time-series and the physical driver, the first integration, and the second integration

year = driver(:,1);
phys = driver(:,2);

ind1 = find((year>min(year)+timespan));  %These are the indices that work for the integration (i.e., more than 1 time span past the beginning of the dataset)

for i=ind1'
    int1(i,1) = year(i);
    ind2 = find(year>year(i)-timespan & year<=year(i)); %indices that are no more than 1 timespan before the measurement date
    int1(i,2) = mean(phys(ind2));
end
int1(find(int1(:,1)==0),:)=[];

%Repeating for a second integration

year = int1(:,1);
phys = int1(:,2);

ind1 = find((year>min(year)+timespan));

for i=ind1'
    int2(i,1) = year(i);
    ind2 = find(year>year(i)-timespan & year<=year(i));
    int2(i,2) = mean(phys(ind2));
end
int2(find(int2(:,1)==0),:)=[];

%Now matching to the biological time-series
%for i=1:height(bio)
%    bio(i,3) = interp1(driver(:,1),driver(:,2),bio(i,1));
%    bio(i,4) = interp1(int1(:,1),int1(:,2),bio(i,1));
%    bio(i,5) = interp1(int2(:,1),int2(:,2),bio(i,1));
%end
%heidi - streamline interpolation w/o loop
    bio(:,3) = interp1(driver(:,1),driver(:,2),bio(:,1));
    bio(:,4) = interp1(int1(:,1),int1(:,2),bio(:,1));
    bio(:,5) = interp1(int2(:,1),int2(:,2),bio(:,1));

%And computing correlations
% heidi - handle NaN entries on edges of interpolated cols
ii = ~isnan(bio(:,3)); correlations(1) = corr(bio(ii,2),bio(ii,3));
ii = ~isnan(bio(:,4)); correlations(2) = corr(bio(ii,2),bio(ii,4));
ii = ~isnan(bio(:,5)); correlations(3) = corr(bio(ii,2),bio(ii,5));
