%marmap data
%https://www.nodc.noaa.gov/archive/arc0001/9800126/1.3/data/0-data/
chl_marmap = readtable("C:\Users\heidi\Downloads\c_4_nodc.txt");
chl_marmap.datetime = datetime(chl_marmap.YEAR, chl_marmap.MON, chl_marmap.DAY, chl_marmap.HR, chl_marmap.MIN,0);
chl_marmap(1,:)
chl_marmap(1,:) = [];
chl_marmap_cruise = chl_marmap.CRUISE;
chl_marmap.CRUISE = [];
%chl_marmap_sne_rows = chl_marmap.LATD<lat_max&chl_marmap.LATD>lat_min&chl_marmap.LOND<lon_max&chl_marmap.LOND>lon_min&chl_marmap.MON<=5&chl_marmap.MON>=3;
chl_marmap_sne_rows = chl_marmap.LATD<lat_max&chl_marmap.LATD>lat_min&chl_marmap.LOND<lon_max&chl_marmap.LOND>lon_min&chl_marmap.MON<=5&chl_marmap.MON>=3 & chl_marmap.DEPTH<=5;

chl_marmap_sne_spring_annual = retime(table2timetable(chl_marmap(chl_marmap_sne_rows,:)), 'yearly', 'mean');
std(chl_marmap_sne_spring_annual.CHLA_log10)

%satellite data
chl = readtable("C:\Users\heidi\Downloads\weekly_chl_ecomon.csv");
chl.Value_log10 = log10(chl.Value);
chl_group = groupsummary(chl, {'EPU' 'season' 'year'}, 'movmean', 'Value_log10');
%chl_group.
chl_group_std = groupsummary(chl_group,{'season' 'EPU'}, 'std', 'mean_Value_log10');

regionlist = unique(chl_group.EPU);
for ii = 1:length(regionlist)
    it = strcmp(regionlist{ii},chl_group.EPU);
    %presumes table sorted by year
    chl_group.runmean(it) = movmean(chl_group.mean_Value_log10(it),5);
end
chl_group_std = groupsummary(chl_group,{'season' 'EPU'}, 'std', 'runmean');
%%
Ffish_std = readtable("C:\Users\heidi\Downloads\SD_trophamp_FFISHv2.csv")
seasonlist = unique(chl_group_std.season);
figure
for ii = 1:length(seasonlist)
   subplot(2,2,ii)
   it = strcmp(seasonlist{ii},chl_group_std.season);
   gscatter(ones(sum(it),1),chl_group_std.std_runmean(it),chl_group_std.EPU(it))
   title(seasonlist{ii})
   it2 = strcmp(seasonlist{ii},Zvol_std.season);
   hold on
   gscatter(2*ones(sum(it2),1),Zvol_std.SD_ts(it2),Zvol_std.Region(it2))
   it3 = strcmp(seasonlist{ii},Zvol_std.season);
   hold on
   gscatter(3*ones(sum(it3),1),Ffish_std.SD_ts(it3),Ffish_std.Region(it3))
   ylim([0 1])
   set(gca, 'xtick', [1 2 3], 'xTickLabel', {'Chl' 'Zoop Vol' 'Forage Fish Larvae'})
   xlim([.9 3.1])
   ylabel('Std dev, 5-yr smoothed log10 annual means')
end
%%
TT_SNEspring.volume_1m2
TT_SNEspring.volume_1m2_log10 = log10(TT_SNEspring.volume_1m2);
Zvol_group = groupsummary(chl, {'EPU' 'season' 'year'}, 'movmean', 'Value_log10');

