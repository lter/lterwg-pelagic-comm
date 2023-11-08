datapath = 'C:\work\LTER\PelagicSynthesisWG\data\';

%marmap data
%https://www.nodc.noaa.gov/archive/arc0001/9800126/1.3/data/0-data/
chl_marmap = readtable([datapath 'c_4_nodc.txt']);
chl_marmap.datetime = datetime(chl_marmap.YEAR, chl_marmap.MON, chl_marmap.DAY, chl_marmap.HR, chl_marmap.MIN,0);
chl_marmap(1,:) = [];
chl_marmap_cruise = chl_marmap.CRUISE;
chl_marmap.CRUISE = [];
%chl_marmap_sne_rows = chl_marmap.LATD<lat_max&chl_marmap.LATD>lat_min&chl_marmap.LOND<lon_max&chl_marmap.LOND>lon_min&chl_marmap.MON<=5&chl_marmap.MON>=3;

% quick and dirty approximate SNE in marmap chl
%ecomon = readtable("C:\work\LTER\PelagicSynthesisWG\data\EcoMon_v3_8_wDateStrata.csv");
%lat_max = max(ecomon.lat(ecomon.region==2)); lat_min = min(ecomon.lat(ecomon.region==2));
%lon_max = max(ecomon.lon(ecomon.region==2)); lon_min = min(ecomon.lon(ecomon.region==2));
lat_max = 41.6850; lat_min = 39.1467;
lon_max = -68.8000; lon_min = -73.9667;
chl_marmap_sne_rows = chl_marmap.LATD<lat_max&chl_marmap.LATD>lat_min&chl_marmap.LOND<lon_max&chl_marmap.LOND>lon_min&chl_marmap.MON<=5&chl_marmap.MON>=3 & chl_marmap.DEPTH<=5;
chl_marmap.CHLA_log10 = log10(chl_marmap.CHLA);
chl_marmap_sne_spring_annual = retime(table2timetable(chl_marmap(chl_marmap_sne_rows,:)), 'yearly', 'mean');
%chl_marmap_sne_spring_annual.CHLA_log10 = log10(chl_marmap_sne_spring_annual.CHLA);
std(chl_marmap_sne_spring_annual.CHLA_log10)

%satellite data
chl = readtable([datapath 'weekly_chl_ecomon.csv']);
chl.Value_log10 = log10(chl.Value);
chl_group = groupsummary(chl, {'EPU' 'season' 'year'}, 'mean', 'Value_log10');
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
%read std tables from Alex's analyses in R
Zvol_std = readtable([datapath 'SD_trophamp_ZP.csv']);
Ffish_std = readtable([datapath 'SD_trophamp_FFISHv2.csv']);
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
print(gcf, [regexprep(datapath, 'data', 'figures') 'std trophic levels v1'], '-dpng')

return

%%
TT_SNEspring.volume_1m2
TT_SNEspring.volume_1m2_log10 = log10(TT_SNEspring.volume_1m2);
Zvol_group = groupsummary(chl, {'EPU' 'season' 'year'}, 'movmean', 'Value_log10');

