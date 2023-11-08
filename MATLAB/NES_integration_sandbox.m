datapath = 'C:\work\LTER\PelagicSynthesisWG\data\';

nao = readtable([datapath 'nao_month_l.csv']);
nao.datetime = datetime(nao.Date, 'Format','MMM uuuu');
nao.driver = nao.nao;
wsw = readtable([datapath 'NES_Wslopewater.csv']);
wsw(~strcmp(wsw.Var,'WSW proportion ne channel'),:) = [];
wsw.datetime = datetime(wsw.Time,1,1);
wsw.driver = wsw.Value;
wsw.driver = (wsw.driver-nanmean(wsw.driver))./std(wsw.driver,'omitnan');
amo = readtable([datapath 'AMO_index.csv']);
amo.datetime = datetime(amo.year, amo.month,1);
amo.driver = amo.amo_index;
%let's get rid of pre-1970 for now
amo(amo.year<1970,:) = [];
%then redo the anomaly for this shorter period
amo.driver = (amo.driver-nanmean(amo.driver))./std(amo.driver,'omitnan');

%%
ecomon = readtable([datapath 'EcoMon_v3_8_wDateStrata.csv']);
ecomon.datetime = ecomon.date + timeofday(datetime(ecomon.time, 'format', 'HH:mm'));
ecomon = movevars(ecomon, 'datetime', 'before', 1);
ecomon.datetime.Format = "dd-MMM-uuuu HH:mm:ss";
%select only spring months (Mar-May) and SNE region (#2)
%SNE = 2, GOM = 4, GB = 1, MAB = 3
rows = ecomon.region==2 & ecomon.month>=3 & ecomon.month<=5;

ec_sne_spring = ecomon(rows,{'datetime' 'region' 'month' 'year'});
cvar = ecomon.Properties.VariableNames;
cvar = cvar(contains(cvar,'m2'));
temp = ecomon{rows,cvar};
temp2 = temp; temp2(temp2==0) = NaN;
ec_sne_spring(:,cvar) = array2table(log10(temp+min(temp2,[],'omitnan')));

sne_spring_group = groupsummary(ec_sne_spring, "year", "mean", cvar);

%%
clear bio driver

%top 10
top10 = {'ctyp_10m2' 'calfin_10m2' 'pseudo_10m2' 'penilia_10m2' 'tlong_10m2' 'cham_10m2'... 
    'echino_10m2' 'larvaceans_10m2' 'para_10m2' 'gas_10m2'};
evar = top10{8};  timespan = 365*2;%days
dvar = 'nao';
%evar = 'euph_10m2'; timespan = 365*2;%days
%evar = 'calfin_10m2'; timespan = 365*2;%days
%evar = 'chaeto_10m2'; timespan = 365*2;%days
bio(:,1) = datenum(sne_spring_group.year,1,1);
temp = sne_spring_group.(['mean_' evar]);
bio(:,2) = log10(temp+min(temp(temp>0))/2);
bio(isnan(bio(:,2)),:) = [];
%anomaly
bio(:,2) = (bio(:,2)-mean(bio(:,2)))./std(bio(:,2));

eval(['temp = ' dvar ';'])
driver(:,1) = datenum(temp.datetime);
driver(:,2) = temp.driver;
driver(isnan(driver(:,2)),:) = [];

[int1,int2,correlations] = doubleintegration(driver,bio,timespan);
correlations

figure
plot(driver(:,1), driver(:,2))
hold on
plot(int1(:,1), int1(:,2), '-', 'linewidth', 2)
plot(int2(:,1), int2(:,2), '-', 'linewidth',2)
plot(bio(:,1), bio(:,2), '*')
title({[evar ' ' dvar ' SNE spring'] ; ['correlations: ' num2str(correlations,2)]}, 'Interpreter','none')
datetick
ylabel('Anomaly log10 values')
print(gcf, [regexprep(datapath, 'data', 'figures') evar ' ' dvar ' SNE spring.png'], '-dpng')