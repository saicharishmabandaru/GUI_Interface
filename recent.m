global bi bs C ca ci le lf mui mus B di Ta rh Rs U Pp e Rl da q G_q t ds xsnow xice A T_iter lv dw

% Read Constants from Excel
    constantsData = xlsread('data_new.xlsx', 'Constants', 'B2:B14');
    %constants = table2struct(constantsData);

    % Read Observed Data from Excel
    observedData = xlsread('data_new.xlsx', 'Observed Data','B2:B12');
    %observed = table2struct(observedData);

% Extract constants
bi = constantsData(1,1);
bs = constantsData(2,1);
C = constantsData(3,1);
ca = constantsData(4,1);
ci = constantsData(5,1);
le = constantsData(6,1);
lf = constantsData(7,1);
lv = constantsData(8,1);
mui = constantsData(9,1);
mus = constantsData(10,1);
B = constantsData(11,1);
di = constantsData(12,1);
dw = constantsData(13,1);


% Extract observed data
Ta = observedData(1,1);
rh = observedData(2,1);
Rs = observedData(3,1);
U = observedData(4,1);
Pp = observedData(5,1);
e = observedData(6,1);
Rl = observedData(7,1);
da = observedData(8,1);
q = observedData(9,1);
G_q = observedData(10,1);
ds = observedData(11,1);

t     = linspace(0, 86400, 3600);     % duration of the observation, taken as 1 day = 86400s, divided into 3600 divisions
xsnow = linspace(0, 0.75, 50);        % snow layer thickness (m)
xice  = linspace(0, 19.25, 1300);     % ice layer thickness (m)
x = linspace(0, 0.75, 50); %snow thickness
s = linspace(0, 19.25, 1300); %ice thickness
t = linspace(0, 86400, 3600);
xyz = Snow_Temperature();
pqr = Ice_Temperature();

yData_endtime = xyz(end, :);
data_enddepth = xyz(:, end);
x;
t;
%a = yData;
disp(size(xyz,1))
disp("displaying yData");
disp(size(pqr,1))
disp("displaying yData");

figure;

plot(x, yData_endtime);


% Add labels and a title
xlabel('Depth');
ylabel('temp with depth for snow');
title('Plot of Last Column of a 2D Matrix');

figure;
plot(t,data_enddepth);

csvwrite('t.csv', t)
csvwrite('data_enddepth.csv',data_enddepth)

% Add labels and a title
xlabel('time');
ylabel('temp with time for snow');
title('Plot of Last Column of a 2D Matrix');

% Display the grid
grid on;

%yData = abc(end, :);
%l = yData;


%disp(yData);

dipData_endtime = pqr(end, :);
dipdata_enddepth = pqr(:, end);
s;
t;
figure;
plot(s,dipData_endtime)

csvwrite('s.csv', s)
csvwrite('dipData_endtime.csv',dipData_endtime)

% Add labels and a title
xlabel('Depth');
ylabel('temp with depth for ice');
title('Plot of Last Column of a 2D Matrix');

figure;
plot(t,dipdata_enddepth)

% Add labels and a title
xlabel('time');
ylabel('temp with time for ice');
title('Plot of Last Column of a 2D Matrix');

%______using data from website_____________

source = 'AWindSpeed1992.nc';%wind speed data
RandPpt = 'Apptandradiation.nc';
time_in_hours = ncread(source,'time');
ncdisp(source);
ncdisp(RandPpt);
t0 = datenum('1900-01-01 00:00:00', 'yyyy-mm-dd HH:MM:SS');
time = time_in_hours/24 + t0;
time = datetime(time,'ConvertFrom','datenum');
lat = ncread(source,'latitude');
long = ncread(source,'longitude');
u = ncread(source,'si10');
rs = ncread(RandPpt,'uvb');
ppt = ncread(RandPpt,'tp');
rs = rs/3600;
ppt = ppt*100000;

PD = [80.007,30.301];%(1,6)
GG = [79.149,30.8];%(18,23)
SP = [77.355,32.421];%(4,6)
ZM = [88.238,27.743];%(91,43)

long_lat = [1 6; 18 23; 4 6; 91 43];

u_PD = u(1,6,:);
u_GG = u(18,23,:);
u_SP = u(4,6,:);
u_ZM = u(91,43,:);


rs_PD = rs(1,6,:);
rs_GG = rs(18,23,:);
rs_SP = rs(4,6,:);
rs_ZM = rs(91,43,:);


ppt_PD = ppt(1,6,:);
ppt_GG = ppt(18,23,:);
ppt_SP = ppt(4,6,:);
ppt_ZM = ppt(91,43,:);


slope_angle = 13;
slope_correction_factor = cosd(13); % calculate the slope correction factor

% main input- Albedo = A
%-----------------------------------Note----------------------------------%
% For sensitivity analysis, the value of albedo is varied from 1 to 0.5 with
% a diffference of 0.05 to obtain surface temperature of glacier and runoff
% variation withchange in albedo

%---------------------------------Literature------------------------------%

% mi, mi1 = Mass of water due to melting of ice (kg) (both mi and mi1)
% mi0 = Initial Mass of water due to melting of ice (kg)

% Hg  = Heat of into glacier (W/m^2)
% Hg0 = Initial Heat of into glacier (W/m^2)

% Ts, Ts1 = Surface temperature of glacier (degree celcius)
% Ts0 = Initial Surface temperature of glacier (degree celcius)

% Fs = Mass of superimposed ice on ice surface below snow (mm.w.e)
% Fc = Mass of water that is refrozen within snow layer (mm.w.e)
% Fc_capillary = Mass of capillary water that is refrozen within wet snow(mm.w.e)

% Pr = Rainfall (mm.w.e.d^-1)
% Ps = Snowfall (mm.w.e.d^-1)
% DS = Runoff   (mm.w.e.d^-1)

%------------------Initialisation------------------------------------------%

albedo = [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8];
%A PD = [0.59883487, 0.61353904, 0.6244669, 0.6248883, 0.63040435, 0.63220996, 0.63295043,0.63626844,0.63753074,0.6416505,0.64568317,0.64894283,0.64984846,0.65084326,0.65176386,0.65237373,0.65237463,0.65250212,0.65271711,0.65293831,0.65332711,0.65332836,0.65351206,0.65355659,0.65355694,0.65356952,0.65358955,0.65359694, 0.65367079, 0.65369797, 0.65370345, 0.65372044, 0.65372288, 0.65375221, 0.65375376, 0.65380442, 0.65381986, 0.6538409, 0.65416652]
%A GG = [0.59152275, 0.59505624, 0.6049614, 0.61089849, 0.61387306, 0.6139589, 0.61444283, 0.61635548, 0.61647862, 0.621333, 0.62281996, 0.62461764, 0.62545967, 0.62639415, 0.62697369, 0.62712008, 0.62738931, 0.62739009, 0.62767321, 0.62781888, 0.62787408, 0.62801087, 0.62802565, 0.62807804, 0.62808359, 0.62812847, 0.62817222, 0.62818319, 0.62819773, 0.62820554, 0.62822956, 0.62823278, 0.6282357, 0.62823606, 0.62824273, 0.62824744, 0.62828493, 0.62829089, 0.62831116, 0.63222039]
%A ZM = [0.46760306, 0.4729085, 0.47917458, 0.48350039, 0.4876276, 0.49561962, 0.49595994, 0.49634272, 0.50642443, 0.50827199, 0.50849235, 0.51746631,0.51831311, 0.51839966, 0.51858681, 0.5205496, 0.52082467, 0.52147043, 0.52237344, 0.52243561, 0.5231697, 0.52362847, 0.52407891, 0.52449125, 0.52464348, 0.52491993, 0.52518564, 0.52538353, 0.52603596, 0.52614641, 0.52615494, 0.52675587, 0.52764356, 0.52776802, 0.53069478, 0.53195834, 0.53268677, 0.53986633, 0.5403527, 0.55143756]
%A SP = [0.52159935, 0.52887732, 0.53576267, 0.54226583, 0.54583514, 0.55148113, 0.5537923, 0.55759686, 0.55956674, 0.56186664, 0.56471324 0.56729442, 0.56908906, 0.57112736, 0.57360232, 0.57434726, 0.57514292, 0.57588124, 0.57715142, 0.57778811,0.57809073, 0.57840586, 0.57855898, 0.57860255, 0.57875502, 0.57875645, 0.57879126, 0.57885867, 0.57888168, 0.57891196, 0.57895666, 0.57896644, 0.578987, 0.5791136, 0.57916594, 0.57921261, 0.57943833, 0.57944179, 0.57947719, 0.58219308 ] 
AR_of_PD = [];
AR_of_GG = [];
AR_of_SP = [];
AR_of_ZM = [];

mi0 = 0; %initially there won't be any melt water ie, in the initial stage, ice won't melt

Hg0 = 0; % Initial Heat transfer into the glacier (W/m^2)

%---------------- Hypsometry data--------------------------------------------------

%Precipitation
filename = 'Precipitation.xls';
sheet = 1;
xlRange = 'B44:B74';

Precipitation = xlsread(filename,sheet,xlRange)
Precipitation = Precipitation';
size(Precipitation);

%Temperature
filename = 'temp.xls';
sheet = 1;
xlRange = 'B42:B72';

temperature = xlsread(filename,sheet,xlRange)
temperature = temperature';
size(temperature);





%ZM

filename = 'ZM_hp.xls';
sheet = 1;
xlRange = 'A1:A31';

ZM_elevation = xlsread(filename,sheet,xlRange)
ZM_elevation = ZM_elevation';
size(ZM_elevation);

filename = 'ZM_hp.xls';
sheet = 1;
xlRange = 'B1:B31';

ZM_Area = xlsread(filename,sheet,xlRange)
ZM_Area = ZM_Area';
size(ZM_Area);


%PD

PD_filename = 'PD_hp.xls';
sheet = 1;
xlRange = 'A1:A27';

PD_elevation = xlsread(filename,sheet,xlRange)
PD_elevation = PD_elevation';
size(PD_elevation);

filename = 'PD_hp.xls';
sheet = 1;
xlRange = 'B1:B27';

PD_Area = xlsread(filename,sheet,xlRange)
PD_Area = PD_Area';
size(PD_Area);


%GG

filename = 'GG_hp.xls';
sheet = 1;
xlRange = 'A1:A26';

GG_elevation = xlsread(filename,sheet,xlRange)
GG_elevation = GG_elevation';
size(GG_elevation);

filename = 'GG_hp.xls';
sheet = 1;
xlRange = 'B1:B26';

GG_Area = xlsread(filename,sheet,xlRange)
GG_Area = GG_Area';
size(GG_Area);


%SP
filename = 'SP_hp.xls';
sheet = 1;
xlRange = 'A1:A17';

SP_elevation = xlsread(filename,sheet,xlRange)
SP_elevation = SP_elevation';
size(SP_elevation);

filename = 'SP_hp.xls';
sheet = 1;
xlRange = 'B1:B17';

SP_Area = xlsread(filename,sheet,xlRange)
SP_Area = SP_Area';
size(SP_Area);

% Create a GUI to prompt the user for the glacier name
glacier_name = inputdlg('Enter the name of the glacier (e.g., PD, GG, SP, ZM):', 'Glacier Name', [1, 50]);

% Extract the glacier name from the cell array returned by inputdlg
glacier_name = glacier_name{1};

% Use the provided glacier name
fprintf('Glacier Name: %s\n', glacier_name);

% Define the mapping of glacier names to elevation, area, u, rs, and ppt data
glacier_elevation_map = containers.Map({'PD', 'GG', 'SP', 'ZM'}, ...
    {PD_elevation, GG_elevation, SP_elevation, ZM_elevation});
glacier_area_map = containers.Map({'PD', 'GG', 'SP', 'ZM'}, ...
    {PD_Area, GG_Area, SP_Area, ZM_Area});
glacier_u_map = containers.Map({'PD', 'GG', 'SP', 'ZM'}, ...
    {u_PD, u_GG, u_SP, u_ZM});
glacier_rs_map = containers.Map({'PD', 'GG', 'SP', 'ZM'}, ...
    {rs_PD, rs_GG, rs_SP, rs_ZM});
glacier_ppt_map = containers.Map({'PD', 'GG', 'SP', 'ZM'}, ...
    {ppt_PD, ppt_GG, ppt_SP, ppt_ZM});

% Check if the provided glacier name is valid
if isKey(glacier_elevation_map, glacier_name)
    % Obtain the data for the specified glacier
    glacier_elevation = glacier_elevation_map(glacier_name);
    glacier_area = glacier_area_map(glacier_name);
    u_glacier = glacier_u_map(glacier_name);
    rs_glacier = glacier_rs_map(glacier_name);
    ppt_glacier = glacier_ppt_map(glacier_name);
    
else
    % Display an error message for an invalid glacier name
    disp('Invalid glacier name. Please enter one of the following: PD, GG, SP, ZM');
end


%----------------------------Program---------------------
Runoff_cumulative_ZM =[];
Runoff_cumulative_PD =[];
Runoff_cumulative_GG =[];
Runoff_cumulative_SP =[];
Runoff_cumulative = [];

 for cnt=1:length(u_PD)
     U=u_PD(cnt);
     Rs = rs_PD(cnt);
     Pp = ppt_PD(cnt);

     
     p =0;
   for A = 0.3:0.05:0.8   % Sensitivity analysis of albedo on runoff
       p=p+1;
            for j = 1:length(PD_elevation) % for each 100m elevation mass balance and runoff components are calculated
              
                    %Ta = temperature(j);
                    %Pp = Precipitation(j);
                    Ta = -10.27 - j*100*(0.0098);
                    Pp    = 670*(1+(PD_elevation(j+1-1)-5000)/10000);
                 
                    
                
                        mi0 = 0; %initially there won't be any melt water ie, in the initial stage, ice won't melt
                        
                        Hg0 = 0; % Initial Heat transfer into the glacier (W/m^2)
                        
                         % Calculation of Ts = Surface Temperature (degree Celcius)
                            
                        Ts0 = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg0); %with Hg0 = 0
                            
                        mi1 = Melted_ice(Ts0, A, Rs, Rl, e, da, U, Ta, rh, q, Hg0);
                            
                        Hg = Heat_Into_Glacier(mi1, Ts0, t, ds, xsnow, xice); %with initial conditions
                            
                        Ts1 = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg); %with obtained Hg
                        
                        
                        while(abs(Ts1 - Ts0) >= 0.1)
                                mi1 = Melted_ice(Ts1, A, Rs, Rl, e, da, U, Ta, rh, q, Hg);
                                Hg = Heat_Into_Glacier(mi1, Ts1, t, ds, xsnow, xice);
                                Ts0 = Ts1; %1st loop Ts
                                Ts1 = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg); %2nd loop Ts
                                    
                        end
                        mi = mi1;
                        Ts = Ts1;
                        Surfacetemp(p,j) = Ts;
                        disp('p')
                        disp(p)
                        disp('j')
                        disp(j)
                        disp('ts')
                        disp(Ts)
                        Fc_capillary=0;
                        Fc=0;
                        Fs=0;
                        if Ts>=0 %ice melts but some water might be refrozen in 3 forms
                                    %mi
                            Fs = Superimposed_ice(ds, xsnow);
                            Fc_capillary = Mass_of_Refrozen_Capillary_Water(ds, Ts);
                            Fc = Mass_of_Snow_Percolated_Water(xice, t, ds);
                                %else % Ts<0 => no meltwater   
                        end
                                
                        F = [Fs,Fc,Fc_capillary];
                            
                                %Precipitations
                        [Pr,Ps] = Precipitations(Pp,Ta); %rainfall and snowfall
                            
                        %Sublimation
                        sub = Sublimation(da, U, rh, q, Ta, Ts,lv );
                                %Final Output of the code:-
                        DS = Pr - (Fc_capillary + Fs + Fc) + mi;
                        if size(DS) ~= [0 0]

                              Runoff(p,j) = DS*1000;    % Output Variable
                        end

                        MB = (Ps +Pr - DS -mi - sub)/dw; 
                        if size(MB) ~= [0 0]
                                Mass_bal(p, j) = MB*1000;  % Output Variable
                        end

                        j = j+1;
                    %end
                    
                
                end
                
                
                Area_cum =0; % runoff is averaged for the glacier
                Runoff_cum =0.0;
               
                for k = 1:length(PD_Area)
                   
                    check = Runoff(p,k)*PD_Area(1,k);
                    if(~isnan(check))
                       Runoff_cum = Runoff_cum + round(check,1);
                    
                    end
                    
                    disp('runoff cum inside loop')
                    disp(Runoff_cum)
                    Area_cum = Area_cum + PD_Area(1,k);
                end
                Runoff_cumulative(p)= (1000*Runoff_cum)/Area_cum;

              
             
   end
 end
 disp('runoff cum after loop')
 disp(Runoff_cumulative)

 Runoff_cumulative_data = Runoff_cumulative';

% Concatenate the data
data = [albedo', Runoff_cumulative_data];

% Write to CSV
filename = 'Runoff_Albedo.csv';
csvwrite(filename, data);

%{

  Mass_balnew = Mass_bal;
  Mass_balnew(isnan(Mass_balnew))=0;
  figure('visible', 'on');
  scatter(Mass_balnew(1,:),ZM_elevation)
  disp('mass balance')
  disp(Mass_balnew)
  xlabel('MASS BALANCE');
  ylabel('ALTITUDE');
  legend('Albedo 0.3');
  
  figure('visible', 'on');
  scatter(Mass_balnew(2,:),ZM_elevation)
  xlabel('MASS BALANCE');
  ylabel('ALTITUDE');
  legend('Albedo 0.35');
  
  figure('visible', 'on');
  scatter(Mass_balnew(3,:),ZM_elevation)
  xlabel('MASS BALANCE');
  ylabel('ALTITUDE');
  legend('Albedo 0.4');
  figure('visible', 'on');
%}

  Runoffnew = Runoff;
  Runoffnew(isnan(Runoffnew))=0;
  rowMeans = mean(Runoffnew, 2); % Compute mean along the rows (dimension 2)

% Create a new matrix with the mean values
%{
  avgrunoff = repmat(rowMeans, 1, size(Runoffnew, 2)); % Replicate means to form a matrix

  disp('runoff, 1st columns is the average of all values, rows are albedo')
  disp(avgrunoff);
  
  %scatter(avgrunoff(:,1),albedo);
  %xlabel('Runoff');
  %ylabel('albedo');

  %scatter(Runoff_cumulative(p), albedo)
  
  %xlabel("Runoff cum")
  %ylabel('albedo')
%}

  stemp = Surfacetemp;
  stemp(isnan(stemp))=0;
  rowMeans = mean(stemp, 2); % Compute mean along the rows (dimension 2)

% Create a new matrix with the mean values

  avgstemp = repmat(rowMeans, 1, size(stemp, 2)); % Replicate means to form a matrix

  disp('surface temp, 1st columns is the average of all values, rows are albedo')
  disp(avgstemp);
  
scatter(albedo,avgstemp);
xlabel('Albedo');
xticks(0.3:0.05:0.8); % Set precise y-axis ticks
xlim([0.3, 0.8]); % Set y-axis limits
ylabel('Surface Temperature (°C)');
title('Surface Temperature vs. Albedo');
saveas(gcf, 'Surface_Temp_vs_Albedo.png');

  csvwrite('avgstemp.csv',Surfacetemp);

%-------------------------EXPORTING data------------------
%Outputting the Mass balance data which vary over the elevation for
%different values of albedo
filename = 'Mass_balance.txt';

% Write the matrix to the csv file
csvwrite(filename, Mass_bal);  
                
%Outputting the Mass balance data which vary over the elevation for
%different values of albedo
filename = 'Runoff.txt';

csvwrite(filename, Runoff_cumulative);

%-------------------------------------------------------------------------%
%Functions
%-------------------------------------------------------------------------%

%--------------------------------Functions--------------------------------%

%1. Calculation of Albedo from Surface sni density
function Acal= Albedo_cal(ssd) % ssd is surface snow density
    global di
    ri = 0.018; % reflectivity of snow
    k1 = 10;  %k is the absorption coefficient of ice (assumed to be 10 m-1)
    S = 10^((-15.32*(10^-9))*(ssd)^3 + (16.65*(10^-6))*(ssd)^2 - (7.3*(10^-3))*ssd + 2.23);
    li = 2/(S*di);%ice layer thickness
    Ri = ri + (((1-ri)^2)*ri*exp(-1*2*k1*li))/(1-(ri^2)*exp(-1*2*k1*li));
    Ti = (((1-ri)^2)*exp(-1*k1*li))/(1-(ri^2)*exp(-1*2*k1*li));
    tou = ((1-Ti)-sqrt((1-Ti)^2 - Ri^2))/Ri;
    Acal = ri + (((1-ri)^2)*tou)/(1-ri*tou);
    
end

%2. Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg)
% Returns the value of surface temperature in degree Celcius
function Ts = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg)
    global C ca le B
    Ts = Ta + ((1-A)*Rs + e*Rl - e* B *(Ta + 273.2)^4 - le*da*C*U*(1-rh)*q*Ta + Hg)/(4*e*B*(Ta + 273.2)^3 + (G_q*le + ca)*da*C*U);
    if Ts>0
        Ts = 0;
    end
end

%3. Melted_ice(Ts, A, Rs, Rl, e, da, U, Ta, rh, q, Hg)
% Return the value of ice that is melted in the glacier in mm.w.e.
function mi = Melted_ice(Ts, A, Rs, Rl, e, da, U, Ta, rh, q, Hg)
   Rn = Net_Radiation(A, Rs, Rl, e, Ts);
   Hs = Sensible_heat(da, U, Ta, Ts);% see the function comments for actual parameters
   Hl = Latent_Heat_Flux(da, U, rh, q, Ta, Ts);
   Hm = Rn+(Hs+Hl)*cosd(13)+Hg; 
   L = 3.36e+5;  %latent heat of fusion of ice (J/kg)
   mi = Hm/(L);
   if mi<0
       mi = 0;
   end
end

%4. Heat_Into_Glacier(mi, miprev, Ts, t, ds, xsnow, xice)
% Returns the value of Hg = heat into the glacier in W/m^2
function Hg = Heat_Into_Glacier(mi, Ts, t, ds, xsnow, xice)
    global bs ci T_iter
    if all(mi == 0) || all(Ts<0) % if no ice is melted ie, when surface temperature is negative
        T_iter = Ts;
        %new line above______________
        [temp_s, t] = Snow_Temperature();
        [temp_ice, t] = Ice_Temperature();
        snow_trap = trapz(temp_s(end,:)); %integration of the temperature variation of snow with limits of depth of snowlayer- from snow-ice interface to glacier surface
        ice_trap = trapz(temp_ice(end,:)); %integration of the temperature variation of snow with limits of depth of snowlayer- from critical depth - zc(from Fujita paper) to snow-ice interface
        Hg = -1*ci*((ice_trap + snow_trap))/t(end) ;
    elseif all(mi == 0)|| all(Ts == 0) % if no ice is melted but surface temperature = 0 degree celcius
        Hg = 0;
    elseif all(Ts<0) && all(mi>0)
        Ks = Snow_Thermal_Conductivity(ds);
        Hg = -Ks*Ts/bs;
    end
 
end

%5. Net_Radiation(A, Rs, Rl, e, Ts)
% Returns net radiation (Rn) in W/m^2 to the function - Melted_ice
function Rn = Net_Radiation(A, Rs, Rl, e, Ts)
    global B
    Rn = ((1-A)*Rs*cosd(13))+e*Rl - e*B*(Ts + 273.2).^4;
end


% 6. Sublimation (Hl, lv)
% Returns the Sublimation or re-sublimation in mm.w.e
function sub = Sublimation(da, U, rh, q, Ta, Ts,lv )
    %global lv
    Hl = Latent_Heat_Flux(da, U, rh, q, Ta, Ts);
    sub = Hl/lv;
end


%7. Sensible_heat(da, U, Ta, Ts)
% Returns Sensible heat (Hs) in W/m^2 to the function - Melted_ice
function Hs = Sensible_heat(da, U, Ta, Ts)
    global ca C
    Hs = ca*da*C*U*(Ta-Ts);
end

%8.Latent_Heat_Flux(da, U, rh, q, Ta, Ts)
% Returns Latent heat (Hl) in W/m^2 to the function - Melted_ice
function Hl = Latent_Heat_Flux(da, U, rh, q, Ta, Ts)
    global le C
    Hl = le*da*C*U*((rh*q*Ta)-(q*Ts));
end

%9. Superimposed_ice(ds, xice)
% Returns Fs = amount of superimposed ice on ice layer below snow layer in mm.w.e
function Fs = Superimposed_ice(ds, xice)
    global di ci lf
    [temp_ice, t] = Ice_Temperature();
    ice_trap = trapz(temp_ice(end,:)); %integration of the temperature variation of ice with limits of depth of icelayer- from critical depth - zc(from Fujita paper) to snow-ice interface
    Fs = di * ci *(ice_trap)/lf ; 
end

% 10. Mass_of_Refrozen_Capillary_Water(ds, Ts)
% Returns Fc_capillary =  Mass of capillary water that is refrozen within wet snow(mm.w.e)
function Fc_capillary = Mass_of_Refrozen_Capillary_Water(ds, Ts)
global lf bs    
    Ks = Snow_Thermal_Conductivity(ds); % Themal conductivity of snow in W/mK
    Fc_capillary = (-1*Ks*Ts)/(lf*bs);
end

%11. Mass_of_Snow_Percolated_Water(xsnow, t, ds)
% Returns Fc = Mass of water that is refrozen within snow layer (mm.w.e)
function Fc = Mass_of_Snow_Percolated_Water(xsnow, t, ds)
    global ci lf    
    [temp_s, ~] = Snow_Temperature();%should use for this function- xice and t
    snow_trap = trapz(temp_s(end,:));%integration of the temperature variation of snow with limits of depth of snowlayer- from snow-ice interface to glacier surface
    Fc = ds*ci*(snow_trap)/lf; 
end


%12. Snow_Temperature(ds, xsnow, t, Ts)
% Returns the variation of temperature of snow layer with time and depth
function [temp_s,t] = Snow_Temperature()
    m = 0;    
    x = linspace(0, 0.75, 50); %snow thickness
    t = linspace(0, 86400, 3600);

    sol_snow = pdepe(m, @temp_var_snow, @snow_ic, @snow_bc, x, t); % partial differential equation

    temp_s = sol_snow(:,:,1);
end

%13. temp_var_snow(x,t,u,dudx) - heat variation in snow
function [c,f,s] = temp_var_snow(x,t,u,dudx)
    global ds ci A Rs mus
    c = 840000;
    f = 0.493*dudx;
    s = 42*exp(-40*x);
end

%14. snow_ic() - initial conditions of snow
% Returns initial conditions of the differential question to the function
% Snow_Temperature(ds, xsnow, t, Ts) for pde
function u0 = snow_ic(x)
    u0 = 0;
end

%15. snow_bc(xl,ul,xr,ur,t,Ts) - boundary conditions of snow
% Returns left and right boundary conditions of the differential question
% to the function Snow_Temperature(ds, xsnow, t, Ts) for pde
function [pl,ql,pr,qr] = snow_bc(~,ul,~,ur,~)
    global T_iter
    pl = ul ; % write (ul - Ts0) here
    ql = 0;
    pr = ur + 5; %-5 C from Fujita Paper
    qr = 0;
end

%16. Ice_Temperature(xice, t)
% Returns the variation of temperature of ice layer with time and depth
% upto critical depth (zc)
function [temp_ice, t] = Ice_Temperature()
    m = 0;
    x = linspace(0, 19.25, 1300); %ice thickness
    t = linspace(0, 86400, 3600);

    sol = pdepe(m, @temp_var_ice, @ice_ic, @ice_bc, x, t);% partial differential equation

    temp_ice = sol(:,:,1);
end


%17. temp_var_ice(x,t,u,dudx) - heat variation in ice
function [c,f,s] = temp_var_ice(x,t,u,dudx)
    global di ci
    c = 1890000;
    f = log(273.2+u)*dudx;
    s = exp(-10*x);
end

%18. ice_ic() - initial conditions of ice
% Returns initial conditions of the differential question to the function
% Ice_Temperature(ds, xice, t, Ts) for pde
function u0 = ice_ic(x)
    u0 = -0.008*x^2-5;
end

%19. ice_bc(xl,ul,xr,ur,t) - boundary conditions of ice
% Returns left and right boundary conditions of the differential question
% to the function Ice_Temperature(ds, xsnow, t, Ts) for pde
function [pl,ql,pr,qr] = ice_bc(~,ul,~,ur,~)
    pl = ul + 5;
    ql = 0;
    pr = ur + 8;
    qr = 0;
end

%20 .Snow_Thermal_Conductivity(ds)
% Returns the thermal conductivity of snow in (W/mK)
function Ks = Snow_Thermal_Conductivity(ds)
    Ks = 0.029*(1+(0.001*ds*ds));
end

% 21. Precipitations(Pp, Ta)
% Returns the rainfall and snowfall in mm.w.e/d
% Pr = Rainfall, Ps = Snowfall
function [Pr, Ps] = Precipitations(Pp, Ta) %units- mm.w.e/d
    if Ta<=0
        Ps = Pp;
    elseif Ta<6
        Ps = Pp * (1-(Ta/6));
    else
        Ps = 0;
    end
    Pr = Pp - Ps;
end