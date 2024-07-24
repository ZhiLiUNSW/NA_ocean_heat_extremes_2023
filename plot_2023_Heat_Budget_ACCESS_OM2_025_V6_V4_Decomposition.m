%% Mixed-Layer Heat Budget Anomaly in 2023, ACCESS-OM2 0.25 JRA55-do IAF Forced Simulation

%% ########################################################################
%% 0. Decomposing SWR penatration into SWR_clim and MLD_clim

%% 0.1 Using monthly SWR and MLD
clc;clear
time_ann=(1981:2010)';

% #########################################################################
% Basin Mask from ACCESS OM2 0.25
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST_and_MLD/')
     lon_025=ncread('full_mixed_layer_heat_budget_year_2022.nc','xt_ocean'); % 100w-20E
     lon_025=lon_025+360;
     lat_025=ncread('full_mixed_layer_heat_budget_year_2022.nc','yt_ocean'); % 0-70N

     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     clear basin* Sxy*
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid(lon_025, lat_025);
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');% figure;imagesc(basin_mask_NA)
       clear lo la lon lat
       [Sxy,~,~]=function_Cgrid_Area_Distance(lon_025,lat_025);
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0 month
% #########################################################################
      

% #########################################################################
% #########################################################################
% MLT budget for 1981-2010
    % #####################################################################
    % Monthly MLD from ACCESS, 1981-2010
    disp(['   MLD from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     lon=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_1980-2010_clim.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_1980-2010_clim.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','mld');
     mld=mld0(721:1200,499:891,:); clear mld0

    % #####################################################################
    % Monthly SWR from ACCESS, 1981-2010
    disp(['   SWR from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_surf');
     SW_surf=SW_surf0(721:1200,499:891,:); clear SW_surf0
     SW_surf(SW_surf==0)=NaN;
     SW_base0=ncread('SW_base_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_base');
     SW_base=SW_base0(721:1200,499:891,:); clear SW_base0
     SW_base(SW_base==0)=NaN;
     % Remove the shortwave penatration
     SW_net=SW_surf+SW_base; clear SW_surf SW_base
     
     
     count_yr=1;
     for year=1981:2010
         disp(['   MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         SWR_monthly=SW_net(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);

         % ##################################################################
         % Mixed layer warming due to Qnet
            rho0 = 1035.0;
            Cp = 3992.10322329649;
            MLT_Qswr_mon(:,:,:,count_yr)  = SWR_monthly./mld_monthly./Cp./rho0.*(second_2_month); % K/month
            clear SWR_monthly
            clear mld_monthly
            % figure;imagesc(MLT_Qswr_mon(:,:,5)); caxis([-3 3])

             % ################################################################
             % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
             dMLT_Qswr_mon0=MLT_Qswr_mon(:,1:302,:,count_yr).*Sxy_NA(:,1:302);
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA(:,:,1)); 
             Sxy_NA0=Sxy_NA(:,1:302,1);
             for month=2:12
                Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
             end
             clear month
             Sxy_NA0(isnan(dMLT_Qswr_mon0))=NaN; 
             dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
             % figure;imagesc(dMLT_Qswr_mon0(:,:,1)); 
             % figure;imagesc(Sxy_NA0(:,:,1)); 

             dMLT_Qswr_mon(1:12,count_yr)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
             clear dMLT_*_mon0 Sxy_NA0
          % ###################################################################
          count_yr=count_yr+1;
     end
     clear year mld SW_net count_yr
     
     MLT_Qswr_clim=squeeze(nanmean(MLT_Qswr_mon(:,:,:,1:30),4));
     clear MLT_Qswr_mon
     dMLT_Qswr_clim=squeeze(nanmean(dMLT_Qswr_mon(:,1:30),2));
     clear dMLT_Qswr_mon
% #########################################################################
% #########################################################################  



% #########################################################################
% #########################################################################
% MLT budget for 2023
    % #####################################################################
    % Monthly MLD from ACCESS, 2023
    disp(['   MLD from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     lon=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','mld');
     mld(:,:,5:8)=mld0(721:1200,499:891,1:4); clear mld0
     mld(:,:,1:4)=NaN; mld(:,:,9:12)=NaN; 

    % #####################################################################
    % Monthly SWR from ACCESS, 2023
    disp(['   SWR from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','SW_surf');
     SW_surf(:,:,5:8)=SW_surf0(721:1200,499:891,1:4); clear SW_surf0
     SW_surf(:,:,1:4)=NaN;SW_surf(:,:,9:12)=NaN;SW_surf(SW_surf==0)=NaN;
     SW_base0=ncread('SW_base_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','SW_base');
     SW_base(:,:,5:8)=SW_base0(721:1200,499:891,1:4); clear SW_base0
     SW_base(:,:,1:4)=NaN;SW_base(:,:,9:12)=NaN;SW_base(SW_base==0)=NaN;
     % Remove the shortwave penatration
     SW_net=SW_surf+SW_base; clear SW_surf SW_base
     
     
     count_yr=1;
     for year=2023
         disp([' MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         SWR_monthly=SW_net(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);

         % ##################################################################
         % Mixed layer warming due to Qnet
            rho0 = 1035.0;
            Cp = 3992.10322329649;
            MLT_Qswr_mon(:,:,:,count_yr)  = SWR_monthly./mld_monthly./Cp./rho0.*(second_2_month); % K/month
            clear SWR_monthly
            clear mld_monthly
            % figure;imagesc(MLT_Qswr_mon(:,:,5)); caxis([-3 3])

             % ################################################################
             % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
             dMLT_Qswr_mon0=MLT_Qswr_mon(:,1:302,:,count_yr).*Sxy_NA(:,1:302);
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA(:,:,1)); 
             Sxy_NA0=Sxy_NA(:,1:302,1);
             for month=2:12
                Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
             end
             clear month
             Sxy_NA0(isnan(dMLT_Qswr_mon0))=NaN; 
             dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA0(:,:,1)); 

             dMLT_Qswr_mon(1:12,count_yr)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
             clear dMLT_*_mon0 Sxy_NA0
          % ###################################################################
          count_yr=count_yr+1;
     end
     clear year mld SW_net count_yr
     
     MLT_Qswr_2023=squeeze(nanmean(MLT_Qswr_mon(:,:,:,1),4));
     clear MLT_Qswr_mon
     dMLT_Qswr_2023=squeeze(nanmean(dMLT_Qswr_mon(:,1),2));
     clear dMLT_Qswr_mon
% #########################################################################
% #########################################################################  

MLT_Qswr_2023_ano  = MLT_Qswr_2023  - MLT_Qswr_clim; 
dMLT_Qswr_2023_ano = dMLT_Qswr_2023 - dMLT_Qswr_clim; 

cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
save('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_1.mat',...
     'MLT_Qswr_2023','dMLT_Qswr_2023','MLT_Qswr_clim','dMLT_Qswr_clim',...
     'MLT_Qswr_2023_ano','dMLT_Qswr_2023_ano','lon_025','lat_025')



%% 0.2 Using monthly SWR and MLD_clim - V2: MLD Climatology
clc;clear
time_ann=(1981:2010)';

% #########################################################################
% Basin Mask from ACCESS OM2 0.25
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST_and_MLD/')
     lon_025=ncread('full_mixed_layer_heat_budget_year_2022.nc','xt_ocean'); % 100w-20E
     lon_025=lon_025+360;
     lat_025=ncread('full_mixed_layer_heat_budget_year_2022.nc','yt_ocean'); % 0-70N

     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     clear basin* Sxy*
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid(lon_025, lat_025);
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');% figure;imagesc(basin_mask_NA)
       clear lo la lon lat
       [Sxy,~,~]=function_Cgrid_Area_Distance(lon_025,lat_025);
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0 month
% #########################################################################
      

% #########################################################################
% #########################################################################
% MLT budget for 1981-2010
    % #####################################################################
    % Monthly MLD from ACCESS, 1981-2010
    disp(['   MLD from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     lon=ncread('MLD_discrete_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_discrete_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_discrete_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','MLD_discrete');
     mld=mld0(721:1200,499:891,:); clear mld0
     
     count_yr=1;
     for year=1981:2010
         mld0(:,:,1:12,count_yr)=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         count_yr=count_yr+1;
     end
     clear mld
     mld=nanmean(mld0,4); clear mld0
     

    % #####################################################################
    % Monthly SWR from ACCESS, 1981-2010
    disp(['   SWR from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_surf');
     SW_surf=SW_surf0(721:1200,499:891,:); clear SW_surf0
     SW_surf(SW_surf==0)=NaN;
     SW_base0=ncread('SW_base_clim_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_base_clim');
     SW_base=SW_base0(721:1200,499:891,:); clear SW_base0
     SW_base(SW_base==0)=NaN;
     % Remove the shortwave penatration
     SW_net=SW_surf+SW_base; clear SW_surf SW_base
     
     
     count_yr=1;
     for year=1981:2010
         disp(['   MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,1:12);
         SWR_monthly=SW_net(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);

         % ##################################################################
         % Mixed layer warming due to Qnet
            rho0 = 1035.0;
            Cp = 3992.10322329649;
            MLT_Qswr_mon(:,:,:,count_yr)  = SWR_monthly./mld_monthly./Cp./rho0.*(second_2_month); % K/month
            clear SWR_monthly
            clear mld_monthly
            % figure;imagesc(MLT_Qswr_mon(:,:,5)); caxis([-3 3])

             % ################################################################
             % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
             dMLT_Qswr_mon0=MLT_Qswr_mon(:,1:302,:,count_yr).*Sxy_NA(:,1:302);
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA(:,:,1)); 
             Sxy_NA0=Sxy_NA(:,1:302,1);
             for month=2:12
                Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
             end
             clear month
             Sxy_NA0(isnan(dMLT_Qswr_mon0))=NaN; 
             dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
             % figure;imagesc(dMLT_Qswr_mon0(:,:,1)); 
             % figure;imagesc(Sxy_NA0(:,:,1)); 

             dMLT_Qswr_mon(1:12,count_yr)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
             clear dMLT_*_mon0 Sxy_NA0
          % ###################################################################
          count_yr=count_yr+1;
     end
     clear year mld SW_net count_yr
     
     MLT_Qswr_clim=squeeze(nanmean(MLT_Qswr_mon(:,:,:,1:30),4));
     clear MLT_Qswr_mon
     dMLT_Qswr_clim=squeeze(nanmean(dMLT_Qswr_mon(:,1:30),2));
     clear dMLT_Qswr_mon
% #########################################################################
% #########################################################################  



% #########################################################################
% #########################################################################
% MLT budget for 2023
    % #####################################################################
    % Monthly MLD from ACCESS, 1981-2010
    disp(['   MLD from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     lon=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
%      mld0=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','mld');
%      mld(:,:,5:8)=mld0(721:1200,499:891,1:4); clear mld0
%      mld(:,:,1:4)=NaN; mld(:,:,9:12)=NaN; 
     mld0=ncread('MLD_discrete_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','MLD_discrete');
     mld=mld0(721:1200,499:891,:); clear mld0
     
     count_yr=1;
     for year=1981:2010
         mld0(:,:,1:12,count_yr)=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         count_yr=count_yr+1;
     end
     clear mld
     mld=nanmean(mld0,4); clear mld0
     

    % #####################################################################
    % Monthly SWR from ACCESS, 1981-2010
    disp(['   SWR from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','SW_surf');
     SW_surf(:,:,5:8)=SW_surf0(721:1200,499:891,1:4); clear SW_surf0
     SW_surf(:,:,1:4)=NaN;SW_surf(:,:,9:12)=NaN;SW_surf(SW_surf==0)=NaN;
     SW_base0=ncread('SW_base_clim_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','SW_base_clim');
     SW_base(:,:,5:8)=SW_base0(721:1200,499:891,1:4); clear SW_base0
     SW_base(:,:,1:4)=NaN;SW_base(:,:,9:12)=NaN;SW_base(SW_base==0)=NaN;
     % Remove the shortwave penatration
     SW_net=SW_surf+SW_base; clear SW_surf SW_base
     
     
     count_yr=1;
     for year=2023
         disp([' MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,1:12);
         SWR_monthly=SW_net(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);

         % ##################################################################
         % Mixed layer warming due to Qnet
            rho0 = 1035.0;
            Cp = 3992.10322329649;
            MLT_Qswr_mon(:,:,:,count_yr)  = SWR_monthly./mld_monthly./Cp./rho0.*(second_2_month); % K/month
            clear SWR_monthly
            clear mld_monthly
            % figure;imagesc(MLT_Qswr_mon(:,:,5)); caxis([-3 3])

             % ################################################################
             % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
             dMLT_Qswr_mon0=MLT_Qswr_mon(:,1:302,:,count_yr).*Sxy_NA(:,1:302);
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA(:,:,1)); 
             Sxy_NA0=Sxy_NA(:,1:302,1);
             for month=2:12
                Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
             end
             clear month
             Sxy_NA0(isnan(dMLT_Qswr_mon0))=NaN; 
             dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA0(:,:,1)); 

             dMLT_Qswr_mon(1:12,count_yr)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
             clear dMLT_*_mon0 Sxy_NA0
          % ###################################################################
          count_yr=count_yr+1;
     end
     clear year mld SW_net count_yr
     
     MLT_Qswr_2023=squeeze(nanmean(MLT_Qswr_mon(:,:,:,1),4));
     clear MLT_Qswr_mon
     dMLT_Qswr_2023=squeeze(nanmean(dMLT_Qswr_mon(:,1),2));
     clear dMLT_Qswr_mon
% #########################################################################
% #########################################################################  

MLT_Qswr_2023_ano  = MLT_Qswr_2023  - MLT_Qswr_clim; 
dMLT_Qswr_2023_ano = dMLT_Qswr_2023 - dMLT_Qswr_clim; 

cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
save('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_2_MLD_clim.mat',...
     'MLT_Qswr_2023','dMLT_Qswr_2023','MLT_Qswr_clim','dMLT_Qswr_clim',...
     'MLT_Qswr_2023_ano','dMLT_Qswr_2023_ano','lon_025','lat_025')




%% 0.3 Using monthly SWR_clim and MLD - V3: Qnet climatology
clc;clear
time_ann=(1981:2010)';

% #########################################################################
% Basin Mask from ACCESS OM2 0.25
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST_and_MLD/')
     lon_025=ncread('full_mixed_layer_heat_budget_year_2022.nc','xt_ocean'); % 100w-20E
     lon_025=lon_025+360;
     lat_025=ncread('full_mixed_layer_heat_budget_year_2022.nc','yt_ocean'); % 0-70N

     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     clear basin* Sxy*
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid(lon_025, lat_025);
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');% figure;imagesc(basin_mask_NA)
       clear lo la lon lat
       [Sxy,~,~]=function_Cgrid_Area_Distance(lon_025,lat_025);
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0 month
% #########################################################################
      

% #########################################################################
% #########################################################################
% MLT budget for 1981-2010
    % #####################################################################
    % Monthly MLD from ACCESS, 1981-2010
    disp(['   MLD from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     lon=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_1980-2010_clim.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_1980-2010_clim.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','mld');
     mld=mld0(721:1200,499:891,:); clear mld0

    % #####################################################################
    % Monthly SWR from ACCESS, 1981-2010
    disp(['   SWR from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_surf');
     SW_surf=SW_surf0(721:1200,499:891,:); clear SW_surf0
     SW_surf(SW_surf==0)=NaN;
     SW_base0=ncread('SW_base_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_base');
     SW_base=SW_base0(721:1200,499:891,:); clear SW_base0
     SW_base(SW_base==0)=NaN;
     % Remove the shortwave penatration
     SW_net=SW_surf+SW_base; clear SW_surf SW_base
     
     count_yr=1;
     for year=1981:2010
         SW_net0(:,:,1:12,count_yr)=SW_net(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         count_yr=count_yr+1;
     end
     clear SW_net
     SW_net=nanmean(SW_net0,4); clear SW_net0
     
     
     count_yr=1;
     for year=1981:2010
         disp(['   MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         SWR_monthly=SW_net(:,:,1:12);

         % ##################################################################
         % Mixed layer warming due to Qnet
            rho0 = 1035.0;
            Cp = 3992.10322329649;
            MLT_Qswr_mon(:,:,:,count_yr)  = SWR_monthly./mld_monthly./Cp./rho0.*(second_2_month); % K/month
            clear SWR_monthly
            clear mld_monthly
            % figure;imagesc(MLT_Qswr_mon(:,:,5)); caxis([-3 3])

             % ################################################################
             % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
             dMLT_Qswr_mon0=MLT_Qswr_mon(:,1:302,:,count_yr).*Sxy_NA(:,1:302);
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA(:,:,1)); 
             Sxy_NA0=Sxy_NA(:,1:302,1);
             for month=2:12
                Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
             end
             clear month
             Sxy_NA0(isnan(dMLT_Qswr_mon0))=NaN; 
             dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
             % figure;imagesc(dMLT_Qswr_mon0(:,:,1)); 
             % figure;imagesc(Sxy_NA0(:,:,1)); 

             dMLT_Qswr_mon(1:12,count_yr)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
             clear dMLT_*_mon0 Sxy_NA0
          % ###################################################################
          count_yr=count_yr+1;
     end
     clear year mld SW_net count_yr
     
     MLT_Qswr_clim=squeeze(nanmean(MLT_Qswr_mon(:,:,:,1:30),4));
     clear MLT_Qswr_mon
     dMLT_Qswr_clim=squeeze(nanmean(dMLT_Qswr_mon(:,1:30),2));
     clear dMLT_Qswr_mon
% #########################################################################
% #########################################################################  



% #########################################################################
% #########################################################################
% MLT budget for 2023
    % #####################################################################
    % Monthly MLD from ACCESS, 2023
    disp(['   MLD from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     lon=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','mld');
     mld(:,:,5:8)=mld0(721:1200,499:891,1:4); clear mld0
     mld(:,:,1:4)=NaN; mld(:,:,9:12)=NaN; 

    % #####################################################################
%     % Monthly SWR from ACCESS, 2023
%     disp(['   SWR from IAPv3, Year#1981-2010'])
%      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
%      SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','SW_surf');
%      SW_surf(:,:,5:8)=SW_surf0(721:1200,499:891,1:4); clear SW_surf0
%      SW_surf(:,:,1:4)=NaN;SW_surf(:,:,9:12)=NaN;SW_surf(SW_surf==0)=NaN;
%      SW_base0=ncread('SW_base_025deg_jra55_iaf_omip2_cycle6_jra55v150_extension_May-Aug2023.nc','SW_base');
%      SW_base(:,:,5:8)=SW_base0(721:1200,499:891,1:4); clear SW_base0
%      SW_base(:,:,1:4)=NaN;SW_base(:,:,9:12)=NaN;SW_base(SW_base==0)=NaN;
%      % Remove the shortwave penatration
%      SW_net=SW_surf+SW_base; clear SW_surf SW_base  

    % Monthly SWR from ACCESS, 1981-2010
    disp(['   SWR from IAPv3, Year#1981-2010'])
     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_SWR_MLD_decomposition')
     SW_surf0=ncread('SW_surf_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_surf');
     SW_surf=SW_surf0(721:1200,499:891,:); clear SW_surf0
     SW_surf(SW_surf==0)=NaN;
     SW_base0=ncread('SW_base_025deg_jra55_iaf_omip2_cycle6_1980-2010.nc','SW_base');
     SW_base=SW_base0(721:1200,499:891,:); clear SW_base0
     SW_base(SW_base==0)=NaN;
     % Remove the shortwave penatration
     SW_net=SW_surf+SW_base; clear SW_surf SW_base
     
     count_yr=1;
     for year=1981:2010
         SW_net0(:,:,1:12,count_yr)=SW_net(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         count_yr=count_yr+1;
     end
     clear SW_net
     SW_net=nanmean(SW_net0,4); clear SW_net0
     
     
     count_yr=1;
     for year=2023
         disp([' MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         SWR_monthly=SW_net(:,:,1:12);

         % ##################################################################
         % Mixed layer warming due to Qnet
            rho0 = 1035.0;
            Cp = 3992.10322329649;
            MLT_Qswr_mon(:,:,:,count_yr)  = SWR_monthly./mld_monthly./Cp./rho0.*(second_2_month); % K/month
            clear SWR_monthly
            clear mld_monthly
            % figure;imagesc(MLT_Qswr_mon(:,:,5)); caxis([-3 3])

             % ################################################################
             % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
             dMLT_Qswr_mon0=MLT_Qswr_mon(:,1:302,:,count_yr).*Sxy_NA(:,1:302);
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA(:,:,1)); 
             Sxy_NA0=Sxy_NA(:,1:302,1);
             for month=2:12
                Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
             end
             clear month
             Sxy_NA0(isnan(dMLT_Qswr_mon0))=NaN; 
             dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
             % figure;imagesc(dMLT_Qswr_mon0(:,:,5)); 
             % figure;imagesc(Sxy_NA0(:,:,1)); 

             dMLT_Qswr_mon(1:12,count_yr)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
             clear dMLT_*_mon0 Sxy_NA0
          % ###################################################################
          count_yr=count_yr+1;
     end
     clear year mld SW_net count_yr
     
     MLT_Qswr_2023=squeeze(nanmean(MLT_Qswr_mon(:,:,:,1),4));
     clear MLT_Qswr_mon
     dMLT_Qswr_2023=squeeze(nanmean(dMLT_Qswr_mon(:,1),2));
     clear dMLT_Qswr_mon
% #########################################################################
% #########################################################################  

MLT_Qswr_2023_ano  = MLT_Qswr_2023  - MLT_Qswr_clim; 
dMLT_Qswr_2023_ano = dMLT_Qswr_2023 - dMLT_Qswr_clim; 

cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
save('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_3_SWR_ciim.mat',...
     'MLT_Qswr_2023','dMLT_Qswr_2023','MLT_Qswr_clim','dMLT_Qswr_clim',...
     'MLT_Qswr_2023_ano','dMLT_Qswr_2023_ano','lon_025','lat_025')




%% ########################################################################
%% 1.1 Plotting: ML Temperature Budget Anomaly in May-August 2023
%  Bar Charts and Maps ####################################################
% #########################################################################
% #########################################################################
% V1
clc; clear
figure('Color',[1 1 1]);  %create a new figure of white color background
% ixs = 0.060; ixe = 0.100;  ixd = 0.045; ixw = (1-ixs-ixe-2*ixd)/3;
% iys = 0.100; iye = 0.100;  iyd = 0.050; iyw = (1-iys-iye-1*iyd)/2;

ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.050; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
ixw./iyw

%          [left            bottom      width height]
pos{11}  = [ixs+0*ixw+0*ixd    iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{12}  = [ixs+1*ixw+1*ixd    iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+2*ixw+2*ixd    iys+2*iyw+2*iyd   ixw 1.0*iyw];

pos{21}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


% #########################################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);
% color0(8,:)=(color0(7,:)+color0(8,:))./2;   
color0(6,:)=(color0(6,:)+color0(5,:))./2;   
% #########################################################################


% #########################################################################
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1981-2010
      count=1;
      for year=1981:2010
          disp(['   Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;


          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Shortwave penetration at MLB
          Qswr(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;


          % (f) Total mixing
          resd0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_submeso');
          resd0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_k33');
          resd0(:,:,:,3)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_diffusion_temp');
          resd0(:,:,:,4)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_gm_temp');
          resd0(:,:,:,5)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'residual');% residual=mixdownslope_temp + temp_sigma_diff  
          resd0(:,:,:,6)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_eta_smooth');
          
          resd0(:,:,:,7)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          resd0(:,:,:,8)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          
          resd0(:,:,:,9)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_rivermix');
          resd0(:,:,:,10)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sfc_hflux_pme');
          resd0(:,:,:,11)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'frazil_3d');
          
          resd(:,:,:,count)=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
          resd(resd==0)=NaN;
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      Qswr_clim=nanmean(Qswr,4); clear Qswr
      resd_clim=nanmean(resd,4); clear resd
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month;
      
      % (b) Entrainment
      entrainment=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month;
      
      % (c) Total advection
      adve0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
      Qnet(Qnet==0)=NaN;
      
      
      % (e) Shortwave penetration at MLB
      Qswr(:,:,:)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sw_heat').*second_2_month;


      % (f) Total mixing
      resd0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_submeso');
      resd0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_k33');
      resd0(:,:,:,3)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_diffusion_temp');
      resd0(:,:,:,4)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_gm_temp');
      resd0(:,:,:,5)=ncread('full_mixed_layer_heat_budget_year_2023.nc','residual');% residual=mixdownslope_temp + temp_sigma_diff  
      resd0(:,:,:,6)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_eta_smooth');
      
      resd0(:,:,:,7)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      resd0(:,:,:,8)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      
      resd0(:,:,:,9)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_rivermix');
      resd0(:,:,:,10)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sfc_hflux_pme');
      resd0(:,:,:,11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','frazil_3d');
      
      resd=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
      resd(resd==0)=NaN;
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      Qswr=Qswr-Qswr_clim; 
      resd=resd-resd_clim; 
      % ###################################################################

     

% #########################################################################       
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('a. MLT tendency (June 2023)','fontsize',20,'FontWeight','bold')
        
        
         
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,6)'+Qswr(:,:,6)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. Surface flux term','fontsize',20,'FontWeight','bold')

        

subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,6)'+resd(:,:,6)'+entrainment(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. Advection + mixing + entrainment','fontsize',19.5,'FontWeight','bold')
        
        
        


% #########################################################################       
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. MLT tendency (July 2023)','fontsize',20,'FontWeight','bold')
        
        
         
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,7)'+Qswr(:,:,7)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. Surface flux term','fontsize',20,'FontWeight','bold')

        

subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,7)'+resd(:,:,7)'+entrainment(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. Advection + mixing + entrainment','fontsize',19.5,'FontWeight','bold')

        
        hBar1 = colorbar('EastOutside','vertical');
        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+3*ixw+2*ixd+0.016 iys+1.5*iyw+1*iyd 0.012 1.0*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.060);
        % ylabel(hBar1, '[ ^oC/month ]','rotation',90);
        m_text(382,95, '[\circC/month]','fontsize',18,'FontWeight','normal')
% #########################################################################
% #########################################################################


    
% #########################################################################    
% #########################################################################
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1981-2010
      count=1;
      for year=1981:2010
          disp(['Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;

          
          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Vertical mixing
          mixz0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          mixz0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          mixz(:,:,:,count)=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
          mixz(mixz==0)=NaN;
          
          
          % ###############################################################
          % heat flux decomposition
          % (a) SW term 
          Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'swflx');
          % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
          Qshortwave(:,:,:,count)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
          Qshortwave(Qshortwave==0)=NaN;
          
          % (b) LW term
          Qlongwave(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'lw_heat').*second_2_month;

          % (c) Sensible term
          Qsensible(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sens_heat').*second_2_month;

          % (d) Latent term
          Qlatent(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'evap_heat').*second_2_month;
          
          % (e) Shortwave penetration at MLB
          Qswr_mld(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;
          % ###############################################################
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      mixz_clim=nanmean(mixz,4); clear mixz
      
      Qshortwave_clim=nanmean(Qshortwave,4); clear Qshortwave
      Qlongwave_clim=nanmean(Qlongwave,4);   clear Qlongwave
      Qsensible_clim=nanmean(Qsensible,4);   clear Qsensible
      Qlatent_clim=nanmean(Qlatent,4);       clear Qlatent
      Qswr_mld_clim=nanmean(Qswr_mld,4);     clear Qswr_mld
       
      % Removing the shortwave penetration from Qnet and SWR
      Qnet_clim=Qnet_clim+Qswr_mld_clim;
      Qshortwave_clim=Qshortwave_clim+Qswr_mld_clim;
      clear Qswr_mld_clim
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month;
      
      % (b) Entrainment
      entrainment=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month;
      
      % (c) Total advection
      adve0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
      Qnet(Qnet==0)=NaN;

      % (e) Vertical mixing
      mixz0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      mixz0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      mixz=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
      mixz(mixz==0)=NaN;
      
                
      % ###############################################################
      % heat flux decomposition
      % (a) SW term 
      Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'swflx');
      % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
      Qshortwave(:,:,:)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
      Qshortwave(Qshortwave==0)=NaN;

      % (b) LW term
      Qlongwave(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'lw_heat').*second_2_month;

      % (c) Sensible term
      Qsensible(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sens_heat').*second_2_month;

      % (d) Latent term
      Qlatent(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'evap_heat').*second_2_month;

      % (e) Shortwave penetration at MLB
      Qswr_mld(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sw_heat').*second_2_month;
      % ###############################################################
      % Removing the shortwave penetration from Qnet and SWR
      Qnet=Qnet+Qswr_mld;
      Qshortwave=Qshortwave+Qswr_mld;
      clear Qswr_mld
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      mixz=mixz-mixz_clim; 

      Qshortwave=Qshortwave-Qshortwave_clim; 
      Qlongwave=Qlongwave-Qlongwave_clim; 
      Qsensible=Qsensible-Qsensible_clim; 
      Qlatent=Qlatent-Qlatent_clim; 
      % ###################################################################

      
    % #########################################################################
      % Basin Mask from ACCESS OM2 0.25 forced by ERA5
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_NA_SST_fields/')
      load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
      [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
      Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
      for month=1:12
          basin_mask_NA0(:,:,month)=basin_mask_NA;
      end
      clear month
      
      dmlt(isnan(basin_mask_NA0))=NaN; 
      entrainment(isnan(basin_mask_NA0))=NaN;
      adve(isnan(basin_mask_NA0))=NaN; 
      Qnet(isnan(basin_mask_NA0))=NaN; 
      mixz(isnan(basin_mask_NA0))=NaN; 

      Qshortwave(isnan(basin_mask_NA0))=NaN; 
      Qlongwave(isnan(basin_mask_NA0))=NaN; 
      Qsensible(isnan(basin_mask_NA0))=NaN; 
      Qlatent(isnan(basin_mask_NA0))=NaN; 
      
      dmlt=dmlt.*Sxy;
      entrainment=entrainment.*Sxy;
      adve=adve.*Sxy;
      Qnet=Qnet.*Sxy;
      mixz=mixz.*Sxy;
      
      Qshortwave=Qshortwave.*Sxy;
      Qlongwave=Qlongwave.*Sxy;
      Qsensible=Qsensible.*Sxy;
      Qlatent=Qlatent.*Sxy;
          
      % NA 0-60N
      dmlt_2023(1:12,1)=squeeze(nansum(nansum(dmlt(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      entrainment_2023(1:12,1)=squeeze(nansum(nansum(entrainment(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      adve_2023(1:12,1)=squeeze(nansum(nansum(adve(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qnet_2023(1:12,1)=squeeze(nansum(nansum(Qnet(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      mixz_2023(1:12,1)=squeeze(nansum(nansum(mixz(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      
      Qshortwave_2023(1:12,1)=squeeze(nansum(nansum(Qshortwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlongwave_2023(1:12,1)=squeeze(nansum(nansum(Qlongwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qsensible_2023(1:12,1)=squeeze(nansum(nansum(Qsensible(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlatent_2023(1:12,1)=squeeze(nansum(nansum(Qlatent(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
    % #########################################################################
      
   

ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.050; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
ixw./iyw

%          [left            bottom      width height]
pos{3}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   3*ixw+2*ixd 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);


clear bar*
bar_dmlt(:,1)=dmlt_2023(1:12,1);       % MLT tendency
bar_dmlt(:,2)=Qnet_2023(1:12,1);       % Qnet
bar_dmlt(:,3)=Qshortwave_2023(1:12,1); % Qnet - SWR
bar_dmlt(:,4)=Qlatent_2023(1:12,1);    % Qnet - Latent
bar_dmlt(:,5)=Qlongwave_2023(1:12,1);  % Qnet - LWR
bar_dmlt(:,6)=Qsensible_2023(1:12,1);  % Qnet - Sensible
bar_dmlt(:,7)=entrainment_2023(1:12,1)+mixz_2023(1:12,1);% Entrainment + Vertical mixing
bar_dmlt(:,8)=adve_2023+(dmlt_2023-Qnet_2023-entrainment_2023-mixz_2023-adve_2023); % Temperature advection + Residual


% SWR decomposition
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_1.mat',...
     'dMLT_Qswr_2023_ano')
bar_dmlt(:,3)=dMLT_Qswr_2023_ano(1:12,1); % Qnet - SWR
clear dMLT_Qswr_2023_ano

subplot('position',pos{3})
    h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.95); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.366, 0.574, 0.188],'EdgeColor',[0.366, 0.574, 0.188])
       hold on
       set(h0(8),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

       
       % ##################################################################
       % Add SWR decomposition
        hold on
        bar_dmlt_Clim(1:9,1:16)=NaN;
        % SWR cimatology
        cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_3_SWR_ciim.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,5)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano
                
        % MLD cimatology
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_2_MLD_clim.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,6)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano
        
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.36); 
           hold on
           set(h0_Clim(5),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(6),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           
           hold on
           plot((6.00:0.001:6.13),1.3*ones(length((5.05:0.001:5.18)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.00:0.001:6.13),1.05*ones(length((5.05:0.001:5.18)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.00:0.001:6.13),1.3*ones(length((5.05:0.001:5.18)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.00:0.001:6.13),1.05*ones(length((5.05:0.001:5.18)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)

           text(6.15, 1.3,'$\mathbf{\overline{Q_{sw}} / MLD^\prime}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
           text(6.15, 1.05,'$\mathbf{{Q_{sw}}^\prime / \overline{MLD}}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
       % ##################################################################



    hold on
%     % set 3 display names for the 3 handles
%     set(h0, {'DisplayName'}, {'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Vertical mixing + Entrainment','Advection + Residual'}')
       
    legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7) h0(8)],...
           'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Vertical mixing + entrainment','Advection + other minor terms',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',17)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
   
    % Legend will show names for each color
    legend() 
    set(gca,'Ylim',[-0.4 1.63],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.6)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2','1.6'},'fontsize',20)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:1:8.5)
    set(gca,'XTickLabel',{[],[],[],[]},'fontsize',20)
%     set(gca,'XTick',5:1:8)
%     set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',20)
    text(4.915,-0.54,'May','fontsize',20,'color','k','FontWeight','normal')
    text(5.930,-0.54,'Jun','fontsize',20,'color','k','FontWeight','normal')
    text(6.945,-0.54,'Jul','fontsize',20,'color','k','FontWeight','normal')
    text(7.920,-0.54,'Aug','fontsize',20,'color','k','FontWeight','normal')


    grid on
%     ax = gca;
%     ax.XGrid = 'off';
%     ax.YGrid = 'on';
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ \circC/month ]'],'fontsize',20,'color','k','FontWeight','normal')

    text(4.6,1.5,'g. MLT budget anomalies','fontsize',20,'color','k','FontWeight','bold')
    
    

% #########################################################################
% % V2 for panel g
% % Reordering bar charts following on Alex's comments:
% % I wonder if the bars could be reordered. Start with the bars that add up to give the MLT tendancy, then afterwards show the breakdown of the heat flux terms.
% % I would also get rid of the vertical lines centered on May, june, July August
% clear bar*
% bar_dmlt(:,1)=dmlt_2023(1:12,1);       % MLT tendency
% bar_dmlt(:,2)=Qnet_2023(1:12,1);       % Qnet
% bar_dmlt(:,3)=entrainment_2023(1:12,1)+mixz_2023(1:12,1);% Entrainment + Vertical mixing
% bar_dmlt(:,4)=adve_2023+(dmlt_2023-Qnet_2023-entrainment_2023-mixz_2023-adve_2023); % Temperature advection + Residual
% bar_dmlt(1:12,5)=NaN;
% bar_dmlt(1:12,6)=NaN;
% bar_dmlt(1:12,7)=NaN;
% bar_dmlt(1:12,8)=NaN;
% 
% bar_dmlt1(1:12,1)=NaN;
% bar_dmlt1(1:12,2)=Qnet_2023(1:12,1);       % Qnet
% 
% 
% bar_dmlt2(1:12,1)=NaN;
% bar_dmlt2(1:12,2)=NaN;
% bar_dmlt2(1:12,3)=NaN;
% bar_dmlt2(1:12,4)=NaN;
% bar_dmlt2(:,5)=Qshortwave_2023(1:12,1); % Qnet - SWR
% bar_dmlt2(:,6)=Qlatent_2023(1:12,1);    % Qnet - Latent
% bar_dmlt2(:,7)=Qlongwave_2023(1:12,1);  % Qnet - LWR
% bar_dmlt2(:,8)=Qsensible_2023(1:12,1);  % Qnet - Sensible
% 
% 
% subplot('position',pos{3})
%     h0=bar(5:8,bar_dmlt(5:8,:));
%        set(h0,'BarWidth',0.96); 
%        set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
%        hold on
%        set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
%        hold on
%        set(h0(3),'FaceColor',[0.366, 0.574, 0.188],'EdgeColor',[0.366, 0.574, 0.188])
%        hold on
%        set(h0(4),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])
%        
% 
% %     hold on
% %     h1=bar(5:8,bar_dmlt1(5:8,:));
% %        set(h1,'BarWidth',1.0); 
% %        hold on
% %        set(h1(2),'FaceColor','none','EdgeColor',[0.959, 0.494, 0.225])
%        
% %     hold on
% %     plot((5.0:0.01:5.4)',0*ones(length((5.0:0.01:5.4)),1),'-','color',[0.959, 0.494, 0.225],'linewidth',2)       
%     hold on
%     plot((5.0:0.01:5.4)',Qnet_2023(5,1)*ones(length((5.0:0.01:5.4)),1),':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot((6.0:0.01:6.4)',Qnet_2023(6,1)*ones(length((6.0:0.01:6.4)),1),':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot((7.0:0.01:7.4)',Qnet_2023(7,1)*ones(length((7.0:0.01:7.4)),1),':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot((8.0:0.01:8.4)',Qnet_2023(8,1)*ones(length((8.0:0.01:8.4)),1),':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%          
%     hold on
%     plot(5.0*ones(length((0:0.01:Qnet_2023(5,1))),1),(0:0.01:Qnet_2023(5,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(5.4*ones(length((0:0.01:Qnet_2023(5,1))),1),(0:0.01:Qnet_2023(5,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(6.0*ones(length((0:0.01:Qnet_2023(6,1))),1),(0:0.01:Qnet_2023(6,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(6.4*ones(length((0:0.01:Qnet_2023(6,1))),1),(0:0.01:Qnet_2023(6,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(7.0*ones(length((0:0.01:Qnet_2023(7,1))),1),(0:0.01:Qnet_2023(7,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(7.4*ones(length((0:0.01:Qnet_2023(7,1))),1),(0:0.01:Qnet_2023(7,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(8.0*ones(length((0:0.01:Qnet_2023(8,1))),1),(0:0.01:Qnet_2023(8,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(8.4*ones(length((0:0.01:Qnet_2023(8,1))),1),(0:0.01:Qnet_2023(8,1))',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%             
%     
%     hold on
%     plot((7.92:0.01:8.46)',0.55*ones(length((7.92:0.01:8.46)),1),':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot((7.92:0.01:8.46)',1.25*ones(length((7.92:0.01:8.46)),1),':','color',[0.959, 0.494, 0.225],'linewidth',1.2)          
%     hold on
%     plot(7.92*ones(length(0.55:0.01:1.25),1),(0.55:0.01:1.25)',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)       
%     hold on
%     plot(8.46*ones(length(0.55:0.01:1.25),1),(0.55:0.01:1.25)',':','color',[0.959, 0.494, 0.225],'linewidth',1.2)  
%     
%     
%     hold on
%     h2=bar(5:8,bar_dmlt2(5:8,:));
%        set(h2,'BarWidth',0.75); 
%        hold on
%        set(h2(5),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
%        hold on
%        set(h2(6),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
%        hold on
%        set(h2(7),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
%        hold on
%        set(h2(8),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933]) 
%        
% 
%        
%     hold on
% %     % set 3 display names for the 3 handles
% %     set(h0, {'DisplayName'}, {'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Vertical mixing + Entrainment','Advection + Residual'}')
%        
%     legend([h0(1) h0(2) h0(3) h0(4) h2(5) h2(6) h2(7) h2(8)],...
%            'MLT tendency','Surface flux term','Vertical mixing + entrainment','Advection + other minor terms','Shortwave','Latent','Longwave','Sensible',...
%            'Location','northeast','Orientation','vertical','NumColumns',2)
%     hold on
%     set(legend,'fontsize',18)
%     hold on
%     % title(leg4,'Monthly SST Anomaly','fontsize',20')
%     legend('boxoff')
%     
%     % Legend will show names for each color
%     legend() 
%     set(gca,'Ylim',[-0.4 1.30],'ycolor','k') 
%     set(gca,'YTick',-0.4:0.4:1.2)
%     set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',20)
%     set(gca,'Xlim',[4.5 8.5]) 
%     set(gca,'XTick',5:1:8)
%     set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',20)
% 
%     grid off
%     set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
%     ylabel(['[ \circC/month ]'],'fontsize',20,'color','k','FontWeight','normal')
% 
%     text(4.6,1.15,'g. MLT budget anomalies','fontsize',20,'color','k','FontWeight','bold')
% #########################################################################
    
    
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% 1.2 Plotting: ML Temperature Budget Anomaly in May-August 2023
%  Bar Charts and Maps ####################################################
% #########################################################################
% #########################################################################
% V2: moving panel (g) to (a)
clc; clear
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.050; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
ixw./iyw

%          [left            bottom      width height]
pos{3}  = [ixs+0*ixw+0*ixd    iys+2*iyw+2.5*iyd   3*ixw+2*ixd 1.05*iyw]; 


% ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
% iys = 0.050; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
% ixw./iyw
% 
% %          [left            bottom      width height]
% pos{11}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
% pos{12}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
% pos{13}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw];
% 
% pos{21}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
% pos{22}  = [ixs+1*ixw+1*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
% pos{23}  = [ixs+2*ixw+2*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);


% #########################################################################    
% #########################################################################
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1981-2010
      count=1;
      for year=1981:2010
          disp(['Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;

          
          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Vertical mixing
          mixz0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          mixz0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          mixz(:,:,:,count)=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
          mixz(mixz==0)=NaN;
          
          
          % ###############################################################
          % heat flux decomposition
          % (a) SW term 
          Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'swflx');
          % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
          Qshortwave(:,:,:,count)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
          Qshortwave(Qshortwave==0)=NaN;
          
          % (b) LW term
          Qlongwave(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'lw_heat').*second_2_month;

          % (c) Sensible term
          Qsensible(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sens_heat').*second_2_month;

          % (d) Latent term
          Qlatent(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'evap_heat').*second_2_month;
          
          % (e) Shortwave penetration at MLB
          Qswr_mld(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;
          % ###############################################################
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      mixz_clim=nanmean(mixz,4); clear mixz
      
      Qshortwave_clim=nanmean(Qshortwave,4); clear Qshortwave
      Qlongwave_clim=nanmean(Qlongwave,4);   clear Qlongwave
      Qsensible_clim=nanmean(Qsensible,4);   clear Qsensible
      Qlatent_clim=nanmean(Qlatent,4);       clear Qlatent
      Qswr_mld_clim=nanmean(Qswr_mld,4);     clear Qswr_mld
       
      % Removing the shortwave penetration from Qnet and SWR
      Qnet_clim=Qnet_clim+Qswr_mld_clim;
      Qshortwave_clim=Qshortwave_clim+Qswr_mld_clim;
      clear Qswr_mld_clim
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month;
      
      % (b) Entrainment
      entrainment=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month;
      
      % (c) Total advection
      adve0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
      Qnet(Qnet==0)=NaN;

      % (e) Vertical mixing
      mixz0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      mixz0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      mixz=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
      mixz(mixz==0)=NaN;
      
                
      % ###############################################################
      % heat flux decomposition
      % (a) SW term 
      Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'swflx');
      % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
      Qshortwave(:,:,:)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
      Qshortwave(Qshortwave==0)=NaN;

      % (b) LW term
      Qlongwave(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'lw_heat').*second_2_month;

      % (c) Sensible term
      Qsensible(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sens_heat').*second_2_month;

      % (d) Latent term
      Qlatent(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'evap_heat').*second_2_month;

      % (e) Shortwave penetration at MLB
      Qswr_mld(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sw_heat').*second_2_month;
      % ###############################################################
      % Removing the shortwave penetration from Qnet and SWR
      Qnet=Qnet+Qswr_mld;
      Qshortwave=Qshortwave+Qswr_mld;
      clear Qswr_mld
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      mixz=mixz-mixz_clim; 

      Qshortwave=Qshortwave-Qshortwave_clim; 
      Qlongwave=Qlongwave-Qlongwave_clim; 
      Qsensible=Qsensible-Qsensible_clim; 
      Qlatent=Qlatent-Qlatent_clim; 
      % ###################################################################

      
    % #########################################################################
      % Basin Mask from ACCESS OM2 0.25 forced by ERA5
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_NA_SST_fields/')
      load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
      [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
      Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
      for month=1:12
          basin_mask_NA0(:,:,month)=basin_mask_NA;
      end
      clear month
      
      dmlt(isnan(basin_mask_NA0))=NaN; 
      entrainment(isnan(basin_mask_NA0))=NaN;
      adve(isnan(basin_mask_NA0))=NaN; 
      Qnet(isnan(basin_mask_NA0))=NaN; 
      mixz(isnan(basin_mask_NA0))=NaN; 

      Qshortwave(isnan(basin_mask_NA0))=NaN; 
      Qlongwave(isnan(basin_mask_NA0))=NaN; 
      Qsensible(isnan(basin_mask_NA0))=NaN; 
      Qlatent(isnan(basin_mask_NA0))=NaN; 
      
      dmlt=dmlt.*Sxy;
      entrainment=entrainment.*Sxy;
      adve=adve.*Sxy;
      Qnet=Qnet.*Sxy;
      mixz=mixz.*Sxy;
      
      Qshortwave=Qshortwave.*Sxy;
      Qlongwave=Qlongwave.*Sxy;
      Qsensible=Qsensible.*Sxy;
      Qlatent=Qlatent.*Sxy;
          
      % NA 0-60N
      dmlt_2023(1:12,1)=squeeze(nansum(nansum(dmlt(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      entrainment_2023(1:12,1)=squeeze(nansum(nansum(entrainment(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      adve_2023(1:12,1)=squeeze(nansum(nansum(adve(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qnet_2023(1:12,1)=squeeze(nansum(nansum(Qnet(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      mixz_2023(1:12,1)=squeeze(nansum(nansum(mixz(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      
      Qshortwave_2023(1:12,1)=squeeze(nansum(nansum(Qshortwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlongwave_2023(1:12,1)=squeeze(nansum(nansum(Qlongwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qsensible_2023(1:12,1)=squeeze(nansum(nansum(Qsensible(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlatent_2023(1:12,1)=squeeze(nansum(nansum(Qlatent(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
    % #########################################################################
      
    

    % #########################################################################
      clear bar*
      bar_dmlt(:,1)=dmlt_2023(1:12,1);       % MLT tendency
      bar_dmlt(:,2)=Qnet_2023(1:12,1);       % Qnet
      bar_dmlt(:,3)=Qshortwave_2023(1:12,1); % Qnet - SWR
      bar_dmlt(:,4)=Qlatent_2023(1:12,1);    % Qnet - Latent
      bar_dmlt(:,5)=Qlongwave_2023(1:12,1);  % Qnet - LWR
      bar_dmlt(:,6)=Qsensible_2023(1:12,1);  % Qnet - Sensible
      bar_dmlt(:,7)=entrainment_2023(1:12,1)+mixz_2023(1:12,1);% Entrainment + Vertical mixing
      bar_dmlt(:,8)=adve_2023+(dmlt_2023-Qnet_2023-entrainment_2023-mixz_2023-adve_2023); % Temperature advection + Residual

    % SWR decomposition
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
      load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_1.mat',...
           'dMLT_Qswr_2023_ano')
      bar_dmlt(:,3)=dMLT_Qswr_2023_ano(1:12,1); % Qnet - SWR
      clear dMLT_Qswr_2023_ano
      
      subplot('position',pos{3})
        h0=bar(5:8,bar_dmlt(5:8,:));
           set(h0,'BarWidth',0.95); 
           set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
           hold on
           set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
           hold on
           set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
           hold on
           set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
           hold on
           set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
           hold on
           set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
           hold on
           set(h0(7),'FaceColor',[0.366, 0.574, 0.188],'EdgeColor',[0.366, 0.574, 0.188])
           hold on
           set(h0(8),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

       
       % ##################################################################
       % Add SWR decomposition
        hold on
        bar_dmlt_Clim(1:9,1:16)=NaN;
        % SWR cimatology
        cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_3_SWR_ciim.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,5)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano
                
        % MLD cimatology
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_2_MLD_clim.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,6)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano
        
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.36); 
           hold on
           set(h0_Clim(5),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(6),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           
           hold on
           plot((6.00:0.001:6.13),1.3*ones(length((5.05:0.001:5.18)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.00:0.001:6.13),1.05*ones(length((5.05:0.001:5.18)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.00:0.001:6.13),1.3*ones(length((5.05:0.001:5.18)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.00:0.001:6.13),1.05*ones(length((5.05:0.001:5.18)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)

           text(6.15, 1.3,'$\mathbf{\overline{Q_{sw}} / MLD^\prime}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
           text(6.15, 1.05,'$\mathbf{{Q_{sw}}^\prime / \overline{MLD}}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
       % ##################################################################


       
        hold on
        legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7) h0(8)],...
               'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Vertical mixing + entrainment','Advection + other minor terms',...
               'Location','northeast','Orientation','vertical','NumColumns',2)
        hold on
        set(legend,'fontsize',17)
        hold on
        % title(leg4,'Monthly SST Anomaly','fontsize',20')
        legend('boxoff')

        % Legend will show names for each color
        legend() 
        set(gca,'Ylim',[-0.4 1.63],'ycolor','k') 
        set(gca,'YTick',-0.4:0.4:1.6)
        set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2','1.6'},'fontsize',20)
        set(gca,'Xlim',[4.5 8.5]) 
        set(gca,'XTick',4.5:1:8.5)
        set(gca,'XTickLabel',{[],[],[],[]},'fontsize',20)
    %     set(gca,'XTick',5:1:8)
    %     set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',20)
        text(4.915,-0.54,'May','fontsize',20,'color','k','FontWeight','normal')
        text(5.930,-0.54,'Jun','fontsize',20,'color','k','FontWeight','normal')
        text(6.945,-0.54,'Jul','fontsize',20,'color','k','FontWeight','normal')
        text(7.920,-0.54,'Aug','fontsize',20,'color','k','FontWeight','normal')


        grid on
    %     ax = gca;
    %     ax.XGrid = 'off';
    %     ax.YGrid = 'on';
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
        ylabel(['[ \circC/month ]'],'fontsize',20,'color','k','FontWeight','normal')

        text(4.6,1.5,'a. MLT budget anomalies','fontsize',20,'color','k','FontWeight','bold')
    
% #########################################################################    
% #########################################################################
clear adv* basin* bar* dmlt* ent* Q* mix* second*



% #########################################################################    
% #########################################################################
ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.050; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
ixw./iyw

%          [left            bottom      width height]
pos{11}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{12}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw];

pos{21}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+2*ixw+2*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


% #########################################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);
color0(6,:)=(color0(6,:)+color0(5,:))./2;   
% #########################################################################


% #########################################################################
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1981-2010
      count=1;
      for year=1981:2010
          disp(['Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;


          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Shortwave penetration at MLB
          Qswr(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;


          % (f) Total mixing
          resd0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_submeso');
          resd0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_k33');
          resd0(:,:,:,3)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_diffusion_temp');
          resd0(:,:,:,4)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_gm_temp');
          resd0(:,:,:,5)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'residual');% residual=mixdownslope_temp + temp_sigma_diff  
          resd0(:,:,:,6)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_eta_smooth');
          
          resd0(:,:,:,7)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          resd0(:,:,:,8)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          
          resd0(:,:,:,9)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_rivermix');
          resd0(:,:,:,10)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sfc_hflux_pme');
          resd0(:,:,:,11)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'frazil_3d');
          
          resd(:,:,:,count)=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
          resd(resd==0)=NaN;
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      Qswr_clim=nanmean(Qswr,4); clear Qswr
      resd_clim=nanmean(resd,4); clear resd
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month;
      
      % (b) Entrainment
      entrainment=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month;
      
      % (c) Total advection
      adve0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
      Qnet(Qnet==0)=NaN;
      
      
      % (e) Shortwave penetration at MLB
      Qswr(:,:,:)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sw_heat').*second_2_month;


      % (f) Total mixing
      resd0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_submeso');
      resd0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_k33');
      resd0(:,:,:,3)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_diffusion_temp');
      resd0(:,:,:,4)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_gm_temp');
      resd0(:,:,:,5)=ncread('full_mixed_layer_heat_budget_year_2023.nc','residual');% residual=mixdownslope_temp + temp_sigma_diff  
      resd0(:,:,:,6)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_eta_smooth');
      
      resd0(:,:,:,7)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      resd0(:,:,:,8)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      
      resd0(:,:,:,9)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_rivermix');
      resd0(:,:,:,10)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sfc_hflux_pme');
      resd0(:,:,:,11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','frazil_3d');
      
      resd=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
      resd(resd==0)=NaN;
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      Qswr=Qswr-Qswr_clim; 
      resd=resd-resd_clim; 
      % ###################################################################

     

% #########################################################################       
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. MLT tendency (June 2023)','fontsize',20,'FontWeight','bold')
        
        
         
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,6)'+Qswr(:,:,6)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. Surface flux term','fontsize',20,'FontWeight','bold')

        

subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,6)'+resd(:,:,6)'+entrainment(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. Advection + mixing + entrainment','fontsize',19.5,'FontWeight','bold')
        
        
        


% #########################################################################       
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. MLT tendency (July 2023)','fontsize',20,'FontWeight','bold')
        
        
         
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,7)'+Qswr(:,:,7)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. Surface flux term','fontsize',20,'FontWeight','bold')

        

subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,7)'+resd(:,:,7)'+entrainment(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('g. Advection + mixing + entrainment','fontsize',19.5,'FontWeight','bold')

        
        hBar1 = colorbar('EastOutside','vertical');
        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+3*ixw+2*ixd+0.016 iys+0.5*iyw+0*iyd 0.012 1.0*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.060);
        % ylabel(hBar1, '[ ^oC/month ]','rotation',90);
        m_text(382,95, '[\circC/month]','fontsize',18,'FontWeight','normal')
% #########################################################################
% #########################################################################


    
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
% #########################################################################
    




%% 1.3 Plotting: ML Temperature Budget Anomaly in May-August 2023
%  Bar Charts and Maps ####################################################
% #########################################################################
% #########################################################################
% V3: moving panel (g) to (a), and add decomposition as panel (b)
clc; clear
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.200; ixe = 0.200;  ixd = 0.025; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.040; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3.3;
ixw./iyw

%          [left            bottom      width height]
pos{3}  = [ixs+0*ixw+0*ixd    iys+2.78*iyw+2.8*iyd   3*ixw+2*ixd 0.53*iyw]; 
pos{4}  = [ixs+0*ixw+0*ixd    iys+2.00*iyw+2.3*iyd   3*ixw+2*ixd 0.79*iyw]; 

% ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
% iys = 0.050; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
% ixw./iyw
% 
% %          [left            bottom      width height]
% pos{11}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
% pos{12}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
% pos{13}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw];
% 
% pos{21}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
% pos{22}  = [ixs+1*ixw+1*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
% pos{23}  = [ixs+2*ixw+2*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);


% #########################################################################    
% #########################################################################
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1981-2010
      count=1;
      for year=1981:2010
          disp(['Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;

          
          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Vertical mixing
          mixz0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          mixz0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          mixz(:,:,:,count)=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
          mixz(mixz==0)=NaN;
          
          
          % ###############################################################
          % heat flux decomposition
          % (a) SW term 
          Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'swflx');
          % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
          Qshortwave(:,:,:,count)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
          Qshortwave(Qshortwave==0)=NaN;
          
          % (b) LW term
          Qlongwave(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'lw_heat').*second_2_month;

          % (c) Sensible term
          Qsensible(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sens_heat').*second_2_month;

          % (d) Latent term
          Qlatent(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'evap_heat').*second_2_month;
          
          % (e) Shortwave penetration at MLB
          Qswr_mld(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;
          % ###############################################################
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      mixz_clim=nanmean(mixz,4); clear mixz
      
      Qshortwave_clim=nanmean(Qshortwave,4); clear Qshortwave
      Qlongwave_clim=nanmean(Qlongwave,4);   clear Qlongwave
      Qsensible_clim=nanmean(Qsensible,4);   clear Qsensible
      Qlatent_clim=nanmean(Qlatent,4);       clear Qlatent
      Qswr_mld_clim=nanmean(Qswr_mld,4);     clear Qswr_mld
       
      % Removing the shortwave penetration from Qnet and SWR
      Qnet_clim=Qnet_clim+Qswr_mld_clim;
      Qshortwave_clim=Qshortwave_clim+Qswr_mld_clim;
      clear Qswr_mld_clim
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month;
      
      % (b) Entrainment
      entrainment=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month;
      
      % (c) Total advection
      adve0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
      Qnet(Qnet==0)=NaN;

      % (e) Vertical mixing
      mixz0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      mixz0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      mixz=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
      mixz(mixz==0)=NaN;
      
                
      % ###############################################################
      % heat flux decomposition
      % (a) SW term 
      Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'swflx');
      % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
      Qshortwave(:,:,:)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
      Qshortwave(Qshortwave==0)=NaN;

      % (b) LW term
      Qlongwave(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'lw_heat').*second_2_month;

      % (c) Sensible term
      Qsensible(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sens_heat').*second_2_month;

      % (d) Latent term
      Qlatent(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'evap_heat').*second_2_month;

      % (e) Shortwave penetration at MLB
      Qswr_mld(:,:,:)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sw_heat').*second_2_month;
      % ###############################################################
      % Removing the shortwave penetration from Qnet and SWR
      Qnet=Qnet+Qswr_mld;
      Qshortwave=Qshortwave+Qswr_mld;
      clear Qswr_mld
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      mixz=mixz-mixz_clim; 

      Qshortwave=Qshortwave-Qshortwave_clim; 
      Qlongwave=Qlongwave-Qlongwave_clim; 
      Qsensible=Qsensible-Qsensible_clim; 
      Qlatent=Qlatent-Qlatent_clim; 
      % ###################################################################

      
    % #########################################################################
      % Basin Mask from ACCESS OM2 0.25 forced by ERA5
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_NA_SST_fields/')
      load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
      [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
      Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
      for month=1:12
          basin_mask_NA0(:,:,month)=basin_mask_NA;
      end
      clear month
      
      dmlt(isnan(basin_mask_NA0))=NaN; 
      entrainment(isnan(basin_mask_NA0))=NaN;
      adve(isnan(basin_mask_NA0))=NaN; 
      Qnet(isnan(basin_mask_NA0))=NaN; 
      mixz(isnan(basin_mask_NA0))=NaN; 

      Qshortwave(isnan(basin_mask_NA0))=NaN; 
      Qlongwave(isnan(basin_mask_NA0))=NaN; 
      Qsensible(isnan(basin_mask_NA0))=NaN; 
      Qlatent(isnan(basin_mask_NA0))=NaN; 
      
      dmlt=dmlt.*Sxy;
      entrainment=entrainment.*Sxy;
      adve=adve.*Sxy;
      Qnet=Qnet.*Sxy;
      mixz=mixz.*Sxy;
      
      Qshortwave=Qshortwave.*Sxy;
      Qlongwave=Qlongwave.*Sxy;
      Qsensible=Qsensible.*Sxy;
      Qlatent=Qlatent.*Sxy;
          
      % NA 0-60N
      dmlt_2023(1:12,1)=squeeze(nansum(nansum(dmlt(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      entrainment_2023(1:12,1)=squeeze(nansum(nansum(entrainment(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      adve_2023(1:12,1)=squeeze(nansum(nansum(adve(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qnet_2023(1:12,1)=squeeze(nansum(nansum(Qnet(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      mixz_2023(1:12,1)=squeeze(nansum(nansum(mixz(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      
      Qshortwave_2023(1:12,1)=squeeze(nansum(nansum(Qshortwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlongwave_2023(1:12,1)=squeeze(nansum(nansum(Qlongwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qsensible_2023(1:12,1)=squeeze(nansum(nansum(Qsensible(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlatent_2023(1:12,1)=squeeze(nansum(nansum(Qlatent(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
    % #########################################################################
      
    

    % #########################################################################
      clear bar*
      bar_dmlt(:,1)=dmlt_2023(1:12,1);       % MLT tendency
      bar_dmlt(:,2)=Qnet_2023(1:12,1);       % Qnet
%       bar_dmlt(:,3)=Qshortwave_2023(1:12,1); % Qnet - SWR
%       bar_dmlt(:,4)=Qlatent_2023(1:12,1);    % Qnet - Latent
%       bar_dmlt(:,5)=Qlongwave_2023(1:12,1);  % Qnet - LWR
%       bar_dmlt(:,6)=Qsensible_2023(1:12,1);  % Qnet - Sensible
      bar_dmlt(:,3)=entrainment_2023(1:12,1)+mixz_2023(1:12,1);% Entrainment + Vertical mixing
      bar_dmlt(:,4)=adve_2023+(dmlt_2023-Qnet_2023-entrainment_2023-mixz_2023-adve_2023); % Temperature advection + Residual


      subplot('position',pos{3})
        h0=bar(5:8,bar_dmlt(5:8,:));
           set(h0,'BarWidth',0.95); 
           set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
           hold on
           set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
           hold on
           set(h0(3),'FaceColor',[0.366, 0.574, 0.188],'EdgeColor',[0.366, 0.574, 0.188])
           hold on
           set(h0(4),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

        hold on
        legend([h0(1) h0(2) h0(3) h0(4)],...
               'MLT tendency','Surface flux term','Vertical mixing + entrainment','Advection + other minor terms',...
               'Location','northeast','Orientation','vertical','NumColumns',3)
        hold on
        set(legend,'fontsize',15)
        hold on
        % title(leg4,'Monthly SST Anomaly','fontsize',20')
        legend('boxoff')

        % Legend will show names for each color
        legend() 
        set(gca,'Ylim',[-0.4 0.8],'ycolor','k') 
        set(gca,'YTick',-0.4:0.4:1.2)
        set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',20)
        set(gca,'Xlim',[4.5 8.5]) 
        set(gca,'XTick',4.5:1:8.5)
        set(gca,'XTickLabel',{[],[],[],[]},'fontsize',20)
%         text(4.915,-0.54,'May','fontsize',20,'color','k','FontWeight','normal')
%         text(5.930,-0.54,'Jun','fontsize',20,'color','k','FontWeight','normal')
%         text(6.945,-0.54,'Jul','fontsize',20,'color','k','FontWeight','normal')
%         text(7.920,-0.54,'Aug','fontsize',20,'color','k','FontWeight','normal')

        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
        ylabel(['[ \circC/month ]'],'fontsize',20,'color','k','FontWeight','normal')

        text(4.6,0.65,'a. MLT budget anomalies','fontsize',20,'color','k','FontWeight','bold')
    
        
        
        
    % #########################################################################
      clear bar*
%       bar_dmlt(:,1)=dmlt_2023(1:12,1);       % MLT tendency
      bar_dmlt(:,1)=Qnet_2023(1:12,1);       % Qnet
      bar_dmlt(:,2)=Qshortwave_2023(1:12,1); % Qnet - SWR
      bar_dmlt(:,3)=Qlatent_2023(1:12,1);    % Qnet - Latent
      bar_dmlt(:,4)=Qlongwave_2023(1:12,1);  % Qnet - LWR
      bar_dmlt(:,5)=Qsensible_2023(1:12,1);  % Qnet - Sensible
%       bar_dmlt(:,7)=entrainment_2023(1:12,1)+mixz_2023(1:12,1);% Entrainment + Vertical mixing
%       bar_dmlt(:,8)=adve_2023+(dmlt_2023-Qnet_2023-entrainment_2023-mixz_2023-adve_2023); % Temperature advection + Residual

      % SWR decomposition
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
      load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_1.mat',...
           'dMLT_Qswr_2023_ano')
      bar_dmlt(:,2)=dMLT_Qswr_2023_ano(1:12,1); % Qnet - SWR
      clear dMLT_Qswr_2023_ano

      subplot('position',pos{4})
        h0=bar(5:8,bar_dmlt(5:8,:));
           set(h0,'BarWidth',0.90); 
           hold on
           set(h0(1),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
           hold on
           set(h0(2),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
           hold on
           set(h0(3),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
           hold on
           set(h0(4),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
           hold on
           set(h0(5),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])

                 
       % ##################################################################
       % Add SWR decomposition
        hold on
        bar_dmlt_Clim(1:9,1:10)=NaN;
        % SWR cimatology
        cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_3_SWR_ciim.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,3)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano
                
        % MLD cimatology
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V6_V4_Decomposition_2_MLD_clim.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,4)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano
        
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.40); 
           hold on
           set(h0_Clim(3),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(4),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           
           hold on
           plot((6.70:0.001:6.83),1.35*ones(length((5.05:0.001:5.18)')),'-','color',[0.929, 0.694, 0.125],'linewidth',16)
           hold on
           plot((6.70:0.001:6.83),1.05*ones(length((5.05:0.001:5.18)')),'-','color',[0.929, 0.694, 0.125],'linewidth',16)
           hold on
           plot((6.70:0.001:6.83),1.35*ones(length((5.05:0.001:5.18)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.70:0.001:6.83),1.05*ones(length((5.05:0.001:5.18)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)

           text(6.85, 1.35,'$\mathbf{\overline{Q_{sw}} / MLD^\prime}$', 'Interpreter', 'latex','fontsize',16,'FontName', 'Aptos')
           text(6.85, 1.05,'$\mathbf{{Q_{sw}}^\prime / \overline{MLD}}$', 'Interpreter', 'latex','fontsize',16,'FontName', 'Aptos')
       % ##################################################################
           
        hold on
        legend([h0(1) h0(2) h0(3) h0(4) h0(5)],...
               'Surface flux term','Shortwave','Latent','Longwave','Sensible',...
               'Location','northeast','Orientation','vertical','NumColumns',2)
        hold on
        set(legend,'fontsize',15)
        hold on
        % title(leg4,'Monthly SST Anomaly','fontsize',20')
        legend('boxoff')

%       clear bar*
%       bar_dmlt0(1:12,1)=NaN;
%       bar_dmlt0(1:12,2)=Qnet_2023(1:12,1);       % Qnet
%       bar_dmlt0(1:12,3)=NaN;
%       bar_dmlt0(1:12,4)=NaN;
%       
% 
%       bar_dmlt(1:12,1:4)=NaN;
%       bar_dmlt(:,5)=Qshortwave_2023(1:12,1); % Qnet - SWR
%       bar_dmlt(:,6)=Qlatent_2023(1:12,1);    % Qnet - Latent
%       bar_dmlt(:,7)=Qlongwave_2023(1:12,1);  % Qnet - LWR
%       bar_dmlt(:,8)=Qsensible_2023(1:12,1);  % Qnet - Sensible
% 
%       
%       subplot('position',pos{4})
%         h0=bar(5:8,bar_dmlt0(5:8,:));
%         set(h0,'BarWidth',0.95); 
%         hold on
%         h1=bar(5:8,bar_dmlt(5:8,:));
%         set(h1,'BarWidth',0.95); 
%            
%            hold on
%            set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
%            hold on
%            set(h1(5),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
%            hold on
%            set(h1(6),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
%            hold on
%            set(h1(7),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
%            hold on
%            set(h1(8),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
% 
%         hold on
%         legend([h0(2) h1(5) h1(6) h1(7) h1(8)],...
%                'Surface flux term','Shortwave','Latent','Longwave','Sensible',...
%                'Location','northeast','Orientation','vertical','NumColumns',2)
%         hold on
%         set(legend,'fontsize',18)
%         hold on
%         % title(leg4,'Monthly SST Anomaly','fontsize',20')
%         legend('boxoff')      
        
        
        % Legend will show names for each color
        legend() 
        set(gca,'Ylim',[-0.4 1.63],'ycolor','k') 
        set(gca,'YTick',-0.4:0.4:1.6)
        set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2','1.6'},'fontsize',20)
        set(gca,'Xlim',[4.5 8.5]) 
        set(gca,'XTick',4.5:1:8.5)
        set(gca,'XTickLabel',{[],[],[],[]},'fontsize',20)
        text(4.915,-0.54,'May','fontsize',20,'color','k','FontWeight','normal')
        text(5.930,-0.54,'Jun','fontsize',20,'color','k','FontWeight','normal')
        text(6.945,-0.54,'Jul','fontsize',20,'color','k','FontWeight','normal')
        text(7.920,-0.54,'Aug','fontsize',20,'color','k','FontWeight','normal')

        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
        ylabel(['[ \circC/month ]'],'fontsize',20,'color','k','FontWeight','normal')

        text(4.6,1.5,'b. Surface flux term','fontsize',20,'color','k','FontWeight','bold')
    
% #########################################################################    
% #########################################################################
clear adv* basin* bar* dmlt* ent* Q* mix* second*




% #########################################################################    
% #########################################################################
ixs = 0.200; ixe = 0.200;  ixd = 0.025; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.040; iye = 0.050;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3.3;
ixw./iyw

%          [left            bottom      width height]
pos{11}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{12}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw];

pos{21}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+2*ixw+2*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


% #########################################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);
color0(6,:)=(color0(6,:)+color0(5,:))./2;   
% #########################################################################


% #########################################################################
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1981-2010
      count=1;
      for year=1981:2010
          disp(['Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;


          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Shortwave penetration at MLB
          Qswr(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;


          % (f) Total mixing
          resd0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_submeso');
          resd0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_k33');
          resd0(:,:,:,3)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_diffusion_temp');
          resd0(:,:,:,4)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_gm_temp');
          resd0(:,:,:,5)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'residual');% residual=mixdownslope_temp + temp_sigma_diff  
          resd0(:,:,:,6)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_eta_smooth');
          
          resd0(:,:,:,7)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          resd0(:,:,:,8)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          
          resd0(:,:,:,9)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_rivermix');
          resd0(:,:,:,10)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sfc_hflux_pme');
          resd0(:,:,:,11)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'frazil_3d');
          
          resd(:,:,:,count)=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
          resd(resd==0)=NaN;
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      Qswr_clim=nanmean(Qswr,4); clear Qswr
      resd_clim=nanmean(resd,4); clear resd
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month;
      
      % (b) Entrainment
      entrainment=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month;
      
      % (c) Total advection
      adve0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
      Qnet(Qnet==0)=NaN;
      
      
      % (e) Shortwave penetration at MLB
      Qswr(:,:,:)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sw_heat').*second_2_month;


      % (f) Total mixing
      resd0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_submeso');
      resd0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_k33');
      resd0(:,:,:,3)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_diffusion_temp');
      resd0(:,:,:,4)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_gm_temp');
      resd0(:,:,:,5)=ncread('full_mixed_layer_heat_budget_year_2023.nc','residual');% residual=mixdownslope_temp + temp_sigma_diff  
      resd0(:,:,:,6)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_eta_smooth');
      
      resd0(:,:,:,7)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      resd0(:,:,:,8)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      
      resd0(:,:,:,9)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_rivermix');
      resd0(:,:,:,10)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sfc_hflux_pme');
      resd0(:,:,:,11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','frazil_3d');
      
      resd=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
      resd(resd==0)=NaN;
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      Qswr=Qswr-Qswr_clim; 
      resd=resd-resd_clim; 
      % ###################################################################

     

% #########################################################################       
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. MLT tendency (June 2023)','fontsize',20,'FontWeight','bold')
        
        
         
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,6)'+Qswr(:,:,6)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. Surface flux term','fontsize',20,'FontWeight','bold')

        

subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,6)'+resd(:,:,6)'+entrainment(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. Advection + mixing + entrainment','fontsize',19,'FontWeight','bold')
        
        
        


% #########################################################################       
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. MLT tendency (July 2023)','fontsize',20,'FontWeight','bold')
        
        
         
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,7)'+Qswr(:,:,7)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('g. Surface flux term','fontsize',20,'FontWeight','bold')

        

subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,7)'+resd(:,:,7)'+entrainment(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('h. Advection + mixing + entrainment','fontsize',19,'FontWeight','bold')

        
        hBar1 = colorbar('EastOutside','vertical');
        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+3*ixw+2*ixd+0.016 iys+0.5*iyw+0*iyd 0.012 1.0*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.065);
        m_text(382,95, '[\circC/month]','fontsize',18,'FontWeight','normal')
% #########################################################################
% #########################################################################


    
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
% #########################################################################
    



