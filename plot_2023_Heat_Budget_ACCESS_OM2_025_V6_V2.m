%% Mixed-Layer Heat Budget Anomaly in 2023, ACCESS-OM2 0.25 JRA55-do IAF Forced Simulation

%% ########################################################################
%% Plotting 1: ML Heat Budget Anomaly in June July 2023
clc;clear
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the extended right screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      [Sxy,~,~]=function_Cgrid_Area_Distance(lon_025,lat_025);
      second_2_month=60*60*24*30;
      
      
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
%           mixz0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_submeso');
%           mixz0(:,:,:,3)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_k33');
%           mixz0(:,:,:,4)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_diffusion_temp');
%           mixz0(:,:,:,4)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_gm_temp');
%           mixz0(:,:,:,5)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'residual');% residual=mixdownslope_temp + temp_sigma_diff  
%           mixz0(:,:,:,6)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_eta_smooth');

%           mixz0(:,:,:,9)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_rivermix');
%           mixz0(:,:,:,10)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sfc_hflux_pme');
%           mixz0(:,:,:,11)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'frazil_3d');
          
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
      
%       mixz0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_submeso');
%       mixz0(:,:,:,3)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_k33');
%       mixz0(:,:,:,4)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_diffusion_temp');
%       mixz0(:,:,:,4)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_gm_temp');
%       mixz0(:,:,:,5)=ncread('full_mixed_layer_heat_budget_year_2023.nc','residual');% residual=mixdownslope_temp + temp_sigma_diff  
%       mixz0(:,:,:,6)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_eta_smooth');
      
%       mixz0(:,:,:,9)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_rivermix');
%       mixz0(:,:,:,10)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sfc_hflux_pme');
%       mixz0(:,:,:,11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','frazil_3d');
      
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
      
  
   
%% ########################################################################
% V1 - Colorful Bar Charts ################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.220; ixe = 0.220;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 

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


subplot('position',pos{101})
    h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',1.0); 
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

    hold on
%     % set 3 display names for the 3 handles
%     set(h0, {'DisplayName'}, {'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Vertical mixing + Entrainment','Advection + Residual'}')
       
    legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7) h0(8)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Vertical mixing + entrainment','Advection + residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',20)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Legend will show names for each color
    legend() 
    set(gca,'Ylim',[-0.4 1.30],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',24)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:0.5:8.5)
    set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',24)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',24,'color','k','FontWeight','normal')

    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')



%% ########################################################################
% V2 - Bar Charts and Maps ################################################
% #########################################################################
% #########################################################################
% Plotting 2.1: ML Heat Budget Anomaly in June July 2023 3*2
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


clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);
% color0(8,:)=(color0(7,:)+color0(8,:))./2;   
color0(6,:)=(color0(6,:)+color0(5,:))./2;   


   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the left screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      second_2_month=60*60*24*30;
      
      
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
      second_2_month=60*60*24*30;
      
      
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

    hold on
%     % set 3 display names for the 3 handles
%     set(h0, {'DisplayName'}, {'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Vertical mixing + Entrainment','Advection + Residual'}')
       
    legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7) h0(8)],...
           'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Vertical mixing + entrainment','Advection + other minor terms',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',18)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
%     % Legend will show names for each color
%     legend() 
%     set(gca,'Ylim',[-0.4 1.30],'ycolor','k') 
%     set(gca,'YTick',-0.4:0.4:1.2)
%     set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',20)
%     set(gca,'Xlim',[4.5 8.5]) 
%     set(gca,'XTick',4.5:0.5:8.5)
%     set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',20)
% %     set(gca,'XTick',5:1:8)
% %     set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',20)

    % Legend will show names for each color
    legend() 
    set(gca,'Ylim',[-0.4 1.30],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',20)
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

    text(4.6,1.15,'g. MLT budget anomalies','fontsize',20,'color','k','FontWeight','bold')
    
    

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





%% ########################################################################
%% ########################################################################
%% Plotting 3: ML Heat Budget Anomaly in June July 2023 - Time series of integrated heat budget
clc;clear
   disp(['MLD and MLT from ACCESS OM2 025...'])
   disp(['Save to the extended right screen...'])
      cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/JRA55-do_IAF_mixed_layer_heat_budget_V4_with_monthly_SST/')
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      [Sxy,~,~]=function_Cgrid_Area_Distance(lon_025,lat_025);
      second_2_month=60*60*24*30;
      
      
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
          
          
          % Monthly SST
          sst0=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sst');
          sst_mon(:,:,count)=sst0(:,:,1); clear sst0
          count=count+1;
      end 
      clear year count
      
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
      % Predicted MLT = accummulated / intgrated heat flux anomalies
        project_mlt=nan(size(dmlt_clim));
        project_mlt_Qnet=nan(size(dmlt_clim));
        % project_mlt_Qent=nan(size(dmlt_clim));
        project_mlt_Qadv=nan(size(dmlt_clim));
        project_mlt_Qmix=nan(size(dmlt_clim));
        
        project_mlt_Qswr=nan(size(dmlt_clim));
        project_mlt_Qlat=nan(size(dmlt_clim));
        project_mlt_Qlon=nan(size(dmlt_clim));
        project_mlt_Qsen=nan(size(dmlt_clim));
        
        % Jan-15
        sst_clim=nanmean(sst_mon,3); clear sst_mon
        project_mlt(:,:,1)=sst_clim(:,:,1);
        project_mlt_Qnet(:,:,1)=sst_clim(:,:,1);
        % project_mlt_Qent(:,:,1)=sst_clim(:,:,1);
        project_mlt_Qadv(:,:,1)=sst_clim(:,:,1);
        project_mlt_Qmix(:,:,1)=sst_clim(:,:,1);
        
        project_mlt_Qswr(:,:,1)=sst_clim(:,:,1);
        project_mlt_Qlat(:,:,1)=sst_clim(:,:,1);
        project_mlt_Qlon(:,:,1)=sst_clim(:,:,1);
        project_mlt_Qsen(:,:,1)=sst_clim(:,:,1);
        
        % Projected mlt from Feb-15
        for month=2:12
            disp(['   Projeted monthly MLT mon#',num2str(month)])
            project_mlt(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(dmlt_clim(:,:,2:month),3));
            project_mlt_Qnet(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(Qnet_clim(:,:,2:month),3));
            % project_mlt_Qent(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(entrainment_clim(:,:,2:month),3));
            project_mlt_Qadv(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(adve_clim(:,:,2:month),3));
            project_mlt_Qmix(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(mixz_clim(:,:,2:month),3))+squeeze(nansum(entrainment_clim(:,:,2:month),3));
            
            project_mlt_Qswr(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(Qshortwave_clim(:,:,2:month),3));
            project_mlt_Qlat(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(Qlatent_clim(:,:,2:month),3));
            project_mlt_Qlon(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(Qlongwave_clim(:,:,2:month),3));
            project_mlt_Qsen(:,:,month)=sst_clim(:,:,1)+squeeze(nansum(Qsensible_clim(:,:,2:month),3));
        end
        clear sst_daily *_clim
        % figure;imagesc(project_sst_Qnet(240:360,91:150,1)); caxis([-3 3])
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
      % Predicted MLT = accummulated / intgrated heat flux anomalies
        project_mlt_2023=nan(size(dmlt));
        project_mlt_Qnet_2023=nan(size(dmlt));
        project_mlt_Qent_2023=nan(size(dmlt));
        project_mlt_Qadv_2023=nan(size(dmlt));
        project_mlt_Qmix_2023=nan(size(dmlt));
        
        project_mlt_Qswr_2023=nan(size(dmlt));
        project_mlt_Qlat_2023=nan(size(dmlt));
        project_mlt_Qlon_2023=nan(size(dmlt));
        project_mlt_Qsen_2023=nan(size(dmlt));
        
        % Jan-15
        % Monthly SST
        sst_2023=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sst');
        project_mlt_2023(:,:,1)=sst_2023(:,:,1);
        project_mlt_Qnet_2023(:,:,1)=sst_2023(:,:,1);
        %project_mlt_Qent_2023(:,:,1)=sst_2023(:,:,1);
        project_mlt_Qadv_2023(:,:,1)=sst_2023(:,:,1);
        project_mlt_Qmix_2023(:,:,1)=sst_2023(:,:,1);
        
        project_mlt_Qswr_2023(:,:,1)=sst_2023(:,:,1);
        project_mlt_Qlat_2023(:,:,1)=sst_2023(:,:,1);
        project_mlt_Qlon_2023(:,:,1)=sst_2023(:,:,1);
        project_mlt_Qsen_2023(:,:,1)=sst_2023(:,:,1);
        
        % Projected mlt from Feb-15
        for month=2:12
            disp(['   Projeted monthly MLT mon#',num2str(month)])
            project_mlt_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(dmlt(:,:,2:month),3));
            project_mlt_Qnet_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(Qnet(:,:,2:month),3));
            %project_mlt_Qent_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(entrainment(:,:,2:month),3));
            project_mlt_Qadv_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(adve(:,:,2:month),3));
            project_mlt_Qmix_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(mixz(:,:,2:month),3))+squeeze(nansum(entrainment(:,:,2:month),3));
            
            project_mlt_Qswr_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(Qshortwave(:,:,2:month),3));
            project_mlt_Qlat_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(Qlatent(:,:,2:month),3));
            project_mlt_Qlon_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(Qlongwave(:,:,2:month),3));
            project_mlt_Qsen_2023(:,:,month)=sst_2023(:,:,1)+squeeze(nansum(Qsensible(:,:,2:month),3));
        end
        clear sst_2023
        % figure;imagesc(project_mlt_2023(:,:,1)); caxis([270 310])
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
%       dmlt=dmlt-dmlt_clim; 
%       entrainment=entrainment-entrainment_clim;
%       adve=adve-adve_clim; 
%       Qnet=Qnet-Qnet_clim; 
%       mixz=mixz-mixz_clim; 
% 
%       Qshortwave=Qshortwave-Qshortwave_clim; 
%       Qlongwave=Qlongwave-Qlongwave_clim; 
%       Qsensible=Qsensible-Qsensible_clim; 
%       Qlatent=Qlatent-Qlatent_clim; 
      
      dmlt=project_mlt_2023-project_mlt; 
      Qnet=project_mlt_Qnet_2023-project_mlt_Qnet; 
      %entrainment=project_mlt_Qent_2023-project_mlt_Qent;
      adve=project_mlt_Qadv_2023-project_mlt_Qadv; 
      mixz=project_mlt_Qmix_2023-project_mlt_Qmix; 

      Qshortwave=project_mlt_Qswr_2023-project_mlt_Qswr; 
      Qlatent=project_mlt_Qlat_2023-project_mlt_Qlat; 
      Qlongwave=project_mlt_Qlon_2023-project_mlt_Qlon; 
      Qsensible=project_mlt_Qsen_2023-project_mlt_Qsen;      
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
      %entrainment(isnan(basin_mask_NA0))=NaN;
      adve(isnan(basin_mask_NA0))=NaN; 
      Qnet(isnan(basin_mask_NA0))=NaN; 
      mixz(isnan(basin_mask_NA0))=NaN; 

      Qshortwave(isnan(basin_mask_NA0))=NaN; 
      Qlongwave(isnan(basin_mask_NA0))=NaN; 
      Qsensible(isnan(basin_mask_NA0))=NaN; 
      Qlatent(isnan(basin_mask_NA0))=NaN; 
      
      dmlt=dmlt.*Sxy;
      %entrainment=entrainment.*Sxy;
      adve=adve.*Sxy;
      Qnet=Qnet.*Sxy;
      mixz=mixz.*Sxy;
      
      Qshortwave=Qshortwave.*Sxy;
      Qlongwave=Qlongwave.*Sxy;
      Qsensible=Qsensible.*Sxy;
      Qlatent=Qlatent.*Sxy;
          
      % NA 0-60N
      dmlt_2023(1:12,1)=squeeze(nansum(nansum(dmlt(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      %entrainment_2023(1:12,1)=squeeze(nansum(nansum(entrainment(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      adve_2023(1:12,1)=squeeze(nansum(nansum(adve(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qnet_2023(1:12,1)=squeeze(nansum(nansum(Qnet(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      mixz_2023(1:12,1)=squeeze(nansum(nansum(mixz(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      
      Qshortwave_2023(1:12,1)=squeeze(nansum(nansum(Qshortwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlongwave_2023(1:12,1)=squeeze(nansum(nansum(Qlongwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qsensible_2023(1:12,1)=squeeze(nansum(nansum(Qsensible(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlatent_2023(1:12,1)=squeeze(nansum(nansum(Qlatent(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
    % #########################################################################
      
    

% #########################################################################
% Simulated SSTa from ACCESS OM2 0.25 forced by ERA5 and JRA55
% % a. SSTa from ERA5 IAF run
% cd('/Users/z5195509/Documents/6_NA_MHW_MLD/7_ML_Heat_Budget_ACCESS-OM2')
% load('plot_2023_SST_Anomaly_ACCESS_OM2_025_V1_1_ERA5_SSTa.mat',...
%      'sst_daily_2023_ano')
% sst_daily_2023_ano_ERA5=sst_daily_2023_ano; 
% clear sst_daily_2023_ano
%  
% b. SSTa from JRA5-do IAF run
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/7_ML_Heat_Budget_ACCESS-OM2')
load('plot_2023_SST_Anomaly_ACCESS_OM2_025_V1_2_JRA55_SSTa.mat',...
     'sst_daily_2023_ano')
sst_daily_2023_ano_JRA55=sst_daily_2023_ano; 
clear sst_daily_2023_ano
% #########################################################################

    

    
%% Time series of integrated heat budget ##############################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.220; ixe = 0.220;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+0*iyw+1*iyd   ixw 2*iyw]; 

% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.

    
subplot('position',pos{101})
    line2023=plot(1:366,sst_daily_2023_ano_JRA55(:,1));
    set(line2023,'color',[0.2,0.2,0.2],'LineWidth',3,'linestyle','-'); 
    
    hold on
    line2023MLT6=plot(16:30:330,adve_2023(1:11,1));
    set(line2023MLT6,'color',[0.594, 0.284, 0.556],'LineWidth',3,'linestyle','-'); 
    line2023MLT3=plot(16:30:330,Qlatent_2023(1:11,1));
    set(line2023MLT3,'color',[0    , 0.447, 0.641],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT4=plot(16:30:330,Qlongwave_2023(1:11,1));
    set(line2023MLT4,'color',[0.101, 0.645, 0.833],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT5=plot(16:30:330,Qsensible_2023(1:11,1));
    set(line2023MLT5,'color',[0.400, 0.850, 0.933],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT2=plot(16:30:330,Qshortwave_2023(1:11,1));
    set(line2023MLT2,'color',[0.929, 0.694, 0.125],'LineWidth',4,'linestyle','-'); 
    
    hold on
    line2023MLT7=plot(16:30:330,mixz_2023(1:11,1));
    set(line2023MLT7,'color',[0.366, 0.574, 0.188],'LineWidth',4,'linestyle','-'); 

    hold on
    line2023MLT1=plot(16:30:330,Qnet_2023(1:11,1));
    set(line2023MLT1,'color',[0.959, 0.494, 0.225],'LineWidth',6,'linestyle','-'); 
    
    hold on
    line2023MLT0=plot(16:30:330,dmlt_2023(1:11,1));
    set(line2023MLT0,'color',[0.850, 0.325, 0.098],'LineWidth',6,'linestyle','-'); 

    
    leg101=legend([line2023 line2023MLT0 line2023MLT1 line2023MLT2 line2023MLT5 line2023MLT4 line2023MLT3 line2023MLT7 line2023MLT6],...
               '2023 Daily SSTA','2023 MLTA','Qnet','Shortwave','Sensible','Longwave','Latent heat','Mix + Ent','Adv + Res',...
               'Location','northwest','NumColumns',1);
    set(leg101,'fontsize',20)
    hold on
    title(leg101,'Heat budget from ACCESS OM2 JRA55-do','fontsize',20')
    legend('boxoff')

    
    set(gca,'Ylim',[-1 2.5],'ycolor','k') 
    set(gca,'YTick',-1.5:0.5:4)
    set(gca,'YTickLabel',{'-1.5','-1.0','-0.5','0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'},'fontsize',24)
    set(gca,'Xlim',[1 366]) 
    set(gca,'XTick',1:30:360)
    set(gca,'XTickLabel',{'           Jan','           Feb','           Mar','           Apr','           May','           Jun',...
                          '           Jul','           Aug','           Sep','           Oct','           Nov','            Dec'},'fontsize',24)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['Difference [ ^oC ]'],'fontsize',24,'color','k','FontWeight','normal')

    % title('a. North Atlantic SST and predicted MLT anomaly','fontsize',24,'color','k','FontWeight','bold')



    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')


