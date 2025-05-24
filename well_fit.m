function well_fit(varargin)

%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: Mar 2024
%-------------------------------------------------------------------------%
%
% Description:
% Well fit is designed to take predicted Pecube-D cooling age files, to
% determine how well input rates descrribe the observed cooling ages and
% exhumation depths, if independent constraints are known. well-fit 
% histories are determined by the chi squared caculation. Note that is is
% different from a best-fit, in that this technique requires that each
% predicted chronometer meet the statistical test. this allows for many
% exhumation hisotries to produce accurate cooling ages.
%
% Usage: 
% well_fit(); (uses all default values below)
% well_fit('integration_time',10,'tolerance',3);
% 
% Required input files:
% Pecube.in - this is the input file for the monte carlo run
% ages_XXXX.csv - this is the csv file containing observed cooling ages.
%                 this needs to be formatted for input into Pecube.
% ages_all_XXXX.txt - text file that includes the predicted cooling ages
%                     from all simulations. 
% rates_all_XXXX.txt - text file that includes the tested exhumation rates
%                     from all simulations. 
%
% Optional Inputs:
% integration_time - the time (Ma) to begin the calculation of the exhumed
%                    thickness from exhumation histories. needs to be an
%                    integer
% tolerance - premissable valued of chi squared
% sigma_perc_ages - a percentage uncertainty that can be applied to ages
%                   during the chi squared calculation (%)
% depth - an independent constraint on the total exhumation (km) over the
%         integration time
% depth_sigma - 1 sigma uncertainty on the independent constraint on the 
%               total exhumation (km)
% sigma_perc_depth - a percentage uncertainty that can be applied to depths
%                   during the chi squared calculation (%)
% save_figs - flag to save figures as eps and png. 'y' for yes.
% err_disp - input to determine how the uncertainty on observed ages are
%            represented in figures. Options: 'pdf' for showing normal
%            distributions based on the mean and 1sigma uncertainty, or
%            'sigma' for showing rectangles based on 2sigma uncertainty
% cull_data - flag to remove data based on fit to independent constraints
%             on exhumed thicknesses. 'y' for yes.
%
% Outputs:
% figures and statistics
% 
%-------------------------------------------------------------------------%
% tashi delek!
%-------------------------------------------------------------------------%

% parse inputs
p = inputParser;         
p.FunctionName = 'well_fit';
addParameter(p,'integration_time',0,@(x) isscalar(x));
addParameter(p,'tolerance',3,@(x) isscalar(x));
addParameter(p,'sigma_perc_ages',0,@(x) isscalar(x));
addParameter(p,'depth',0,@(x) isscalar(x));
addParameter(p,'depth_sigma',0,@(x) isscalar(x));
addParameter(p,'sigma_perc_depth',0,@(x) isscalar(x));
addParameter(p,'save_figs','n',@(x) ischar(x));
addParameter(p,'err_disp','pdf',@(x) ischar(x));
addParameter(p,'cull_data','n',@(x) ischar(x));
    
parse(p,varargin{:});
integration_time = p.Results.integration_time;
tolerance = p.Results.tolerance;
sigma_perc_ages = p.Results.sigma_perc_ages;
sigma_perc_ages = sigma_perc_ages/100;
depth = p.Results.depth;
depth_sigma = p.Results.depth_sigma;
sigma_perc_depth = p.Results.sigma_perc_depth;
sigma_perc_depth = sigma_perc_depth/100;
save_figs = p.Results.save_figs;
err_disp = p.Results.err_disp;
cull_data = p.Results.cull_data;

% open the input file and pull out the time slices
file_id = fopen('Pecube.in');
line = fgets(file_id);
while ischar(line)
    if contains(line,'mc_time_slices:')
        t_string = line(17:end);
        break
    end
    line = fgets(file_id);
end
fclose(file_id);
t_string = strrep(t_string,',',' ');
t_string = strrep(t_string,'-',' ');
times = str2num(t_string); %#ok<ST2NM>
times = sort(times,'descend');
simpler_times = fliplr(unique(times));
num_slices = length(simpler_times) - 1;

% read in observed data
fid = dir('*csv');
model_name = fid.name(6:length(fid.name)-4);
obs = dlmread(fid.name,';',1,1);
obs_ages = obs(:,1);
obs_sigma = obs(:,2);
all = textread(fid.name,'%s'); %#ok<DTXTRD>
for i = 1:length(all)-3
    values = strsplit(all{i+3},';');
    chron_names{i} = values{:,1};
end

% open all txt files and pull out the important information (i.e., 
% velocities and predicted ages)
fid = dir('ages_all*');
ages = dlmread(fid.name,' ');
[m,~] = find(ages <= 0);
m = unique(m);
ages(m,:) = [];

fid = dir('rates_all*');
rates = dlmread(fid.name,' ');
rates(m,:) = [];

% calculate chi squared metric
for i = 1:length(obs_ages)
    if sigma_perc_ages > 0
        chi_sim(:,i) = (obs_ages(i) - ages(:,i)).^2/(sigma_perc_ages*obs_ages(i)).^2;
    else
        chi_sim(:,i) = (obs_ages(i) - ages(:,i)).^2/obs_sigma(i).^2;
    end
end

% filter rates by chi squared
[m,~] = find(chi_sim > tolerance);
m = unique(m);
rates(m,:) = [];
ages(m,:) = [];

% calculate exhumed thicknesses based on the acceptable exhumation
% histories
if integration_time > 0
    start = find(simpler_times <= integration_time,1);
    for i = 1:length(rates(:,1))
        for j = start:num_slices
            E_hist(j) = rates(i,j)*(simpler_times(j) - simpler_times(j+1));
        end

        if integration_time > simpler_times(start)
            E_hist(start-1) = rates(i,start-1)*(integration_time - simpler_times(start));
        end
        exhumed_thickness(i) = sum(E_hist(:));
    end
    
    % filter simulations to fit observations of exhumed thickness
    if cull_data == 'y'
        m = find(exhumed_thickness < (depth + depth_sigma*2) & exhumed_thickness > (depth - depth_sigma*2));
        m = unique(m);
        rates_no_good = rates;
        rates_no_good(m,:) = [];
        ages_no_good = ages;
        ages_no_good(m,:) = [];
        rates = rates(m,:);
        ages = ages(m,:);
        exhumed_thickness = exhumed_thickness(m);
    end
end

% calculate mean and standard deviation of the predicted thicknesses
h_mean = mean(exhumed_thickness);
h_std = std(exhumed_thickness);

% plot observed and predicted age data and stats
for i = 1:length(chron_names)
    synth_ages = (randn(1,1e8).*obs_sigma(i))+obs_ages(i);
    figure(i)
    h1 = histogram(synth_ages,1e3);
    bin_centers = h1.BinEdges-(h1.BinWidth/2);
    bin_centers = bin_centers(1:end-1);
    bin_tops = h1.BinCounts;
    centers_save(i,:) = bin_centers; %#ok<*AGROW>
    tops_save(i,:) = bin_tops;
    
    hold off
    yyaxis left
    h2 = histogram(ages(:,i),round(length(ages)/2),'FaceColor','none','EdgeColor',[0 0.4470 0.7410],'LineWidth',1);
    hold on
    yyaxis left
    h = plot(bin_centers,bin_tops/max(bin_tops)*max(h2.BinCounts),'-k','LineWidth',2);
    xline(obs_ages(i)+obs_sigma(i),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
    xline(obs_ages(i)-obs_sigma(i),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
    xline(obs_ages(i)+2*obs_sigma(i),'LineWidth',0.5,'Color',[0.5 0.5 0.5],'LineStyle','--');
    xline(obs_ages(i)-2*obs_sigma(i),'LineWidth',0.5,'Color',[0.5 0.5 0.5],'LineStyle','--');
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
    txt = {[chron_names{i} ' (Ma)'] [num2str(obs_ages(i)) ' ' char(177) ' ' num2str(obs_sigma(i)*2) ' ' num2str(2) '\sigma (obs)'] [num2str(round(mean(ages(:,i)),1)) ' ' char(177) ' ' num2str(round(std(ages(:,1)),1)*2) ' ' num2str(2) '\sigma (pred)']};
    annotation('textbox',[0.125,0.725,0.22,0.2],'String',txt,'LineStyle','none','FontSize',8,'FontWeight','bold','FontName','Arial');
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
    xlabel('Cooling Age (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel('Count','FontSize',18,'FontWeight','bold','FontName','Arial')
    
    yyaxis right
    test_ages = h.Parent.XLim(1):(h.Parent.XLim(2)-h.Parent.XLim(1))/1e5:h.Parent.XLim(2);
    if sigma_perc_ages > 0 && obs_ages(i)*sigma_perc_ages > obs_sigma(i)
        chi_sigma = obs_ages(i)*sigma_perc_ages;
        chi_squared = (test_ages - obs_ages(i)).^2/chi_sigma;
        h = plot(test_ages,chi_squared,'.r','MarkerSize',3);
        ylabel(['\chi^2' ' (based on model uncertainty)'],'FontSize',18,'FontWeight','bold','FontName','Arial')
    else
        chi_squared = (test_ages - obs_ages(i)).^2/obs_sigma(i);
        h = plot(test_ages,chi_squared,'.r','MarkerSize',3);
        ylabel(['\chi^2' ' (based on observation uncertainty)'],'FontSize',18,'FontWeight','bold','FontName','Arial')
    end
    box on
    hold on
    location = find(chi_squared - tolerance < 0);
    plot([test_ages(location(1)),test_ages(location(end))],[chi_squared(location(1)),chi_squared(location(end))],'.r','MarkerSize',15)
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
    h.Parent.YColor = [0 0 0];
end

% calculate and plot exhumed thicknesses based on the acceptable exhumation
% histories
if integration_time > 0
    figure(length(obs_ages)+1)
    synth_depths = (randn(1,1e8).*depth_sigma)+depth;
    h1 = histogram(synth_depths,1e3);
    bin_centers = h1.BinEdges-(h1.BinWidth/2);
    bin_centers = bin_centers(1:end-1);
    bin_tops = h1.BinCounts;
    
    hold off
    yyaxis left
    h2 = histogram(exhumed_thickness,round(length(exhumed_thickness)/2),'FaceColor','none','EdgeColor',[0 0.4470 0.7410],'LineWidth',1);
    hold on
    yyaxis left
    h = plot(bin_centers,bin_tops/max(bin_tops)*max(h2.BinCounts),'-k','LineWidth',2);
    xline(depth+depth_sigma,'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
    xline(depth-depth_sigma,'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
    xline(depth+2*depth_sigma,'LineWidth',0.5,'Color',[0.5 0.5 0.5],'LineStyle','--');
    xline(depth-2*depth_sigma,'LineWidth',0.5,'Color',[0.5 0.5 0.5],'LineStyle','--');
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
    txt = {'Depth (km)' [num2str(depth) ' ' char(177) ' ' num2str(depth_sigma*2) ' ' num2str(2) '\sigma' ' (obs)'] [num2str(round(h_mean,1)) ' ' char(177) ' ' num2str(round(h_std,1)*2) ' ' num2str(2) '\sigma' ' (pred)']};
    annotation('textbox',[0.125,0.725,0.22,0.2],'String',txt,'LineStyle','none','FontSize',8,'FontWeight','bold','FontName','Arial');
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
    xlabel('Exhumed thickness (km)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel('Count','FontSize',18,'FontWeight','bold','FontName','Arial')
    
    yyaxis right
    test_depths = h.Parent.XLim(1):(h.Parent.XLim(2)-h.Parent.XLim(1))/1e5:h.Parent.XLim(2);
    if sigma_perc_depth > 0 && depth*sigma_perc_depth > depth_sigma
        chi_sigma = depth*sigma_perc_depth;
        chi_squared = (test_depths - depth).^2/chi_sigma;
        h = plot(test_depths,chi_squared,'.r','MarkerSize',3);
        ylabel(['\chi^2' ' (based on model uncertainty)'],'FontSize',18,'FontWeight','bold','FontName','Arial')
    else
        chi_squared = (test_depths - depth).^2/depth_sigma;
        h = plot(test_depths,chi_squared,'.r','MarkerSize',3);
        ylabel(['\chi^2' ' (based on observation uncertainty)'],'FontSize',18,'FontWeight','bold','FontName','Arial')
    end
    box on
    hold on
    location = find(chi_squared - tolerance < 0);
    plot([test_depths(location(1)),test_depths(location(end))],[chi_squared(location(1)),chi_squared(location(end))],'.r','MarkerSize',15)
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
    h.Parent.YColor = [0 0 0];
end

%plot exhumation history data
figure(length(obs_ages)+2)
mean_rate = mean(rates);
std_rate = std(rates);
for m = 1:length(mean_rate)
    new_mean(2*m-1) = mean_rate(m);
    new_mean(2*m) =  mean_rate(m);
    new_std(2*m-1) = std_rate(m);
    new_std(2*m) =  std_rate(m);
end

if time_fudge > 0
    new_mean = [new_mean 0 0];
    new_std = [new_std 0 0];
    times = [times times(end) 0];
end

yyaxis left
h = plot(times, new_mean,'-k','LineWidth',2);
hold on
box on
plot(times, new_mean-new_std,'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-')
plot(times, new_mean+new_std,'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
ylabel('Exhumation rate (km/Myr)','FontSize',18,'FontWeight','bold','FontName','Arial')
h.Parent.YColor = [0 0 0];
ylim([0 (max(new_mean+new_std)*1.4)])

yyaxis right
cmap = jet(length(chron_names));
set(gca,'Ydir','reverse')
ylim([0 16])
for i = 1:length(chron_names)
    for j = 1:length(centers_save(i,:))
        fade = tops_save(i,:)/max(tops_save(i,:));
        h = plot([centers_save(i,j) centers_save(i,j)],[i i+0.4],'-','LineWidth',1,'Color',fade(j)*cmap(i,:)+(1-fade(j))*[1 1 1]);
        hold on
        drawnow
    end
    
    hold on
    y_plot = ones(1,length(ages(:,i))).*(i+0.2);
    plot(ages(:,i),y_plot,'o','Color',[0 0 0],'MarkerSize',10,'LineWidth',1)
    set(gca,'ytick',[])
end
h.Parent.YColor = [0 0 0];
xlim([0 max(times)])

if cull_data == 'y'
    figure(length(obs_ages)+3)
    [dim_1,dim_2] = size(rates_no_good);
    for l = 1:dim_1
        clear new_rates
        for m = 1:dim_2
            new_rates(2*m-1) = rates_no_good(l,m);
            new_rates(2*m) =  rates_no_good(l,m);
        end
            
        if time_fudge > 0
            new_rates = [new_rates 0 0];
        end
        
        yyaxis left
        h = plot(times,new_rates,'-','LineWidth',1);
        h.Color = [0 0 0 0.1];
        hold on
        box on
        set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
        xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
        ylabel('Exhumation rate (km/Myr)','FontSize',18,'FontWeight','bold','FontName','Arial')
        h.Parent.YColor = [0 0 0];
        ylim([0 inf])
        xlim([0 inf])
    end
    
    figure(length(obs_ages)+3)
    [dim_1,dim_2] = size(rates);
    for l = 1:dim_1
        clear new_rates
        for m = 1:dim_2
            new_rates(2*m-1) = rates(l,m);
            new_rates(2*m) =  rates(l,m);
        end
            
        if time_fudge > 0
            new_rates = [new_rates 0 0];
        end
        
        yyaxis left 
        h = plot(times, new_rates,'-','LineWidth',1);
        h.Color = [1 0 0 0.5];
        hold on
        box on
        set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
        xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
        ylabel('Exhumation rate (km/Myr)','FontSize',18,'FontWeight','bold','FontName','Arial')
        h.Parent.YColor = [0 0 0];
        ylim([0 inf])
        xlim([0 inf])
    end
    
    yyaxis right
    hp = [];
    cmap = jet(length(chron_names));
    [dim,~] = size(ages_no_good);
    for i = 1:length(chron_names)
        if strcmp(err_disp,'pdf') == 1
            h = plot(centers_save(i,:),tops_save(i,:)/max(tops_save(i,:)),'LineWidth',1,'LineStyle','-','Color',cmap(i,:));
        else
            xline(obs_ages(i)+obs_sigma(i)*2,'LineWidth',1,'Color',cmap(i,:));
            xline(obs_ages(i)-obs_sigma(i)*2,'LineWidth',1,'Color',cmap(i,:));
        end
        hold on
        my_field = strcat('p',num2str(i));
        variable.(my_field) = plot(ages_no_good(:,i),(1:1:dim)./dim./40+0.1,'o','MarkerSize',3,'Color',cmap(i,:),'DisplayName',chron_names{i});
        hp = [hp variable.(my_field)];
        set(gca,'ytick',[])
    end
    
    yyaxis right
    hp = [];
    cmap = jet(length(chron_names));
    [dim,~] = size(ages);
    for i = 1:length(chron_names)
        hold on
        my_field = strcat('p',num2str(i));
        variable.(my_field) = plot(ages(:,i),(1:1:dim)./dim./40+0.01,'.','MarkerSize',10,'Color',cmap(i,:),'DisplayName',chron_names{i});
        hp = [hp variable.(my_field)];
        set(gca,'ytick',[])
    end
    h.Parent.YColor = [0 0 0];
    ylim([0 1])
    legend(hp,'Location','northwest')
end

% save figures
if save_figs == 'y'
	for i = 1:length(obs_ages)+2
        figure(i)
        exportgraphics(gcf,['fig_' num2str(i) '.eps'],'BackgroundColor','none','ContentType','vector')
        saveas(i,['fig_' num2str(i)],'png')
	end
end