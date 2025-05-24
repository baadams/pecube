function bring_the_heat(varargin)

%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: Mar 2024
%-------------------------------------------------------------------------%
%
% Description:
% Function extracts and plots data from Pecube-D 1D output files (see  
% required files below). This function only works with output data from a 
% single Pecube simulation. It is not designed to work with amalgamated 
% monte carlo data. The script expects the current directory to be the main
% simulation directory containing the Pecube.in file and output directory.
%
% Usage: 
% bring_the_heat(); (uses all default values below)
% bring_the_heat('target_depth',10,'iso_temp',120,'rate',2,'fig_num',4);
% 
% Required input files:
% Pecube.in - this is the input file for the monte carlo run
% ages_XXXX.csv - this is the csv file containing observed cooling ages.
%                 this needs to be formatted for input into Pecube.
% Ages_tecXXXX.dat (one file for each time card)
% Temps_tecXXXX.dat (one file for each time card)
% time_temperature_history_XXXX.txt (one for each time card. **NOTE** these
%                                   cannot be the .bin files!)
%
% Optional Inputs:
% target_depth - the depth used to calculate the shallow geothermal
%                gradient (km)
% iso_temp - temperature of an isotherm of interest, the depth to this
%            isotherm will be monitored over time.
% rate - the frame rate of the output movie, if one is to be written (fps)
% fig_num - the figure number of figure that is to be saved as a movie. an
%           entry of 0 will not save a movie.
% save_figs - flag to save figures as eps and png. 'y' for yes.
% video_format - the format of the output video. Options: 'avi' or 'mp4'
%
% Outputs:
% figures and a movie file if one is desired
% 
%-------------------------------------------------------------------------%
% tashi delek!
%-------------------------------------------------------------------------%

% parse inputs
p = inputParser;         
p.FunctionName = 'bring_the_heat';
addParameter(p,'target_depth',2,@(x) isscalar(x));
addParameter(p,'iso_temp',70,@(x) isscalar(x));
addParameter(p,'rate',10,@(x) isscalar(x));
addParameter(p,'save_figs','y',@(x) ischar(x));
addParameter(p,'video_format','mp4',@(x) ischar(x));
    
parse(p,varargin{:});
target_depth = p.Results.target_depth;
iso_temp = p.Results.iso_temp;
rate = p.Results.rate;
save_figs = p.Results.save_figs;
video_format = p.Results.video_format;

% create video file and set frame rate
if strcmp(video_format,'avi')
    writerObj = VideoWriter('pecube_movie.avi');
elseif strcmp(video_format,'mp4')
    writerObj = VideoWriter('pecube_movie.mp4','MPEG-4');
end
writerObj.FrameRate = rate;
open(writerObj);
 
% open the input file and pull out the important information (i.e., model 
% depth, number of z nodes, times, temp at the surface)   
file_id = fopen('Pecube.in');
line = fgets(file_id);
data = [];
while ischar(line)
    if contains(line,'values below should agree with what is specified in Input 5.')
        line = fgets(file_id); %#ok<NASGU>
        line = fgets(file_id);
        stuff = split(line(1:end));
        stuff = str2double(stuff);
        dimx = stuff(1);
        dimy = stuff(2);
    elseif contains(line,'final model condition')
        line = fgets(file_id); %#ok<NASGU>
        line = fgets(file_id);
        num_t_steps = str2double(line(1:end))+1;
	elseif contains(line,'i j k')
        line = fgets(file_id);
        for i = 1:num_t_steps
            stuff = split(line(1:end));
            test  = str2double(stuff);
            data = [data;transpose(test(1:6))];
            line = fgets(file_id);
        end
	elseif contains(line,'No brittle shear heating')
        line = fgets(file_id); %#ok<NASGU>
        line = fgets(file_id);
        stuff = regexp(line(1:end),' ','split');
        stuff = str2double(stuff);
        model_depth = stuff(1);
        z_nodes = stuff(2);
        conductivity = stuff(3);
        line = fgets(file_id);
        stuff = regexp(line(1:end),' ','split');
        stuff = str2double(stuff);
        T_surf = stuff(2);
    end
    line = fgets(file_id);
end
fclose(file_id);

% calculate node depths
time = transpose(data(:,1));
z = 0:(model_depth/(z_nodes-1)):model_depth;
dz = abs(z(1) - z(2));

% pull the exhumation rates from the erates files
cd output
E_files = dir('time*.txt');
for i = 1:num_t_steps-1
	E_names{i} = E_files(i,1).name; %#ok<*SAGROW>
	rates = dlmread(E_names{i},',',1,1); %#ok<*DLMRD>
    E(i) = rates(5,11);
end
E = [E(1) E];

% calculate the sample depth history
step_dur(1) = 0;
delta_z(1) = 0;
sam_depth(1) = 0;
for k = 2:length(time)
    step_dur(k) = time(k-1) - time(k);
    delta_z(k) = fliplr(E(k).*step_dur(k));
    sam_depth(k) = sam_depth(k-1) + delta_z(k);
end

% pull the temperature histories from the samples at the surface at the
% at the end of the model
hist_files = dir('time*');
hist_files.name;
filename = ans; %#ok<NOANS>
file_id = fopen(filename);
line = fgets(file_id); %#ok<NASGU>
line = fgets(file_id);
data = [];
while ischar(line)
    stuff = split(line(1:end));
    data = [data;transpose(str2double(stuff))]; %#ok<*AGROW>
    line = fgets(file_id);
end
data = sortrows(data,4,'ascend');
keep_rows = 1:4:length(data);
keep_cols = [4,6,9,12];
data = data(keep_rows,keep_cols);
sam_time = data(:,1);
delta_s_t = round(sam_time(2) - sam_time(1),1);
sam_time = [0;(sam_time + delta_s_t)];
sam_temp = data(:,2);
sam_temp = [T_surf;sam_temp];

% pull the geothermal gradients from each time step
T_files = dir('Temps*');
for i = 1:num_t_steps
	T_names{i} = T_files(i,1).name; %#ok<*SAGROW>
	temps = dlmread(T_names{i},'',4,0);
	for j = 1:z_nodes
    	T(j,i) = temps(j,5);
	end
end
T = flipud(T);

% open txt files and pull out the predicted ages
A_files = dir('Ages_tec*');
fid = A_files(num_t_steps-1).name;
ages = dlmread(fid,'',7,7);
pred_ages = ages(1,:);
cd ../

% read in observed data
fid = dir('*csv');
model_name = fid.name(6:length(fid.name)-4);
obs = dlmread(fid.name,';',1,1);
obs_ages = obs(:,1);
obs_sigma = obs(:,2);

% for each predicted chronometer calculate the predicted closure temp
for i = 1:length(pred_ages)
    location = find(sam_time < pred_ages(i),1,'last');
    diff_t = pred_ages(i) - sam_time(location);
    delta_t = abs(sam_time(location) - sam_time(location+1));
    rat = diff_t/delta_t;
    delta_T = abs(sam_temp(location) - sam_temp(location+1));
    pred_Tcb(i) = rat*delta_T + sam_temp(location);
end

% calculate the shallow geothermal gradient and the depth to the isotherm
% of interest
diff = abs(z - target_depth);
G_base = find(diff == min(diff));
for k = 1:num_t_steps
    if rem(time(k),1) == 0 || time(k) == 0
	    location = find(T(:,k) < iso_temp,1,'last');
	    diff_T = iso_temp - T(location,k);
	    delta_T = abs(T(location,k) - T(location+1,k));
	    rat = diff_T/delta_T;
	    iso_depth(k) = rat*dz + z(location);
        G(k) = T(G_base,k)/target_depth; %#ok<FNDSB>
        time_new(k) = time(k);
    end
end
G(G==0) = [];
iso_depth(iso_depth==0) = [];
time_new(time_new==0) = [];
time_new = [time_new 0];
heat_flux = conductivity*G*dimx*dimy*1e6*31556952000000;
for k = 1:length(G)
    if k == 1
        cum_flux(k) = heat_flux(k);
    else
        cum_flux(k) = cum_flux(k-1) + heat_flux(k);
    end
end
    
% ***MODEL HISTORY DATA***
for j = 1:num_t_steps 
    if rem(time(j),0.1) == 0 || time(j) == 0
	    figure(1) % temperature vs depth
	    h = plot(T(:,j),z,'-k','LineWidth',1);
	    box on
	    set(gca,'Ydir','reverse')
	    set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
	    xlabel(['Temperature (' char(176) 'C)'],'FontSize',18,'FontWeight','bold','FontName','Arial')
	    ylabel('Depth below surface (km)','FontSize',18,'FontWeight','bold','FontName','Arial')
        title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
        txt = [num2str(round(time(j),3)) ' Ma'];    
        an = annotation('textbox',[0.125,0.025,0.22,0.2],'String',txt,'LineStyle','none','FontSize',12,'FontWeight','bold','FontName','Arial');
        y_tick = h.Parent.YTick(2) - h.Parent.YTick(1);
        z_max = ceil(model_depth/y_tick)*y_tick;
        x_tick = h.Parent.XTick(2) - h.Parent.XTick(1);
        T_max = ceil(max(max(T))/x_tick)*x_tick;
 	    ylim([0 z_max])
	    xlim([0 T_max])
                
	    % capture the frame to make a movie
        pecube_movie(j) = getframe(figure(1));
        writeVideo(writerObj,pecube_movie(j))
        delete(an)
    end
end
close(writerObj);

figure(2) % exhumation rate vs time
	hold on
	h = plot(time,E,'-k','LineWidth',1.5);
	box on
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
	set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
	xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
 	ylabel('Exhumation rate (km/Myr)','FontSize',18,'FontWeight','bold','FontName','Arial')
    y_tick = h.Parent.YTick(2) - h.Parent.YTick(1);
    E_max = ceil(max(E)/y_tick)*y_tick;
    x_tick = h.Parent.XTick(2) - h.Parent.XTick(1);
    t_max = ceil(max(time)/x_tick)*x_tick;
 	ylim([0 E_max])
	xlim([0 t_max])

figure(3) % isotherm depth vs time
	hold on
	h = plot(time_new,iso_depth,'-k','LineWidth',1.5);
	box on
	set(gca,'Ydir','reverse')
	set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
	xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel('Depth below surface (km)','FontSize',18,'FontWeight','bold','FontName','Arial')
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
    txt = [num2str(iso_temp) char(176) 'C isotherm'];    
    annotation('textbox',[0.125,0.725,0.22,0.2],'String',txt,'LineStyle','none','FontSize',12,'FontWeight','bold','FontName','Arial');
    y_tick = h.Parent.YTick(2) - h.Parent.YTick(1);
    iso_max = ceil(max(iso_depth)/y_tick)*y_tick;
    x_tick = h.Parent.XTick(2) - h.Parent.XTick(1);
    t_max = ceil(max(time_new)/x_tick)*x_tick;
	ylim([0 iso_max])
	xlim([0 t_max])

figure(4) % shallow geothermal gradient vs time
	hold on
	h = plot(time_new,G,'-k','LineWidth',1.5);
	box on
	set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
	xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel(['Geothermal gradient (' char(176) 'C/km)'],'FontSize',18,'FontWeight','bold','FontName','Arial')
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
    txt = [num2str(target_depth) '-km geotherm'];    
    annotation('textbox',[0.125,0.725,0.3,0.2],'String',txt,'LineStyle','none','FontSize',12,'FontWeight','bold','FontName','Arial');
	y_tick = h.Parent.YTick(2) - h.Parent.YTick(1);
    G_max = ceil(max(G)/y_tick)*y_tick;
    x_tick = h.Parent.XTick(2) - h.Parent.XTick(1);
    t_max = ceil(max(time_new)/x_tick)*x_tick;
	ylim([0 G_max])
	xlim([0 t_max])

figure(5) % cummulative energy vs time
	hold on
	h = plot(time_new,cum_flux,'-k','LineWidth',1.5);
	box on
	set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
	xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel('Cummulative advected energy (J)','FontSize',18,'FontWeight','bold','FontName','Arial')
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
    power = ceil(log10(cum_flux(end))-1);
    txt = {'Total advected energy' [num2str(round(cum_flux(end)/1^power),4) ' (J)']};    
    annotation('textbox',[0.125,0.1,0.4,0.2],'String',txt,'LineStyle','none','FontSize',12,'FontWeight','bold','FontName','Arial');
	y_tick = h.Parent.YTick(2) - h.Parent.YTick(1);
    J_max = ceil(max(cum_flux)/y_tick)*y_tick;
    x_tick = h.Parent.XTick(2) - h.Parent.XTick(1);
    t_max = ceil(max(time_new)/x_tick)*x_tick;
	ylim([0 J_max])
	xlim([0 t_max])

figure(6) % sample depth and temperature vs time
	hold on    
	yyaxis left
	h1 = plot(sam_time,sam_temp,'-b','Linewidth',1.5,'DisplayName','Temperature');
    plot(pred_ages,pred_Tcb,'.','MarkerSize',15,'MarkerEdgeColor','b','MarkerFaceColor','b','DisplayName','Predicted')
    e = errorbar(obs_ages,pred_Tcb,obs_sigma*2,'horizontal','.','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Observed');
    e.Color = 'k';
    e.CapSize = 5;
    e.LineWidth = 1;
	set(gca,'Ydir','reverse')
	box on
	xlabel('Time (Ma)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel(['Temperature (' char(176) 'C)'],'FontSize',15,'FontWeight','bold','FontName','Arial')
    title(model_name,'FontSize',18,'FontWeight','bold','FontName','Arial')
	y_tick = h1.Parent.YTick(2) - h1.Parent.YTick(1);
    temp_max = ceil(max(sam_temp)/y_tick)*y_tick;
    x_tick = h1.Parent.XTick(2) - h1.Parent.XTick(1);
    t_max = ceil(max(sam_time)/x_tick)*x_tick;
    h1.Parent.YColor = [0 0 0];
	ylim([0 temp_max])
	xlim([0 t_max])
    
	yyaxis right
	h2 = plot(time,abs(sam_depth-max(sam_depth)),'-k','Linewidth',1.5,'DisplayName','Depth');
	set(gca,'Ydir','reverse')
	set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial')
	ylabel('Depth below surface (km)','FontSize',18,'FontWeight','bold','FontName','Arial')
	y_tick = h2.Parent.YTick(2) - h2.Parent.YTick(1);
    depth_max = ceil(max(sam_depth)/y_tick)*y_tick;
    x_tick = h2.Parent.XTick(2) - h2.Parent.XTick(1);
    t_max = ceil(max(time)/x_tick)*x_tick;
    h2.Parent.YColor = [0 0 0];
	ylim([0 depth_max])
	xlim([0 t_max])

% save figures
if save_figs == 'y'
	for i = 2:6
        saveas(i,['fig_' num2str(i)],'epsc')
        saveas(i,['fig_' num2str(i)],'png')
        saveas(i,['fig_' num2str(i)],'fig')
	end
end