function [frequencies,r,theta] = main_ZI_sweep(sigout_ch,amplitude,num_points,TC_settle,settle_time,start_f,stop_f,samp_avg, ...
    scan_order,aux1_out,aux2_out,aux3_out,settings_filename,filter_order,filter_BW,samp_rate,diff_sense,sense_range)
%Run a frequency sweep using 1 excitation signal and 2 input demodulators 

%connect to device
device = 'dev447';
%device = 'dev1418';
clear ziDAQ;

%{
if ~exist('device_id', 'var')
    error(['No value for device_id specified. The first argument to the ' ...
           'example should be the device ID on which to run the example, ' ...
           'e.g. ''dev2006'' or ''uhf-dev2006''.'])
end
%}

% Check the ziDAQ MEX (DLL) and Utility functions can be found in Matlab's path.
if ~(exist('ziDAQ') == 3) && ~(exist('ziDevices', 'file') == 2)
    fprintf('Failed to either find the ziDAQ mex file or ziDevices() utility.\n')
    fprintf('Please configure your path using the ziDAQ function ziAddPath().\n')
    fprintf('This can be found in the API subfolder of your LabOne installation.\n');
    fprintf('On Windows this is typically:\n');
    fprintf('C:\\Program Files\\Zurich Instruments\\LabOne\\API\\MATLAB2012\\\n');
    return
end

%{

% Determine the device identifier from it's ID.
device = lower(ziDAQ('discoveryFind', device_id));
% Get the device's default connectivity properties.
props = ziDAQ('discoveryGet', device);
% The maximum API level supported by this example.
apilevel_example = 5;
% The maximum API level supported by the device class, e.g., MF.
apilevel_device = props.apilevel;
% Ensure we run the example using a supported API level.
apilevel = min(apilevel_device, apilevel_example);
% See the LabOne Programming Manual for an explanation of API levels.

%}
% Create a connection to a Zurich Instruments Data Server (a API session)
% using the device's default connectivity properties.


%serverHostname = 'boyd-labpc.eecs.umich.edu';
serverHostname = 'localhost';
serverPort = 8005;
apilevel = 1; 
%ziDAQ('connect', props.serveraddress, props.serverport, apilevel);
ziDAQ('connect',serverHostname,serverPort,apilevel);

% Check that the device is visible to the Data Server.
if ~ismember(device, ziDevices())
    message = ['The specified device `', device, '` is not visible to the Data ', ...
               'Server. Please ensure the device is connected by using the LabOne ', ...
               'User Interface or ziControl (HF2 Instruments).'];
    error(message);
end

% Get the device type and its options (in order to set correct device-specific
% configuration).
devtype = ziDAQ('getByte', ['/' device '/features/devtype']);
options = ziDAQ('getByte', ['/' device '/features/options']);

fprintf('Will run frequency sweep on `%s`, an `%s` with options `%s`.\n', device, ...
        devtype, regexprep(options, '\n' , '|'));

settings = struct('sigout_ch',sigout_ch,'amplitude',amplitude,'num_points',num_points, ...
    'TC_settle',TC_settle,'settle_time',settle_time,'start_f',start_f,'stop_f',stop_f,'samp_avg',samp_avg,...
    'scan_order',scan_order,'aux1_out',aux1_out,'aux2_out',aux2_out,'aux3_out',aux3_out,'settings_filename',settings_filename ...
    ,'filter_order',filter_order,'filter_BW',filter_BW,'samp_rate',samp_rate,'diff_sense',diff_sense,'sense_range',sense_range);      
    
[frequencies,r,theta] = ZI_sweep(device,settings);

end

function [frequencies,r,theta] = ZI_sweep(device,settings)

%isposscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
%isone_or_two = @(x) isposscalar(x) && ((x==1)||(x==2));
%p = inputParser;
%{
p.addParameter('sigout_ch',1,isposscalar);
p.addParameter('sweep_samplecount', 512, isposscalar); %default 512 points
p.addParameter('amplitude', 100e-3, isposscalar); %default 100 mV amplitude
disp(varargin{:});
p.parse(varargin{:});
%}


if(settings.sigout_ch == 1)  
    demod_idx = [1, 2]; % 1-based indexing
    out_c = '0'; % signal output channel
    out_mixer_c = '6';
    osc_c = '0'; % oscillator
else
    demod_idx = [4, 5]; % 1-based indexing
    out_c = '1'; % signal output channel
    out_mixer_c = '7';
    osc_c = '1'; % oscillator
end

demod_c = cell(1,2);
for d=1:length(demod_idx)
  % demod_c, demod channel, 0-based indexing for paths on the device
  demod_c{1,d} = num2str((demod_idx(d)-1));
end

ziLoadSettings(device, settings.settings_filename);
%ziLoadSettings(device, 'base_sweep_settings_12_22_16.xml');

for d=1:length(demod_idx)
   ziDAQ('setInt',['/' device '/demods/' demod_c{1,d} '/order'],ceil(settings.filter_order/6));
   ziDAQ('setDouble',['/' device '/demods/' demod_c{1,d} '/timeconstant'], compute_TC_per_stage(ceil(settings.filter_order/6),settings.filter_BW));
   ziDAQ('setDouble',['/' device '/demods/' demod_c{1,d} '/rate'],settings.samp_rate);                                                   
end

ziDAQ('setDouble', ['/' device '/auxouts/0/offset'], settings.aux1_out); %quadrature tuning voltage setting (in-axis)
ziDAQ('setDouble', ['/' device '/auxouts/1/offset'], settings.aux2_out); %(cross-axis)
ziDAQ('setDouble', ['/' device '/auxouts/2/offset'], settings.aux3_out); %(quadrature tuning)

%ziLoadSettings(device, '1_24_17_base_sweep_quad_n3p2V.xml');

%set input sensing parameters
ziDAQ('setInt', ['/' device '/sigins/0/diff'], settings.diff_sense);
ziDAQ('setInt', ['/' device '/sigins/1/diff'], settings.diff_sense);
ziDAQ('setDouble', ['/' device '/sigins/0/range'], settings.sense_range);
ziDAQ('setDouble', ['/' device '/sigins/1/range'], settings.sense_range);

%set amplitude and enable the output 
ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/amplitudes/' out_mixer_c], settings.amplitude);
ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/enables/' out_mixer_c], 1);
ziDAQ('setInt', ['/' device '/sigouts/' out_c '/on'], 1);

%create sweeper object
%% Sweeper settings
%============================================================
% Create a thread for the sweeper 
timeout = 500; % milliseconds
h = ziDAQ('sweep');
% Device on which sweeping will be performed
ziDAQ('set', h, 'sweep/device', device);
% Sweeping setting is the frequency of the output signal
ziDAQ('set', h, 'sweep/gridnode', ['oscs/' osc_c '/freq']);
% Start frequency 
ziDAQ('set', h, 'sweep/start', settings.start_f);
% Stop frequency
ziDAQ('set', h, 'sweep/stop', settings.stop_f);
%number of sweep points
ziDAQ('set', h, 'sweep/samplecount', settings.num_points);

%===================================================================
%other sweeper settings
% Single sweep 
ziDAQ('set', h, 'sweep/loopcount', 1);
% linear sweep mode
ziDAQ('set', h, 'sweep/xmapping', 0);
% sequential scan type/binary/bidirectional/reverse 
ziDAQ('set', h, 'sweep/scan', settings.scan_order);

ziDAQ('set', h, 'sweep/bandwidthcontrol', 0);
ziDAQ('set', h, 'sweep/bandwidth', 11.5); %this value should be ignored in manual BWcontrol mode

%=================================================================
%settling and averaging time
ziDAQ('set', h, 'sweep/settling/time', settings.settle_time);
ziDAQ('set', h, 'sweep/settling/tc',settings.TC_settle);

ziDAQ('set', h, 'sweep/averaging/tc',0); %no averaging
ziDAQ('set', h, 'sweep/averaging/sample',settings.samp_avg) %no averaging 
%====================================================================

ziDAQ('subscribe', h, ['/' device '/demods/' demod_c{1,1} '/sample']);
ziDAQ('subscribe', h, ['/' device '/demods/' demod_c{1,2} '/sample']);

% Start sweeping
ziDAQ('execute', h);

data = [];
if(settings.scan_order == 2)
    frequencies = nan(1, (1*settings.num_points));
    r = nan(length(demod_c),(1*settings.num_points));
    theta = nan(length(demod_c), (1*settings.num_points));
else
    frequencies = nan(1, settings.num_points);
    r = nan(length(demod_c),settings.num_points);
    theta = nan(length(demod_c), settings.num_points);    
end
%clf;

figure(settings.sigout_ch); clf; %changed code here! from figure(1)
timeout = 7200; %timeout (1200)
t0 = tic;

% Read and plot intermediate data until the sweep has finished.
while ~ziDAQ('finished', h)
    pause(10);
    tmp = ziDAQ('read', h);
    fprintf('Sweep progress %0.0f%%\n', ziDAQ('progress', h) * 100);
    fprintf('Remaining Time: %g minutes %g seconds\n', floor(tmp.remainingtime/60) , floor(mod(tmp.remainingtime,60)));
    
    % Using intermediate reads we can plot a continuous refinement of the ongoing
    % measurement. If not required it can be removed.
    if ziCheckPathInData(tmp, ['/' device '/demods/' demod_c{1,1} '/sample'])
        if ~isempty(tmp.(device).demods(demod_idx(1)).sample{1}) && ...
                ~isempty(tmp.(device).demods(demod_idx(2)).sample{1})
            data = tmp;

            % Get the magnitude and phase of demodulator from the sweeper result.
            sample = data.(device).demods(demod_idx(1)).sample{1};
            r(1, :) = sample.r;
            theta(1, :) = sample.phase;

            sample = data.(device).demods(demod_idx(2)).sample{1};
            r(2, :) = sample.r;
            theta(2, :) = sample.phase;

            % Frequency values at which measurement points were taken
            frequencies = sample.grid;
            valid = ~isnan(frequencies);
            plot_data(frequencies(valid), r(:, valid), theta(:, valid),settings.amplitude, '-',settings.sigout_ch,0)
            drawnow;
        end
    end
    if toc(t0) > timeout
        error('Timeout: Sweeper failed to finish after %f seconds.', timeout)
    end
end

fprintf('Sweep time %f seconds\n', toc(t0));

% now read the data. This command can also be executed during the waiting.
tmp = ziDAQ('read', h);

% unsubscribe from the node; stop filling the data from that node to the
% internal buffer in the server
ziDAQ('unsubscribe', h, ['/' device '/demods/*/sample']);

ziDAQ('setInt', ['/' device '/sigouts/' out_c '/on'], 0);

% Process any remainging data returned by read().
if ziCheckPathInData(tmp, ['/' device '/demods/' demod_c{1,1} '/sample'])
    if ~isempty(tmp.(device).demods(demod_idx(1)).sample{1}) && ...
            ~isempty(tmp.(device).demods(demod_idx(2)).sample{1})
        data = tmp;

        % Get the magnitude and phase of demodulator from the sweeper result.
        sample = data.(device).demods(demod_idx(1)).sample{1};
        r(1, :) = sample.r;
        theta(1, :) = sample.phase;

        sample = data.(device).demods(demod_idx(2)).sample{1};
        r(2, :) = sample.r;
        theta(2, :) = sample.phase;

        % Frequency values at which measurement points were taken
        frequencies = sample.grid;
        valid = ~isnan(frequencies);
        plot_data(frequencies(valid), r(:, valid), theta(:, valid), settings.amplitude, '-',settings.sigout_ch,1)
        drawnow;
    
    end
end

end

function plot_data(frequencies, r, theta, amplitude, style,sigout_ch,legend_on)

if(sigout_ch == 2)
   colors = ['g' 'm']; 
else
   colors = ['b' 'r'];
end

% Plot data
figure(sigout_ch); 
clf
subplot(2, 1, 1)
s = plot(frequencies, 20*log10(r(1, :)), style);
hold on;
set(s, 'LineWidth', 2)
set(s, 'Color', colors(1));
s = plot(frequencies, 20*log10(r(2, :)), style);
set(s, 'LineWidth', 2)
set(s, 'Color', colors(2));
grid on

if(legend_on)
legend('in axis', 'cross axis');
end

xlabel('Frequency [Hz]')
ylabel('Amplitude [dBVrms]')

subplot(2, 1, 2)
s = plot(frequencies, theta(1, :)*180/pi, style);
hold on;
set(s, 'LineWidth', 2)
set(s, 'Color', colors(1));
s = plot(frequencies, theta(2, :)*180/pi, style);
set(s, 'LineWidth', 2)
set(s, 'Color', colors(2));
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
end