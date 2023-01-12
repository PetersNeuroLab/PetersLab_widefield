function [U,Vrec,im_avg_color,frame_info] = AP_preprocess_widefield_hc(im_path)
% [U,Vrec,im_color_avg,frame_info] = AP_preprocess_widefield_pco(im_path)
%
% SVD-compress widefield imaging from PCO Edge 5.5 camera
% Assumes: 
% - alternating 2-color imaging that resets order on recording start
% - recordings are defined by timestamp gaps of >2s
% - binary and ASCII timestamps are turned on


verbose = false; % Turn off messages by default


%% Get image filenames

im_files = dir(fullfile(im_path,'*.tif'));


%% Get information from all frames
% PCO timestamps in binary coded decimal (BCD) as:
% Pixels code: 
% 1:4          5 6 7 8 9 10 11 12    13    14
% framenum(x4) Y Y M D H M  S  10kus 100us us
% BCD converts to decimal as: a decimal number that codes for one byte, 
% with two nibbles that are one decimal each
% (e.g. 10dec = 16bcd: 16 -> 0001 0000 = "1,0" = 10)

tic;
if verbose; disp('Getting image headers...'); end

% % Loop through all files, get frame number and timestamp
% frame_info = struct('frame_num',cell(size(im_files)),'time_from_last',cell(size(im_files)), 'time_from_start',cell(size(im_files)));

for curr_im=1:length(im_files)
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    file_info = imfinfo(curr_im_fn);
    
    % get file metadata
    file_md = {file_info.ImageDescription};

    % using the first one, save out the information for the whole recording
    frame_md = file_md{1};

    % find indexes for the title and the content
    all_title_beg_idx = strfind(frame_md, '[');
    all_title_end_idx = strfind(frame_md, ']');
    all_section_beg_idx = all_title_end_idx+1;
    all_section_end_idx = [all_section_beg_idx(2:end)-1, length(frame_md)];   

    % go through each section and pull content
    for n_section=1:length(all_section_beg_idx)-1 % exclude the last section
        
        % get section title and get rid of extra characters
        title = frame_md(all_title_beg_idx(n_section):all_title_end_idx(n_section));
        title = lower(strrep(erase(title,{'[ ', ' ]'}),' ', '_'));
    
        % get section contents
        section_str = frame_md(all_section_beg_idx(n_section):all_section_end_idx(n_section));
    
        % split by .. = ..
        name_pattern = '(?<name>(\w+).*?)';
        values_pattern = '(?<value>.*?(\w+).*?\n)';
        expr = [name_pattern ' = ' values_pattern];
        test = regexp(section_str, expr, 'names');
    
        % save in frame_info 
        for val_idx=1:length(test)
            name = strrep(test(val_idx).name,' ', '_');
            value = test(val_idx).value;
    %         value = erase(value, '\b\n')
            frame_info(curr_im).(title).(name) = value(1:end-2); % line above doesn't work so had to do this
        end
    end

    % get timestamp
    timestamp_section = frame_md(1:all_title_beg_idx(1)-1);
    expr = '\d*.[a-z]+.\d*.\d*:\d*:\d*';
    timestamp_str = regexpi(timestamp_section, expr, 'match');
    frame_info(curr_im).timestamp = repmat(datetime(timestamp_str,...
        'InputFormat','dd MMM yyyy HH:mm:ss', ...
        'Format', 'yyyy-MMM-dd HH:mm:ss.SSSS'), ...
        1, length(file_info));

    % get times for all frames
    name_pattern = '(?<name>Time_From_+(Start|Last))';
    time_pattern = '(?<time>\d\d:\d\d:\d\d.\d\d\d\d)';
    expr = [name_pattern ' = ' time_pattern];
    capture_times = regexp(file_md, expr, 'names');

    % convert both time from start and time from last to seconds
    F = 'hh:mm:ss.SSSS';
    time_from_last = cellfun(@(X) ...
        seconds(duration(X(strcmp({X.name},'Time_From_Last')).time, ...
        'InputFormat', F, 'Format', F)), capture_times);
    time_from_start = cellfun(@(X) ...
        seconds(duration(X(strcmp({X.name},'Time_From_Start')).time, ...
        'InputFormat', F, 'Format', F)), capture_times);

    % Store info for file
    frame_info(curr_im).frame_num(:) = 1:length(file_info);
    frame_info(curr_im).time_from_last(:) = time_from_last;
    frame_info(curr_im).time_from_start(:) = time_from_start;
    frame_info(curr_im).timestamp.Second = frame_info(curr_im).timestamp.Second + time_from_start;

    % --------------- temp took out because I don't have it
    % --------------------------------------------------------------------------------------------------------------------------------
    % AP_print_progress_fraction(curr_im,length(im_files));
end


if verbose; toc; end


%% Get illumination color for each frame
% Assume: blue/violet alternating, with blue starting whenever there is a >
% 2s gap (= new recording, which starts on blue)

% Sanity check: make sure frame numbers are consecutive
if any(diff(vertcat(frame_info.frame_num)) ~= 1)
    error('Frame numbers not consecutive');
end

% Get difference between frames in miliseconds - without the first one (because it's 0)
timestamp_diff = frame_info.time_from_last(2:end);

% Find recording boundaries (time between frames > threshold)
recording_boundary_thresh = 20000; % miliseconds between frames to define recording
recording_frame_boundary = [0,find(timestamp_diff > recording_boundary_thresh),frame_info(end).frame_num(end)]+1;

% Get recording for each frame
im_rec_idx = mat2cell( ....
    interp1(recording_frame_boundary,1:length(recording_frame_boundary), ...
    [1:frame_info(end).frame_num(end)],'previous'), ...
    1,cellfun(@length,{frame_info.frame_num}));
[frame_info.rec_idx] = im_rec_idx{:};

% Get illumination color of frame (alternate starting at each recording)
n_colors = 2;
im_frame_color = mat2cell( ...
    1+mod([1:frame_info(end).frame_num(end)] - ...
    interp1(recording_frame_boundary,recording_frame_boundary,...
    1:frame_info(end).frame_num(end),'previous'),n_colors), ...
    1,cellfun(@length,{frame_info.frame_num}));
[frame_info.color] = im_frame_color{:};


%% Set pixels to keep 
% % (no timestamp to remove)
% 
im_info = imfinfo(fullfile(im_path,im_files(1).name));
im_size = [im_info(1).Height,im_info(1).Width];

im_px_loc = {[1,im_size(1)],[1,im_size(2)],[1,Inf]};
im_grab_size = cellfun(@(x) diff(x)+1,im_px_loc(1:2));


%% Create moving- and total-average images

% Set moving average number of frames
n_frame_avg = 15;

tic;
if verbose; disp('Building image moving average by color...'); end

% Loop through images, cumulatively build averages by illumination color
im_avg_color = zeros(im_grab_size(1),im_grab_size(2),n_colors);
im_color_mov_avg = cell(length(im_files),n_colors);
for curr_im = 1:length(im_files)
   
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    im = single(tiffreadVolume(curr_im_fn,'PixelRegion',im_px_loc));

    % Loop through illumination colors
    for curr_color = 1:n_colors
        % Get color index forifif frames in current image
        curr_frame_color_idx = frame_info(curr_im).color == curr_color;

        % Cumulatively add average image
        curr_color_partial = ...
            sum(im(:,:,curr_frame_color_idx)./sum([frame_info.color] == curr_color),3);
        im_avg_color(:,:,curr_color) = im_avg_color(:,:,curr_color) + ...
            curr_color_partial;

        % Get moving average (truncate based on moving avg modulus)
        curr_n_frames = sum(curr_frame_color_idx);
        curr_frame_avg_idx = find(curr_frame_color_idx, ...
            curr_n_frames - mod(curr_n_frames,n_frame_avg));
        im_color_mov_avg{curr_im,curr_color} = ...
            permute(mean(reshape(im(:,:,curr_frame_avg_idx), ...
            size(im,1),size(im,2),n_frame_avg,[]),3),[1,2,4,3]);
    end

    % temp took out --------------------------------------------------------
    % -----------------------------------------------------------------------------------------------------------------------------
%     AP_print_progress_fraction(curr_im,length(im_files));

end

if verbose; toc; end


%% Do SVD on moving-average images
% (keep U and S, don't keep V since U will be applied to full dataset next)

tic;
if verbose; disp('Running SVD on moving average images by color...'); end

[U,~,~] = arrayfun(@(color) ...
    svd(reshape(cat(3,im_color_mov_avg{:,color}),prod(im_grab_size),[]),'econ'), ...
    1:n_colors,'uni',false);

% Reshape U into pixels row x pixels column x components
U = cellfun(@(x) reshape(x,im_grab_size(1),im_grab_size(2),[]),U,'uni',false);

toc;


%% Apply spatial components (U's) from SVD to full data
% (note: spatial components are applied as U' * mean-subtracted data, so
% the resulting matrix is equivalent to S*V but just called 'V')

tic;
if verbose; disp('Applying SVD spatial components to full data...'); end

V = cell(length(im_files),n_colors);
for curr_im = 1:length(im_files)
   
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    im = single(tiffreadVolume(curr_im_fn,'PixelRegion',im_px_loc));

    % Loop through illumination colors
    for curr_color = 1:n_colors
        % Get color index for frames in current image
        curr_frame_color_idx = frame_info(curr_im).color == curr_color;

        % Apply spatial components to mean-subtracted data
        V{curr_im,curr_color} = ...
            reshape(U{curr_color},[],size(U{curr_color},3))' * ...
            (reshape(im(:,:,curr_frame_color_idx),[],sum(curr_frame_color_idx)) - ...
            reshape(im_avg_color(:,:,curr_color),[],1));
    end
% temp took out ----------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------
%     AP_print_progress_fraction(curr_im,length(im_files));

end

if verbose; toc; end


%% Split V's by recording (instead of by file)

if verbose; disp('Applying SVD spatial components to full data...');end

% Store V's as recordings x color
Vrec = cell(length(recording_frame_boundary)-1,n_colors);

frame_num_cat = vertcat(frame_info.frame_num);

frame_color_cat = horzcat(frame_info.color)';
frame_rec_idx_cat = horzcat(frame_info.rec_idx)';
for curr_color = 1:n_colors

    % Split concatenated V by recording index
    curr_V_cat = horzcat(V{:,curr_color});
    Vrec(:,curr_color) = ...
        mat2cell(curr_V_cat,size(curr_V_cat,1), ...
        accumarray(frame_rec_idx_cat(frame_color_cat == curr_color),1));
    
end

if verbose; disp('Finished SVD.'); end















































