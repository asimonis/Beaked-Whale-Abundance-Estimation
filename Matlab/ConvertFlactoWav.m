%% Configuration
input_dir  = "D:/sg679_CalCurCEAS_Aug2024/audio-FLAC/20240920";  % single subfolder to test
output_dir = 'D:/sg679_CalCurCEAS_Aug2024/audio-WAV/20240920';
dry_run    = false;                             % set to false to actually convert

%% Find all .flac files
files = dir(fullfile(input_dir, '*.flac'));
fprintf('Found %d .flac files in: %s\n', length(files), input_dir);

if isempty(files)
    error('No .flac files found — check your input_dir path.');
end

%% Create output directory if needed
if ~dry_run && ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

if dry_run
    fprintf('\n*** DRY RUN — no files converted. Set dry_run = false to execute. ***\n');
    fprintf('Would convert %d files to: %s\n', length(files), output_dir);
    return
end

%% Convert
converted  = 0;
errors     = 0;
error_log  = {};

fprintf('\nConverting %d files...\n', length(files));

for i = 1:length(files)
    src_path  = fullfile(input_dir, files(i).name);
    wav_name  = strrep(files(i).name, '.flac', '.wav');
    dest_path = fullfile(output_dir, wav_name);

    try
        [data, fs] = audioread(src_path);   % read flac (preserves sample rate)
        audiowrite(dest_path, data, fs);    % write wav
        converted = converted + 1;

        if mod(i, 10) == 0
            fprintf('  %d / %d done...\n', i, length(files));
        end

    catch e
        fprintf('ERROR converting: %s\n  %s\n', files(i).name, e.message);
        error_log{end+1} = files(i).name; %#ok<AGROW>
        errors = errors + 1;
    end
end

fprintf('\nDone. Converted: %d  |  Errors: %d\n', converted, errors);

%% Write error log if needed
if ~isempty(error_log)
    log_path = fullfile(output_dir, 'conversion_errors.txt');
    fid = fopen(log_path, 'w');
    for i = 1:length(error_log)
        fprintf(fid, '%s\n', error_log{i});
    end
    fclose(fid);
    fprintf('Error log written to: %s\n', log_path);
end