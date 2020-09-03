function output_file = convertPlexonToRawBinary(file, obj)
% function output_file = convertPlexonToRawBinary(file, obj)
% convert plexon file to raw binary file for kilosort2

if ~exist('mexPlex', 'file')
    error('Cannot find plexon library, download and add to matlab path')
end

if nargin<1
    file = find_files('*.pl2');
    file=file{1};
end

if nargin<2
    obj = [];
end

%% check output file
output_file = [file(1:end-3), 'bin'];
if exist(output_file,'file')
    log(obj, 'Found existing binary file, no conversion needed');
    log(obj, sprintf('Slected file %s', output_file));
    return;
end

%% find channel index
chan_name_str_groups = {'WB', 'SPKC'};
chan_index_groups = cell(size(chan_name_str_groups));
n_chan_groups = length(chan_name_str_groups);
[n_adchans, adchan_names] = plx_adchan_names(file);

for i=1:n_chan_groups
    chan_name_str = chan_name_str_groups{i};
    ind = false(n_adchans,1);
    for j=1:n_adchans
        ind(j) = is_str_start_with(adchan_names(j,:), chan_name_str);
    end
    chan_index_groups{i} = find(ind);
end

% get the number of total sample
[~, adchan_samplecounts] = plx_adchan_samplecounts(file);

% find channel group with samples > 0
chan_found = false;
for i=1:n_chan_groups
    nsamples = adchan_samplecounts(chan_index_groups{i}(1));
    if nsamples > 0
        chan_found = true;
        chan_index = chan_index_groups{i}';
        break;
    end
end

if ~chan_found
    error('Cannot find continuous channel data, please check plexon file')
end

%% get sampling frequency
if isobject(obj)
    [~, freqs] = plx_adchan_freqs(file);
    sampling_frequency = freqs(chan_index(1));
    obj.H.settings.setFsEdt.String = num2str(sampling_frequency);
end


%% read by blocks
nch = length(chan_index);
BLOCK_SIZE = 1024^3/(nch*2); % 1G block
n_block = ceil(nsamples/BLOCK_SIZE);

log(obj, 'Converting 0%...');

fid = fopen(output_file, 'w');
data_block = zeros(nch, BLOCK_SIZE, 'int16');
for i = 1:n_block
    startCount = (i-1)*BLOCK_SIZE+1;
    endCount = i*BLOCK_SIZE;
    
    if i == n_block % last block
        endCount = nsamples;
        nsample = endCount-startCount+1;
        data_block = zeros(nch, nsample, 'int16');
    end
    
    for ich = 1:nch
        [~, ~, ad] = plx_ad_span(file, chan_index(ich)-1, startCount, endCount);
        data_block(ich,:) = int16(ad);
    end
    
    % write file
    fwrite(fid, data_block, 'int16');
    log(obj, sprintf('Converting %.0f%%...', i/n_block*100))
end
fclose(fid);

log(obj, 'Conversion completed');
log(obj, sprintf('Output file %s', output_file));

end

function log(obj, message)

if isobject(obj) && ismethod(obj, 'log')
    obj.log(message)
else
    fprintf('%s\n', message);
end

end

function TF = is_str_start_with(str, start_str)
% TF = is_str_start_with(str, start_str)

n = length(start_str);
TF = length(str) >= n && strcmp(str(1:n), start_str);

end