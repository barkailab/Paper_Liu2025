% this function is for reading meme files for motifs, to get PWM in mat
% format
function motifs = load_meme_pwm(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file.');
    end

    motifs = struct('name', {}, 'matrix', {});
    motifIdx = 0;

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        % Look for MOTIF definition
        if startsWith(line, 'MOTIF ')
            motifIdx = motifIdx + 1;
            tokens = split(line);
            motifs(motifIdx).name = tokens{2};
            motifs(motifIdx).matrix = [];

        % Look for letter-probability matrix
        elseif startsWith(line, 'letter-probability matrix')
            % Start reading the matrix rows
            while true
                pos = ftell(fid); % Remember position in case we need to step back
                nextLine = strtrim(fgetl(fid));
                if isempty(nextLine) || startsWith(nextLine, 'MOTIF') || startsWith(nextLine, 'URL') || startsWith(nextLine, 'letter-probability matrix')
                    fseek(fid, pos, 'bof'); % Go back one line
                    break;
                end

                row = str2num(nextLine); %#ok<ST2NM> % Convert to numeric
                if ~isempty(row)
                    motifs(motifIdx).matrix(end+1, :) = row;
                end
            end
        end
    end

    fclose(fid);
end
