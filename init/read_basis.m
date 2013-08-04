% reads basis set definition in Turbomole format

function basisdef = read_basis(file)
	f = fopen(file, 'r');
	basisdef = struct('Z', {}, 'l', {}, 'zeta', {}, 'd', {});
    while true
        l = fgets(f);
        if ~isempty(regexp(l, '^\$basis', 'once'))
            fgets(f);
            break
        end
    end
    moments = get_moments();
    while true
        l = fgets(f);
        if ~isempty(regexp(l, '^\$end', 'once'))
            return
        end
        Z = regexp(l, '^(\w+)', 'tokens');
        Z = get_element(Z{1}{1});
        fgets(f);
        while true
            l = fgets(f);
            if regexp(l, '^\*')
                break
            end
            l = regexp(strtrim(l), '\s+', 'split');
            n = str2double(l{1});
            moment = moments.(l{2});
            [zeta, d] = deal(zeros(n, 1));
            for j = 1:n
                l = fgets(f);
                l = sscanf(strrep(l, 'D', 'e'), '%f', 2);
                zeta(j) = l(1);
                d(j) = l(2);
            end
            basisdef(end+1) = struct('Z', Z, 'l', moment,...
                                     'zeta', zeta, 'd', d);
        end
    end
	fclose(f);
end