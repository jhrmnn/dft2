% reads xyz geometry

function atoms = read_xyz(file)
	bohr = 1.889725989;
	f = fopen(file, 'r');
    l = fgets(f);
	n = sscanf(l, '%i');
	fgets(f);
    atoms = struct('Z', cell(1, n), 'R', []);
	for i = 1:n
        l = fgets(f);
        l = regexp(strtrim(l), '\s+', 'split');
        R = bohr*cellfun(@str2double, l(2:4));
        atoms(i) = struct('Z', get_element(l{1}), 'R', R);
	end
	fclose(f);
end

