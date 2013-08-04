% reads geometry, basis and orbitals from Molden-formatted file

function [atoms, basis, orbitals] = read_molden(file)
    f = fopen(file);
    atoms = struct('Z', {}, 'R', {});
    basis = struct('R', {}, 'l', {}, 'zeta', {}, 'd', {});
    orbitals = struct('sym', {}, 'ene', {}, 'spin', {},...
                      'occ', {}, 'coeff', {});
    rheader = '\[([^\[\]]*)\]\s*(\w.*)?';
    l = fgets(f);
    while ~feof(f)
        while isempty(strtrim(l))
            l = fgets(f);
        end
        tok = regexp(l, rheader, 'tokens');
        [section, info] = deal(tok{1}{:});
        switch section
            case 'Atoms'
                bohr = 1.889725989;
                while true
                    l = fgets(f);
                    if ~isempty(regexp(l, rheader, 'once'))
                        break
                    end
                    l = regexp(strtrim(l), '\s+', 'split');
                    R = cellfun(@str2double, l(4:6));
                    if strcmpi(info, 'angs')
                        R = bohr*R;
                    end 
                    atoms(end+1) = struct('Z', sscanf(l{3}, '%i'), 'R', R);
                end
            case 'GTO'
                moments = get_moments();
                for i = 1:length(atoms)
                    l = fgets(f);
                    iatom = sscanf(l, '%i', 1);
                    while true
                        l = fgets(f);
                        if isempty(strtrim(l))
                            break
                        end
                        tok = regexp(l, '(\w)\s+(\d+)\s+1.00', 'tokens');
                        [moment, n] = deal(moments.(tok{1}{1}),...
                                           sscanf(tok{1}{2}, '%i'));
                        [zeta, d] = deal(zeros(n, 1));
                        for j = 1:n
                            l = fgets(f);
                            l = sscanf(strrep(l, 'D', 'e'), '%f', 2);
                            if l(2) == 0
                                continue
                            end
                            zeta(j) = l(1);
                            d(j) = l(2);
                        end
                        zeta(d==0) = [];
                        d(d==0) = [];
                        basis(end+1) = struct('R', atoms(iatom).R,...
                                              'l', moment, 'zeta', zeta,...
                                              'd', d);
                    end
                end
            case 'MO'
                l = [basis.l];
                nbas = sum((l+1).*(l+2)/2);
                while ~feof(f)
                    for j = 1:4
                        l = fgets(f);
                        l = regexp(strtrim(l), '=\s+', 'split');
                        data.(l{1}) = l{2};
                    end
                    coeff = fscanf(f, '%f', [2 nbas])';
                    fgets(f);
                    if size(coeff, 1) < nbas
                        warning('Orbital %i set to zero', length(orbitals)+1);
                        for j = 1:nbas-size(coeff, 1)
                            fgets(f);
                        end
                        coeff = zeros(nbas, 2);
                    end
                    orbitals(end+1) = struct('sym', data.Sym,...
                                             'ene', str2double(data.Ene),...
                                             'spin', strcmp(data.Spin, 'Alpha'),...
                                             'occ', str2double(data.Occup),...
                                             'coeff', coeff(:, 2));
                end
            otherwise
                while true
                    l = fgets(f);
                    if ~isempty(regexp(l, rheader, 'once'))
                        break
                    end
                end
        end
    end
    fclose(f);
end

