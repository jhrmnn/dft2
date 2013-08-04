function E = nuclear(atoms)
    E = 0;
    for i = 1:length(atoms)
        for j = 1:i-1
            E = E+atoms(i).Z*atoms(j).Z/norm(atoms(i).R-atoms(j).R);
        end
    end
end