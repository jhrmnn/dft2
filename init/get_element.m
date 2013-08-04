% converts element symbol on atomic number or vice versa. 12/04/07

function n = get_element(s)
	table = regexp(fileread('elements.txt'), '\n', 'split');
	if isnumeric(s)
		n = table{s};
    else
        n = find(strcmpi(s, table), 1);
	end
end

