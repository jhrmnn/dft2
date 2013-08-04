function moments = get_moments()
    alphabet = 'spdfghi';
    for i = 1:length(alphabet)
        moments.(alphabet(i)) = i-1;
    end
end