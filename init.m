function init()
    dirs = {'init' 'ints' 'misc' ['misc' filesep 'tprod'] 'data' 'dft'};
    for i = 1:length(dirs)
        dirs{i} = [pwd filesep dirs{i}];
    end
    addpath(dirs{:});
end