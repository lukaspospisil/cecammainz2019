function [ out ] = signalclustering( in )

% CPU sequential code
if or(strcmp(in.options.type,'FEMH1'),strcmp(in.options.type,'SPAH1'))
    out = signalclustering_seqcpu(in);
end

end
