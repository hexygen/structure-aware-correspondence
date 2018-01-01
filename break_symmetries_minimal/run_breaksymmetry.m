% clear all;
close all;


resultsdir = '../results/shrec/';

pairs = dir([resultsdir '*.mat']);

n = length(pairs);

TB = [];
for i=1:n
    name = pairs(i).name(1:end-4);
    sf = strfind(name, 'nosym');
    
    if (isempty(sf))

        load_name = [resultsdir name];
        save_name = [resultsdir name '_nosym'];

        dowork = 0;
        if (exist([save_name '.mat'], 'file'))
            display([save_name ' - file exists.']);

        elseif (~exist([load_name '.mat'], 'file'))
            display([load_name ' - file does not exist!']);
        else
            dowork = 1;
        end

        if (dowork)

            % Load symmetric matching:
            x = load(load_name);

            R1 = x.R;

            if (isfield(R1.M1.shape, 'PCD'))
                R1.M1 = CloudToTris(R1.M1);
            end

            if (isfield(R1.M2.shape, 'PCD'))
                R1.M2 = CloudToTris(R1.M2);
            end

            t1 = tic;

            R1 = BreakSymmetries(R1);

            t = toc(t1);
            TB(end + 1) = t;
            
            R1.bs_time = t;
            
            % Visualize results:
            VisualizeMatching(R1, save_name);   

            close all;

            display(['Saved ' save_name ', time = ' num2str(t) ' seconds.']);
        end

    end
end

