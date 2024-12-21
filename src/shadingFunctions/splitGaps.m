function [x_split,y_split,shaded_split] = splitGaps(x,y,shaded)
% Stupid little function to split up (x,y) position vectors to be plotted as segments
% for shading +plant motion plots

% Find where logical vector changes between zero and one
d_shaded = [diff(shaded); 0];

ind_start = find(shaded==0 & d_shaded~=0) + 1;  % where 1's start
ind_end= find(shaded==1 & d_shaded~=0);        % where 1's end

% if the array is entirely 1's...
if isempty(ind_start) || isempty(ind_end)
    x_split = x;
    y_split = y;
    shaded_split = shaded;
    return
end

% if the array starts with 1's...
if ind_start(1)>ind_end(1)
    ind_start = [1; ind_start];
end

% Number of subset vectors to split off
N = max(length(ind_start),length(ind_end));

% Initialize cell arrays
shaded_split = cell(1,N);
x_split = cell(1,N);
y_split = cell(1,N);

% Create each subset vector
for i=1:N

    % if the array ends with 1's...
    if i <= length(ind_end)
        shaded_split{i} = shaded(ind_start(i):ind_end(i));
        x_split{i} = x(ind_start(i):ind_end(i));
        y_split{i} = y(ind_start(i):ind_end(i));
        % keyboard
    elseif i > length(ind_end)
        shaded_split{i} = shaded(ind_start(i):end);
        x_split{i} = x(ind_start(i):end);
        y_split{i} = y(ind_start(i):end);
        % keyboard
    end

    
end

% end function
end