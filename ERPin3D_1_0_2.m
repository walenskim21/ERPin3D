% Reset workspace
%clear, close all

%change constants here
voxelsize_x = 75;
voxelsize_y = 75;
time_range=[-200 1200];
electrode_xyplane_extent=[-500 500];
voxeltransparency = .3; % 0 = transparent, 1 = opaque

% Create a custom dialog box using questdlg
dlg_title = 'File Selection';
prompt = 'Select an option:';
options = {'Choose Electrode Coordinate File', 'Choose Data File', 'Cancel'};
default_option = 'Cancel';

% Initialize variables
location_data = [];
values = [];

while true
    % Display the custom dialog box and get the user's choice
    choice = questdlg(prompt, dlg_title, options{:}, default_option);

    switch choice
        case 'Choose Electrode Coordinate File'
            % Importing Electrode Coordinate Datafiles
            disp('Select Electrode Coordinate Data File:');
            [filename, pathname] = uigetfile('*.txt', 'Choose Electrode Coordinate File');
            delimiterIn = '\t';

            % Check if the user canceled the selection
            if isequal(filename,0) || isequal(pathname,0)
                disp('File selection canceled. Exiting...')
                return;
            end

            % Read electrode datafile
            location_data = readcell(fullfile(pathname, filename), 'Delimiter', delimiterIn);

            % Get variable names from column 1, 1 row down
            location_names = location_data(2:end, 1); 

            % Get data values from rest of data file
            location_values = location_data(2:end, 4:5);

            %convert data values to a matrix
            electrodes = cell2table(location_values,"RowNames", location_names, "VariableNames", {'x','y'});

        case 'Choose Data File'
            % Importing EEG Data
            disp('Select EEG Data File:');
            [filename, pathname] = uigetfile('*.txt', 'Choose Data File');
            delimiterIn = '\t';

            % Check if the user canceled the selection
            if isequal(filename,0) || isequal(pathname,0)
                disp('File selection canceled. Exiting...')
                return;
            end

            % Read in data into cell array
            values = readcell(fullfile(pathname, filename), 'Delimiter', delimiterIn);

            % Get variable names from column 1
            variable_names = values(:,1);

            % Get data values from rest of data file and transpose data
            values = values(:,2:end).';

            % Convert cell array to table with data organized in columns with
            % Appropriate variable names. To access data, use dot notation with variable
            % name. (ex. EEGdata.time would return the time values.
            values = cell2table(values, 'VariableNames', variable_names);

        case 'Cancel'
            disp('File selection canceled. Exiting...')
            return;
    end
    
    % Check if both files have been selected
    if ~isempty(location_data) && ~isempty(values)
        break; % Exit the loop when both files have been selected
    end
end

% Variables
voxelsize_z = abs(values{1,3} - values{1,2});

% Define colormap choices
colormap_choices = {'jet', 'hsv', 'hot', 'cool', 'gray', 'turbo', 'parula'};

% Create a custom dialog box using questdlg
dlg_title = 'Settings Selection';
prompt = {
    'Select a colormap (eg. jet, hsv, hot, cool, gray, turbo, parula):'
    'Enter the lower bound of data values to exclude:'
    'Enter the upper bound of data values to exclude:'
};
num_lines = [1, 50];
default_answers = {
    colormap_choices{1}
    ''
    ''
};

% Display the custom dialog box
answers = inputdlg(prompt, dlg_title, num_lines, default_answers);

% Check if the user canceled the selection
if isempty(answers)
    disp('Settings Selection. Exiting...')
    return;
end

% Get the input values
colormap_choice = answers{1};

%initialize variables
lower_range = 0;
upper_range = 0;
lower_range = str2double(answers{2});
upper_range = str2double(answers{3});

% Calculate the overall highest value and overall lowest value
overall_min = min(values{:, 2:end}, [], 'all');
overall_max = max(values{:, 2:end}, [], 'all');

% Creates the empty plot.
clc
close all
fig=figure('Visible','off');
ax = axes(fig);
hold(ax, 'on');
xlim(ax, electrode_xyplane_extent);
ylim(ax, electrode_xyplane_extent);
zlim(ax, time_range);
view(ax, 3)
hold(ax, 'on')
grid on
colormap(ax, colormap_choice); % Set colormap
% Add axis labels
labelPosition = [-600, 300, -200]; % x, y, z coordinates of the label
labelString = 'Anterior';
text(labelPosition(1), labelPosition(2), labelPosition(3), labelString, 'HorizontalAlignment', 'center');
labelPosition = [-600, -175, -200]; % x, y, z coordinates of the label
labelString = 'Posterior';
text(labelPosition(1), labelPosition(2), labelPosition(3), labelString, 'HorizontalAlignment', 'center');
labelPosition = [-260 -600, -200]; % x, y, z coordinates of the label
labelString = 'Left';
text(labelPosition(1), labelPosition(2), labelPosition(3), labelString, 'HorizontalAlignment', 'center');
labelPosition = [275, -600, -200]; % x, y, z coordinates of the label
labelString = 'Right';
text(labelPosition(1), labelPosition(2), labelPosition(3), labelString, 'HorizontalAlignment', 'center');
zlabel('Time (ms)')

% Adds the figure into the plot.
I = imread("blank_head_ERPin3D.jpg","jpg");
z_off = 0;
X = [electrode_xyplane_extent(1) size(I,2)];
Y = [size(I,1) electrode_xyplane_extent(1)];
Z = z_off + zeros(2);
surf(ax, X,Y,Z,I,'FaceColor','texture','EdgeColor','none')
plot3(ax, rand(1,10), rand(1,10), rand(1,10), 'Color', [.5 .5 .5], 'LineWidth', 4)
view(3)


% Create a color scale using the selected colormap
num_colors = 64; % Number of colors in the colormap
color_scale = colormap(colormap_choice); % Selected colormap

% Voxel properties
% Voxel size in x, y, z directions
voxel_size = [voxelsize_x, voxelsize_y, voxelsize_z];
% Transparency of the voxels
voxel_transparency = voxeltransparency;

% Define voxel faces
faces = [
    1, 2, 3, 4;
    2, 6, 7, 3;
    6, 5, 8, 7;
    5, 1, 4, 8;
    1, 2, 6, 5;
    4, 3, 7, 8
    ];

% Voxel vertices
x = [0; voxel_size(1); voxel_size(1); 0; 0; voxel_size(1); voxel_size(1); 0];
y = [-voxel_size(2)/2; -voxel_size(2)/2; voxel_size(2)/2; voxel_size(2)/2; -voxel_size(2)/2; -voxel_size(2)/2; voxel_size(2)/2; voxel_size(2)/2];
z = [0; 0; 0; 0; 0 + voxel_size(3); 0 + voxel_size(3); 0 + voxel_size(3); 0 + voxel_size(3)];


%start voxel loops
count_draw_calls=1;
count_faces=0;
count_vertices=0;
all_vertices={};
all_faces={};
all_colors={};
       
%start the new loop
for i = 1:length(location_names)

    if ~ismember(location_names{i}, values.Properties.VariableNames)
        continue;
    end

    %variable for current node
    current = location_names{i};

    %get the current electrode
    electrode_name = current;

    %get x and y  and z cords
    x_cord = electrodes{current,"x"};
    y_cord = electrodes{current,"y"};
    z_cord = time_range;

    %add text labels
    text(ax,x_cord, y_cord, time_range(2), electrode_name, 'BackgroundColor', [0.7 0.9 0.7], 'FontWeight', 'bold', 'FontSize', 8, 'HorizontalAlignment', 'center');

    % Plot vertical line
    plot3(ax,[x_cord x_cord], [y_cord y_cord], z_cord, 'k', 'LineWidth', 0.75);

    % Get data values for the current voxel
    data_value = values{:, current};

    % Define the range of the data values
    data_min = min(data_value);
    data_max = max(data_value);

    % Create an interpolated color scale based on the range of data values
    color_map_range = linspace(overall_min, overall_max, num_colors);
    interpolated_colors = interp1(linspace(overall_min, overall_max, size(color_scale, 1)), color_scale, color_map_range, 'linear');
    z_values = table2array(values(1:end,1));
    % Loop through the desired range and plot the voxels
        for i = 1:length(z_values) 
        z_val=z_values(i);
        % Update the z-coordinate of voxel vertices
        vertices = [x + x_cord, y + y_cord, z + z_val];

        % Get the corresponding values for the current range
        if z_val <= time_range(2)
            range_values = data_value(i, :);
        else
            break; % Stop the loop if the index exceeds the number of elements in 'values'
        end

        % Skip voxel creation if values are between lower range and upper range
        if all(range_values >= lower_range) && all(range_values <= upper_range) 
            continue;
        end

        % Compute the color index based on the range values
        color_index = round((range_values - overall_min) / (overall_max - overall_min) * (num_colors - 1)) + 1;
        color_index = max(color_index, 1); % Ensure indices are within the valid range
        color_index = min(color_index, num_colors);
      
        % Retrieve the corresponding color for each voxel
        color = interpolated_colors(color_index, :);

        % Repeat the color for each voxel
        num_voxels = size(vertices, 1) / size(faces, 1);
        color_repeated = repmat(color, num_voxels * size(faces, 1), 1);

        % Plot voxel using patch function with specified colors and transparency
        %patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', color_repeated, 'FaceColor', 'interp', 'FaceAlpha', voxel_transparency, 'EdgeColor', 'none');
         all_vertices{count_draw_calls}=vertices;
         all_faces{count_draw_calls}=faces+(count_draw_calls-1)*max(faces(:));
         all_colors{count_draw_calls}=color_repeated;
         count_draw_calls=count_draw_calls+1;
         count_faces=count_faces+size(faces,1);
         count_vertices=count_vertices+size(count_vertices,1);
        end
end

vertices=zeros(count_vertices,3);
colors=zeros(count_vertices,3);
faces=zeros(count_faces,4);
count1=1; %counting final number of vertices and colors
count2=1; %counting final number of faces
for i=1:length(all_vertices)
  vertices(count1:count1+size(all_vertices{i},1)-1,:)=all_vertices{i};
  colors(count1:count1+size(all_colors{i},1)-1,:)=all_colors{i};
  count1=count1+size(all_vertices{i},1);
  faces(count2:count2+size(all_faces{i},1)-1,:)=all_faces{i};
  count2=count2+size(all_faces{i},1);
end

patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'FaceAlpha', voxel_transparency, 'EdgeColor', 'none');

% Create the colormap scale
tf = isMATLABReleaseOlderThan("R2022a");
if (tf ==1)
    caxis(ax,[overall_min, overall_max]);
else
    clim(ax,[overall_min, overall_max]);
end
colorbar(ax);
fig.Visible='on';
drawnow;

