function[] = clear_saved_data(force)
% clear_saved_data
%
% Removes all saved data for "fast" routines. This deletes all .mat files in
% the "data" subdirectory.

if (nargin < 1) || (force ~= 'f')

  response = input('Delete all .mat files in the "data" subdirectory? (y/n) ', 's');
  if ~strcmp(response, 'y')
    fprintf('Returning without deleting data.\n')
    return
  end

end

mfiledir = fileparts(mfilename('fullpath'));

cmd = fullfile(mfiledir, 'data', '*.mat');
disp(cmd)

delete(cmd);
fprintf('Cleared .mat files in "data" subdirectory\n');
