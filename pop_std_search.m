function [STUDY ALLEEG ] = pop_std_search( ALLEEG,EEG,CURRENTSET,varargin)

% [STUDY ALLEEG ] = pop_std_search( ALLEEG,EEG,CURRENTSET,varargin)
% 
% search for datasets to include in a STUDY.
% 
% optional key-value pairs of inputs: 
%           'rootdir', string path to root directory to search files in
%                       (recursively) 
%           'whatfiles', string search file pattern (regular expression)
%           'studyname', string study name.
%   

g = vararg2struct(varargin);
struct2ws(g);

if not(isfield(g,'rootdir'))
    rootdir = uigetdir(cd,'Select root directory of study data.');
    if isequal(rootdir, 0)
        STUDY  = [];
        ALLEEG = [];
        return
    end
end
if ~isfield(g,'whatfiles') || ~isfield(g,'studyname')
    try def{1} = g.whatfiles;catch; def{1} = '.*\.set';end
    try def{2} = g.studyname;catch; def{2} = 'New STUDY';end
    [rep] = inputdlg({'Search pattern' 'STUDY name'},'Search for files to create STUDY...',1,def);
    if isempty(rep)
        STUDY  = [];
        ALLEEG = [];
        return
    end
    whatfiles = rep{1};
    studyname = rep{2};
end

fs = flister(whatfiles,'dir',rootdir);
if isempty(fs)
    errordlg('No files found...')
    STUDY  = [];
    ALLEEG = [];
    return
end

cmds = {};
for i = 1:numel(fs)
    ncmd = {'index',i, 'load',fs(i).name,'condition','AllConds'};
    if isfield(fs,'suj')
        ncmd = [ncmd {'subject',fs(i).suj}];
    elseif  isfield(fs,'subject')
        ncmd = [ncmd {'subject',fs(i).subject}];
    end
    cmds = [cmds ncmd];
end
    
STUDY=[];
if isempty(CURRENTSET)
    CURRENTSET = 1;
end
%creating the study
[STUDY ALLEEG]= std_editset(STUDY, ALLEEG,'name', studyname, 'commands', cmds);
[EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',[1:numel(ALLEEG)] ,'study',1);
CURRENTSTUDY = 1;




