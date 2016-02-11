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




function [varargout] = flister(re,varargin)

% [filenames, factnames,factidx,factlevelnames] = flister(re)
% [fstruct] = flister(re,'key',val, ...)
% 
% This function lists all files matching re in the current directory and
% subdirectories
% re is of the type 'AP(?<suj>.*)_YN1x2_(?<Noise>Yes|No)Noise_(?<SDT>Hit|Miss)\.set$'
% (see regular expressions, named tokens)
%
% optional input key val pairs : 
%           'exclude', str  : file names matching regexp str will be
%                            excluded from the list 
%           'dir',str       : root directory where to search for the files 
%           'isdir',boolean : whether or not listed files should be
%                           directories (defaults to NaN to ignore this filter)
%           'list',cellstr  : any list to process with flister (instead of
%                           listing files)
%           'rename',str    : renaming pattern, using tokens defined in re,
%                           passed to regexprep
% outputs:
%           filenames: cell of filenames with absolute path
%           factnames: cell: all the factor names of the named tokens of re
%           factlevelnames cell: the level names of each factor
%           factidx:   cell: for each factor, the level to which each file belongs
% output method 2:
%           f: 	structure with fields name (name of the file), and each of
%               the factor names with the corresponding level to which the file belongs.

% v 1 Max: Basic functionality
% v 1.1 Max: added list input method 01/06/2015

if numel(varargin) == 1 && isstruct(varargin{1})
    varargin = struct2vararg(varargin{1});
end
g = finputcheck( varargin, ...
    {
    'exclude'   { 'string';'cell' }    []   ''
    'dir' { 'string'}    []   cd
    'isdir' {'integer'} [NaN 0 1] NaN
    'list' {'cell'} {''} {}
    'rename' {'string'} [] ''
    });
if ischar(g)
    error(g);
end
rootdir = g.dir;
% we'll first list
% everything and then filter out.
if isempty(g.list)
    [dum,dum,all] = dirr(rootdir,'name');
    if isempty(all)
        disp('No files found')
        varargout{1} = [];
        return
    end
else
    all = g.list;
end
filenames = all(regexpcell(all,re));
if isempty(filenames)
    disp('All files filtered out');
    varargout{1} = [];
    return
end
if not(isempty(g.exclude))
    filenames(regexpcell(filenames,g.exclude)) = [];
end
if not(isnan(g.isdir))
    aredir = cellfun(@isdir,filenames);
    if g.isdir
        filenames = filenames(aredir);
    else
        filenames = filenames(~aredir);
    end
end
    

names = regexp(filenames,re,'names');
if isempty(names) || numel(fieldnames(names{1})) == 0
    % disp('Warning: Could not find tokens in regular expression')
    factnames = {};
    f = struct('name',filenames);
else
    for i = 1:numel(names)
        f(i) = names{i};
    end
    factnames = fieldnames(f);
    [f.name] = filenames{:};
end

for i = 1:numel(factnames)
    factlevelnames{i} = unique({f.(factnames{i})});
    finder = regexptranslate('escape',{f.(factnames{i})});
    factidx{i} = regexpcell(factlevelnames{i},finder,'exact');
end
if not(isempty(factnames))
    for i = 1:numel(filenames)
        for j = 1:numel(factnames)
            f(i).([factnames{j} 'idx']) = factidx{j}(i);
        end
    end
end

if not(isempty(g.rename))
    nf = regexprep({f.name},re,g.rename);
    cellfun(@(x,y)fprintf('%s --> %s\n',x,y),{f.name},nf)
    r = input('Is this ok? (y/n)','s');
    if strcmpi(r,'y')
        for i_f = 1:numel(f)
            movefile(f(i_f).name,nf{i_f});
        end
    end
end

if nargout > 1
    varargout{1} = filenames;
    varargout{2} = factnames;
    varargout{3} = factidx;
    varargout{4} = factlevelnames;
elseif nargout == 1
    varargout{1} = f;
else
    for i = 1:numel(f)
        if isdir(f(i).name)
            continue
        end
        [p fn e] = fileparts(f(i).name);
        fl = strrep(p,'''','''''');
        switch e
            case '.m'
                out = sprintf(['<a href="matlab: %s">run</a> %s%s<a href="matlab:edit(''%s'')">%s%s</a>'], fn, p, filesep, f(i).name, fn, e);
            case '.mat'
                out = sprintf(['    %s%s<a href="matlab:load(''%s'')">%s%s</a>'], p, filesep, f(i).name, fn, e);
            otherwise
                out = sprintf(['    %s%s%s%s'], p, filesep, fn, e);
        end
        out = strrep(out,'\','\\');
        out = [out '\n'];
        fprintf(char(out));
    end% for
end
    
 
    
% finputcheck() - check Matlab function {'key','value'} input argument pairs
%
% Usage: >> result = finputcheck( varargin, fieldlist );
%        >> [result varargin] = finputcheck( varargin, fieldlist, ... 
%                                              callingfunc, mode, verbose );
% Input:
%   varargin  - Cell array 'varargin' argument from a function call using 'key', 
%               'value' argument pairs. See Matlab function 'varargin'.
%               May also be a structure such as struct(varargin{:})
%   fieldlist - A 4-column cell array, one row per 'key'. The first
%               column contains the key string, the second its type(s), 
%               the third the accepted value range, and the fourth the 
%               default value.  Allowed types are 'boolean', 'integer', 
%               'real', 'string', 'cell' or 'struct'.  For example,
%                       {'key1' 'string' { 'string1' 'string2' } 'defaultval_key1'}
%                       {'key2' {'real' 'integer'} { minint maxint } 'defaultval_key2'} 
%  callingfunc - Calling function name for error messages. {default: none}.
%  mode        - ['ignore'|'error'] ignore keywords that are either not specified 
%                in the fieldlist cell array or generate an error. 
%                {default: 'error'}.
%  verbose     - ['verbose', 'quiet'] print information. Default: 'verbose'.
%
% Outputs:
%   result     - If no error, structure with 'key' as fields and 'value' as 
%                content. If error this output contain the string error.
%   varargin   - residual varagin containing unrecognized input arguments.
%                Requires mode 'ignore' above.
%
% Note: In case of error, a string is returned containing the error message
%       instead of a structure.
%
% Example (insert the following at the beginning of your function):
%	result = finputcheck(varargin, ...
%               { 'title'         'string'   []       ''; ...
%                 'percent'       'real'     [0 1]    1 ; ...
%                 'elecamp'       'integer'  [1:10]   [] });
%   if isstr(result)
%       error(result);
%   end
%
% Note: 
%   The 'title' argument should be a string. {no default value}
%   The 'percent' argument should be a real number between 0 and 1. {default: 1}
%   The 'elecamp' argument should be an integer between 1 and 10 (inclusive).
%
%   Now 'g.title' will contain the title arg (if any, else the default ''), etc.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode, verbose )

	if nargin < 2
		help finputcheck;
		return;
	end;
	if nargin < 3
		callfunc = '';
	else 
		callfunc = [callfunc ' ' ];
	end;
    if nargin < 4
        mode = 'do not ignore';
    end;
    if nargin < 5
        verbose = 'verbose';
    end;
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
	
	varargnew = {};
	% create structure
	% ----------------
	if ~isempty(vararg)
        if isstruct(vararg)
            g = vararg;
        else
            for index=1:length(vararg)
                if iscell(vararg{index})
                    vararg{index} = {vararg{index}};
                end;
            end;
            try
                g = struct(vararg{:});
            catch
                vararg = removedup(vararg, verbose);
                try
                    g = struct(vararg{:});
                catch
                    g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
                end;
            end;
        end;
	else 
		g = [];
	end;
	
	for index = 1:size(fieldlist,NAME)
		% check if present
		% ----------------
		if ~isfield(g, fieldlist{index, NAME})
			g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
		end;
		tmpval = getfield( g, {1}, fieldlist{index, NAME});
		
		% check type
		% ----------
        if ~iscell( fieldlist{index, TYPE} )
            res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
            if isstr(res), g = res; return; end;
        else 
            testres = 0;
            tmplist = fieldlist;
            for it = 1:length( fieldlist{index, TYPE} )
                if ~iscell(fieldlist{index, VALS})
                     res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}, tmpval, callfunc );
                else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}{it}, tmpval, callfunc );
                end;
                if ~isstr(res{it}), testres = 1; end;
            end;
            if testres == 0,
                g = res{1};
                for tmpi = 2:length(res)
                    g = [ g 10 'or ' res{tmpi} ];
                end;
                return; 
            end;
        end;
	end;
    
    % check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	for index=1:length(allfields)
		if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
			if ~strcmpi(mode, 'ignore')
				g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
			end;
			varargnew{end+1} = allfields{index};
			varargnew{end+1} = getfield(g, {1}, allfields{index});
		end;
	end;


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
    g = [];
    
    switch fieldtype
     case { 'integer' 'real' 'boolean' 'float' }, 
      if ~isnumeric(tmpval) && ~islogical(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
      end;
      if strcmpi(fieldtype, 'boolean')
          if tmpval ~=0 && tmpval ~= 1
              g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
          end;  
      else 
          if strcmpi(fieldtype, 'integer')
              if ~isempty(fieldval)
                  if (any(isnan(tmpval(:))) && ~any(isnan(fieldval))) ...
                          && (~ismember(tmpval, fieldval))
                      g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                  end;
              end;
          else % real or float
              if ~isempty(fieldval) && ~isempty(tmpval)
                  if any(tmpval < fieldval(1)) || any(tmpval > fieldval(2))
                      g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                  end;
              end;
          end;
      end;  
      
      
     case 'string'
      if ~isstr(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
      end;
      if ~isempty(fieldval)
          if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
              g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
          end;
      end;

      
     case 'cell'
      if ~iscell(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
      end;
      
      
     case 'struct'
      if ~isstruct(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
      end;
      
      
     case '';
     otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
    end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella, verbose)
% make sure if all the values passed to unique() are strings, if not, exist
%try
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        myfprintf(verbose,'Note: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
%catch
    % some elements of cella were not string
%    error('some ''key'' values are not string.');
%end;    

function myfprintf(verbose, varargin)

if strcmpi(verbose, 'verbose')
    fprintf(varargin{:});
end;
function idx = regexpcell(c,pat, cmds)

% idx = regexpcell(c,pat, cmds)
%
% Return indices idx of cells in c that match pattern(s) pat (regular expression).
% Pattern pat can be char or cellstr. In the later case regexpcell returns
% indexes of cells that match any pattern in pat.
%
% cmds is a string that can contain one or several of these commands:
% 'inv' return indexes that do not match the pattern.
% 'ignorecase' will use regexpi instead of regexp
% 'exact' performs an exact match (regular expression should match the whole strings in c).
% 'all' (default) returns all indices, including repeats (if several pat match a single cell in c).
% 'unique' will return unique sorted indices.
% 'intersect' will return only indices in c that match ALL the patterns in pat.
% 
% v1 Maximilien Chaumon 01/05/09
% v1.1 Maximilien Chaumon 24/05/09 - added ignorecase
% v2 Maximilien Chaumon 02/03/2010 changed input method.
%       inv,ignorecase,exact,combine are replaced by cmds

error(nargchk(2,3,nargin))
if not(iscellstr(c))
    error('input c must be a cell array of strings');
end
if nargin == 2
    cmds = '';
end
if not(isempty(regexpi(cmds,'inv', 'once' )))
    inv = true;
else
    inv = false;
end
if not(isempty(regexpi(cmds,'ignorecase', 'once' )))
    ignorecase = true;
else
    ignorecase = false;
end
if not(isempty(regexpi(cmds,'exact', 'once' )))
    exact = true;
else
    exact = false;
end
if not(isempty(regexpi(cmds,'unique', 'once' )))
    combine = 2;
elseif not(isempty(regexpi(cmds,'intersect', 'once' )))
    combine = 3;
else
    combine = 1;
end

if ischar(pat)
    pat = cellstr(pat);
end

if exact
    for i_pat = 1:numel(pat)
        pat{i_pat} = ['^' pat{i_pat} '$'];
    end
end
for i_pat = 1:length(pat)
    if ignorecase
        trouv = regexpi(c,pat{i_pat}); % apply regexp on each pattern
    else
        trouv = regexp(c,pat{i_pat}); % apply regexp on each pattern
    end
    idx{i_pat} = find(not(cellfun('isempty',trouv)));
end
if isempty(pat)
    idx = {};
end
switch combine
    case 1
        idx = [idx{:}];
    case 2
        idx = unique([idx{:}]);
    case 3
        for i_pat = 2:length(pat)
            idx{1} = intersect(idx{1},idx{i_pat});
        end
        idx = idx{1};
end
if inv % if we want to invert result, then do so.
    others = 1:numel(trouv);
    others(idx) = [];
    idx = others;
end

function s = vararg2struct(v,tag)

% s = vararg2struct(v,tag)
% 
% translate a sequence of varargin 'name', value pairs into a structure.
% substructure fields can be defined by using underscores (if tag is
% provided, another character can be used)
% 
% ex:   v = {'name','toto','size',55,'hair_style','cool','hair_color','blue'}
%       s = vararg2struct(v)
% s = 
%     name: 'toto'
%     size: 55
%     hair: [1x1 struct]
% s.hair
% ans = 
%     style: 'cool'
%     color: 'blue'
% 

if not(exist('tag','var'))
    tag = '_';
end
s = struct;
f = regexp(v(1:2:end),['[^'  regexptranslate('escape',tag) ']*'],'match');
for i_f = 1:numel(f)
    str = 's';
    for i_ff = 1:numel(f{i_f})
        str = [str '.' f{i_f}{i_ff}];
    end
    str = [str ' = v{i_f*2};'];
    eval(str);
end

function struct2ws(s,varargin)

% struct2ws(s,varargin)
%
% Description : This function returns fields of scalar structure s in the
% current workspace
% __________________________________
% Inputs : 
%   s (scalar structure array) :    a structure that you want to throw in
%                                   your current workspace.
%   re (string optional) :          a regular expression. Only fields
%                                   matching re will be returned
% Outputs :
%   No output : variables are thrown directly in the caller workspace.
%
% Examples :
% 
%   Example 1:
%     >> who
% 
%     Your variables are:
% 
%     params  
% 
%      >>struct2ws(params)
%      >> who
% 
%     Your variables are:
% 
%     blanc         grille        ratio_ecran   unitX         
%     c_map         gris          rect          unitY         
%     centre        magni_jitt    taille        window        
%     dim_grille    noir          taille_cr     zoomzoom      
%     epais_cr      params        taille_items  
%
%   Example 2:
%     >> struct2ws(params,'unit')
%     >> who
% 
%     Your variables are:
% 
%     params  unitX   unitY   
% 
%
% _____________________________________
% See also : ws2struct ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if length(s) > 1
    error('Structure should be scalar.');
end
if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = fieldnames(s);
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i})
        continue
    end
    varsmatch(end+1) = i;
end
for i = varsmatch
    assignin('caller',vars{i},s.(vars{i}));
end

function [list,sumbytes,varargout] = DIRR(chemin,varargin)

%
%       
% DIRR
% Lists all files in the current directory and sub directories
% recursively.
% 
% [LIST] = DIRR(PATH)
% Returns a structure LIST with the same fieldnames as returned 
% by LIST = DIR(PATH)
% PATH can contain wildcards * and ? after the last \ or / (filename
% filter)
% The content of each directory in PATH is listed inside its 'isdir'
% field with the same format. The 'bytes' field is NOT zero but the
% sum of all filesizes inside the directory.
% 
% [LIST,BYTES] = DIRR(PATH)
% BYTES is a structure with fields 'total' and 'dir'. 'total' is the total
% size of PATH. 'dir' is a recursive substructure that contains the
% same fields ('total' and 'dir') for the subdirectories.
% 
% [...] = DIRR(PATH,FILTER)
% Lists only files matching the string FILTER (non case sensitive
% regular expression).
% N.B.: FILTER is optional and must not be equal to a fieldname
% ('name' or 'date' ... will never be interpreted as filters)
% 
% [LIST,BYTES,FIELDOUT] = DIRR(PATH,FIELDIN, ...)
% FIELDIN is a string specifying a field (of the structure LIST) that
% will be listed in a separate cell array of strings in FIELDOUT for
% every file with absolute path at the begining of the string.
% Multiple fields can be specified.
% 
% [LIST,BYTES,FIELDOUT] = DIRR(PATH,FIELDIN,FILTER, ...)
% Only files for which FIELDIN matches FILTER will be returned.
% Multiple [FIELDIN, FILTER] couples may be specified.
% Recursion can be avoided here by setting 'isdir' filter to '0'.
% For bytes, numeric comparison will be performed.
% 
% 
% EXAMPLES :
% 
% DIRR
% Lists all files (including path) in the current directory and it's
% subdirectories recursively.
% 
% DIRR('c:\matlab6p5\work\*.m')
% Lists all M-files in the c:\matlab6p5\work directory and it's
% subdirectories recursively.
% 
% Music = DIRR('G:\Ma musique\&Styles\Reggae\Alpha Blondy')
% Returns a structure Music very similar to what DIR returns 
% but containing the information on the files stored in
% subdirectories of 'G:\Ma musique\&Styles\Reggae\Alpha Blondy'.
% The structure Music is a bit difficult to explore though.
% See next examples.
% 
% [Files,Bytes,Names] = DIRR('c:\matlab6p5\toolbox','\.mex\>','name')
% Lists all MEX-files in the c:\matlab6p5\toolbox directory in the cell
% array of strings Names (including path).
% Note the regexp syntax of the filter string.
% Bytes is a structure with fields "total" and "dir". total is the
% total size of the directory, dir is a recursive substructure with
% the same fields as bytes for the subdirectories. 
% 
% [Files,Bytes,Names] = DIRR('c:\toto'...
%       ,'name','bytes','>50000','isdir','0')
% Lists all files larger than 50000 bytes NOT recursively.
% 
% [Files,Bytes,Dates] = DIRR('c:\matlab6p5\work','date','2005')
% Lists all dates of files from year 2005. (With path in front of
% date in the cell array of strings Dates)
% 
% 
% 
%       v1.03 % ignore files starting with .
%       Maximilien Chaumon
%       maximilien.chaumon@gmail.com
%       2012 02 22
%       

verbose = 0;
% set to 1 to get folders list in command window

if nargin == 0
    chemin = cd;
end
if nargout == 0
    dum = varargin;
    varargin{1} = 'name';
    varargin = [varargin(1) dum];
end

fields = {'name' 'date' 'bytes' 'isdir'};

if regexp(chemin,'[\*\?]') % if chemin contains any ? or *
    filt = regexprep(chemin,'.*[\\/](.*\>)','$1');% get filter
    filt = regexprep(filt,'\.','\.');% in regexp format
    filt = regexprep(filt,'\*','.*');
    filt = regexprep(filt,'\?','.');
    filt = regexprep(filt,'(.*)','\\<$1');
    chemin = regexprep(chemin,'(.*)[\\/].*\>','$1');% and chemin
end

if not(isempty(varargin)) % if additional fields were provided after chemin
    for i = 1:length(fields)
        if strcmp(varargin{1},fields{i})% if first varargin matches a fieldname,
            % assume no filter was provided,
            
            if not(exist('filt','var'))% or it was in chemin and was set just before
                filt = '.*';% set it to wildcard
                break
                
            end
        end
    end
    if not(exist('filt','var'))% else
        filt = varargin{1};% first varargin is the filter
        varargin(1) = [];
    end
else% if no additional fields were provided and filter was not in chemin
    if not(exist('filt','var'))
        filt = '.*';
    end
end
% determine which varargin are fieldnames
whicharefields = zeros(1,length(varargin));
for i = 1:length(varargin)
    for j = 1:length(fields)
        if strcmp(varargin{i},fields{j})
            whicharefields(i) = 1;
            break
        end
    end
end
% set f2out and f2outfilt
f2out = {}; f2outfilt = {};
idx = 0;
if not(isempty(varargin))
    for i = 1:length(varargin)
        if whicharefields(i)
            idx = idx + 1;
            f2out{idx} = varargin{i};
            f2outfilt{idx} = '';
        else % if nargin{i} is not a fieldname, assume it's a filter
            f2outfilt{idx} = varargin{i}; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%% START
if verbose
    disp(chemin);
end

list = dir(chemin);
if isempty(list)
    disp([chemin ' not found']);
    if nargout == 0
        clear list
    else
        for i = 1:nargout - 2
            varargout{i} = [];
        end
        sumbytes = 0;
    end
    return
end
% remove . and ..
i_file = 1;
while i_file <= length(list)
    if strcmp(list(i_file).name(1),'.')
        list(i_file) = [];
    else
        i_file = i_file + 1;
    end
end

% set sumbytes
sumbytes = struct('total',0,'dir',{});
sumbytes(1).total = 0;
i_dir = 0;
% and all output fields
for i = 1:size(f2out,2)
    f2out{2,i} = {};
end
filenames = {};
todel = 0;
r = 1;
for i_out = 1:size(f2out,2)
    if strcmp(f2out{1,i_out},'isdir')
        if strcmp(f2outfilt{i_out},'0') % check if no recursion is wanted
            r = 0;
        end
    end
end

% for each item in list
for i_file = 1:length(list)
    for i_out = 1:size(f2out,2) % for every output field
        if not(isempty(f2outfilt{i_out}))% if there is a filter
            if strcmp(f2out{1,i_out},'bytes') % if field is 'bytes'
                line = [num2str(list(i_file).(f2out{1,i_out})) f2outfilt{i_out} ';']; % compare with filter numerically
                if eval(line)% if passes the filter
                    continue % continue to next field
                else
                    todel(end+1) = i_file; % else set to be deleted
                end
            elseif not(strcmp(f2out{1,i_out},'isdir'))% if field is 'name' or 'date'
                if regexpi(list(i_file).(f2out{1,i_out}),f2outfilt{i_out}) % apply filter
                    continue % continue to next field
                else
                    todel(end+1) = i_file; % else set to be deleted
                end
            end
        end
    end
    % once checked for every field's filter
    if todel(end) == i_file % if one didn't pass,
        if not(list(i_file).isdir) % and it's not a directory
            continue % skip this file and continue
        end
    else
        if regexpi(list(i_file).name,filt) % else, check for general filter on filename
            sumbytes(1).total = sumbytes(1).total + list(i_file).bytes; % sum bytes of that level
            for i_out = 1:size(f2out,2)% and assign all output fields with the values of that file
                f2out{2,i_out}{end+1} = [chemin filesep num2str(list(i_file).(f2out{1,i_out}))];
            end
        else
            todel(end+1) = i_file; % else the file will be removed from the list structure
        end
    end
    if list(i_file).isdir % if it's a directory
        if not(r)
            continue
        end
        i_dir = i_dir + 1;
        cheminext = strcat(chemin,filesep,list(i_file).name);
        % get it's content by recursion
        % write the line to enter eval
        line = '[list(i_file).isdir,sumbytes.dir(i_dir)';
        for i_out = 1:size(f2out,2)% with all the requested fields as temporary variables
            line = [line ',f2outtemp{' num2str(i_out) '}'];
        end
        line = [line '] = dirr(cheminext,filt'];
        for i_out = 1:size(f2out,2)
            line = [line ',f2out{1,' num2str(i_out) '}'];
            if f2outfilt{i_out}
                line = [line ',f2outfilt{' num2str(i_out) '}'];
            end
        end
        line = [line ');'];
        eval(line);
        
        for i_out = 1:size(f2out,2)
            f2out{2,i_out} = [f2out{2,i_out} f2outtemp{i_out}]; % catenate temporary variables with f2out
        end
        % sum bytes 
        sumbytes(1).total = sumbytes(1).total + sumbytes(1).dir(i_dir).total; % that level + the next one
        list(i_file).bytes = sumbytes(1).dir(i_dir).total; % and set list(i_file).bytes to that value
        if list(i_file).bytes & todel(end) == i_file
            todel(end) = [];
        end
    end
end
todel(1) = [];
list(todel) = [];


for i_out = 1:size(f2out,2)
    varargout{i_out} = f2out{2,i_out};
end
if nargout == 0
    clear list
    disp(char(f2out{2,1}));
end
