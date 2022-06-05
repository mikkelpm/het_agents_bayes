function varargout = plot_depfun(foo,plot_name,varargin)
%plot_depfun(foo,varargin)
%plots a tree of the dependencies of function foo
%ignores built-in function, any any function names given in varargin
%
%see also    : depfun, plot_subfun (file exchange ID 46070)



opt.ignore  = varargin;
opt.me      = which(foo);
if isempty(opt.me); error('this file could not be found'); end
opt.me      = sub_fileparts(opt.me);
opt         = sub_deps(opt);
opt.us      = sub_fileparts2(opt.us,opt.me);
opt.us.color = ones(size(opt.us.short));
file_subfolder = opt.us.subfolder;
all_subfolder = ...
    {'firm_model\','firm_model\auxiliary_functions\likelihood\','firm_model\auxiliary_functions\sim\','firm_model\dynare\',...
    'functions\','functions\likelihood\','functions\mcmc\','functions\sim\',...
    'hh_model\','hh_model\auxiliary_functions\likelihood\','hh_model\auxiliary_functions\sim\','hh_model\dynare\'};
n_subfolder = length(all_subfolder);
for i_subfolder = 1:n_subfolder
    opt.us.color(strcmp(file_subfolder,all_subfolder{i_subfolder})) = i_subfolder+1;
end
opt.us.color(~contains(opt.us.fol,opt.me.fol)) = n_subfolder+2;
sub_plot(opt,plot_name);

if nargout;
    varargout{1} = opt;
end
end

function out = sub_fileparts(foos)
if iscell(foos)
    for i=1:numel(foos)
        if ~exist(foos{i},'file')
            disp(foos{i});
            error('this file does not exist');
        end
    end
end

fols = regexprep(foos,'[^/|\\]+$','');

%short name of each function, and the folder it occurs
short = regexp(foos,'[^/|\\]+$','match','once');
short = regexprep(short,'\..+','');

out.full  = foos;
out.fol   = fols;
out.short = short;
end

function out = sub_fileparts2(foos,me)
if iscell(foos)
    for i=1:numel(foos)
        if ~exist(foos{i},'file')
            disp(foos{i});
            error('this file does not exist');
        end
    end
end

fols = regexprep(foos,'[^/|\\]+$','');

%short name of each function, and the folder it occurs
short = regexp(foos,'[^/|\\]+$','match','once');
short = regexprep(short,'\..+','');
is_subfolder = contains(fols,me.fol) & ~strcmp(fols,me.fol);
subfolder = repmat({''},size(short));
subfolder(is_subfolder) = erase(fols(is_subfolder),me.fol);

out.full  = foos;
out.fol   = fols;
out.short = short;
out.subfolder = subfolder;
end

function opt = sub_deps(opt)
%find dependencies of opt.me
%uses recursive calls to depfun with -toponly
%culls builtin functions for speed

names = {opt.me.full}; %list of all files found (will grow)
done = false;          %which files have been examined
from = [];             %dependency - parent
to   = [];             %dependency - child
t = now;
while any(~done)
    for i=find(~done)
        if verLessThan('matlab','8.3');
            new = depfun(names{i},'-toponly','-quiet')';
            %remove any that are built in
            keep = cellfun('isempty',strfind(new,matlabroot));
        else
            new = matlab.codetools.requiredFilesAndProducts(names{i},'toponly');
            %built-in already removed
            keep = true(size(new));
        end
%         keep = ~cellfun('isempty',strfind(new,'C:\Users\Laura\Dropbox\projects\sandbox\version_track\winberry_20180719_assemble'));
        %catch any strange return sizes from other os/versions
        if size(new,1)>1; new = new'; end
        %remove self
        keep(ismember(new,names{i})) = false;
        %remove ignored : full filename
        keep(ismember(new,opt.ignore)) = false;
        %remove ignored : short filename
        short = regexp(new,'[^/|\\]+$','match','once');
        keep(ismember(short,opt.ignore)) = false;
        %remove ignored : short filename no extension
        short = regexprep(short,'\..+','');
        keep(ismember(short,opt.ignore)) = false;
        %reduce the set of new
        new = new(keep);
                
        %add to list of names any new that are not already in it
        [~,~,I] = setxor(names,new);
        names = [names new(I)]; %#ok<AGROW>
        %rearrange I because apparently mac and pc do things differently
        if size(I,1)>1; I = I'; end
        
        %add to list of done
        done_new = ~contains(new(I),opt.me.fol);
        done = [done done_new]; %#ok<AGROW>
        done(i) = true;
        
        %new dependencies
        [~,newkid] = ismember(new,names);
        newdad = repmat(i,size(newkid));
        from = [from newdad]; %#ok<AGROW>
        to   = [to   newkid]; %#ok<AGROW>

        %report every 10 seconds
        if now-t>10;
            t = now;
            fprintf(1,'%d dependencies found so far',numel(names));
        end
    end
end

%sort names
[names,order] = sort(names);

%sort from/to references to match new names order
[~,rev] = sort(order);
from = rev(from);
to   = rev(to);

%export results
opt.us = names;
opt.from = from;
opt.to   = to;
end

function sub_plot(opt,plot_name)
if isempty(opt.from); disp('this function has no dependencies'); return; end
colors_cool = colormap(cool(8));
colors_autumn = colormap(autumn(12));
colors_summer = colormap(summer(6));
plot_graph(opt.us.short,opt.us.subfolder,opt.from,opt.to,'-colour',opt.us.color,...
    [ones(1,3); colors_cool(4:-1:1,:); colors_autumn(5:2:11,:); colors_summer(1:4,:); ones(1,3)*.7]);
saveas(gcf,['../doc/fig/' plot_name],'png')
end

%% DEVNOTES
%140329 added handling of functions with no dependencies
%150311 bugfix : was not ignoring function calls correctly if they didn't
%                have an extension.
