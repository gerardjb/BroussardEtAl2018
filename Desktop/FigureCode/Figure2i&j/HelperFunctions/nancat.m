function C=nancat(DIM, varargin)
%   X = nancat( DIM, x1, x2, x3... )
% Concatenate, with nan padding.
%
% Works just like cat. Doesn't complain if x2 is longer or shorter than x1,
% it simply pads the one that is shorter with nans. Works with arbitrary
% number of dimensions.
%
% 'padvalue', VAL  (defaults to nan) 
%        VAL specifies a single value that is used to fill unused space. 
%        It should be the same class as the items being concatenated, 
%        and it can't be a vector or matrix.
% 'alignend'  
%        aligns each item at the ends of its dimensions than
%        beginnings. (default = 0, align at beginning of each dimension)
% 2013 - now deals with cell arrays and character arrays.
%        the default padding for cells = { []  [] ... }, and chars = ' '.
%        If any concatenands are cells, the output is a cell.
%        If any of the concatenands are numeric, they are turned into cell
%        arrays with num2cell. 
%        -- if the first operand isn't a cell, but later ones are, we
%        default to padding with { [nan] }, rather than { [] }.
%
% (c) sanjay manohar 2007

if(nargin<2)
    error('syntax: X=nancat(DIMENSION, X1, X2, ...)');
end
C=varargin{1};
if(nargin==2) return ; end; % nothing to concatenate?  return

% default padding values
if      isnumeric(C) || islogical(C),   padvalue = nan;
elseif  iscell(C)                       padvalue = { [] };
elseif  ischar(C)                       padvalue = ' ';
else    error('nancat: unknown type, %s', class(C));
end

% check for options
alignend=0; remove=[];
for(i=1:length(varargin))
    if(strcmpi(varargin{i},'alignend'))
        alignend=1;
        remove=[remove i];
    end
    if strcmpi(varargin{i},'padvalue')
      padvalue=varargin{i+1};
      remove=[remove i i+1];
    end
end;  varargin(remove)=[];

% create padding of a given size
makepadding = @(sz) repmat( padvalue, sz );

for(i=2:length(varargin)) % for each operand
    v = varargin{i}; % concatenate v onto C
    if iscell(C) && isnumeric(v), v=num2cell(v); end % convert v to cell
    if iscell(v) && isnumeric(C), C=num2cell(C); padvalue=cell(padvalue); end % switch to cell mode
    sizeC=size(C); sizeV=size(v);
    if all(sizeC==0), C=v; continue; end; % earlier operands were empty?
    if    (length(sizeC)>length(sizeV)) sizeV=[sizeV ones(1,length(sizeC)-length(sizeV))]; 
    elseif(length(sizeV)>length(sizeC)) sizeC=[sizeC ones(1,length(sizeV)-length(sizeC))]; 
    end % ensure same number of dimensions.
    uneqdim = find(~(sizeV==sizeC)); % for each of the unequal dimensions
    for(j=1:length(uneqdim))
        if uneqdim(j)==DIM; continue; end % (except for the one you're catting along)
        d=uneqdim(j);
        if(sizeC(d)>sizeV(d)) % V is too small: add nans to d if it's shorter
            diff= sizeC(d)-sizeV(d); % difference in size along dimension d
            addsize=size(v); addsize(d)=diff;
            if(alignend)   v=cat(d, makepadding( addsize ), v);
            else           v=cat(d, v, makepadding( addsize ));
            end
        elseif(sizeV(d)>sizeC(d)) % C is too small
            diff= sizeV(d)-sizeC(d); % difference in size along dimension d
            addsize=size(C); addsize(d)=diff;
            if(alignend)    C=cat(d, makepadding( addsize ), C);
            else            C=cat(d, C, makepadding( addsize ));
            end
        end
    end
    C=cat(DIM,C,v);
end

    