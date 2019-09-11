function dat_sf = proc_multiBandFunction(dat,varargin)
%PROC_MULTIBANDLINEARDERIVATION - Apply linear derivation to multi-band
% data
%
%Synopsis:
% DAT_SF = proc_multiBandLinearDerivation(DAT, func)
%
%Arguments:
% DAT - data structure of epoched and multi-band bandpass filtered data
% func- function pointer
%
%Returns:
% DAT_SF - updated data structure, containing the projected data of all
%          bands, appended as channels
%
%Description:
% The input data structure is assumed to contain pre-filtered channels as
% returned by the function proc_filterbank. proc_multiBandLinearDerivation
% can be used either alone or as a processing step in the crossvalidation
% function.
%
% See also: proc_multiBandSpatialFilter crossvalidation proc_filterbank

% 07-2019: miklody@tu-berlin.de

props = {'func'   @proc_linearDerivation          '!FUNC|CELL'
    };

if nargin==0;
    dat_sf= props; return
end

if misc_isproplist(varargin{1}),
    opt= opt_proplistToStruct(varargin{:});
else
    if nargin>2&&iscell(varargin{2})
        if nargin>3&&iscell(varargin{3})
            if nargin>4&&iscell(varargin{4})
                if nargin>5&&iscell(varargin{5})
                    opt= opt_proplistToStruct(varargin{6:end});
                    opt.arg = varargin{2};
                    opt.arg2 = varargin{3};
                    opt.arg3 = varargin{4};
                    opt.arg4 = varargin{5};
                else
                    opt= opt_proplistToStruct(varargin{5:end});
                    opt.arg = varargin{2};
                    opt.arg2 = varargin{3};
                    opt.arg3 = varargin{4};
                end
            else
                opt= opt_proplistToStruct(varargin{4:end});
                opt.arg = varargin{2};
                opt.arg2 = varargin{3};
            end
        else
            opt= opt_proplistToStruct(varargin{3:end});
            opt.arg = varargin{2};
        end
    else
        opt= opt_proplistToStruct(varargin{2:end});
    end
    opt.func= varargin{1};
    
end

[Fcn, Par]= misc_getFuncParam(opt.func);

misc_checkType(dat,'STRUCT(x clab y)');


% get number of frequency bands
band_ix = zeros(1,length(dat.clab));
flt_ix = cellfun(@(x) strfind(x,'flt'),dat.clab);
for ii = 1:length(dat.clab)
    band_ix(ii) = str2double(dat.clab{ii}(flt_ix(ii)+3:end));
end
n_bands = max(band_ix);

% band-wise apply spatial filters to data
dat_sf = [];
for bi = 1:n_bands
    dat2 = proc_selectChannels(dat,sprintf('*flt%d',bi));
    if isfield(opt,'arg')
        if isfield(opt,'arg2')
            if isfield(opt,'arg3')
                if isfield(opt,'arg4')
                    dat2 = Fcn(dat2, opt.arg{bi},opt.arg2{bi},opt.arg3{bi},opt.arg4{bi}, Par{:});
                else
                    dat2 = Fcn(dat2, opt.arg{bi},opt.arg2{bi},opt.arg3{bi}, Par{:});
                end
            else
                dat2 = Fcn(dat2, opt.arg{bi},opt.arg2{bi}, Par{:});
            end
            else
                dat2 = Fcn(dat2, opt.arg{bi}, Par{:});
            end
        else
            dat2 = Fcn(dat2, Par{:});
        end
        dat_sf = proc_appendChannels(dat_sf,dat2);
    end