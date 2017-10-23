function [dat,state] = online_BeerLambert(dat,state,varargin)

if isempty(state)
      state= opt_proplistToStruct(varargin{:});
%     params = rmfield(state,'init');
    if ~isfield(dat,'wavelengths') && isfield(state,'wavelengths')
        dat.wavelengths = state.wavelengths;
    end
    [dat, state] = proc_BeerLambert(dat,state,varargin{:});
    state.dat = dat;
%     state.init = 0;
    
else
    
%     opt = opt_structToProplist(state);
%     opt = rmfield(opt,{'data','init'});
    dat = proc_BeerLambert(dat,state);
    state.dat.x = cat(1,dat.x,state.dat.x);
    
end



end