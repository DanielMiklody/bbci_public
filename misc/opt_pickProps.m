function [opt]= opt_pickProps(opt, props)
%OPT_PICKPROPS picks the props according to a proplist
%
%Synopsis:
%  [OPT, ISDEFAULT]= opt_setDefaults(OPT, PROPSPEC)
%
%Arguments:
%  OPT:      STRUCT of optional properties
%  PROPSPEC: PROPSPECLIST - Property specification list, i.e., CELL of size
%      [N 2] or [N 3], with the first column all being strings.
%
%Returns:
%  OPT: STRUCT with added properties from the property specification list
%          PROPSPEC
% 09-2019 Daniel Miklody

misc_checkType(opt, 'STRUCT|CHAR');
misc_checkType(props, 'PROPSPEC');

if ~isempty(opt),
  for Fld=fieldnames(opt)',
    fld= Fld{1};
    idx= find(strcmpi(fld, props(:,1)));
    if length(idx)>1,
      error('Case ambiguity in propertylist');
    elseif isempty(idx),
        %remove field if not found
      opt= rmfield(opt, fld);
    end
  end
end