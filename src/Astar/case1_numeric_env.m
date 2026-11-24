function values=case1_numeric_env(name,requiredCount,defaultValues)
% Strict numeric parsing; never evaluate environment text as MATLAB code.
raw=strtrim(getenv(name));
if isempty(raw)
    assert(nargin>=3,'Multidiameter:MissingEnvironment','Missing environment variable %s',name);
    values=defaultValues;return;
end
if startsWith(raw,'[') && endsWith(raw,']'),raw=strtrim(raw(2:end-1));end
parts=regexp(raw,'[,;\s]+','split');
values=str2double(parts);
assert(~isempty(values) && all(isfinite(values)) && isreal(values),...
    'Multidiameter:NumericEnvironment','%s must contain finite numbers only.',name);
values=reshape(values,1,[]);
if requiredCount>0
    assert(numel(values)==requiredCount,'Multidiameter:NumericEnvironment','%s must have %d values.',name,requiredCount);
end
end
