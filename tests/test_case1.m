function test_case1
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
out=fullfile(tempdir,['case1_test_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
r=run_case1('GFDHS',.05,1,out);
assert(strcmp(r.status,'COMPLETE')&&r.candidate_count==18&&numel(r.routes)==3);
assert(all(cellfun(@(v)v.endpoints_verified,r.routes)));
assert(isequal(sort(cellfun(@(v)v.input_request_id,r.routes)),[1 2 3]));
rejected=false;
try,run_case1('GFDHS',.05,1,out);catch ex,rejected=strcmp(ex.identifier,'Case1:OutputExists');end
assert(rejected,'Case1:Test','Existing outputs must be protected.');
rejected=false;
try,run_case1('GFDHS',.03,1,[out '_invalid']);catch ex,rejected=strcmp(ex.identifier,'Case1:Resolution');end
assert(rejected,'Case1:Test','Unsupported resolution must be rejected.');
fprintf('CASE1_SMOKE_TEST_PASSED\n');
end
