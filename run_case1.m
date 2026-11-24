function result=run_case1(method,resolution_m,seed,output_directory)
%RUN_CASE1 Run the three Case 1 routing requests (six candidates each).
% Use a fresh MATLAB R2024a session. Example:
%   result=run_case1('GFDHS',0.05,1,'results/gfdhs_h005_s1');
if nargin<1,method='GFDHS';end
if nargin<2,resolution_m=0.05;end
if nargin<3,seed=1;end
root=fileparts(mfilename('fullpath'));
method=validatestring(method,{'GFDHS','GBH','Astar'});
validateattributes(resolution_m,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
h=round(resolution_m/0.01);
assert(ismember(h,[1 5 10 15 20 25]) && abs(h*.01-resolution_m)<1e-12,...
    'Case1:Resolution','Use 0.01, 0.05, 0.10, 0.15, 0.20 or 0.25 m.');
if nargin<4,output_directory=fullfile(root,'results',sprintf('%s_h%g_seed%d',method,h,seed));end
output_directory=char(output_directory);
if ~is_absolute(output_directory),output_directory=fullfile(pwd,output_directory);end
assert(~isfolder(output_directory)&&~isfile(output_directory),'Case1:OutputExists','Choose a new output directory; existing runs are never overwritten.');
mkdir(output_directory);
oldpath=path;oldrng=rng;oldtimers=timerfindall;
names={'MULTIRES_ALG','MULTIRES_H','GFDHS_STEP_MAX','CASE1_TRIAL_SEED','CASE1_PIPE_DIAMETERS',...
    'MULTIDIAM_CONFIG','MULTIDIAM_PHASE','BATCH_OUTPUT','BATCH_WARMUP'};
oldenv=cellfun(@getenv,names,'UniformOutput',false);
cleanup=onCleanup(@()restore_state(oldpath,oldrng,oldtimers,names,oldenv));
values={method,num2str(h),num2str(h),num2str(seed),'10,10,10',...
    sprintf('%s_h%g_seed%d',method,h,seed),'example',output_directory,'0'};
for k=1:numel(names),setenv(names{k},values{k});end
% Exactly one implementation directory: avoid collisions between native helpers.
restoredefaultpath;addpath(root,fullfile(root,'src',method));
clear case1_native_entry native_search_observe aimFcn_1 aimFcncell aimFcn_astar
config=struct('method',method,'resolution_m',resolution_m,'seed',seed,...
    'coordinate_unit_m',0.01,'candidate_budget',18,'outer_optimization_generations',0,...
    'additional_physical_clearance_m',0,'energy_dsafe_m',0.03,'version','Case1-source-v0.1.0');
write_json(fullfile(output_directory,'config.json'),config);
clock=tic;
try
    result=case1_native_entry();
    result.config=config;result.wall_time_s=toc(clock);result.routes={};
    result.runtime=struct('version',version,'computer',computer);
    result.candidate_count=sum(cellfun(@(c)c.is_search,result.calls));
    assert(result.candidate_count==18,'Case1:CandidateBudget','Expected exactly 18 constructor calls.');
    assert(numel(result.physical_routes)==3 && numel(result.interfaces)==3,'Case1:Incomplete','The three requests were not completed.');
    counts=zeros(1,3);first_calls={};seen=[];
    for k=1:numel(result.calls)
        c=result.calls{k};if c.is_search
            counts(c.input_request_id)=counts(c.input_request_id)+1;
            if ~ismember(c.input_request_id,seen),first_calls{end+1}=c;seen(end+1)=c.input_request_id;end
        end
    end
    assert(isequal(counts,[6 6 6]),'Case1:CandidateBudget','Expected six calls per input request.');
    for k=1:3
        v=first_calls{k};
        route=full_port_route(result.physical_routes{k},v.physical_start,v.physical_goal,v.diameter);
        [id,~]=case1_request_identity(v.physical_start,v.physical_goal,[10 10 10]);
        route.input_request_id=id;route.points_m=route.points*0.01;
        route.polyline_length_m=route.length*0.01;
        result.routes{k}=route;
        writematrix(route.points_m,fullfile(output_directory,sprintf('request_%d_points_m.csv',id)));
    end
    result.status='COMPLETE';
    result.metric_note='Exported lengths measure the sharp-vertex polylines, including endpoint connections.';
    save(fullfile(output_directory,'result.mat'),'result');
    summary=rmfield(result,{'selected_search_points','selected_search_codes','physical_routes'});
    write_json(fullfile(output_directory,'result.json'),summary);
catch ex
    write_json(fullfile(output_directory,'FAILED.json'),struct('status','FAILED','config',config,...
        'identifier',ex.identifier,'message',ex.message,'report',getReport(ex,'extended','hyperlinks','off')));
    rethrow(ex)
end
fprintf('\nCASE1 COMPLETE: %s, h=%.2f m, seed=%d, 3 requests, 18 candidates\n',method,resolution_m,seed);
end
function restore_state(oldpath,oldrng,oldtimers,names,oldenv)
        ts=timerfindall;
        for j=1:numel(ts)
            if ~any(ts(j)==oldtimers),stop(ts(j));delete(ts(j));end
        end
        path(oldpath);rng(oldrng);
        for j=1:numel(names),setenv(names{j},oldenv{j});end
end
function yes=is_absolute(p)
yes=startsWith(p,filesep)||~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'));
end
function write_json(p,v)
f=fopen(p,'w','n','UTF-8');assert(f>0,'Case1:Output','Cannot open %s',p);
c=onCleanup(@()fclose(f));fprintf(f,'%s',jsonencode(v,'PrettyPrint',true));
end
