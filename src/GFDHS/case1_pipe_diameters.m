function diameters=case1_pipe_diameters()
diameters=case1_numeric_env('CASE1_PIPE_DIAMETERS',3);
assert(all(diameters>0),'Multidiameter:Diameter','All physical diameters must be positive.');
end
