function batch_ram_event(event)
global BATCH_RAM_EVENTS BATCH_RAM_CLOCK
BATCH_RAM_EVENTS{end+1}=struct('event',event,'configuration',getenv('MULTIDIAM_CONFIG'),'pipe_diameters_input_order',case1_pipe_diameters(),'elapsed_s',toc(BATCH_RAM_CLOCK));
end
