function output = td500_EPSP_latency(pointer);
% function output = td500_EPSP_latency(spden,sem,stim_ts,baseline,numtrials);
% Latency = timepoint where spden exceeds baseline+2sd for at least 20ms
global lsnconfig NeuronData

output = zeros(length(pointer),5);
for ff=1:length(pointer),
    newname=char(NeuronData.plxname(pointer(ff))); newunit=char(NeuronData.unitname(pointer(ff)));
    load([lsnconfig.rsvp500spks,filesep,newname(1:12),'-',newunit,'-500_NeuronData.mat']);
    
    baseline=mean([spikestruct.faces_avg_epsp(1:100) spikestruct.fruit_avg_epsp(1:100) ...
        spikestruct.places_avg_epsp(1:100) spikestruct.bodyp_avg_epsp(1:100) spikestruct.objct_avg_epsp(1:100)]);
    baselineSTD=std([spikestruct.faces_avg_epsp(1:100) spikestruct.fruit_avg_epsp(1:100) ...
        spikestruct.places_avg_epsp(1:100) spikestruct.bodyp_avg_epsp(1:100) spikestruct.objct_avg_epsp(1:100)]);    
    cutoff=baseline+(3*baselineSTD);
    output(ff,1)=spinLoop(spikestruct.faces_avg_epsp,cutoff);
    output(ff,2)=spinLoop(spikestruct.fruit_avg_epsp,cutoff);
    output(ff,3)=spinLoop(spikestruct.bodyp_avg_epsp,cutoff);
    output(ff,4)=spinLoop(spikestruct.places_avg_epsp,cutoff);
    output(ff,5)=spinLoop(spikestruct.objct_avg_epsp,cutoff);
    
    clear baseline baselineSEM cutoff counter spikestruct respstructsingle graphstructsingle
end
return

function latency = spinLoop(data,cutoff)
counter=0; latency=nan;
for tt=150:500,
    if data(tt)>cutoff & counter==0,
        counter=1;
    end
    if data(tt)>cutoff & counter>0,
        counter=counter+1;
    end
    if data(tt)>cutoff & counter==20,
        latency = tt-100; % subtract 100 (to align to trace and another 20 to account for the counter)
        counter=0;
        break
    end
    if data(tt)<cutoff & counter>0,
        counter=0;
    end
end
return