function [avgfunc semfunc]=td500_avgSpden(pointer,spden_choice,lsnconfig)
global NeuronData lsnconfig 

% Initialize structures
graphdata=struct(...
    'faces_avg',[],'fruit_avg',[],'bodyp_avg',[],'places_avg',[],'objct_avg',[],...
    'NORM_faces_avg',[],'NORM_fruit_avg',[],'NORM_bodyp_avg',[],'NORM_places_avg',[],'NORM_objct_avg',[]);
avgfunc=[];

% Load graphdata
for ff=1:length(pointer),
    newname=char(NeuronData.plxname(pointer(ff))); newunit=char(NeuronData.unitname(pointer(ff)));
    load([lsnconfig.rsvp500spks,filesep,newname(1:12),'-',newunit,'-500_NeuronData.mat']);
    if spden_choice==1,
        face_field='faces_avg';
        fruit_field='fruit_avg';
        bodyp_field='bodyp_avg';
        places_field='places_avg';
        objct_field='objct_avg';
    elseif spden_choice==2,
        face_field='faces_avg_smlK';
        fruit_field='fruit_avg_smlK';
        bodyp_field='bodyp_avg_smlK';
        places_field='places_avg_smlK';
        objct_field='objct_avg_smlK';
    elseif spden_choice==3,
        face_field='faces_avg_epsp';
        fruit_field='fruit_avg_epsp';
        bodyp_field='bodyp_avg_epsp';
        places_field='places_avg_epsp';
        objct_field='objct_avg_epsp';
    end
    graphdata.faces_avg=[graphdata.faces_avg;getfield(spikestruct,face_field)];
    graphdata.fruit_avg=[graphdata.fruit_avg;getfield(spikestruct,fruit_field)];
    graphdata.bodyp_avg=[graphdata.bodyp_avg;getfield(spikestruct,bodyp_field)];
    graphdata.places_avg=[graphdata.places_avg;getfield(spikestruct,places_field)];
    graphdata.objct_avg=[graphdata.objct_avg;getfield(spikestruct,objct_field)];
       
    normalizer=max([max(graphdata.faces_avg) max(graphdata.fruit_avg) max(graphdata.bodyp_avg)...
        max(graphdata.bodyp_avg) max(graphdata.bodyp_avg)]);
      
    graphdata.NORM_faces_avg=[graphdata.NORM_faces_avg;(getfield(spikestruct,face_field)/normalizer)];
    graphdata.NORM_fruit_avg=[graphdata.NORM_fruit_avg;(getfield(spikestruct,fruit_field)/normalizer)];
    graphdata.NORM_bodyp_avg=[graphdata.NORM_bodyp_avg;(getfield(spikestruct,bodyp_field)/normalizer)];
    graphdata.NORM_places_avg=[graphdata.NORM_places_avg;(getfield(spikestruct,places_field)/normalizer)];
    graphdata.NORM_objct_avg=[graphdata.NORM_objct_avg;(getfield(spikestruct,objct_field)/normalizer)];

    clear spikestruct
end

avgfunc(1,:)= mean(graphdata.faces_avg,1);
avgfunc(2,:)= mean(graphdata.fruit_avg,1);
avgfunc(3,:)= mean(graphdata.bodyp_avg,1);
avgfunc(4,:)= mean(graphdata.places_avg,1);
avgfunc(5,:)= mean(graphdata.objct_avg,1);

avgfunc(6,:)= mean(graphdata.NORM_faces_avg,1);
avgfunc(7,:)= mean(graphdata.NORM_fruit_avg,1);
avgfunc(8,:)= mean(graphdata.NORM_bodyp_avg,1);
avgfunc(9,:)= mean(graphdata.NORM_places_avg,1);
avgfunc(10,:)=mean(graphdata.NORM_objct_avg,1);

semfunc(1,:)= sem(graphdata.faces_avg);
semfunc(2,:)= sem(graphdata.fruit_avg);
semfunc(3,:)= sem(graphdata.bodyp_avg);
semfunc(4,:)= sem(graphdata.places_avg);
semfunc(5,:)= sem(graphdata.objct_avg);

semfunc(6,:)= sem(graphdata.NORM_faces_avg);
semfunc(7,:)= sem(graphdata.NORM_fruit_avg);
semfunc(8,:)= sem(graphdata.NORM_bodyp_avg);
semfunc(9,:)= sem(graphdata.NORM_places_avg);
semfunc(10,:)=sem(graphdata.NORM_objct_avg);

return

