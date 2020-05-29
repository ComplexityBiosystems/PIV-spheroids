accum_vel_prof_e1_CDH6=zeros(length(av_vel_profile_moreCDH_6_stripes_thr10)...
    ,length(av_vel_profile_moreCDH_6_stripes_thr10{1}));

accum_vort_prof_e1_CDH6=zeros(length(av_curlz_profile_moreCDH_6_stripes_thr10)...
    ,length(av_curlz_profile_moreCDH_6_stripes_thr10{1}));


for i=1:5
    
    accum_vel_prof_e1_CDH6(i,:)=...
        av_vel_profile_moreCDH_6_stripes_thr10{i};
    accum_vort_prof_e1_CDH6(i,:)=...
        av_curlz_profile_moreCDH_6_stripes_thr10{i};    
end


subplot(1,2,1)
errorbar((1:num_prof_points)*dR, nanmean(accum_vel_prof_e1_CDH6),...
    nanstd(accum_vel_prof_e1_CDH6))

subplot(1,2,2)

errorbar((1:num_prof_points)*dR, nanmean(accum_vort_prof_e1_CDH6),...
    nanstd(accum_vort_prof_e1_CDH6))
