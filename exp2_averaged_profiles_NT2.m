
R_border_2=470;


num_prof_points=50;
dR=R_border_2/num_prof_points;


accum_vel_prof_e2_NT2=zeros(length(av_vel_profile_bluredNT_2_stripes_thr10_exp2)...
    ,length(av_vel_profile_bluredNT_2_stripes_thr10_exp2{1}));

accum_vort_prof_e2_NT2=zeros(length(av_curlz_profile_bluredNT_2_stripes_thr10_exp2)...
    ,length(av_curlz_profile_bluredNT_2_stripes_thr10_exp2{1}));



for i=1:5
    
    accum_vel_prof_e2_NT2(i,:)=...
        av_vel_profile_bluredNT_2_stripes_thr10_exp2{i};
    accum_vort_prof_e2_NT2(i,:)=...
        av_curlz_profile_bluredNT_2_stripes_thr10_exp2{i};    
end


subplot(1,2,1)
errorbar((1:num_prof_points)*dR, nanmean(accum_vel_prof_e2_NT2),...
    nanstd(accum_vel_prof_e2_NT2))



subplot(1,2,2)

errorbar((1:num_prof_points)*dR, nanmean(accum_vort_prof_e2_NT2),...
    nanstd(accum_vort_prof_e2_NT2))

write_3_column_table('av_vel_profile_exp2_NT2.dat',(1:num_prof_points)*dR, nanmean(accum_vel_prof_e2_NT2),...
nanstd(accum_vel_prof_e2_NT2))
write_3_column_table('av_vort_profile_exp2_NT2.dat',(1:num_prof_points)*dR, nanmean(accum_vort_prof_e2_NT2),...
nanstd(accum_vort_prof_e2_NT2))
