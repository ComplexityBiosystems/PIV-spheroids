subplot(1,2,1)

errorbar((1:num_prof_points)*dR, nanmean([accum_vel_prof_e1_CDH6;...
accum_vel_prof_e2_CDH6]), nanstd([accum_vel_prof_e1_CDH6;...
accum_vel_prof_e2_CDH6]))

subplot(1,2,2)

errorbar((1:num_prof_points)*dR, nanmean([accum_vort_prof_e1_CDH6;...
accum_vort_prof_e2_CDH6]), nanstd([accum_vort_prof_e1_CDH6;...
accum_vort_prof_e2_CDH6]))



write_3_column_table('av_vel_profile_exp1_and_exp2_CDH6.dat',...
    (1:num_prof_points)*dR,...
    nanmean([accum_vel_prof_e1_CDH6;...
accum_vel_prof_e2_CDH6]), nanstd([accum_vel_prof_e1_CDH6;...
accum_vel_prof_e2_CDH6]))
write_3_column_table('av_vort_profile_exp1_and_exp2_CDH6.dat',...
    (1:num_prof_points)*dR, nanmean([accum_vort_prof_e1_CDH6;...
accum_vort_prof_e2_CDH6]), ...
nanstd([accum_vort_prof_e1_CDH6; accum_vort_prof_e2_CDH6]))