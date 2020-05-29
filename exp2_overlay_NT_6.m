

master_folder_circ_cells_exp2=...
    '/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/Another_DATA_set/';

fig=figure; %('Position', [10 10 1200 600]);

master_fldr_out_exp2=...
    '/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/Another_DATA_set/';
    %'/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/the_THRSH_8_all_files/NT_6_num_';
%'/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/THRSH_8_Fields_Movies/';

for i=1:5
    
    %R_border=470;
    %R_bulk=3/4*R_border;
    %R_core=R_border/4;
    
    th=0:0.01:2*pi;
    x_c=exp2_NT_6_mg_center.xs(i);
    y_c=exp2_NT_6_mg_center.ys(i);
    
    loc_fldr=['exp2_BF_NT_6_' num2str(i) '/'];
    
    the_folder=[master_folder_circ_cells_exp2 loc_fldr];
    
    the_tif_files_h=dir([the_folder '*tif']);
    
    
    
    
%     [~, ~,...
%         ~, ~,...
%         ~, pivDIVall_NT_6mg_blured_thr10,...
%         pivCURLz_NT_6mg_blured_thr10, pivCURLav_NT_6mg_blured_thr10,...
%         pivCURLz_all_NT_6mg_blured_thr10, pivCURLav_all_NT_6mg_blured_thr10] = ...
%         piv_obj_to_matrix(exp2_pivobj_NT_6_mg_blured_stripes_thr10{i});
    
    
    for t=70:89
        
        im_NT_6_mg_h=imread([the_folder the_tif_files_h(t).name]);
        %subplot(1,2,1)
        imshow(im_NT_6_mg_h)
        
        hold on
        
        quiver(exp2_pivobj_NT_6_mg_blured_stripes_thr10{i}.xs(:,t-70+1),...
            exp2_pivobj_NT_6_mg_blured_stripes_thr10{i}.ys(:,t-70+1),...
            exp2_pivobj_NT_6_mg_blured_stripes_thr10{i}.vxs(:,t-70+1),...
            exp2_pivobj_NT_6_mg_blured_stripes_thr10{i}.vys(:,t-70+1))
        
        
        %plot(x_c+R_bulk*cos(th),...
        %    y_c+R_bulk*sin(th),'y')
        %plot(x_c+R_border*cos(th),...
        %    y_c+R_border*sin(th),'r')
        %plot(x_c+R_core*cos(th),...
         %   y_c+R_core*sin(th),'c')
        
        
        drawnow
        
        saveas(fig,...
            [the_folder ...
            'exp2_field_velocity_thr10_bluredNT_6_mg_'...
            num2str(i) '_t_' num2str(t,'%04d') '.tif'])
        
        hold off
        
        
    end
end