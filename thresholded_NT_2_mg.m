R_border=470;
R_bulk=3/4*R_border;
R_core=R_border/4;
num_prof_points=50;
dR=R_border/num_prof_points;
pivobj_NT_2_mg_better_thr=cell(1,5);


fldr_NT_2_mg_master_better=...
    '/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/Better_Pictures/NT_2_mg_';
%NT_2_mg_1_OI/';

av_vel_profile_betterNT_2_thr=cell(1,5);

av_curlz_profile_betterNT_2_thr=cell(1,5);
av_div_profile_betterNT_2_thr=cell(1,5);

vxs_core_all_t_NT_2_mg_better_thr=cell(1,5);
vys_core_all_t_NT_2_mg_better_thr=cell(1,5);


vxs_bulk_all_t_NT_2_mg_better_thr=cell(1,5);
vys_bulk_all_t_NT_2_mg_better_thr=cell(1,5);

vxs_border_all_t_NT_2_mg_better_thr=cell(1,5);
vys_border_all_t_NT_2_mg_better_thr=cell(1,5);


vort_core_NT_2_mg_better_thr=cell(1,5);
vort_bulk_NT_2_mg_better_thr=cell(1,5);
vort_border_NT_2_mg_better_rep=cell(1,5);


for i=1:5
    
    av_vel_profile_betterNT_2_thr{i}=zeros(1,num_prof_points);

    av_curlz_profile_betterNT_2_thr{i}=zeros(1,num_prof_points);
    av_div_profile_betterNT_2_thr{i}=zeros(1,num_prof_points);
    
    fldr_NT_2_mg_h=[fldr_NT_2_mg_master_better num2str(i) '_OI' '/'];
    
    
    fls_NT_2mg_h=dir([fldr_NT_2_mg_h 'thr*.txt']);
    
    data_test=importdata(...
        [fldr_NT_2_mg_h fls_NT_2mg_h(1).name]);
    len = length(data_test.data);
    
    
    pivobj_NT_2_mg_better_thr{i}.xs=zeros(len, length(fls_NT_2mg_h));
    pivobj_NT_2_mg_better_thr{i}.ys=zeros(len, length(fls_NT_2mg_h));
    pivobj_NT_2_mg_better_thr{i}.vxs=zeros(len, length(fls_NT_2mg_h));
    pivobj_NT_2_mg_better_thr{i}.vys=zeros(len, length(fls_NT_2mg_h));
    for t=1:length(fls_NT_2mg_h)
        data_test=importdata(...
            [fldr_NT_2_mg_h fls_NT_2mg_h(t).name]);
        pivobj_NT_2_mg_better_thr{i}.xs(:,t) = data_test.data(:,1);
        pivobj_NT_2_mg_better_thr{i}.ys(:,t) = data_test.data(:,2);
        pivobj_NT_2_mg_better_thr{i}.vxs(:,t) = data_test.data(:,3);
        pivobj_NT_2_mg_better_thr{i}.vys(:,t) = data_test.data(:,4);
        
    end
    
    [pivX_NT_2mg_better_thr, pivY_NT_2mg_better_thr,...
        pivVX_NT_2mg_better_thr, pivVY_NT_2mg_better_thr,...
        pivDIV_NT_2mg_better_thr, pivDIVall_NT_2mg_better_thr,...
        pivCURLz_NT_2mg_better_thr, pivCURLav_NT_2mg_better_thr,...
        pivCURLz_all_NT_2mg_better_thr, pivCURLav_all_NT_2mg_better_thr] = ...
        piv_obj_to_matrix(pivobj_NT_2_mg_better_thr{i});
    
    
    x_c=NT_2_mg_center.xs(i);
    y_c=NT_2_mg_center.ys(i);
    
    
    cond_core=(pivobj_NT_2_mg_better_thr{i}.xs-x_c).^2+...
        (pivobj_NT_2_mg_better_thr{i}.ys-y_c).^2<=R_core^2;
    
    vxs_core_all_t_NT_2_mg_better_thr{i}=...
        pivobj_NT_2_mg_better_thr{i}.vxs(...
        cond_core);
    
    vys_core_all_t_NT_2_mg_better_thr{i}=...
        pivobj_NT_2_mg_better_thr{i}.vys(...
        cond_core);
    
    
    
    cond_bulk=(pivobj_NT_2_mg_better_thr{i}.xs-x_c).^2+...
        (pivobj_NT_2_mg_better_thr{i}.ys-y_c).^2<=R_bulk^2 & ...
        (pivobj_NT_2_mg_better_thr{i}.xs-x_c).^2+...
        (pivobj_NT_2_mg_better_thr{i}.ys-y_c).^2>R_core^2;
    
    vxs_bulk_all_t_NT_2_mg_better_thr{i}=...
        pivobj_NT_2_mg_better_thr{i}.vxs(...
        cond_bulk);
    
    vys_bulk_all_t_NT_2_mg_better_thr{i}=...
        pivobj_NT_2_mg_better_thr{i}.vys(...
        cond_bulk);
    
    cond_border=(pivobj_NT_2_mg_better_thr{i}.xs-x_c).^2+...
        (pivobj_NT_2_mg_better_thr{i}.ys-y_c).^2>R_bulk^2 & ...
        (pivobj_NT_2_mg_better_thr{i}.xs-x_c).^2+...
        (pivobj_NT_2_mg_better_thr{i}.ys-y_c).^2<=R_border^2;
    
    vxs_border_all_t_NT_2_mg_better_thr{i}=...
        pivobj_NT_2_mg_better_thr{i}.vxs(...
        cond_border);
    
    vys_border_all_t_NT_2_mg_better_thr{i}=...
        pivobj_NT_2_mg_better_thr{i}.vys(...
        cond_border);
    
    figure
    
    h=histogram(sqrt(vxs_core_all_t_NT_2_mg_better_thr{i}.^2+...
        vys_core_all_t_NT_2_mg_better_thr{i}.^2),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_core_thr_betterNT_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    
    
    
    hold on
    h=histogram(sqrt(vxs_bulk_all_t_NT_2_mg_better_thr{i}.^2+...
        vys_bulk_all_t_NT_2_mg_better_thr{i}.^2),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_bulk_thr_betterNT_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    
    
    hold on
    h=histogram(sqrt(vxs_border_all_t_NT_2_mg_better_thr{i}.^2+...
        vys_border_all_t_NT_2_mg_better_thr{i}.^2),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_border_thr_betterNT_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    
    
    vort_core_NT_2_mg_better_thr{i}=pivCURLz_all_NT_2mg_better_thr(cond_core);  %pivobj_NT_2_mg_better{i}.loaded_vort
    vort_bulk_NT_2_mg_better_thr{i}=pivCURLz_all_NT_2mg_better_thr(cond_bulk); %  pivobj_NT_2_mg_better{i}.loaded_vort
    vort_border_NT_2_mg_better_rep{i}=pivCURLz_all_NT_2mg_better_thr(cond_border); %pivobj_NT_2_mg_better{i}.loaded_vort
    
    figure
    
    h=histogram(vort_core_NT_2_mg_better_thr{i},...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_core_thr_betterNT_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    
    hold on
    
    h=histogram(vort_bulk_NT_2_mg_better_thr{i},...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_bulk_thr_betterNT_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    
    hold on
    h=histogram(vort_border_NT_2_mg_better_rep{i},...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_border_thr_betterNT_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    R_positions=sqrt((pivobj_NT_2_mg_better_thr{i}.xs-x_c).^2 + ...
        (pivobj_NT_2_mg_better_thr{i}.ys-y_c).^2);
    absv_all=sqrt(pivobj_NT_2_mg_better_thr{i}.vxs.^2+pivobj_NT_2_mg_better_thr{i}.vys.^2);
    abscurlz_all=pivCURLz_all_NT_2mg_better_thr;
    absdiv_all=pivDIVall_NT_2mg_better_thr;
    
    for j=1:num_prof_points
        
        cond_RdR=(R_positions<(j*dR)) & (R_positions>=((j-1)*dR));
        
        av_vel_profile_betterNT_2_thr{i}(j)=nanmean(absv_all(cond_RdR));
        
        av_curlz_profile_betterNT_2_thr{i}(j)=nanmean(abs(abscurlz_all(cond_RdR)));
        av_div_profile_betterNT_2_thr{i}(j)=nanmean((absdiv_all(cond_RdR)));
    end
    
    
    figure;
    plot((1:num_prof_points)*dR,av_vel_profile_betterNT_2_thr{i})
    hold on
    
    write_2_column_table(['velocity_profile_thr_betterNT_2_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_vel_profile_betterNT_2_thr{i})
    
    
    
    
    xlabel('r')
    ylabel('mean v')
    legend('Location','best')
    
    title('NT 2 mg')
    
    set(gca,'FontSize', 14)
    
    figure
    plot((1:num_prof_points)*dR,av_curlz_profile_betterNT_2_thr{i},...
        'DisplayName',['set ' num2str(i)])
    hold on
    
    write_2_column_table(['better_absvort_profile_thr_betterNT_2_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_curlz_profile_betterNT_2_thr{i})
    
    figure
    
    plot((1:num_prof_points)*dR,av_div_profile_betterNT_2_thr{i},...
        'DisplayName',['set ' num2str(i)])
    hold on
    
    write_2_column_table(['better_div_profile_thr_betterNT_2_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_div_profile_betterNT_2_thr{i})
    
    
    
    
end
