%R_border=470;
%R_bulk=3/4*R_border;
%R_core=R_border/4;

R_core=50;

R_bulk_1=200;
R_bulk_2=250;

R_border_1=420;
R_border_2=470;


num_prof_points=50;
dR=R_border_2/num_prof_points;
pivobj_CDH_2_mg_better_stripes_thr8=cell(1,5);


fldr_CDH_2_mg_master_thrsh8=...
    '/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/the_THRSH_8_all_files/CDH_2_num_';
    


av_vel_profile_betterCDH_2_stripes_thr8=cell(1,5);

av_curlz_profile_betterCDH_2_stripes_thr8=cell(1,5);
av_div_profile_betterCDH_2_stripes_thr8=cell(1,5);

vxs_core_all_t_CDH_2_mg_better_stripes_thr8=cell(1,5);
vys_core_all_t_CDH_2_mg_better_stripes_thr8=cell(1,5);


vxs_bulk_all_t_CDH_2_mg_better_stripes_thr8=cell(1,5);
vys_bulk_all_t_CDH_2_mg_better_stripes_thr8=cell(1,5);

vxs_border_all_t_CDH_2_mg_better_stripes_thr8=cell(1,5);
vys_border_all_t_CDH_2_mg_better_stripes_thr8=cell(1,5);


vort_core_CDH_2_mg_better_stripes_thr8=cell(1,5);
vort_bulk_CDH_2_mg_better_stripes_thr8=cell(1,5);
vort_border_CDH_2_mg_better_stripes_thr8=cell(1,5);


for i=1:5
    
    av_vel_profile_betterCDH_2_stripes_thr8{i}=zeros(1,num_prof_points);

    av_curlz_profile_betterCDH_2_stripes_thr8{i}=zeros(1,num_prof_points);
    av_div_profile_betterCDH_2_stripes_thr8{i}=zeros(1,num_prof_points);
    
    fldr_CDH_2_mg_h=[fldr_CDH_2_mg_master_thrsh8 num2str(i) '/'];
    
    
    fls_CDH_2mg_h=dir([fldr_CDH_2_mg_h '*.txt']);
    
    data_test=importdata(...
        [fldr_CDH_2_mg_h fls_CDH_2mg_h(1).name]);
    len = length(data_test.data);
    
    
    pivobj_CDH_2_mg_better_stripes_thr8{i}.xs=zeros(len, length(fls_CDH_2mg_h));
    pivobj_CDH_2_mg_better_stripes_thr8{i}.ys=zeros(len, length(fls_CDH_2mg_h));
    pivobj_CDH_2_mg_better_stripes_thr8{i}.vxs=zeros(len, length(fls_CDH_2mg_h));
    pivobj_CDH_2_mg_better_stripes_thr8{i}.vys=zeros(len, length(fls_CDH_2mg_h));
    for t=1:length(fls_CDH_2mg_h)
        data_test=importdata(...
            [fldr_CDH_2_mg_h fls_CDH_2mg_h(t).name]);
        pivobj_CDH_2_mg_better_stripes_thr8{i}.xs(:,t) = data_test.data(:,1);
        pivobj_CDH_2_mg_better_stripes_thr8{i}.ys(:,t) = data_test.data(:,2);
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vxs(:,t) = data_test.data(:,3);
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vys(:,t) = data_test.data(:,4);
        
    end
    
    [pivX_CDH_2mg_better_thr, pivY_CDH_2mg_better_thr,...
        pivVX_CDH_2mg_better_thr, pivVY_CDH_2mg_better_thr,...
        pivDIV_CDH_2mg_better_thr, pivDIVall_CDH_2mg_better_thr,...
        pivCURLz_CDH_2mg_better_thr, pivCURLav_CDH_2mg_better_thr,...
        pivCURLz_all_CDH_2mg_better_thr, pivCURLav_all_CDH_2mg_better_thr] = ...
        piv_obj_to_matrix(pivobj_CDH_2_mg_better_stripes_thr8{i});
    
    
    x_c=CDH_2_mg_center.xs(i);
    y_c=CDH_2_mg_center.ys(i);
    
    
    cond_core=(pivobj_CDH_2_mg_better_stripes_thr8{i}.xs-x_c).^2+...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.ys-y_c).^2<=R_core^2;
    
    vxs_core_all_t_CDH_2_mg_better_stripes_thr8{i}=...
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vxs(...
        cond_core);
    
    vys_core_all_t_CDH_2_mg_better_stripes_thr8{i}=...
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vys(...
        cond_core);
    
    
    
    cond_bulk=(pivobj_CDH_2_mg_better_stripes_thr8{i}.xs-x_c).^2+...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.ys-y_c).^2<=R_bulk_2^2 & ...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.xs-x_c).^2+...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.ys-y_c).^2>R_bulk_1^2;
    
    vxs_bulk_all_t_CDH_2_mg_better_stripes_thr8{i}=...
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vxs(...
        cond_bulk);
    
    vys_bulk_all_t_CDH_2_mg_better_stripes_thr8{i}=...
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vys(...
        cond_bulk);
    
    cond_border=(pivobj_CDH_2_mg_better_stripes_thr8{i}.xs-x_c).^2+...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.ys-y_c).^2>R_border_1^2 & ...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.xs-x_c).^2+...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.ys-y_c).^2<=R_border_2^2;
    
    vxs_border_all_t_CDH_2_mg_better_stripes_thr8{i}=...
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vxs(...
        cond_border);
    
    vys_border_all_t_CDH_2_mg_better_stripes_thr8{i}=...
        pivobj_CDH_2_mg_better_stripes_thr8{i}.vys(...
        cond_border);
    
    figure
    
    
    fh=sqrt(vxs_core_all_t_CDH_2_mg_better_stripes_thr8{i}.^2+...
        vys_core_all_t_CDH_2_mg_better_stripes_thr8{i}.^2);
    h=histogram(fh(~isnan(fh)),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_core_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vabs_core_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))
    
    
    hold on
    
    fh=sqrt(vxs_bulk_all_t_CDH_2_mg_better_stripes_thr8{i}.^2+...
        vys_bulk_all_t_CDH_2_mg_better_stripes_thr8{i}.^2);
    h=histogram(fh(~isnan(fh)),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_bulk_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vabs_bulk_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))
    
    
    
    hold on
    
    fh=sqrt(vxs_border_all_t_CDH_2_mg_better_stripes_thr8{i}.^2+...
        vys_border_all_t_CDH_2_mg_better_stripes_thr8{i}.^2);
    h=histogram(fh(~isnan(fh)),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_border_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vabs_border_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))
    
    
    vort_core_CDH_2_mg_better_stripes_thr8{i}=pivCURLz_all_CDH_2mg_better_thr(cond_core);  %pivobj_CDH_2_mg_better{i}.loaded_vort
    vort_bulk_CDH_2_mg_better_stripes_thr8{i}=pivCURLz_all_CDH_2mg_better_thr(cond_bulk); %  pivobj_CDH_2_mg_better{i}.loaded_vort
    vort_border_CDH_2_mg_better_stripes_thr8{i}=pivCURLz_all_CDH_2mg_better_thr(cond_border); %pivobj_CDH_2_mg_better{i}.loaded_vort
    
    figure
    
    
    fh=vort_core_CDH_2_mg_better_stripes_thr8{i};
    h=histogram(fh(~isnan(fh)),...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_core_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vort_core_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))

    
    hold on
    
    fh=vort_bulk_CDH_2_mg_better_stripes_thr8{i};
    h=histogram(fh(~isnan(fh)),...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_bulk_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vort_bulk_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))

    
    hold on
    
    fh=vort_border_CDH_2_mg_better_stripes_thr8{i};
    h=histogram(fh(~isnan(fh)),...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_border_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vort_border_stripes_thr8_betterCDH_2_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))

    R_positions=sqrt((pivobj_CDH_2_mg_better_stripes_thr8{i}.xs-x_c).^2 + ...
        (pivobj_CDH_2_mg_better_stripes_thr8{i}.ys-y_c).^2);
    absv_all=sqrt(pivobj_CDH_2_mg_better_stripes_thr8{i}.vxs.^2+pivobj_CDH_2_mg_better_stripes_thr8{i}.vys.^2);
    abscurlz_all=pivCURLz_all_CDH_2mg_better_thr;
    absdiv_all=pivDIVall_CDH_2mg_better_thr;
    
    for j=1:num_prof_points
        
        cond_RdR=(R_positions<(j*dR)) & (R_positions>=((j-1)*dR));
        
        av_vel_profile_betterCDH_2_stripes_thr8{i}(j)=nanmean(absv_all(cond_RdR));
        
        av_curlz_profile_betterCDH_2_stripes_thr8{i}(j)=nanmean(abs(abscurlz_all(cond_RdR)));
        av_div_profile_betterCDH_2_stripes_thr8{i}(j)=nanmean((absdiv_all(cond_RdR)));
    end
    
    
    figure;
    plot((1:num_prof_points)*dR,av_vel_profile_betterCDH_2_stripes_thr8{i})
    hold on
    
    write_2_column_table(['velocity_profile_stripes_thr8_betterCDH_2_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_vel_profile_betterCDH_2_stripes_thr8{i})
    
    
    
    
    xlabel('r')
    ylabel('mean v')
    legend('Location','best')
    
    title('CDH 2 mg')
    
    set(gca,'FontSize', 14)
    
    figure
    plot((1:num_prof_points)*dR,av_curlz_profile_betterCDH_2_stripes_thr8{i},...
        'DisplayName',['set ' num2str(i)])
    hold on
    
    write_2_column_table(['better_absvort_profile_stripes_thr8_betterCDH_2_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_curlz_profile_betterCDH_2_stripes_thr8{i})
    
    figure
    
    plot((1:num_prof_points)*dR,av_div_profile_betterCDH_2_stripes_thr8{i},...
        'DisplayName',['set ' num2str(i)])
    hold on
    
    write_2_column_table(['better_div_profile_stripes_thr8_betterCDH_2_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_div_profile_betterCDH_2_stripes_thr8{i})
    
    
    
    
end
