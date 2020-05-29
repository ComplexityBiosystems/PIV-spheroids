%data extraction from experimental files CDH 6

R_core=50;

R_bulk_1=150;
R_bulk_2=200;

R_border_1=300;
R_border_2=350;




num_prof_points=50;
dR=R_border_2/num_prof_points;
exp2_pivobj_CDH_6_mg_blured_stripes_thr10=cell(1,5);


fldr_CDH_6_mg_master_thrsh10_orig_exp2=...
    '/Volumes/Seagate Expansion Drive/workspace/Circular_Cells/Another_DATA_set/exp2_CDH_6mg_';
    
    


av_vel_profile_bluredCDH_6_stripes_thr10_exp2=cell(1,5);

av_curlz_profile_bluredCDH_6_stripes_thr10_exp2=cell(1,5);
av_div_profile_bluredCDH_6_stripes_thr10_exp2=cell(1,5);

vxs_core_all_t_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);
vys_core_all_t_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);


vxs_bulk_all_t_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);
vys_bulk_all_t_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);

vxs_border_all_t_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);
vys_border_all_t_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);


vort_core_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);
vort_bulk_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);
vort_border_CDH_6_mg_blured_stripes_thr10_exp2=cell(1,5);


for i=1:5
    
    av_vel_profile_bluredCDH_6_stripes_thr10_exp2{i}=zeros(1,num_prof_points);

    av_curlz_profile_bluredCDH_6_stripes_thr10_exp2{i}=zeros(1,num_prof_points);
    av_div_profile_bluredCDH_6_stripes_thr10_exp2{i}=zeros(1,num_prof_points);
    
    fldr_CDH_6_mg_h=[fldr_CDH_6_mg_master_thrsh10_orig_exp2 num2str(i) '/'];
    
    
    fls_CDH_6mg_h=dir([fldr_CDH_6_mg_h 'P*.txt']);
    
    data_test=importdata(...
        [fldr_CDH_6_mg_h fls_CDH_6mg_h(1).name]);
    len = length(data_test.data);
    
    
    exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs=zeros(len, length(fls_CDH_6mg_h));
    exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys=zeros(len, length(fls_CDH_6mg_h));
    exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vxs=zeros(len, length(fls_CDH_6mg_h));
    exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vys=zeros(len, length(fls_CDH_6mg_h));
    for t=1:length(fls_CDH_6mg_h)
        data_test=importdata(...
            [fldr_CDH_6_mg_h fls_CDH_6mg_h(t).name]);
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs(:,t) = data_test.data(:,1);
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys(:,t) = data_test.data(:,2);
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vxs(:,t) = data_test.data(:,3);
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vys(:,t) = data_test.data(:,4);
        
    end
    
    [pivX_CDH_6mg_thr10_exp2, pivY_CDH_6mg_thr10_exp2,...
        pivVX_CDH_6mg_thr10_exp2, pivVY_CDH_6mg_thr10_exp2,...
        pivDIV_CDH_6mg_thr10_exp2, pivDIVall_CDH_6mg_thr10_exp2,...
        pivCURLz_CDH_6mg_thr10_exp2, pivCURLav_CDH_6mg_thr10_exp2,...
        pivCURLz_all_CDH_6mg_thr10_exp2, pivCURLav_all_CDH_6mg_thr10_exp2] = ...
        piv_obj_to_matrix(exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i});
    
    
    x_c=exp2_CDH_6_mg_center.xs(i);
    y_c=exp2_CDH_6_mg_center.ys(i);
    
    
    cond_core=(exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs-x_c).^2+...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys-y_c).^2<=R_core^2;
    
    vxs_core_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}=...
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vxs(...
        cond_core);
    
    vys_core_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}=...
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vys(...
        cond_core);
    
    
    
    cond_bulk=(exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs-x_c).^2+...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys-y_c).^2<=R_bulk_2^2 & ...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs-x_c).^2+...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys-y_c).^2>R_bulk_1^2;
    
    vxs_bulk_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}=...
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vxs(...
        cond_bulk);
    
    vys_bulk_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}=...
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vys(...
        cond_bulk);
    
    cond_border=(exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs-x_c).^2+...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys-y_c).^2>R_border_1^2 & ...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs-x_c).^2+...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys-y_c).^2<=R_border_2^2;
    
    vxs_border_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}=...
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vxs(...
        cond_border);
    
    vys_border_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}=...
        exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vys(...
        cond_border);
    
    figure
    
    
    fh=sqrt(vxs_core_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}.^2+...
        vys_core_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}.^2);
    h=histogram(fh(~isnan(fh)),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_core_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vabs_core_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))
    
    
    hold on
    
    fh=sqrt(vxs_bulk_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}.^2+...
        vys_bulk_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}.^2);
    h=histogram(fh(~isnan(fh)),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_bulk_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vabs_bulk_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))
    
    
    
    hold on
    
    fh=sqrt(vxs_border_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}.^2+...
        vys_border_all_t_CDH_6_mg_blured_stripes_thr10_exp2{i}.^2);
    h=histogram(fh(~isnan(fh)),...
        linspace(0,25,100),'Normalization','pdf');
    write_2_column_table(['pdf_vabs_border_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vabs_border_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))
    
    
    vort_core_CDH_6_mg_blured_stripes_thr10_exp2{i}=pivCURLz_all_CDH_6mg_thr10_exp2(cond_core);  %pivobj_CDH_6_mg_blured{i}.loaded_vort
    vort_bulk_CDH_6_mg_blured_stripes_thr10_exp2{i}=pivCURLz_all_CDH_6mg_thr10_exp2(cond_bulk); %  pivobj_CDH_6_mg_blured{i}.loaded_vort
    vort_border_CDH_6_mg_blured_stripes_thr10_exp2{i}=pivCURLz_all_CDH_6mg_thr10_exp2(cond_border); %pivobj_CDH_6_mg_blured{i}.loaded_vort
    
    figure
    
    
    fh=vort_core_CDH_6_mg_blured_stripes_thr10_exp2{i};
    h=histogram(fh(~isnan(fh)),...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_core_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vort_core_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))

    
    hold on
    
    fh=vort_bulk_CDH_6_mg_blured_stripes_thr10_exp2{i};
    h=histogram(fh(~isnan(fh)),...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_bulk_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vort_bulk_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))

    
    hold on
    
    fh=vort_border_CDH_6_mg_blured_stripes_thr10_exp2{i};
    h=histogram(fh(~isnan(fh)),...
        linspace(-1.5,1.5,100),'Normalization','pdf');
    write_2_column_table(['pdf_vort_border_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        h.BinEdges(1:end-1),h.Values)
    write_1_column_table(['all_vort_border_stripes_thr10_exp2_bluredCDH_6_mg_' num2str(i) '.dat'],...
        fh(~isnan(fh)))

    R_positions=sqrt((exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.xs-x_c).^2 + ...
        (exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.ys-y_c).^2);
    absv_all=sqrt(exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vxs.^2+exp2_pivobj_CDH_6_mg_blured_stripes_thr10{i}.vys.^2);
    abscurlz_all=pivCURLz_all_CDH_6mg_thr10_exp2;
    absdiv_all=pivDIVall_CDH_6mg_thr10_exp2;
    
    for j=1:num_prof_points
        
        cond_RdR=(R_positions<(j*dR)) & (R_positions>=((j-1)*dR));
        
        av_vel_profile_bluredCDH_6_stripes_thr10_exp2{i}(j)=nanmean(absv_all(cond_RdR));
        
        av_curlz_profile_bluredCDH_6_stripes_thr10_exp2{i}(j)=nanmean(abs(abscurlz_all(cond_RdR)));
        av_div_profile_bluredCDH_6_stripes_thr10_exp2{i}(j)=nanmean((absdiv_all(cond_RdR)));
    end
    
    
    figure;
    plot((1:num_prof_points)*dR,av_vel_profile_bluredCDH_6_stripes_thr10_exp2{i})
    hold on
    
    write_2_column_table(['velocity_profile_stripes_thr10_exp2_bluredCDH_6_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_vel_profile_bluredCDH_6_stripes_thr10_exp2{i})
    
    
    
    
    xlabel('r')
    ylabel('mean v')
    legend('Location','best')
    
    title('CDH 2 mg')
    
    set(gca,'FontSize', 14)
    
    figure
    plot((1:num_prof_points)*dR,av_curlz_profile_bluredCDH_6_stripes_thr10_exp2{i},...
        'DisplayName',['set ' num2str(i)])
    hold on
    
    write_2_column_table(['blured_absvort_profile_stripes_thr10_exp2_bluredCDH_6_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_curlz_profile_bluredCDH_6_stripes_thr10_exp2{i})
    
    figure
    
    plot((1:num_prof_points)*dR,av_div_profile_bluredCDH_6_stripes_thr10_exp2{i},...
        'DisplayName',['set ' num2str(i)])
    hold on
    
    write_2_column_table(['blured_div_profile_stripes_thr10_exp2_bluredCDH_6_mg_set_'...
        num2str(i) '.dat' ], (1:num_prof_points)*dR ,...
        av_div_profile_bluredCDH_6_stripes_thr10_exp2{i})
    
    
    
    
end
