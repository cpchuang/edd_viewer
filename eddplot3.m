% eddplot to visdualize 6BM data.
%
%   +1.5 2020/01/16 support different size (for 10-element detector)
%   +1.4 2018/11/27 new detpart format
%   +1.3 2018/04/17 add x scaling factor
%   +1.2 2018/02/03 switch yscale for 1-D plot
%   +1.1 2017/09/20 bug fix
%
% Copyright 2017-2020 Andrew Chuang (chuang.cp@gmail.com)
% $Revision: 1.5 $  $Date: 2020/01/16 $

function eddplot3(da,opt)

fsa = 13;
fst = 18;
cc = lines(20);

if nargin == 1
    opt = '';
end

%%%%%%  define default option %%%%%

if isfield(opt,'datalim')
    datalim = opt.datalim;
else
    % range for (Ch/E/d),(posno),(Intensity)
    datalim = {'auto','auto','auto'}; 
end

if isfield(opt,'normalize')
    normalize = opt.normalize;
else
    normalize = 0;
end

if isfield(opt,'detno')
    detno = opt.detno;
else
    detno = 1;
end

if isfield(opt,'phase')
    phase = opt.phase;
else
    phase = 1;
end

if isfield(opt,'pk')
    pk = opt.pk;
else
    pk = 5;
end

if isfield(opt,'do_export')
    do_export = opt.do_export;
else
    do_export = 0;
end

if isfield(opt,'x_unit')
    x_unit = lower(opt.x_unit);
else
    x_unit = 'ch';
end

if isfield(opt,'yscale')
    switch opt.yscale
        case 1
            y_scale = 'linear';
        otherwise
            y_scale = 'log';
    end
else
    y_scale = 'log';
end

if isfield(opt,'xscaling')
    x_scaling = opt.xscaling;
else
    x_scaling = 1;
end

if isfield(opt,'x_range')
    x_range = opt.x_range;
else
    x_range = ':';
end

% if isfield(opt,'fig_handle')&&ishandle(opt.fig_handle)
%     fig_handle = opt.fig_handle;
% else
%     fig_handle = 0;
% end

if isfield(opt,'avg_1D')
    avg_1D = opt.avg_1D;
else
    avg_1D = 0;
end

if isfield(opt,'title')
    title_text = opt.title;
else
    title_text = '';
end

if isfield(opt,'dzero')
    dzero = opt.dzero;
else
    dzero = zeros(1,length(pk))+1;
end

if isfield(opt,'posno')
    posno = opt.posno;
else
    posno = 1;
end

if isfield(opt,'xoffset')
    xoffset = opt.xoffset;
else
    xoffset = 0;
end

if isfield(opt,'Eoffset')
    Eoffset = opt.Eoffset;
else
    Eoffset = [0 0];
end

if isfield(opt,'type')
    type = opt.type;
else
    type = '2draw';
end

if isfield(opt,'scno')
    scno = opt.scno(1);   % only plot one scan at a time. (May.17)
else
    scno = 1;
end

if isfield(da,'Material')
    label_peak = 1;
else
    label_peak = 0;
end

if ~strcmp(type,'raw')&isfield(da(1),'Material')
    pkname = da(1).Material(phase).hkls;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,'2draw')&&size(da(scno).data{detno},1)==1
    type = 'raw';
end

if isfield(opt,'Inst')
    da(1).Inst = opt.Inst;
elseif ~isfield(da(scno),'Inst')
    fprintf('No Instrument parameters provided!!\n');
    fprintf('Use default TOA  = 3 (deg)\n');
    fprintf('Use default Ch2E = [0 0.03 0];\n');
    da(scno).Inst.detpar = [3 0 0.03 0;...
                            3 0 0.03 0];
end

% calculate E_grid, d_grid
hc = 12.398419057638671;
for i = 1:2
    da(scno).Inst.E_grid(i,:) = polyval(flip(da(1).Inst.detpar(i,2:end)), 1:size(da(scno).data{detno},2));
    da(scno).Inst.d_grid(i,:) = hc./da(scno).Inst.E_grid(i,:)*0.5/sind(da(1).Inst.detpar(i,1)/2);
end

%%%% initiate figure for step plot
if strcmpi(type,'step_plot')
    hfig = findall(0,'Tag','edd_fig_step_plot');
    if ishandle(hfig)
        fig = hfig(1);
        clf(fig,'reset')
        set(fig,'Tag','edd_fig_step_plot');
    else
        fig = figure(163);
        set(fig,'Position',[900 100 800 600],'Tag','edd_fig_step_plot');
    end
    figure(163);
else
    
    
    %%%%%% switch between detectors
    switch detno
        case 1
            if strcmp(type,'2draw')
                hfig = findall(0,'Tag','edd_fig_det1_map');
                if ishandle(hfig)
                    fig = hfig(1);
                    clf(fig,'reset');
                    set(fig,'Tag','edd_fig_det1_map');
                else
                    fig = figure(161);
                    set(fig,'Position',[50 100 800 600],'Tag','edd_fig_det1_map');
                end
                figure(161);
            else
                hfig = findall(0,'Tag','edd_fig_det1_line');
                if ishandle(hfig)
                    fig = hfig(1);
                    clf(fig,'reset');
                    set(fig,'Tag','edd_fig_det1_line');
                else
                    fig = figure(162);
                    set(fig,'Position',[50 700 800 600],'Tag','edd_fig_det1_line');
                end
                figure(162);
            end
        case 2
            if strcmp(type,'2draw')
                hfig = findall(0,'Tag','edd_fig_det2_map');
                if ishandle(hfig)
                    fig = hfig(1);
                    clf(fig,'reset');
                    set(fig,'Tag','edd_fig_det2_map');
                else
                    fig = figure(261);
                    set(fig,'Position',[100 100 800 600],'Tag','edd_fig_det2_map');
                end
                figure(261);
            else
                hfig = findall(0,'Tag','edd_fig_det2_line');
                if ishandle(hfig)
                    fig = hfig(1);
                    clf(fig,'reset');
                    set(fig,'Tag','edd_fig_det2_line');
                else
                    fig = figure(262);
                    set(fig,'Position',[100 700 800 600],'Tag','edd_fig_det2_line');
                end
                figure(262);
            end
            %         for i = 1:length(da)
            %             da(i).data = da(i).data2;
            %         end
            
            %da(scno).data = cat(3,da(scno).data,da(scno).data2);
        otherwise
            fprintf('3 or more detector config is not supported yet!!\n');
            return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start ploting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = lines(100);
ax = axes('parent',hfig,'fontsize',fsa,'box','on');
grid on;

switch lower(type)
    case 'cen'
        xdata = da(scno).motorpos(x_range)+xoffset;
        fprintf('\nAveraged Peak Position (Mean%sStd) for %s\n',char(177),title_text);
        for i = 1:length(pk)
            ydata = da(scno).fit(phase,pk(i),detno).cen(x_range);
            flag  = find(ydata~=0);
            ydata = ydata(flag);
            if normalize == 1 % normalize to mean
                line(xdata(flag),((mean(ydata)./ydata)-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
                fprintf('[%s] = %8.2f %s %8.2f, ',num2str(pkname(pk(i),:)),mean(((ydata./mean(ydata))-1)*1E6),char(177),std(((ydata./mean(ydata))-1)*1E6));
            elseif normalize == 0
                line(xdata(flag),ydata,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Peak Center (keV)');
                fprintf('[%s] = %8.8f %s %6.4f, ',num2str(pkname(pk(i),:)),mean(ydata),char(177),std(ydata));
            elseif normalize == 3  % show absolute diff
                line(xdata(flag),ydata-mean(ydata),'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Peak Center (keV)');
                fprintf('[%s] = %8.8f %s %6.4f, ',num2str(pkname(pk(i),:)),mean(ydata),char(177),std(ydata));
            else % normalize to user defined d_zero
                line(xdata(flag),((dzero(pk(i))./ydata)-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                fprintf('[%s] = %8.2f %s %8.2f, ',num2str(pkname(pk(i),:)),mean((ydata./dzero(pk(i))-1).*1E6),char(177),std((ydata./dzero(pk(i))-1).*1E6));
                %std((ydata./dzero(pk(i))-1).*1E6)
            end
            %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
            
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1});
        ylim(datalim{3});
        title(title_text,'fontsize',fst);
        leg = legend('toggle');

    case 'dspac'    % d-spacing was calculated during fitting
        xdata = (da(scno).motorpos(x_range)+xoffset)*x_scaling;
        fprintf('\nAveraged Peak Position (Mean%sStd) for %s\n',char(177),title_text);
        for i = 1:length(pk)
            %da(scno).fit(phase,pk(i),detno);
            ydata = da(scno).fit(phase,pk(i),detno).dspac(x_range);
            flag  = find(ydata~=0);
            ydata = ydata(flag);
            if normalize == 1
                line(xdata(flag),((ydata./mean(ydata))-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
                fprintf('[%s] = %8.2f %s %8.2f,\n',num2str(pkname(pk(i),:)),mean(((ydata./mean(ydata))-1)*1E6),char(177),std(((ydata./mean(ydata))-1)*1E6));
            elseif normalize == 0
                line(xdata(flag),ydata,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('d-spacing (Angstron)');
                fprintf('[%s] = %8.8f %s %6.4f,\n',num2str(pkname(pk(i),:)),mean(ydata),char(177),std(ydata));
            else
                line(xdata(flag),((ydata./dzero(pk(i)))-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{6})');
                fprintf('[%s] = %8.2f %s %8.2f,\n',num2str(pkname(pk(i),:)),mean((ydata./dzero(pk(i))-1).*1E6),char(177),std((ydata./dzero(pk(i))-1).*1E6));
                %std((ydata./dzero(pk(i))-1).*1E6)
            end
            %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
            
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1});
        ylim(datalim{3});
        title(title_text,'fontsize',fst);
        leg = legend('toggle');        

    case 'dspac2'     % calculate d-spacing based on tth value
        xdata = (da(scno).motorpos(x_range)+xoffset)*x_scaling;
        fprintf('\nAveraged Peak Position (Mean%sStd) for %s\n',char(177),title_text);
        for i = 1:length(pk)
            dspac_cal = hc./(da(scno).fit(phase,pk(i),detno).cen+Eoffset(detno))/2/sind(da(1).Inst.detpar(detno,1)/2);
            ydata = dspac_cal(x_range);
            flag  = find(ydata>0 & ydata<20);
            ydata = ydata(flag);
            if normalize == 1
                line(xdata(flag),((ydata./mean(ydata))-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
                fprintf('[%s] = %8.2f %s %8.2f, \n',num2str(pkname(pk(i),:)),mean(((ydata./mean(ydata))-1)*1E6),char(177),std(((ydata./mean(ydata))-1)*1E6));
            elseif normalize == 0
                line(xdata(flag),ydata,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('d-spacing (Angstron)');
                ppfit = polyfit(xdata(flag),ydata,1);
                fprintf('[%s] = %8.8f %s %6.4f, d_0 = %1.6f\n',num2str(pkname(pk(i),:)),mean(ydata),char(177),std(ydata),ppfit(2));
                
            else
                line(xdata(flag),((ydata./dzero(pk(i)))-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{6})');
                fprintf('[%s] = %8.2f %s %8.2f,\n',num2str(pkname(pk(i),:)),mean((ydata./dzero(pk(i))-1).*1E6),char(177),std((ydata./dzero(pk(i))-1).*1E6));
                %std((ydata./dzero(pk(i))-1).*1E6)
            end
            %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
            
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1});
        ylim(datalim{3});
        title(title_text,'fontsize',fst);
        leg = legend('toggle');              
        
    case 'int'
        xdata = (da(scno).motorpos+xoffset)*x_scaling;
        for i = 1:length(pk)
            ydata = da(scno).fit(phase,pk(i)).int;
            yref  = da(scno).fit(phase,1).int;    % use fitst peak as refernce
            flag  = find(ydata~=0);
            ydata = ydata(flag);
            if normalize == 0
                line(xdata(flag),ydata.*da(1).exp_time,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Intensity (counts)');
            elseif normalize == 1
                line(xdata(flag),ydata/max(ydata),'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Normalized Intensity (cps)');
            else
                line(xdata(flag),ydata./yref(flag),'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Normalized Intensity (cps)');
            end
            text(0.4,0.95,sprintf('Total count time: %d sec',da(1).exp_time),'sc','fontsize',fsa);
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1})
        ylim(datalim{3})
        title(title_text,'fontsize',fst)
        leg = legend('toggle');        
        
    case 'fwhm'
        xdata = da(scno).motorpos+xoffset;
        for i = 1:length(pk)
            ydata = da(scno).fit(phase,pk(i)).fwhm;
            flag  = find(ydata~=0);
            line(xdata(flag),ydata(flag),'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
            ylabel('FWHM (keV)');
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1})
        ylim(datalim{3})
        title(title_text,'fontsize',fst)
        leg = legend('toggle');        
        
    case 'raw'
        yrange = datalim{3};
        switch x_unit
            case 'd'
                %xxdata = da(point).Inst.d_grid;
                %xxdata = da(scno).Inst(detno).d_grid;
                xxdata = da(scno).Inst.d_grid(detno,:);
                xrange = [0 10];
                xlab = sprintf('d (%s)',char(197));
            case 'e'
                %xxdata = da(point).Inst.E_grid;
                %xxdata = da(scno).Inst(detno).E_grid;
                xxdata = da(scno).Inst.E_grid(detno,:);
                xrange = [0 200];
                xlab = sprintf('E (keV)');
            case 'ch'
                xxdata = 1:size(da(scno).data{detno},2);
                xrange = [1 size(da(scno).data{detno},2)];
                xlab = sprintf('Channel');
        end
        if avg_1D
            if normalize == 0
                yydata = sum(da(scno).data{detno}(posno,:),1)./length(posno);
                tx_yax = 'Intensity (counts)';
            else
                yydata = sum(da(scno).data{detno}(posno,:),1)./length(posno)./da(scno).exp_time;
                tx_yax = 'Normalized Intensity (cps)';
            end
            line(xxdata, yydata,'color',cc(1,:),'marker','.','displayname','raw curve')
        else
            if normalize == 0
                yydata = double(da(scno).data{detno}(posno,:));
                tx_yax = 'Intensity (counts)';
            else
                yydata = double(da(scno).data{detno}(posno,:))./da(scno).exp_time;
                %yydata = sum(da(scno).data(posno,:),1)./length(posno);
                tx_yax = 'Normalized Intensity (cps)';
            end
            for k = 1:length(posno)
                line(xxdata, yydata(k,:),'color',cc(k,:),'marker','.','displayname',sprintf('pos-%d',posno(k)))
            end
        end

%         if strcmp(x_unit,'d')
%             xxx = 0.5:0.001:10;
%             yyy = interp1(xxdata,yydata,xxx,'spline');
%             line(xxx,yyy,'color','r','marker','o','markersize',3);
%         end
        
%         if label_peak
%             ccode = {'b' '[0 0.8 0]' '[0.73 0.7 0.7]'};
%             %pk_list = [da(point).pkid_fit{:}];
%             pk_list = [da(1).pkid_fit{:}];
%             for m = 1:length(da(1).Material)
%                 phase_name_list{m} = da(1).Material(m).name;
%                 switch x_unit
%                     case 'd'
%                         xxsample = da(1).Material(m).d_hkl;
%                     case 'E'
%                         xxsample = da(1).Material(m).E_hkl;
%                     case 'Ch'
%                         xxsample = da(1).Material(m).Ch_hkl;
%                 end
%                 line(xxsample, ones(length(xxsample), 1),'marker','^','color',ccode{m},'markersize',9,'linestyle','none','Displayname',sprintf('%s',da(1).Material(m).name))
%                 lnfit = line(xxsample(pk_list([da(1).phase_fit{:}]==m)), ones(length(pk_list([da(1).phase_fit{:}]==m)), 1),'marker','o','color','r','markersize',5,'linestyle','none','displayname','peak fitted');
%                 if m ~= 1; set(get(get(lnfit,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');end
%             end
%         end
        if ~strcmp(datalim{1},'auto'); xrange = datalim{1};end
        grid on
        set(gca,'yscale',y_scale,'fontsize',fsa,'box','on');
        xlabel(xlab)
        ylabel(tx_yax)
        xlim(xrange);
        ylim(yrange);
        legend('toggle');
        if length(posno)>1
            title_text = sprintf('Det-%d, Positon = %d - %d',detno,posno(1),posno(end));
        else
            title_text = sprintf('Det-%d, Positon = %d',detno,posno(1));
        end
        title(title_text,'fontsize',fst)
        text(0.1,0.05,sprintf('Total count time: %d sec',da(1).exp_time*length(posno)),'sc','fontsize',fsa);
        
    case '2draw'
        switch x_unit
            case 'ch'
                xrange = datalim{1};
                yrange = datalim{2};  %yrange = 'auto';
                int_range = datalim{3};
                if ~strcmp(yrange,'auto'); yrange=yrange+[-0.5 0.5];end
                xlab = sprintf('Channel');
                xgrid =1:size(da(scno).data{detno},2);
                ygrid =1:size(da(scno).data{detno}(:,:),1);
            case 'e'
                xrange = datalim{1};
                yrange = datalim{2};  %yrange = 'auto';
                int_range = datalim{3};
                if ~strcmp(yrange,'auto'); yrange=yrange+[-0.5 0.5];end
                xlab = sprintf('Energy (keV)');
                xgrid= da(scno).Inst.E_grid(detno,:);
                ygrid = 1:size(da(scno).data{detno}(:,:),1);
            case 'd'
                xrange = datalim{1};
                yrange = datalim{2};
                int_range = datalim{3};
                if ~strcmp(yrange,'auto'); yrange=yrange+[-0.5 0.5];end
                xlab = sprintf('d (Angstron)');
                if ~isfield(da(scno),'d_data')
                    for k = 1:size(da(scno).data{detno}(:,:),1)
                        da(scno).Inst.d_grid2(detno,:) = 0.4:0.001:10;
                        xx = da(scno).Inst.d_grid(detno,:);
                        yy = double(da(scno).data{detno}(k,:));
                        qy = interp1(xx,yy,da(scno).Inst.d_grid2(detno,:),'spline');
                        qy(qy<0.1|qy>max(yy)*2)=0;
                        da(scno).d_data(k,:,detno) = qy;
                    end
                end
                xgrid= da(scno).Inst.d_grid2(detno,:);
                ygrid = 1:size(da(scno).d_data(:,:,detno),1);
                %assignin('base','da',da);
            otherwise
                fprintf('The desire x_unit: %s is not supported yet!!\n',x_unit);
                return;
        end
        ylab = sprintf('%s (#) steps',da(scno).motorname);
        linesap = find(da(scno).log==char(10));
        wordsap = find(da(scno).log(1:linesap(1)-1)==' ');
        title_text=sprintf('%s',da(scno).log(wordsap(3)+1:linesap(1)-1));
        
        if strcmp(x_unit,'d')
            imagesc(xgrid,ygrid,log10(da(scno).d_data(:,:,detno)));
            %set(gca,'clim',[]);
        else
            imagesc(xgrid,ygrid,log10(double(da(scno).data{detno}(:,:))));
        end
        
        if ~strcmp(xrange,'auto'); xlim(xrange);end
        if ~strcmp(yrange,'auto'); ylim(yrange);end
        if ~strcmp(int_range,'auto'); set(gca,'clim',[log10(int_range(1)) log10(int_range(2))]);end
        xlabel(xlab,'fontsize',13);
        ylabel(ylab,'fontsize',13);
        title(title_text,'fontsize',17);
    case 'step_plot'
        title_text=sprintf('%s, %2.3f to %2.3f, %2.3f/step',da(scno).motorname, da(scno).motorpos(1), da(scno).motorpos(end), da(scno).motorstep);
        xrange = datalim{1};
        yrange = datalim{2};
        switch x_unit
            case 'ch'
                flag = [1:8192]>=xrange(1) & [1:8192] <= xrange(2);
            case 'e'
                % use only det-1 to calculate sum range
                flag = da(scno).Inst.E_grid(1,:) >=xrange(1) & da(scno).Inst.E_grid(1,:) <= xrange(2);
            case 'd'
                % use only det-1 to calculate sum range
                flag = da(scno).Inst.d_grid(1,:) >=xrange(1) & da(scno).Inst.d_grid(1,:) <= xrange(2);
            otherwise
                fprintf('unknown x_unit, not ready yet!!\n');
                return
            
        end
        if sum(flag)==0
            flag = 1;
        end
        if normalize
            yy1data = sum(da(scno).data{1}(:,flag),2)/max(sum(da(scno).data{1}(:,flag),2));
            yy2data = sum(da(scno).data{2}(:,flag),2)/max(sum(da(scno).data{2}(:,flag),2));
        else
            yy1data = sum(da(scno).data{1}(:,flag),2);
            yy2data = sum(da(scno).data{2}(:,flag),2);
        end
        xxdata = da(scno).motorpos;
        line(xxdata,yy1data,'marker','o','color',cc(1,:),'Displayname','Det-1');       
        line(xxdata,yy2data,'marker','o','color',cc(2,:),'Displayname','Det-2');
        xlabel(sprintf('%s position (mm)',da(scno).motorname),'fontsize',15);
        ylabel(sprintf('Sum Intensity'),'fontsize',15);
        title(title_text,'fontsize',17);
       
       
    otherwise
        fprintf('Selected plot type is not supported');
        
end