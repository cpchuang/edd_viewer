function edd_cali_energy(varargin)
% detector channel to energy calibration
%
%   + 1.1.2 2020/07/07
%           - change inital guess by changing detpar
%   + 1.1.1  2019/10/23
%           - select steps to fit
%   + 1.1.0  2018/08/06
%           - fit dual source (Co-57/Cd-109)
%           - fit ch/E with 2-deg polynomial
%   + 1.0.0  2017/11/03
%           - initial release
%
%  To-DO: (1) group peaks to fit
%
% Copyright 2017-2019 Andrew Chuang
% $Revision: 1.1.2 $  $Date: 2021/07/07 $

% Check if the panel exist or not
heddparMain = findall(0,'Tag','eddcaliEnergymain_Fig');

if isempty(heddparMain)
    if nargin == 1
        opt = varargin{1};
    else
        opt.source = 1;
        opt.detno  = 1;
        opt.posno  = 1;
        opt.detpar = [0 0.026 0];
    end
    opt.detpar = [0 0.0264 0];   
    % default options
    opt.table_selection = [];
    if ~isfield(opt,'detno')
        opt.detno = 1;             % detector to be calibrated
    end
    % co-57 (in keV)  http://www.spectrumtechniques.com/products/sources/cobalt-57/
    opt.emission(1).name    = 'Co-57'; 
    opt.emission(1).energy  = [136.47356 122.06065 14.41295 7.058 (6.40384*33.2+6.39084*16.8)/50]; 
    %opt.emission(1).channel = [3925 3508 416 203 184];
    opt.emission(1).channel = (opt.emission(1).energy-opt.detpar(1))/opt.detpar(2);
    opt.emission(1).width   = [  33   30  10  10  10];
    opt.emission(1).int     = [ 0.1    1 0.3 0.02 0.08];
    opt.emission(1).default = [true true true false false];

    % cd-109 (in keV)  http://www.spectrumtechniques.com/products/sources/cadmium-109/
    opt.emission(2).name    = 'Cd-109'; 
    %opt.emission(2).energy  = [88.033610 (25.512 25.4567) (25.146 24.9427 24.9118) (22.16317 21.9906)]; 
    opt.emission(2).energy  = [88.033610 25.4664 24.9332 22.1031];    % [gamma kb2 kb1 ka]
    %opt.emission(2).channel = [2530 732 718 636];
    opt.emission(2).channel = (opt.emission(2).energy-opt.detpar(1))/opt.detpar(2);
    opt.emission(2).width   = [18 10 10 10];
    opt.emission(2).int     = [0.035 0.04 0.2 1];
    opt.emission(2).default = [true true true true];
    
    opt.emission(3).name    = 'dual_57/109'; 
    opt.emission(3).energy  = [136.47356 122.06065 88.033610 25.4664 24.9332 22.1031 14.41295 7.058 6.3995];    % [gamma_Co gamma_Co gamma_Cd kb2 kb1 ka]
    %opt.emission(3).channel = [3925 3508 2530 732 718 636 416 203 184];
    opt.emission(3).channel = (opt.emission(3).energy-opt.detpar(1))/opt.detpar(2);
    opt.emission(3).width   = [27 24 17 8 10 10 9 7 8];
    opt.emission(3).int     = [0.01 0.1 0.025 0.025 0.15 0.7 0.03 0.005 0.008];
    opt.emission(3).default = [true true true true true true false true true];

    
    initFigure(opt);

else
%     if nargin == 1
%         hh = get_handle;
%         hh.opt.data = varargin{1}.data;
%         %hh.opt.detno = varargin{1}.detno;
%         update_handle(hh)
%     end
    figure(heddparMain);
    
end

%================================================================
% --- initialize main figure layout
%================================================================
function initFigure(opt)
% --------------------------------
% --- main figure handle
% --------------------------------
pos = get(0,'PointerLocation');
figureSize = [570 280];   % [550 245]
figurePos  = [pos(1) pos(2)-180 figureSize];
tablesize = 145;
row1 = tablesize+15;
row2 = tablesize+60;  % currently not in-use.
row3 = tablesize+50;
heddsetparmain = figure(...
    'BackingStore','on',...
    'Units','pixels',...
    'DockControls','off',...
    'Resize','off',...
    'PaperOrient','portrait',...
    'PaperPositionMode','auto',...
    'IntegerHandle','off',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'Toolbar','none',...
    'Name','Detector Channel_to_Energy Calibration',...
    'Position',figurePos,...
    'HandleVisibility','callback',...
    'Tag','eddcaliEnergymain_Fig',...
    'CloseRequestFcn',@eddparmain_CloseRequestFcn,...
    'UserData',[]);
hAxes = axes(...
    'Parent',heddsetparmain,...
    'Units','pixels',...
    'Position',[0 0 figureSize],...
    'Xlim',[0 figureSize(1)],...
    'Ylim',[0 figureSize(2)],...
    'Tag','eddgetparmain_Axes');
hPatchMain = patch(...
    'Parent',hAxes,...
    'XData',[0 figureSize(1) figureSize(1) 0],...
    'YData',[0 0 figureSize(2) figureSize(2)],...
    'FaceColor',[1 1 0.55],...
    'EdgeColor',[1 1 0.55]);

% ----------------------------------------
% --- layout of control panel putshbuttons
% ----------------------------------------
% text('Parent',hAxes,...
%     'Position',[475 148],...
%     'String','E = a * Ch + b',...
%     'Fontsize',13,...
%     'Units','pixel',...
%     'HorizontalAlignment','Center',...
%     'Color',[0.4 0.3 0]);

hpop(1) = uicontrol('Style','PopupMenu',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[10 row1 80 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String',{'Det-1','Det-2'},...
    'Value',1,...
    'HorizontalAlignment','left',...
    'callback',@(a,b)fprintf('select Det-%d to calibrate\n',a.Value),...
    'Tag','main_pop_detno');

hpop(2) = uicontrol('Style','PopupMenu',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[100 row1 80 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String',{'Co-57','Cd-109','dual(Cd/Co)'},...
    'Value',1,...
    'HorizontalAlignment','left',...
    'callback',@main_select_emission_Callback,...
    'Tag','main_pop_emission');

hpop(3) = uicontrol('Style','PopupMenu',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[340 row1 50 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String',{'1','2'},...
    'TooltipString','fit with n-deg polynomial',...
    'Value',1,...
    'HorizontalAlignment','center',...
    'callback',@main_select_npoly_Callback,...
    'Tag','main_pop_npoly');

hpush(1) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[185 row1 40 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','Add',...
    'HorizontalAlignment','center',...
    'callback',@main_add_par_Callback,...
    'Enable','off',...
    'Tag','main_push_add');

hpush(2) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[230 row1 40 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','Del',...
    'HorizontalAlignment','center',...
    'callback',@main_del_par_Callback,...
    'Enable','off',...
    'Tag','main_push_del');

hpush(3) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[275 row1 60 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','fit_pk',...
    'HorizontalAlignment','center',...
    'callback',@main_fit_peak_Callback,...
    'Enable','on',...
    'Tag','main_push_fit_peak');

hpush(4) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[395 row1 80 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','fit_det',...
    'HorizontalAlignment','center',...
    'callback',@main_fit_det_Callback,...
    'Enable','on',...
    'Tag','main_push_fit_det');

hpush(5) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[480 row1 80 30],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','Accept',...
    'HorizontalAlignment','center',...
    'callback',@main_accept_Callback,...
    'Enable','on',...
    'Tag','main_push_accept');

% htext(1) = uicontrol('Style','text',...
%     'Parent',heddsetparmain,...
%     'units','pixels',...
%     'Position',[10 row2-5 50 30],...
%     'backgroundcolor',[1 1 0.55],...
%     'fontsize',12,...
%     'String','a =',...
%     'HorizontalAlignment','left',...
%     'Enable','on',...
%     'Tag','main_text_par1');
% 
% htext(2) = uicontrol('Style','text',...
%     'Parent',heddsetparmain,...
%     'units','pixels',...
%     'Position',[235 row2-5 50 30],...
%     'backgroundcolor',[1 1 0.55],...
%     'fontsize',12,...
%     'String','b =',...
%     'HorizontalAlignment','left',...
%     'Enable','on',...
%     'Tag','main_text_par2');
% 
% hedit(1) = uicontrol('Style','edit',...
%     'Parent',heddsetparmain,...
%     'units','pixels',...
%     'Position',[50 row2 175 30],...
%     'backgroundcolor',[1 1 1],...
%     'fontsize',12,...
%     'String',sprintf('%4.8f',0.035),...
%     'TooltipString','Energy (keV) = a*channel + b',...
%     'HorizontalAlignment','center',...
%     'callback',@(a,b)disp('Not in use'),...
%     'Enable','on',...
%     'Tag','main_edit_par1');
% 
% hedit(2) = uicontrol('Style','edit',...
%     'Parent',heddsetparmain,...
%     'units','pixels',...
%     'Position',[275 row2 175 30],...
%     'backgroundcolor',[1 1 1],...
%     'fontsize',12,...
%     'String',sprintf('%4.8f',0),...
%     'TooltipString','Energy (keV) = a*channel + b',...
%     'HorizontalAlignment','center',...
%     'callback',@(a,b)disp('Not in use'),...
%     'Enable','on',...
%     'Tag','main_edit_par2');

% construct table (peak parameters)
data = num2cell([opt.emission(1).energy' opt.emission(1).channel' opt.emission(1).width' opt.emission(1).int']);
data = cat(2,data,{true;true;true;false;false});

htablepar = uitable('Parent',heddsetparmain,...
    'Units','pixels',...
    'Position',[10 10 figureSize(1)-20 tablesize],...
    'ColumnName',{'Emission (keV)',...
                  'Channel',...
                  'Width',...
                  'Intensity',...
                  'Use',...
                  'Error (eV)'},...
    'ColumnFormat',{'numeric','numeric','numeric','numeric','logical','numeric'},...
    'ColumnWidth',{100,100,100},...
    'ColumnEditable',[true,true,true,true,true,false],...
    'Data',data,...
    'Fontsize',11,...
    'RowName',[],...
    'CellSelectionCallback',@uitable_selection_Callback,...
    'Tag','par_table');

% construct table (energy-channel conversion)
data2 = num2cell([0 0.3 0]);

htablepar(2) = uitable('Parent',heddsetparmain,...
    'Units','pixels',...
    'Position',[10 row3 figureSize(1)-20 70],...
    'ColumnName',{'p0','p1','p2'},...
    'ColumnFormat',{'long','long','long'},...
    'ColumnWidth',{round((figureSize(1)-20)/3)-2,round((figureSize(1)-20)/3),round((figureSize(1)-20)/3)},...
    'ColumnEditable',[true,true,true],...
    'Data',data2,...
    'Fontsize',13,...
    'RowName',[],...
    'CellSelectionCallback',@uitable2_selection_Callback,...
    'CellEditCallback',@uitable2_edit_Callback,...
    'Tag','en2ch_table');


allhandle = guihandles(heddsetparmain);
allhandle.opt = opt;
update_handle(allhandle)

%==========================================================================
% --- Callbacks for  
%==========================================================================
function main_select_emission_Callback(~,eventdata,~)
h = get_handle;
mtab = h.par_table;
selection = eventdata.Source.Value;
data = num2cell([h.opt.emission(selection).energy' h.opt.emission(selection).channel' h.opt.emission(selection).width' h.opt.emission(selection).int']);
use  = num2cell(h.opt.emission(selection).default');
data = cat(2,data,use);
set(mtab,'Data',data);
h.opt.source = selection;
update_handle(h)

function main_add_par_Callback(~,~,~)
h = get_handle;
mtab = h.par_table;
data = mtab.Data;
empty={0,0,1,'false'};
if isempty(h.opt.table_selection)
    data=cat(1,data,empty);
else
    select=h.opt.table_selection(1);
    data=cat(1,data(1:select,:),empty,data(select+1:end,:));
end
set(mtab,'Data',data);

function main_del_par_Callback(~,~,~)
h = get_handle;
mtab = h.par_table;
data = mtab.Data;

if isempty(h.opt.table_selection)
    if size(data,1)==1
        data = {0,0,1,'false'};
    else
        data = data(1:end-1,:);
    end
else
    select=h.opt.table_selection(1);
    if size(data,1)==1
        data = {0,0,1,'false'};
    else
        data=cat(1,data(1:select-1,:),data(select+1:end,:));
    end
end
set(mtab,'Data',data);

function main_fit_peak_Callback(~,~,~)
h = get_handle;
posno = sort(str2num(h.opt.xpos_handle.String));
mtab = h.par_table;
source = h.opt.source;
detno = h.main_pop_detno.Value;
h.opt.detno = detno;

if ~isfield(h.opt,'data')
    warndlg('read scan data first!!');
    return;
end
switch source
    case 1
        y_ = mean(h.opt.data{detno}(posno,:),1);
        
        %%% initial guess of peak position
        p0.pkid = {1 2 3 [4 5]};
        p0.cen  = {mtab.Data{1,2}, mtab.Data{2,2}, mtab.Data{3,2}, [mtab.Data{4,2} mtab.Data{5,2}]};
        p0.int  = {max(y_)*mtab.Data{1,4}, max(y_)*mtab.Data{2,4}, max(y_)*mtab.Data{3,4}, max(y_)*[mtab.Data{4,4} mtab.Data{5,4}]};
        p0.fwhm =  {mtab.Data{1,3}, mtab.Data{2,3}, mtab.Data{3,3}, [mtab.Data{4,3} mtab.Data{5,3}]};
        
        opt.vis_fit = 1;
        opt.dup   = [80 80 80 60];
        opt.ddown = [80 80 80 40];
        
        pfit2 = edd_fit_emission(p0,1:length(y_),y_,opt);
        if opt.vis_fit == 1
            close(figure(2))
        end
        for i = 1:length(pfit2.int)
            mtab.Data{i,2} = pfit2.cen(i);
            mtab.Data{i,3} = pfit2.fwhm(i);
            mtab.Data{i,4} = pfit2.int(i)./max(pfit2.int);
        end
        update_handle(h);
    case 2
        y_ = mean(h.opt.data{detno}(posno,:),1);
        
        %%% initial guess of peak position
        p0.pkid = {1 [2 3] 4};
        p0.cen  = {mtab.Data{1,2}, [mtab.Data{2,2}, mtab.Data{3,2}], mtab.Data{4,2}};
        p0.int  = {max(y_)*mtab.Data{1,4}, [max(y_)*mtab.Data{2,4} max(y_)*mtab.Data{3,4}], max(y_)*mtab.Data{4,4}};
        p0.fwhm =  {mtab.Data{1,3}, [mtab.Data{2,3} mtab.Data{3,3}], mtab.Data{4,3}};
        
        opt.vis_fit = 1;
        opt.dup   = [60 40 40 60];
        opt.ddown = [60 40 40 60];
        
        pfit2 = edd_fit_emission(p0,1:length(y_),y_,opt);
        if opt.vis_fit == 1
            close(figure(2))
        end
        for i = 1:length(pfit2.int)
            mtab.Data{i,2} = pfit2.cen(i);
            mtab.Data{i,3} = pfit2.fwhm(i);
            mtab.Data{i,4} = pfit2.int(i)./max(y_);
        end
        update_handle(h);
    case 3
        size(h.opt.data{detno}(posno,:))
        y_ = mean(h.opt.data{detno}(posno,:),1);
        
        %%% initial guess of peak position
        p0.pkid = {1 2 3 [4 5] 6 7 [8 9]};
        p0.cen  = {mtab.Data{1,2}, mtab.Data{2,2}, mtab.Data{3,2}, [mtab.Data{4,2} mtab.Data{5,2}], mtab.Data{6,2}, mtab.Data{7,2}, [mtab.Data{8,2} mtab.Data{9,2}]};
        p0.int  = {max(y_)*mtab.Data{1,4}, max(y_)*mtab.Data{2,4}, max(y_)*mtab.Data{3,4}, [max(y_)*mtab.Data{4,4} max(y_)*mtab.Data{5,4}], max(y_)*mtab.Data{6,4}, max(y_)*mtab.Data{7,4}, [max(y_)*mtab.Data{8,4} max(y_)*mtab.Data{9,4}]};
        p0.fwhm =  {mtab.Data{1,3}, mtab.Data{2,3}, mtab.Data{3,3}, [mtab.Data{4,3} mtab.Data{5,3}], mtab.Data{6,3}, mtab.Data{7,3}, [mtab.Data{8,3} mtab.Data{9,3}]};
        
        opt.vis_fit = 1;
        opt.dup   = [100 100 80 60 60 50 50];
        opt.ddown = [100 100 80 60 80 40 40];
        
        pfit2 = edd_fit_emission(p0,1:length(y_),y_,opt);
        if opt.vis_fit == 1
            close(figure(2))
        end
        for i = 1:length(pfit2.int)
            mtab.Data{i,2} = pfit2.cen(i);
            mtab.Data{i,3} = pfit2.fwhm(i);
            mtab.Data{i,4} = pfit2.int(i)./max(y_);
        end
        update_handle(h);        
    otherwise
        warndlg('Only support Co-57/Cd-109 as calibration source.');        
end

function main_fit_det_Callback(~,~,~)
h = get_handle;
mtab  = h.par_table;
mtab2 = h.en2ch_table;
ideal_emission_energy = cell2mat(mtab.Data(:,1))';
peak_pos_in_channel   = cell2mat(mtab.Data(:,2))';
peak_to_use           = find(cell2mat(mtab.Data(:,5))'*1);

xdata = peak_pos_in_channel(peak_to_use);
ydata = ideal_emission_energy(peak_to_use);

npoly = h.main_pop_npoly.Value;

switch npoly
    case 1
        pfit = polyfit(xdata, ydata, 1);
        yfit = polyval(pfit,peak_pos_in_channel);
        yfit_use = polyval(pfit,xdata);
        mtab2.Data{1} = pfit(2);
        mtab2.Data{2} = pfit(1);
        mtab2.Data{3} = 0;
        title_txt = sprintf('Energy  = %f(ch) + %f',pfit(1),pfit(2));
    case 2
        % method-1
        pfit = polyfit(xdata, ydata, 2);
        yfit = polyval(pfit,peak_pos_in_channel);
        yfit_use = polyval(pfit,xdata);
        mtab2.Data{1} = pfit(3);
        mtab2.Data{2} = pfit(2);
        mtab2.Data{3} = pfit(1);
        title_txt = sprintf('Energy  = %2.4e(ch)^{2} + %f(ch) + %f',pfit(1),pfit(2),pfit(3));

%         % method-2 (generic way of doing fitting, keep for future reference)
%         fn = {'polyval'}; npar = [3];
%         pfix = [nan nan nan]; % fix fitting parameters (currently none are fixed)
%         opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
%         p02 =  [   0 0.04   0];
%         pUB2 = [ inf +0.2  +1];
%         pLB2 = [-inf    0  -1];
%         y0   = sumfun1(p02,pfix,npar,fn,xdata);
%         
%         [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
%         yfit = sumfun1(pfit,pfix,npar,fn,peak_pos_in_channel);
%         yfit_use = sumfun1(pfit,pfix,npar,fn,xdata);
end

error = (ideal_emission_energy-yfit)*1E3;
for i = 1:size(mtab.Data,1)
    mtab.Data{i,6} = round(error(i),1);
end

figure(169); clf(gcf,'reset'); set(gcf,'Tag','eddgetpar_fit_fig');
plot(peak_pos_in_channel,error,'-o','displayname','Co57')
line(xdata,(ydata-yfit_use')*1E3,'marker','^','linestyle','none','color','r','displayname','Co57 (used to fit)')
%title(sprintf('Ch to Energy  = [%f %f]',pfit(1),pfit(2)));
title(title_txt);
xlabel('Channel')
ylabel('Error in Energy (eV)')
ylim([-60 60])

function main_accept_Callback(~,~,~)
h = get_handle;

hroot = findall(0,'Tag','eddviewmain_Fig');
hroothandledata = getappdata(hroot,'allhandle');

switch h.opt.detno
    case 1
        %hroothandledata.expinfo_edit_ch2e1a.String=num2str(h.en2ch_table.Data{2},'%.10g');
        %hroothandledata.expinfo_edit_ch2e1b.String=num2str(h.en2ch_table.Data{3},'%.10g');
        %hroothandledata.config.visopt.Inst(1).Ch2E(1) = str2double(h.main_edit_par1.String);
        %hroothandledata.config.visopt.Inst(1).Ch2E(2) = str2double(h.main_edit_par2.String);
        %hroothandledata.config.visopt.Inst(1).Ch2E(1) = h.en2ch_table.Data{2};
        %hroothandledata.config.visopt.Inst(1).Ch2E(2) = h.en2ch_table.Data{3};
        
        hroothandledata.config.visopt.Inst.detpar(1,2:end) = cell2mat(h.en2ch_table.Data);
        hroothandledata.expinfo_table_detpar.Data(1,2:end) = h.en2ch_table.Data;
        
    case 2
        %hroothandledata.expinfo_edit_ch2e2a.String=num2str(h.en2ch_table.Data{2},'%.10g');
        %hroothandledata.expinfo_edit_ch2e2b.String=num2str(h.en2ch_table.Data{3},'%.10g');
        %hroothandledata.config.visopt.Inst(2).Ch2E(1) = str2double(h.main_edit_par1.String);
        %hroothandledata.config.visopt.Inst(2).Ch2E(2) = str2double(h.main_edit_par2.String);
        %hroothandledata.config.visopt.Inst(2).Ch2E(1) = h.en2ch_table.Data{2};
        %hroothandledata.config.visopt.Inst(2).Ch2E(2) = h.en2ch_table.Data{3};
        
        hroothandledata.config.visopt.Inst.detpar(2,2:end) = cell2mat(h.en2ch_table.Data);
        hroothandledata.expinfo_table_detpar.Data(2,2:end) = h.en2ch_table.Data;

    otherwise
        warndlg('Not supportted yet!!')
end
setappdata(hroot,'allhandle',hroothandledata)

function main_select_npoly_Callback(hObject,~)

switch hObject.Value
    case 1
        fprintf('\nFit with linear\n')
        fprintf('   energy (keV) = a * ch + b\n')
    case 2
        fprintf('\nFit with Trinomial\n')
        fprintf('   energy (keV) = a * ch^2 + b * ch + c\n')
end

%%%% call back for edit
function uitable2_edit_Callback(hhobj,eventdata)
    h = get_handle;
    
    % get updated detpar
    detpar = [hhobj.Data{:}];
    % get theoratical emission energy
    emission_eng = [h.par_table.Data{:,1}]';
    % change expected channel number based on detpar and emission energy
    h.par_table.Data(:,2) = num2cell((emission_eng-detpar(1))./detpar(2));
    
    update_handle(h);


%%%% call back to report selection
function uitable_selection_Callback(~,eventdata)
    h = get_handle; 
    h.opt.table_selection = eventdata.Indices;
    update_handle(h);

function uitable2_selection_Callback(~,eventdata)
    h = get_handle;
    fprintf('select %d\n',eventdata.Indices);
    h.opt.table_selection = eventdata.Indices;
    eventdata.Source
    update_handle(h);
    

%==========================================================================
% --- close request function of eddgetparmain fig 
%==========================================================================
function eddparmain_CloseRequestFcn(hObject,eventdata)                         %#ok<INUSD>
    delete(findall(0,'Tag','eddcaliEnergymain_Fig'))
    delete(findall(0,'-regexp','Tag','eddgetpar*'))  

%==========================================================================
% --- Load all handle object function
%==========================================================================
function h = get_handle
    hmain = findall(0,'Tag','eddcaliEnergymain_Fig');
    h = getappdata(hmain,'allhandle');

%==========================================================================
% --- Save all handle object function
%==========================================================================
function update_handle(handle_to_save)
    hmain = findall(0,'Tag','eddcaliEnergymain_Fig');
    setappdata(hmain,'allhandle',handle_to_save);