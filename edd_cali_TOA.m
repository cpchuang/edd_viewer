function edd_cali_TOA(varargin)
% detector TOA calibration
%
%   + 1.0  2018/08/08
%          - initial release
%
%  To-DO: (1) add other materials
%         (2) add fluorescence line
%
% Copyright 2018 Andrew Chuang
% $Revision: 1.0 $  $Date: 2018/08/08 $

heddparMain = findall(0,'Tag','eddcaliTOAmain_Fig');
if isempty(heddparMain)
    if nargin == 1
        opt = varargin{1};
    else
        opt.mat   = 1;
        opt.detno = 1;
        %opt.Inst(1).par = [5 0 0.3 0];
        %opt.Inst(2).par = [5 0 0.3 0];
        opt.detpar = [5 0 0.3 0;...
                      5 0 0.3 0];
    end
        
    % default options
    opt.table_selection = [];
    if ~isfield(opt,'detno')
        opt.detno = 1;             % detector to be calibrated
    end
    % constants
    % hc = 12.398419057638671;
    
    % Materials
    %%%%% ceria
    opt.material(1).name    = 'CeO2'; 
    opt.material(1).lat     = 5.41165;               % IN Angstrom
    [d_,hkl_] = d0('fcc',opt.material(1).lat,10,0);  % calculate hkl_list and d-spacing
    opt.material(1).hkls    = hkl_;
    opt.material(1).d_hkl   = d_;
    opt.material(1).width   = [ 0.44 0.46 0.56 0.64 0.65 0.77 0.85 0.85 0.92 0.95];
    opt.material(1).int     = [ 0.92 0.33 0.40 0.15 0.02 0.01 0.01 0.01 0.01 0.01];
    opt.material(1).fit_seq = [1 2 3 4 4 5 6 6 7 8];
    opt.material(1).default = logical([1 1 1 1 1 0 0 0 0 0]);

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
figureSize = [550 280];
figurePos  = [pos(1) pos(2)-180 figureSize];
tableh  = 180;
table2h = 42;
row1 = tableh+12 + table2h + 5;
row2 = tableh+12;

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
    'Name','Detector TOA Calibration',...
    'Position',figurePos,...
    'HandleVisibility','callback',...
    'Tag','eddcaliTOAmain_Fig',...
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
%     'Position',[400 148],...
%     'String','E = a * Ch + b',...
%     'Fontsize',13,...
%     'Units','pixel',...
%     'HorizontalAlignment','Center',...
%     'Color',[0.4 0.3 0]);

hpop(1) = uicontrol('Style','PopupMenu',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[10 row1 70 23],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String',{'Det-1','Det-2'},...
    'Value',1,...
    'HorizontalAlignment','left',...
    'callback',@main_select_det_Callback,...
    'Tag','main_pop_detno');

hpop(2) = uicontrol('Style','PopupMenu',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[85 row1 70 23],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String',{'ceria'},...
    'Value',1,...
    'HorizontalAlignment','left',...
    'callback',@main_select_material_Callback,...
    'Tag','main_pop_emission');

hpush(1) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[160 row1 30 23],...
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
    'Position',[195 row1 30 23],...
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
    'Position',[230 row1 40 23],...
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
    'Position',[275 row1 60 23],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','fit_TOA',...
    'HorizontalAlignment','center',...
    'callback',@main_fit_TOA_Callback,...
    'Enable','on',...
    'Tag','main_push_fit_TOA');

hpush(5) = uicontrol('Style','pushbutton',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[340 row1 60 23],...
    'backgroundcolor',[1 1 1],...
    'fontsize',10,...
    'String','Accept',...
    'HorizontalAlignment','center',...
    'callback',@main_accept_Callback,...
    'Enable','on',...
    'Tag','main_push_accept');

hcheck(1) = uicontrol('Style','CheckBox',...
    'Parent',heddsetparmain,...
    'units','pixels',...
    'Position',[410 row1 60 23],...
    'backgroundcolor', [1 1 1],...
    'fontsize', 10,...
    'String',{'guess'},...
    'TooltipString','Guess initial peak position (in keV)',...
    'HorizontalAlignment','center',...
    'Tag','main_check_guess');

% constants
hc = 12.398419057638671;

% construct par table
matno = opt.mat;
detno = opt.detno;

%detpar = opt.Inst(detno).par;
detpar = opt.detpar(detno,:);
d_  = opt.material(matno).d_hkl;
E_  = hc/2./d_/sind(detpar(1)/2);
da_ = zeros(10,6);
da_(:,1) = d_(1:10);
da_(:,2) = E_(1:10);
da_(:,3) = hc/2./d_(1:10)/sind(detpar(1)/2);
da_(:,4) = opt.material(matno).width;
da_(:,5) = opt.material(matno).int;
da_(:,6) = opt.material(matno).fit_seq;
data = cat(2,num2cell(da_),num2cell(opt.material(matno).default'));

htabledet = uitable('Parent',heddsetparmain,...
    'Units','pixels',...
    'Position',[10 row2 figureSize(1)-20 table2h],...
    'ColumnName',{'TOA',...
                  'p0',...
                  'p1',...
                  'p2'},...
    'ColumnFormat',{'long','long','long','long'},...
    'ColumnWidth',{120,120,120,120},...
    'ColumnEditable',[true,true,true,true],...
    'Data',opt.detpar(detno,:),...
    'Fontsize',11,...
    'RowName',[],...
    'CellEditCallback',@uitable_det_edit_Callback,...
    'CellSelectionCallback',@uitable_det_select_Callback,...
    'Tag','det_table');

htablepar = uitable('Parent',heddsetparmain,...
    'Units','pixels',...
    'Position',[10 10 figureSize(1)-20 tableh],...
    'ColumnName',{'d_ideal',...
                  'E_ideal',...
                  'E_fit',...
                  'Width',...
                  'Intensity',...
                  'seq',...
                  'Use',...
                  'Error (eV)'},...
    'ColumnFormat',{'long','long','long','long','numeric',cellstr(strsplit(num2str(0:10))),'logical','numeric'},...
    'ColumnWidth',{80,80,80,80,80,40,30,80},...
    'ColumnEditable',[true,true,true,true,true,true,true,false],...
    'Data',data,...
    'Fontsize',11,...
    'RowName',[],...
    'CellSelectionCallback',@uitable_par_select_Callback,...
    'Tag','par_table');

allhandle = guihandles(heddsetparmain);
allhandle.opt = opt;
update_handle(allhandle)

%==========================================================================
% --- Callbacks for  
%==========================================================================
function main_select_det_Callback(~,eventdata,~)
detno = eventdata.Source.Value;
h = get_handle;
h.opt.detno = detno;
%h.det_table.Data(2:4) = h.opt.Inst(detno).par(2:4);
h.det_table.Data(2:4) = h.opt.detpar(detno,2:4);
update_handle(h)


function main_select_material_Callback(~,eventdata,~)
h = get_handle;
h.opt.mat = eventdata.Source.Value;
%selection = eventdata.Source.Value
%data = num2cell([h.opt.emission(selection).energy' h.opt.emission(selection).channel' h.opt.emission(selection).width' h.opt.emission(selection).int']);
%use  = num2cell(h.opt.emission(selection).default');
%data = cat(2,data,use);
%set(mtab,'Data',data);
%h.opt.source = selection;
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
    % constants
    hc = 12.398419057638671;
    h = get_handle;

    if ~isfield(h.opt,'data')
        warndlg('read scan data first!!');
        return;
    end

    det_par = flip(h.det_table.Data(2:end));
    %det_TOA = h.det_table.Data(1);

    mtab = h.par_table;
    detno = h.opt.detno;
    h_posno = findall(0,'tag','visopt_edit_xpos');
    posno = str2num(h_posno.String);
    exp_time = h.opt.exp_time;

    if exp_time>0
        y_ = mean(h.opt.data{detno}(posno,:),1)/exp_time;
    else
        y_ = mean(h.opt.data{detno}(posno,:),1);
    end

    % create peak fitting sequence
    pk_ind = cell2mat(mtab.Data(:,6)); 
    fit_list = unique(cell2mat(mtab.Data(:,6)));

    opt.vis_fit = 1;
    opt.dup   = [100 100 100 100];
    opt.ddown = [100 100 100 100];

    %assignin('base','xxdata',polyval(det_par,1:length(y_)));

    for i = 1:length(fit_list)
        if fit_list(i) ~= 0   % ignore seq == 0
            fprintf('fit %d-th peak\n',fit_list(i));
            flag  = (pk_ind==fit_list(i));
            %%% initial guess of peak position
            p0.pkid = {find(pk_ind==fit_list(i))};
            p0.cen  = {cell2mat(mtab.Data(flag,3))};
            p0.int = {max(y_)*cell2mat(mtab.Data(flag,5))};
            p0.fwhm = {cell2mat(mtab.Data(flag,4))};
            %%% Do peak fitting (in Energy space)
            x_ = polyval(det_par,1:length(y_)); % x-data in Energy (keV)
            pfit2 = edd_fit_emission(p0,x_,y_,opt);
            %%%% assign parameters
            for j = 1:length(p0.pkid{:})
                pid = p0.pkid{1}(j);
                %mtab.Data{pid,2} = hc/2./pfit2.cen(pid)/sind(det_TOA/2);
                mtab.Data{pid,3} = pfit2.cen(pid);
                mtab.Data{pid,4} = pfit2.fwhm(pid);
                mtab.Data{pid,5} = pfit2.int(pid)./max(y_);
                mtab.Data{pid,8} = (mtab.Data{pid,3} - mtab.Data{pid,2}) * 1E3;
            end
        end
    end

    if opt.vis_fit == 1
        close(figure(2))
    end

    update_handle(h);


function main_fit_TOA_Callback(~,~,~)
    % constants
    hc = 12.398419057638671;
    
    h = get_handle;
    det_TOA = h.det_table.Data(1);
    mtab = h.par_table;
    d_ideal = cell2mat(mtab.Data(:,1))';
    %peak_pos_in_d  = cell2mat(mtab.Data(:,2))';
    peak_pos_in_E  = cell2mat(mtab.Data(:,3))';
    peak_to_use    = find(cell2mat(mtab.Data(:,7))'*1);

    % method-2
    xdata = d_ideal(peak_to_use)';           % ideal d-spacing
    yobserv = peak_pos_in_E(peak_to_use);    % observed E_hkl (keV)

    fn = {'eddBragg2'}; npar = [1];
    pfix = [nan]; % fix fitting parameters (currently none are fixed)
    opt = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-8,'tolf',1e-8,'MaxIter',500);
    p0  = det_TOA;
    pUB = [10.0];
    pLB = [2.0];

    [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p0,pLB,pUB,opt,pfix,npar,fn,xdata,yobserv,ones(size(yobserv)));

    %yfit   = sumfun1(pfit,pfix,npar,fn,xdata);
    %ycalc  = sumfun1(pfit,pfix,npar,fn,ideal_d(peak_to_use));   % ideal Energy (KeV)
    
    %d_observe = hc/2./peak_pos_in_E/sind(pfit/2);
    E_ideal = hc/2./d_ideal/sind(pfit/2);
    E_error = (peak_pos_in_E-E_ideal)*1E3; % (in eV)

    % assign fit value
    mtab.Data(:,8) = num2cell(E_error);
    h.det_table.Data(1) = pfit;
    
    % plot result
    fig = figure(8); clf(fig, 'reset');
    set(fig,'Tag','eddcaliTOA_fig_fitTOA');

    %plot(eddBragg2(pfit,ideal_d),d_observe-ideal_d,'marker','o')
    %line(eddBragg2(pfit,ideal_d((peak_to_use))),d_observe(peak_to_use)-ideal_d(peak_to_use),'marker','x','color','r','linestyle','none')
    %xlabel('Energy (KeV)');
    %ylabel('\Deltad-spacing (A)');

    plot(eddBragg2(pfit,d_ideal), E_error,'marker','o')
    line(eddBragg2(pfit,d_ideal((peak_to_use))),E_error(peak_to_use),'marker','x','color','r','linestyle','none')
    xlabel('Energy (KeV)');
    ylabel('E_{fit} - E_{ideal} (eV)');
    
    title(sprintf('TOA = %2.15f, rn = %f',pfit(1),rn))
    xlim([30 160])
    ylim([min([-50 min(E_error)]) max([50 max(E_error)])])

    
function main_accept_Callback(~,~,~)
h = get_handle;
det_TOA = h.det_table.Data(1);
detno = h.main_pop_detno.Value;

hroot = findall(0,'Tag','eddviewmain_Fig');
hroot = getappdata(hroot,'allhandle');
switch detno
    case 1
        %hroot.expinfo_edit_toa1.String = num2str(det_TOA,'%1.8f');
        %hroot.config.visopt.Inst(1).TOA = det_TOA;
        hroot.config.visopt.Inst.detpar(1,1) = det_TOA;
        hroot.expinfo_table_detpar.Data{1,1} = det_TOA;
    case 2
        %hroot.expinfo_edit_toa2.String = num2str(det_TOA,'%1.8f');
        %hroot.config.visopt.Inst(2).TOA = det_TOA;
        hroot.config.visopt.Inst.detpar(2,1) = det_TOA;
        hroot.expinfo_table_detpar.Data{2,1} = det_TOA;
    otherwise
        warndlg('Not supportted yet!!')
end
    

function uitable_det_edit_Callback(~,eventdata,~)
    % constants
    hc = 12.398419057638671;

    h = get_handle;
    if_guess = h.main_check_guess.Value;
    detpar = h.det_table.Data;
    mtab = h.par_table;
    col_select = eventdata.Indices(2);
    switch col_select
        case 1
            d_ideal = cell2mat(mtab.Data(:,1));
            if if_guess
                mtab.Data(:,3) = num2cell(hc/2./d_ideal/sind(detpar(1)/2));
            else
                mtab.Data(:,2) = num2cell(hc/2./d_ideal/sind(detpar(1)/2));
                mtab.Data(:,8) = num2cell((cell2mat(mtab.Data(:,3))-hc/2./d_ideal/sind(detpar(1)/2))*1E3);
            end
%         case 2
%             fprintf('change gain\n')
%             d_ideal = cell2mat(mtab.Data(:,1));
%             E_guess_old = cell2mat(mtab.Data(:,3));
%           case 3
%             d_ideal = cell2mat(mtab.Data(:,1));
%             E_ideal = hc/2./d_ideal/sind(detpar(1)/2);
%             mtab.Data(:,3) = num2cell(E_ideal+eventdata.NewData);
    end

    update_handle(h);

%%%% call back to report selection
function uitable_par_select_Callback(~,eventdata,~)
    h = get_handle; 
    h.opt.table_selection = eventdata.Indices;
    fprintf('Select row: %d  col: %d\n', eventdata.Indices);
    update_handle(h);

function uitable_det_select_Callback(~,eventdata,~)
    h = get_handle;
    h.opt.table_selection = eventdata.Indices;
    fprintf('Select col: %d\n', eventdata.Indices(2));
    update_handle(h);
    
%==========================================================================
% --- close request function of eddgetparmain fig 
%==========================================================================
function eddparmain_CloseRequestFcn(hObject,eventdata)                         %#ok<INUSD>
    delete(findall(0,'Tag','eddcaliTOAmain_Fig'))
    delete(findall(0,'-regexp','Tag','eddgetpar*'))
    delete(findall(0,'-regexp','Tag','eddcaliTOA_fit*'))

%==========================================================================
% --- Load all handle object function
%==========================================================================
function h = get_handle
    hmain = findall(0,'Tag','eddcaliTOAmain_Fig');
    h = getappdata(hmain,'allhandle');

%==========================================================================
% --- Save all handle object function
%==========================================================================
function update_handle(handle_to_save)
    hmain = findall(0,'Tag','eddcaliTOAmain_Fig');
    setappdata(hmain,'allhandle',handle_to_save);
