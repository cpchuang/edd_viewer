function edd_calibrant(varargin)
% GUI to calibrate detector
%
%   + 1.0.2 2020/01/31
%           - [bug fix]send proper scan number to "edd_cali_TOA"
%   + 1.0.1 2019/10/23
%           - Energy calibration can do individual step
%   + 1.0   2017/11/03
%           - initial version
%
% Copyright 2017-2020 Andrew Chuang (chuang.cp@gmail.com)
% $Revision: 1.0.2 $  $Date: 2020/01/31 $

% choice = questdlg('Which detector do you want to calibrate?', ...
%     'Choose Detector', ...
%     'Det-1 (V)','Det-2 (H)','None','Det-1 (V)');
% % Handle response
% switch choice
%     case 'Det-1 (V)'
%         disp([choice ', will be calibrated!!'])
%         detno = 1;
%     case 'Det-2 (H)'
%         disp([choice ', will be calibrated!!'])
%         detno = 2;
%     case 'None'
%         disp('Okay!! Bye....')
%         return;
% end

choice = questdlg('Which type of calibration do you want to do?', ...
    'Choose Calibration', ...
    'Energy Calibration','TOA Calibration','None','Energy Calibration');
% Handle response
switch choice
    case 'Energy Calibration'
        disp(['Do ' choice '!!'])
        h = varargin{1};
        if ~isfield(h,'data')
            warndlg({'Please read the data first!!'});
            return;
        end
        opt.data = h.data.data;
        opt.detno = 1;
        opt.source = 1;
        opt.xpos_handle = h.visopt_edit_xpos;  % pass xpos handle to "cali_Energy"
        edd_cali_energy(opt)
    case 'TOA Calibration'
        %warndlg({'This is a paid function!!','Please upgrade to the full version.'},'!! Warning !!')
        %return;
        disp(['Do ' choice '!!'])
        h = varargin{1};
        if ~isfield(h,'data')
            warndlg({'Please read the data first!!'});
            return;
        end
        opt.data = h.data(h.config.visopt.scno).data;
        %opt.Inst(1).par = [h.config.visopt.Inst(1).TOA 0 h.config.visopt.Inst(1).Ch2E];
        %opt.Inst(2).par = [h.config.visopt.Inst(2).TOA 0 h.config.visopt.Inst(2).Ch2E];
        opt.detpar = h.config.visopt.Inst.detpar;
        opt.detno = 1;
        opt.source = 1;
        opt.posno = h.config.visopt.posno;
        opt.mat = 1;
        opt.exp_time = h.data.exp_time;

        edd_cali_TOA(opt);
    case 'None'
        disp('Okay!! Bye....')
        return;
end


% %%%% select which type of calibration
% str = {'Co-57','Cd-109'};
% [selection,~] = listdlg('PromptString','Select type',...
%                 'SelectionMode','single',...
%                 'ListString',str);
% if isempty(selection)            
%     disp("No selection was made\n")
%     return
% end
%         
% switch selection
%     case 1
% %         prompt = {'Enter Channel Number','Enter colormap name:'};
% %         dlg_title = 'Energy = a * Channel + b';
% %         num_lines = 1;
% %         defaultans = {'20','20'};
% %         answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%     case 2
%         disp("Only work for Co-57 as of today!!")
%         return
% end




