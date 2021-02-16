%
% Peak fitting function for EDD detector calibration
%
% by Andrew Chuang, 

function p = edd_fit_emission(p0,xxdata,yydata,opt)

if nargin < 3
    fprintf('Not enough input parameters!!\n');
    return;
elseif nargin == 3
    vis_fit = 0;
elseif nargin == 4
    vis_fit= opt.vis_fit;
else
    fprintf('Too many input parameters!!\n');
    return;
end

pk_id  = p0.pkid;
pcen0  = p0.cen;
pint0  = p0.int;
pfwhm0 = p0.fwhm;

if isfield(opt,'dup')
    dup = opt.dup;
else
    dup = repmat(50,1,length(pk_id));
end

if isfield(opt,'dup')
    ddown = opt.ddown;
else
    ddown = repmat(50,1,length(pk_id));
end

for i = 1:length(pk_id)
    npeaks = length(pcen0{i});
    [~, idx_min] = min(abs(xxdata-min(pcen0{i})));
    [~, idx_max] = min(abs(xxdata-max(pcen0{i})));
    %x_ind = round(min(pcen0{i}),0)-ddown(i):round(max(pcen0{i}),0)+dup(i);
    x_ind = max([idx_min-ddown(i) 1]):min([idx_max+dup(i) length(xxdata)]);
    xdata = xxdata(x_ind);
    ydata = yydata(x_ind);
    
    switch npeaks
        case 1
%             fn = {'psv1'}; npar = [4];
%             pfix = [nan nan nan nan]; % fix fitting parameters (currently none are fixed)
%             opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
%             p02 =  [pint0{i}     pcen0{i}    pfwhm0{i} 0.5];
%             pUB2 = [pint0{i}*10  pcen0{i}+10       inf 1.0];
%             pLB2 = [      0      pcen0{i}-10         0   0];
%             y0   = sumfun1(p02,pfix,npar,fn,xdata);
%             
%             [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
%             yfit = sumfun1(pfit,pfix,npar,fn,xdata);
%             txt_int  = sprintf('Int: %2.1f',pfit(1));
%             txt_cen  = sprintf('Cen: %2.1f',pfit(2));
%             txt_fwhm = sprintf('FWHM: %2.1f',pfit(3));
%             txt_gratio = sprintf('Ratio: %2.1f',pfit(4));
%             %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(5),pfit(6));
%             p.int(pk_id{i}) = pfit(1);
%             p.cen(pk_id{i}) = pfit(2);
%             p.fwhm(pk_id{i}) = pfit(3);
%             p.gratio(pk_id{i}) = pfit(4);
            fn = {'gs1' 'backg1'}; npar = [3 1];
            pfix = [nan nan nan nan]; % fix fitting parameters (currently none are fixed)
            opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
            p02 =  [pint0{i}     pcen0{i}    pfwhm0{i} median(yydata([1:5 end-4:end]))];
            pUB2 = [pint0{i}*10  pcen0{i}+10       inf max(yydata)/2];
            pLB2 = [      0      pcen0{i}-10         0 -max(yydata)];
            y0   = sumfun1(p02,pfix,npar,fn,xdata);
            
            [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
            yfit = sumfun1(pfit,pfix,npar,fn,xdata);
            txt_int  = sprintf('Int: %2.1f',pfit(1));
            txt_cen  = sprintf('Cen: %2.1f',pfit(2));
            txt_fwhm = sprintf('FWHM: %2.1f',pfit(3));
            %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(4),pfit(4));
            txt_bg   = sprintf('BG: %2.6f / %2.6f',0,pfit(4));
            p.int(pk_id{i}) = pfit(1);
            p.cen(pk_id{i}) = pfit(2);
            p.fwhm(pk_id{i}) = pfit(3);
            p.bg(pk_id{i}) = pfit(4);
       case 2
%             fn = {'psv1','psv1'}; npar = [4 4];
%             pfix = [nan nan nan nan nan nan nan nan]; % fix fitting parameters (currently none are fixed)
%             opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
%             p02 =  [pint0{i}(1)     pcen0{i}(1)    pfwhm0{i}(1) 0.5 pint0{i}(2)     pcen0{i}(2)    pfwhm0{i}(2) 0.5];
%             pUB2 = [pint0{i}(1)*10  pcen0{i}(1)+10          inf 1.0 pint0{i}(2)*10  pcen0{i}(2)+10       inf    1.0];
%             pLB2 = [      0         pcen0{i}(1)-10            0   0          0      pcen0{i}(2)-10         0      0];
%             y0   = sumfun1(p02,pfix,npar,fn,xdata);
%             
%             [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
%             yfit = sumfun1(pfit,pfix,npar,fn,xdata);
%             txt_int  = sprintf('Int: %4.1f/%4.1f',pfit(1),pfit(5));
%             txt_cen  = sprintf('Cen: %4.1f/%4.1f',pfit(2),pfit(6));
%             txt_fwhm = sprintf('FWHM: %4.1f/%4.1f',pfit(3),pfit(7));
%             txt_gratio = sprintf('Ratio: %2.1f/%2.1f',pfit(4),pfit(8));
%             %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(5),pfit(6));
%             p.int(pk_id{i}(1)) = pfit(1);
%             p.cen(pk_id{i}(1)) = pfit(2);
%             p.fwhm(pk_id{i}(1)) = pfit(3);
%             p.gratio(pk_id{i}(1)) = pfit(4);
%             p.int(pk_id{i}(2)) = pfit(5);
%             p.cen(pk_id{i}(2)) = pfit(6);
%             p.fwhm(pk_id{i}(2)) = pfit(7);
%             p.gratio(pk_id{i}(2)) = pfit(8);
            fn = {'gs1' 'gs1' 'backg1'}; npar = [3 3 1];
            pfix = [nan nan nan nan nan nan nan]; % fix fitting parameters (currently none are fixed)
            opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
            p02 =  [pint0{i}(1)     pcen0{i}(1)    pfwhm0{i}(1) pint0{i}(2)     pcen0{i}(2)    pfwhm0{i}(2) median(yydata([1:5 end-4:end]))];
            pUB2 = [pint0{i}(1)*10  pcen0{i}(1)+10          inf pint0{i}(2)*10  pcen0{i}(2)+10       inf    max(yydata)/2];
            pLB2 = [      0         pcen0{i}(1)-10            0          0      pcen0{i}(2)-10         0   -max(yydata)];
            y0   = sumfun1(p02,pfix,npar,fn,xdata);
            
            [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
            yfit = sumfun1(pfit,pfix,npar,fn,xdata);
            txt_int  = sprintf('Int: %4.1f/%4.1f',pfit(1),pfit(4));
            txt_cen  = sprintf('Cen: %4.1f/%4.1f',pfit(2),pfit(5));
            txt_fwhm = sprintf('FWHM: %4.1f/%4.1f',pfit(3),pfit(6));
            %txt_gratio = sprintf('Ratio: %2.1f/%2.1f',pfit(4),pfit(8));
            %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(7),pfit(7));
            txt_bg   = sprintf('BG: %2.6f / %2.6f',0,pfit(7));
            p.int(pk_id{i}(1)) = pfit(1);
            p.cen(pk_id{i}(1)) = pfit(2);
            p.fwhm(pk_id{i}(1)) = pfit(3);
            p.int(pk_id{i}(2)) = pfit(4);
            p.cen(pk_id{i}(2)) = pfit(5);
            p.fwhm(pk_id{i}(2)) = pfit(6);
            p.bg(pk_id{i}(2)) = pfit(7);
        case 3
            fn = {'gs1' 'gs1' 'gs1' 'backg1'}; npar = [3 3 3 1];
            pfix = [nan nan nan nan nan nan nan nan nan nan]; % fix fitting parameters (currently none are fixed)
            opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
            p02 =  [pint0{i}(1)     pcen0{i}(1)    pfwhm0{i}(1) pint0{i}(2)     pcen0{i}(2)    pfwhm0{i}(2) pint0{i}(3)     pcen0{i}(3)    pfwhm0{i}(3) median(yydata([1:5 end-4:end]))];
            pUB2 = [pint0{i}(1)*10  pcen0{i}(1)+10          inf pint0{i}(2)*10  pcen0{i}(2)+10          inf pint0{i}(3)*10  pcen0{i}(3)+10          inf   max(yydata)/2];
            pLB2 = [      0         pcen0{i}(1)-10            0          0      pcen0{i}(2)-10            0          0      pcen0{i}(3)-10            0  -max(yydata)];
            y0   = sumfun1(p02,pfix,npar,fn,xdata);
            
            [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
            yfit = sumfun1(pfit,pfix,npar,fn,xdata);
            txt_int  = sprintf('Int: %4.1f/%4.1f/%4.1f',pfit(1),pfit(4),pfit(7));
            txt_cen  = sprintf('Cen: %4.1f/%4.1f/%4.1f',pfit(2),pfit(5),pfit(8));
            txt_fwhm = sprintf('FWHM: %4.1f/%4.1f/%4.1f',pfit(3),pfit(6),pfit(9));
            %txt_gratio = sprintf('Ratio: %2.1f/%2.1f',pfit(4),pfit(8));
            %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(7),pfit(7));
            txt_bg   = sprintf('BG: %2.6f / %2.6f / %2.6f',0,pfit(10));
            p.int(pk_id{i}(1))  = pfit(1);
            p.cen(pk_id{i}(1))  = pfit(2);
            p.fwhm(pk_id{i}(1)) = pfit(3);
            p.int(pk_id{i}(2))  = pfit(4);
            p.cen(pk_id{i}(2))  = pfit(5);
            p.fwhm(pk_id{i}(2)) = pfit(6);
            p.int(pk_id{i}(3))  = pfit(7);
            p.cen(pk_id{i}(3))  = pfit(8);
            p.fwhm(pk_id{i}(3)) = pfit(9);            
            p.bg(pk_id{i}(3))   = pfit(10);
        otherwise
            fprintf('%d peak fitting is not supported yet!!',npeaks)
            return
    end
    
    if vis_fit
        figure(2);
        subplot(4,1,1:3)
        plot(xdata,ydata,'o',xdata,y0,'g',xdata,yfit,'b')
        set(gca,'yscale','log');
        ylim([max([min(ydata)/5 1E-4]) max(ydata)*5])
        text(0.02,0.95,txt_int,'sc')
        text(0.02,0.90,txt_cen,'sc')
        text(0.02,0.85,txt_fwhm,'sc')
        %text(0.02,0.80,txt_gratio,'sc')
        text(0.02,0.80,txt_bg,'sc')
        subplot(4,1,4)
        plot(xdata,rs/pfit(1),'-x')
        text(0.05,0.1,sprintf('residual norm: %4.2f',rn),'sc');
        waitforbuttonpress
        %pause(0.5)
    end
end