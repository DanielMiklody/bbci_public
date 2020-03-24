function [ kest, Dest, kfull ] = clsutil_estimate_DF( xTr,yTr,D, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin>3
    shrinkage=varargin{1};
else
    shrinkage=0;
end

switch nargin
    case 1
        if isstruct(xTr)
            yTr=xTr.y;
            xTr=squeeze(xTr.x);
            D=mean(xTr(:,yTr(2,:)==1),2);        
        else
            error('Unsupported arguments.');
        end
    case 2
        if isstruct(xTr)
            D=yTr;
            yTr=xTr.y;
            xTr=squeeze(xTr.x);
        else
            D=mean(xTr(:,yTr(2,:)==1),2);
        end
end


nClasses= size(yTr,1);
nChannels=size(xTr,1);
vars=zeros(nChannels,nClasses);
if ~shrinkage
    for ii=1:nClasses
        vars(:,ii)=var(xTr(:,yTr(ii,:)==1),[],2);
    end
else
    for ii=1:nClasses
        vars(:,ii)=diag(clsutil_shrinkage(xTr(:,yTr(ii,:)==1)));
    end
end

kfull=[2*(1-D).^2./vars(:,1),(2*D.^2./vars(:,2))];
kest=kfull;
kest(D>0.5,1)=kest(D>0.5,2);
kest=kest(:,1);
Dest=D;
end

