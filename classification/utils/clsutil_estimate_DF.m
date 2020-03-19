function [ kest, Dest ] = clsutil_estimate_DF( xTr,yTr,D )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
for ii=1:nClasses
    vars(:,ii)=var(xTr(:,yTr(ii,:)==1),[],2);
end

kest=[2*(1-D)./vars(:,1),(2*D./vars(:,2))];
kest(D>0.5,1)=kest(D>0.5,2);
kest=kest(:,1);
Dest=D;
end

