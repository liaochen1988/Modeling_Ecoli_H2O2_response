function [b,freq,idx]=unique2(a,flag)
%UNIQUE2   Set unique
%   [B,F,I]=UNIQUE2(A) returns the unique values B and occuring frequencies
%   F. B is sorted. 
%
%   [B,F,I]=UNIQUE2(A,'rows') returns the unique rows and corresponding
%   frequencies in A. B is sorted along from the 1st column to the last
%   column.
%
%   See also UNIQUE
%
%   Tony
%   08/07/2004

r=a;
if(nargin==1)
    [b,m,n]=unique(sort(a(:)));
    [tf,idx]=ismember(b,r(:));
elseif(nargin==2)
    if ~strcmpi(flag,'rows')
        error('Unknown flag.');
    end
    [rows,cols]=size(a);
    a=sortrows(a,[1:cols]);
    [b,m,n]=unique(a,flag);
    [tf,idx]=ismember(b,r,flag);
else
    error('Wrong number of input arguments.');
end
if(length(m)==1)
    freq=m;
else
    m2=m(2:end);
    freq=[m(1);m2-m(1:end-1)];
end
    