function colloc = collocate(knots,k,x,nd)
% function to determine spline basis functions based on specified knots, spline order and data sites.
%    Usage: colloc = collocate2(knots,k,x,nd)
%             where the knot series "knots" has appropriate multiplicity for order k at beginiing
%             and end, x are data sites, and nd is the derivative count (0 for no
%             derivatives, 1 for first derivative only, 2 for first and second, etc.
% 
%             colloc is an array of size 
%                 (number of data sites, number of coefficients, 1+number of requested derivatives).  
%                 Thus the array(:,:,1) is the function, array(:,:,2) is the first derivative, etc.
%  JMB 2019 - based on de Boor's book and help from the mathworks toolbox.
%              No dependencies

nd=nd+1; % 
npk=length(knots); %number of knots
n=npk-k;  % number of control points
npts = length(x); %number of data
km1 = k-1;

% Error check the input:
dk=diff(knots);  
dx=diff(x);
id=find(dk>0);

if(min(dk)<0),error('knots do not increase sequentially'),end
if(nd>k),error('requested derivative is too high for spline order)');end   
if(id(1)~=k),error('wrong knot multiplicity at beginning of sequence'),end
if(id(end)~=npk-k),error('wrong knot multiplicity at end of sequence'),end
if(min(dx)<0),error('data series does not increase sequentially'),end
if(x(1)<knots(1)),error('data site lies below first knot'),end
if(x(end)>knots(end)),error('data site lies beyond last knot'),end

%Determine the left side knot index for each data point
naug = length(knots)-k;
pts = x(:); knots = knots(:);
meshsites=knots(1:naug);
[~,index] = sort([meshsites(:).' pts(:).']);
savl = max(find(index>length(meshsites))-(1:length(pts)),k);

% create block array for basis functions-one row for each data point and k wide
b = zeros(npts,k);

if nd==1 % If no derivatives are required do less work:
   % initialize the  b  array.
   b(:,1) = ones(npts,1);
   % run the recurrence simultaneously for all  pts .
   for j=1:km1
      saved = zeros(npts,1);
      for r=1:j
         tr = knots(savl+r)-pts;
         tl = pts-knots(savl+r-j);
         term = b(:,r)./(tr+tl);
         b(:,r) = saved+tr.*term;
         saved = tl.*term;
      end
      b(:,j+1) = saved;
   end

else % if derivatives are required, more calculations are done
% initialize the  bb  array.
   bb = repmat([1 zeros(1,km1)],nd*npts,1);
   lptss = nd*[1:npts];
% run the recurrence simultaneously for all  pts .
   for j=1:k-nd
      saved = zeros(npts,1);
      for r=1:j
         tr = knots(savl+r)-pts;
         tl = pts-knots(savl+r-j);
         term = bb(lptss,r)./(tr+tl);
         bb(lptss,r) = saved+tr.*term;
         saved = tl.*term;
      end
      bb(lptss,j+1) = saved;
   end
% save the B-spline values in successive blocks in  bb .
   for jj=1:nd-1
      j = k-nd+jj; saved = zeros(npts,1); lptsn = lptss-1;
      for r=1:j
         tr = knots(savl+r)-pts;
         tl = pts-knots(savl+r-j);
         term = bb(lptss,r)./(tr+tl);
         bb(lptsn,r) = saved+tr.*term;
         saved = tl.*term;
      end
      bb(lptsn,j+1) = saved; lptss = lptsn;
   end
% now use the fact that derivative values can be obtained by differencing:
   for jj=nd-1:-1:1
      j = k-jj;
      temp = repmat([jj:nd-1].',1,npts)+repmat(lptsn,nd-jj,1); lptss=temp(:);
      for r=j:-1:1
         temp = repmat((knots(savl+r)-knots(savl+r-j)).'/j,nd-jj,1);
         bb(lptss,r) = -bb(lptss,r)./temp(:);
         bb(lptss,r+1) = bb(lptss,r+1) - bb(lptss,r);
      end
   end 
   b=bb;
end

% At this point b contains the block of basis functions (for all derivative
% levels and savl contains the index of the position of each set of basis
% functions.  The trick is to load the data into an array that is n wide by
% npts long by nd deep.  This can be done using indexing and repmat.  Here
% its done with do loops - something that can be fixed.
colloc=zeros(npts,n,nd);
ids=savl(:)-km1+(0:km1);
for i=1:nd
  for j=1:npts
      colloc(j,ids(j,:),i)=b(i+nd*(j-1),:);
  end
end



