function sub=ind2subV(N,ind)

%IND2SUBV Implements a vectorial version of the ind2sub function
%The method inherits from the example in 
%https://stackoverflow.com/questions/8918137/return-subscripts-of-a-variable-dimension-matrix
%   SUB=IND2SUBV(N,IND)
%   * N is the size of the multidimensional array the subscripts will correspond to
%   * IND is a set of indexes to the multidimensional array
%   * SUB is a matrix of subscripts arranged as # indexes times # dimensions
%

ND=numel(N);
ind=double(ind(:));
assert(gather(all(ind>0 & ind<=prod(N))),'Indexes outside the range given by dimensions sizes');
c=cell([1 ND]);
if length(N)>1;[c{:}]=ind2sub(N,ind);else c{:}=ind;end
sub=cat(2,c{:});
