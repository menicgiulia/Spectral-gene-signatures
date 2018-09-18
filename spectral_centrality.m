function sc=spectral_centrality(A,k)
%Inputs:  
% A = square symmetric adjacency matrix
% k = index of the centrality (default = 1)
 


if nargin==1
    k = 1;
end

%Graph Laplacian
L=diag(sum(A))-A;
%Compute its eigendata
[V,E]=eig(L);
[~,idx]=sort(diag(E));
V=V(:,idx);
%relevant eigenvector
vec=V(:,k+1);

%node centralities
sc=sum(A.*((vec*ones(1, size(A,1)))-(ones(size(A,1),1)*(vec'))).^2,2);