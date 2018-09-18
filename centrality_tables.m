function [Tnode, Tnodeselection]=centrality_tables(A, genes, mutatedgenes,nc)

    
    Tnode=table;
    k0=(sum(A)==0);
    A=A(~k0, ~k0);
    G=graph(A);
    n = numnodes(G);
    gnot0=genes(~k0);

    a=conncomp(G)';
    uc=unique(a);
    b=zeros(size(uc));
    for ind=1:length(uc)
        b(ind)=sum(a==uc(ind));
    end

    Tnode.components=a;
    Tnode.K=sum(A,2);

    wbc = centrality(G,'betweenness','Cost',(G.Edges.Weight).^(-1));
    Tnode.Bninv=2*wbc./((n-2)*(n-1));
    
    Tnode.SCN1=zeros(height(Tnode),1);
    for ind=1:length(b)
        cluster=A(a==ind, a==ind); 
        %spectral centrality k=1
        spectral_1=spectral_centrality(cluster);
        Tnode.SCN1(a==ind)=spectral_1;
    end

    Tnode.Properties.RowNames=gnot0;
    Tnode.mutation=ismember(Tnode.Properties.RowNames,mutatedgenes);
    
    [~,ib]=sort(b,'descend');
    compsel=ib(1:nc);
    Tnodeselection=Tnode(ismember(Tnode.components, compsel),:);

end