function [GO_parents, GO_children] = FindGOParent(GeneOntology_relation,GeneAssociation,spectra)
%[GO_PARENTS, GO_CHILDREN] =
%   FINDGOPARENT(GENEONTOLOGY_RELATION,GENEASSOCIATION,SPECTRA) identifies
%   GO terms that are parents (GO_PARENTS) of other GO terms (GO_CHILDREN)
%   but account for the same fraction.
%
%   GENEONTOLOGY_RELATION:  a structure specifying parent-to-child relation
%   between GO terms
%
%   GENEASSOCIATION:    a structure specifying gene-to-GO relation
%
%   SPECTRA:   an array of SPECTRA for which the parent GO terms are to
%   identified.
%

[unigene,freq_gene] = unique2(spectra);
uniGO = unique(GeneAssociation.GO);

%construct the matrix containing abundance info
matrix = false(length(unigene),length(uniGO));
for i = 1:length(unigene)
    ind = strmatch(unigene(i),GeneAssociation.gene);
    zgo = GeneAssociation.GO(ind);
    [tf,loc] = ismember(zgo,uniGO);
    loc(loc==0) = [];
    matrix(i,loc) = true;
end

%construct matrix containing children-parent info
matrix_cp = false(length(uniGO));
for i = 1:length(uniGO)
    tf = strcmp(uniGO{i},GeneOntology_relation.GO1);
    [tf,loc] = ismember(uniGO,GeneOntology_relation.GO2(tf));
    matrix_cp(i,tf) = true;
end

%generate indice for all combinations of two GO terms
ind = combinator(int8(length(uniGO)),2,'c');

ind_pc = [];%indice for parent and children
for i = 1:size(ind,1)
    zm = matrix(:,ind(i,:));
    zs = sum(zm,2);
    if sum(zs==1)==0
        if matrix_cp(ind(i,1),ind(i,2))
            ind_pc = [ind_pc;ind(i,[2 1])];
        elseif matrix_cp(ind(i,2),ind(i,1))
            ind_pc = [ind_pc;ind(i,:)];
        end
    end
end

GO_parents = uniGO(ind_pc(:,1));
GO_children = uniGO(ind_pc(:,2));
return;
    
