function [GO_lists_k,ind_lists_k,fractionOfSector_lists_k,degree_overlap_lists_k] = EnumerateLists_k(GeneAssociation_sector,spectra_sector,k)
%[GO_LISTS_K,IND_LISTS_K,FRACTIONOFSECTOR_LISTS_K,DEGREE_OVERLAP_LISTS_K] =
%   ENUMERATELISTS_K(GENEASSOCIATION_SECTOR,SPECTRA_SECTOR,K) returns
%   GO_LISTS_K, a cell array of GO terms, IND_LISTS_K, a double matrix
%   specifying to a list's indices of GO terms, FRACTIONOFSECTOR_LISTS_K, a
%   double array specifying the fraction of sector for the list,
%   DEGREE_OVERLAP_LISTS_K, a double array specifying the degree of overlap
%   for the lists, and DEGREE_OVERLAP_LISTS_K, a double array specifying
%   the gene coverage for the lists.
%
%   GENEASSOCIATION_SECTOR:    a structure specifying gene-to-GO relation
%
%
%   SPECTRA_SECTOR:   list of spectra for the sector
%
%   K:  a scalar specifying the size of the GO list


[genes_sector,freq_genes] = unique2(spectra_sector);
GO_lists_k = unique(GeneAssociation_sector.GO);

%contruct genes to GO terms relation matrix
matrix = false(length(genes_sector),length(GO_lists_k));
for i = 1:length(genes_sector)
    ind = strmatch(genes_sector(i),GeneAssociation_sector.gene);
    GOTerms_temp = GeneAssociation_sector.GO(ind);
    [tf,loc] = ismember(GOTerms_temp,GO_lists_k);
    loc(loc==0) = [];
    matrix(i,loc) = true;
end

%generate indice for all lists with length of k
ind_lists_k = combinator(int8(length(GO_lists_k)),k,'c');

%
fractionOfSector_lists_k = zeros(size(ind_lists_k,1),1);
degree_overlap_lists_k = zeros(size(ind_lists_k,1),1);
for i = 1:size(ind_lists_k,1)
    submatrix = matrix(:,ind_lists_k(i,:));
    sum_temp = sum(submatrix,2);
    fractionOfSector_lists_k(i) = freq_genes'*(~~sum_temp);
    fraction1 = freq_genes'*submatrix;%fraction of sector that each GO term accounts for
    submatrix(sum_temp>1,:) = 0;
    fraction2 = freq_genes'*submatrix;%fraction of sector that each GO term accounts for uniquely
    degree_overlap_lists_k_GOTerms = 1 - fraction2./fraction1;%degree of overlap for each GO term (i.e., O_{i,l,t} in the text)
    degree_overlap_lists_k(i) = max(degree_overlap_lists_k_GOTerms);%degree of overlap for the list (i.e., O_{i,l} in the text)
end

fractionOfSector_lists_k = fractionOfSector_lists_k/length(spectra_sector);

return;
    
