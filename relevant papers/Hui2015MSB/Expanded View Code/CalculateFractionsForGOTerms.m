function [GO_terms,fractionOfSector_terms,fractionOfProteome_terms] = CalculateFractionsForGOTerms(GeneAssociation_complete,spectra_proteome,spectra_sector)
%[GO_TERMS,FRACTIONOFSECTOR_TERMS,FRACTIONOFPROTEOME_TERMS] =
%   CALCULATEFRACTIONSFORGOTERMS(GENEASSOCIATION_COMPLETE,SPECTRA_PROTEOME,
%   SPE CTRA_SECTOR) returns GO_TERMS,the list of unique GO terms,
%   FRACTIONOFSECTOR_TERMS, fraction of sector for each of the unique GO
%   terms, and FRACTIONOFPROTEOME_TERMS, fraction of proteome for each of
%   the unqiue GO terms.
%
%   GENEASSOCIATION_COMPLETE:   a structure indicating the association
%   between genes and their correponding GO terms. Note that it is a
%   multiple-to-multiple relationship.
%
%   SPECTRA_PROTEOME: cell array of strings; each cell is a gene name,
%   representing a spectrum in the protein abundance data set. The absolute
%   abundance of a protein is calculated as the number of occurences of the
%   protein in the array divided by the length of the cell array.
%
%   SPECTRA_SECTOR:   list of spectra for the sector
%

%Obtain "fraction of proteome" (Omega_t) for a GO term t
[genes_proteome,freq_genes_proteome] = unique2(spectra_proteome);
GO_terms_proteome = [];
for i = 1:length(genes_proteome)
    ind = strmatch(genes_proteome{i},GeneAssociation_complete.gene,'exact');
    GO_terms_proteome = [GO_terms_proteome;repmat(GeneAssociation_complete.GO(ind),freq_genes_proteome(i),1)];
end
[GO_terms_proteome,freq_GO_terms_proteome] = unique2(GO_terms_proteome);
fraction_GO_terms_proteome = freq_GO_terms_proteome/length(spectra_proteome);

%Obtain "fraction of sector" (omega_i,t) for a GO term t and for a sector i
[genes_sector,freq_genes_sector] = unique2(spectra_sector);
GO_terms = [];
for i = 1:length(genes_sector)
    ind = strmatch(genes_sector{i},GeneAssociation_complete.gene,'exact');
    GO_terms = [GO_terms;repmat(GeneAssociation_complete.GO(ind),freq_genes_sector(i),1)];
end
[GO_terms,freq_GO_terms] = unique2(GO_terms);
fractionOfSector_terms = freq_GO_terms/length(spectra_sector);

%focus only on the GO terms that are present in the sector
[tf,loc] = ismember(GO_terms,GO_terms_proteome);
fractionOfProteome_terms = fraction_GO_terms_proteome(loc);


return;