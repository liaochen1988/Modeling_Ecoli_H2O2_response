%load data
load data.mat;
%   It includes the following variables.
%
%   GeneAssociation_complete:   a structure specifying gene-to-GO relation
%
%   GeneOntology_relation:  a structure specifying parent-child relation
%   between GO terms
%
%   spectra_proteome:   list of spectra for the whole proteome, specifying
%   abundance for each protein. Each spectrum is represented as a gene
%   name. The abundance for a protein is calculated by the number of
%   occurence of its gene name divided by the size of the list.
%
%   sectorC:    list of genes belonging to C-sector
%   sectorA:    list of genes belonging to A-sector
%   sectorR:    list of genes belonging to R-sector
%   sectorU:    list of genes belonging to U-sector
%   sectorS:    list of genes belonging to S-sector
%   sectorO:    list of genes belonging to O-sector

%specify the set of genes the GO analysis is carried out for
genes_sector = sectorC;

%obtain the spectra for the sector
[tf,loc] = ismember(spectra_proteome,genes_sector);
spectra_sector = spectra_proteome(tf);%list of spectra for the sector

%==========filter out GO terms with the four criteria==========

[GO_terms,fractionOfSector_terms,fractionOfProteome_terms]...
    = CalculateFractionsForGOTerms(GeneAssociation_complete,spectra_proteome,spectra_sector);

%apply the first three filters
tf1 = fractionOfSector_terms<0.1;%first criterion
tf2 = fractionOfProteome_terms>0.4;%second criterion
tf3 = fractionOfSector_terms<fractionOfProteome_terms;%third criterion
tf = tf1|tf2|tf3;
GO_terms(tf) = [];
fractionOfSector_terms(tf) = [];
fractionOfProteome_terms(tf) = [];
%simplify the structure GeneAssociation_complete, to only focus on the GO
%terms that are the sector and passed the filters
[tf,loc] = ismember(GeneAssociation_complete.GO,GO_terms);
GeneAssociation_sector.gene = GeneAssociation_complete.gene(tf);
GeneAssociation_sector.GO = GeneAssociation_complete.GO(tf);

%apply the fourth filter
[GO_parents,GO_children] = FindGOParent(GeneOntology_relation,GeneAssociation_sector,spectra_sector);
[tf,loc] = ismember(GeneAssociation_sector.GO,GO_parents);
GeneAssociation_sector.GO(tf) = [];
GeneAssociation_sector.gene(tf) = [];
[tf,loc] = ismember(GO_terms,GO_parents);
GO_terms(tf) = [];
fractionOfSector_terms(tf) = [];
fractionOfProteome_terms(tf) = [];

%==========filter out GO term lists with three criteria==========
GO_lists = [];
ind_lists = [];
fractionOfSector_lists = [];
degree_overlap_lists = [];%
gene_coverage_lists = [];
k = 0;
degree_overlap_lists_k = 0;
while min(degree_overlap_lists_k)<0.05
    k = k+1;
    k
    %determine the degree of overlap and fraction of sector for all lists of k
    [GO_lists_k,ind_lists_k,fractionOfSector_lists_k,degree_overlap_lists_k] = ... 
        EnumerateLists_k(GeneAssociation_sector,spectra_sector,k);
    
    tf1 = fractionOfSector_lists_k<0.6;
    tf2 = degree_overlap_lists_k>0.05;
    tf = tf1|tf2;
    fractionOfSector_lists_k(tf) = [];
    degree_overlap_lists_k(tf) = [];
    ind_lists_k(tf,:) = [];
    
    %remove GO terms not present in ind_k, and update ind_k
    GO_lists_k_temp = GO_lists_k(unique(ind_lists_k(:)));
    [tf,loc] = ismember(GO_lists_k,GO_lists_k_temp);
    ind_lists_k = loc(ind_lists_k);
    GO_lists_k = GO_lists_k_temp;
    
    %calculate gene coverage for each list
    GO_proteome = unique(GeneAssociation_complete.GO);
    [tf,GOid_proteome] = ismember(GeneAssociation_complete.GO,GO_proteome);
    [tf,GOid_sector] = ismember(GO_lists_k,GO_proteome);
    gene_coverage_lists_k = zeros(size(ind_lists_k,1),1);
    for i = 1:size(ind_lists_k,1)
        tf = false(length(GOid_proteome),1);
        for j = 1:k
            tf_temp = GOid_sector(ind_lists_k(i,j))== GOid_proteome;
            tf = tf|tf_temp;
        end
        %total number of unique genes associated with the list of GO terms
        n0 = numel(unique(GeneAssociation_complete.gene(tf)));
        n = numel(intersect(GeneAssociation_complete.gene(tf),genes_sector));
        gene_coverage_lists_k(i) = n/n0;
    end
    
    GO_lists = [GO_lists;{GO_lists_k}];
    ind_lists = [ind_lists;{ind_lists_k}];
    fractionOfSector_lists = [fractionOfSector_lists;{fractionOfSector_lists_k}];
    degree_overlap_lists = [degree_overlap_lists;{degree_overlap_lists_k}];
    gene_coverage_lists = [gene_coverage_lists;{gene_coverage_lists_k}];
    
    %keep the while loop going if nothing found so far
    if sum(cellfun('isempty',degree_overlap_lists))==k
        degree_overlap_lists_k = 0;
    end
end

%The main output of the script
lists.GO = GO_lists;
lists.ind = ind_lists;
lists.fractionOfSector = fractionOfSector_lists;
lists.degree_overlap = degree_overlap_lists;
lists.gene_coverage = gene_coverage_lists;

return;

