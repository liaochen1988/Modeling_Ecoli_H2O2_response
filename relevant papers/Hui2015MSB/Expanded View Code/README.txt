GOAnalysis.m
  The main Matlab script that calls various functions below. It is the only script that needs to be run from the command line. One needs to specify the proteome sector for which the GO analysis is done.

data.mat
  Matlab data file that the 'GOAnalysis.m' loads at the beginning. The data file includes the following Matlab variables.
  It includes the following variables:

  GeneAssociation_complete:   a structure specifying gene-to-GO relation

  GeneOntology_relation:  a structure specifying parent-child relation
  between GO terms

  spectra_proteome:   list of spectra for the whole proteome, specifying
  abundance for each protein. Each spectrum is represented as a gene
  name. The abundance for a protein is calculated by the number of
  occurence of its gene name divided by the size of the list.

  sectorC:    list of genes belonging to C-sector
  sectorA:    list of genes belonging to A-sector
  sectorR:    list of genes belonging to R-sector
  sectorU:    list of genes belonging to U-sector
  sectorS:    list of genes belonging to S-sector
  sectorO:    list of genes belonging to O-sector

CalculateFractionsForGOTerms.m
  A function called by the script 'GOAnalysis.m' to calculate fraction of sector and fraction of proteome for each of given GO terms.

FindGOParent.m
  A function called by the script 'GoAnalsis.m' to find parent-child GO term pairs which have the same values of fraction of sector.

EnumerateLists_k.m
  A function called by the script 'GOAnalysis.m' to numerate all combinations of lists with k number of GO terms. It also calculates the fraction of sector and degree of overlap for each of the lists.

combinator
  A folder containing the function 'combinator.m' (and its sub-function) that is called by the function 'NumerateLists_k.m'. The functions are third party and a distribution license is included.

unique2.m
  A routine that returns the unique entries and the corresponding number of occurrence in the original array.