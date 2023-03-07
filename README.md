# OrthoPhy

 A program to construct ortholog datasets using taxonomic information

## Trouble shooting

The program IslandPath-DIMOB, which is used to predict horizontally transferred genes, requires bioperl. However, if you install it using the ubuntu apt command, you may get the such error, "Bio::Perl not installed". In this case, replace the line "use Bio::Perl;" with "use BioPerl;" in the GenomeUtils.pm file in IslandPath-DIMOB.
