require 'yaml'
require_relative 'lib/entrez_gene_info'


file = dir = "#{Dir.home}/Data/Homo_sapiens.gene_info.20121231"

out = File.new("#{Dir.home}/Data/sky-cgh/Hs_gene_band.txt", 'w')
cols = ['EntrezGeneId', 'OfficialSymbol', 'EnsembleId', 'Chromosome', 'Band']
out.write( cols.join("\t") + "\n" )
File.open(file, 'r').each_with_index do |line, index|
  next if line.start_with?"#"
  gi = EntrezGeneInfo.parse_line(line)
  next unless (gi.band and gi.tax_id.eql?'9606')
  next if (gi.chr.eql?'MT')# or !gi.type.eql?'protein-coding')

  puts index

  cols = [gi.gene_id, gi.official_symbol, gi.alt_ids['Ensembl'], gi.chr, "#{gi.chr}#{gi.band}"]
  out.write( cols.join("\t") + "\n")

end
out.close