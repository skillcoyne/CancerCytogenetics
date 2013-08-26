

file = "#{Dir.home}/Data/TSGene/Human_716_TSGs.txt"

genes = []


File.open(file, 'r').each_with_index do |line, i|
  next if i.eql?0
  line.chomp!

  (id, symbol, gene_alias, xref, chr, band, name, type, description, nuc_seq, prot_seq) = line.split("\t")

  #cols = line.split("\t")
  #puts cols
  ensembl_id = xref.match(/(ENSG\d+)\|/)
  ensembl_id = ensembl_id.captures.first unless ensembl_id.nil?

  genes << {'entrez_id' => id, 'official_symbol' => symbol, 'ensembl' => ensembl_id, 'chr' => chr } unless id.nil?

end

genes.sort_by!{|e| e['chr']}


File.open("#{Dir.home}/Data/TSGene/tumor_supressor_genes.txt", 'w') {|f|
  genes.map{|e|
    f.write [e['chr'], e['official_symbol'], e['entrez_id'], e['ensembl']].join("\t") + "\n"
  }
}

