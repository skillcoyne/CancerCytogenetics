require 'fileutils'
require 'biomart'
require 'yaml'

def process_results(filters)
  results = $hsgene.search(
      :filters => filters,
      :attributes => [ 'ensembl_gene_id','external_gene_id'],
      :process_results => true
  )

  unless results.nil?
    results.uniq!

    return results
    #puts results.map{|e| e['external_gene_id']}.sort
  end

end


biomart = Biomart::Server.new("http://www.ensembl.org/biomart")
$hsgene = biomart.datasets['hsapiens_gene_ensembl']


dir = "#{Dir.home}/Data/TCGA"

files = Dir.glob("#{dir}/chr*_variants.txt")


files.each do |file|
  regions = []

  puts file

  current_chr = nil
  mut_genes = []
  File.open(file, 'r').each_with_index do |line, i|
    next if i.eql?0
    line.chomp!
    (center, patient, cancer, chr, start, stop, varclass, type, ref, tumor1, tumor2) = line.split("\t")

    regions << "#{chr}:#{start}-#{chr}:#{stop}"
    if regions.length >= 400
		 	print "." 
      filters = {
          'chromosomal_region' => regions,
          'biotype' => ["protein_coding"],
          'status' => ["KNOWN"]
      }
      mut_genes << process_results(filters)
      regions.clear
    end
    current_chr = chr
  end

  mut_genes.flatten!
  mut_genes.uniq!

  puts "#{current_chr} mutated genes #{mut_genes.size}"

  mut_genes.map!{|e| [e['ensembl_gene_id'], e['external_gene_id']] }

  filters = {
      "chromosome_name" => [current_chr],
      'transcript_status' => ["KNOWN"],
      'biotype' => ["protein_coding"],
      'status' => ["KNOWN"]
  }
  all_genes = process_results(filters)
  all_genes.map!{|e| [e['ensembl_gene_id'], e['external_gene_id']] }

  unmod_genes = all_genes - mut_genes

  puts "#{current_chr} All genes: #{all_genes.length}, cancer genes: #{mut_genes.length}, unmodified genes: #{unmod_genes.length}"

  File.open("#{Dir.home}/Data/CancerUnmodified/chr#{current_chr}-unmod.txt", 'w') {|f|
    f.write "## #{current_chr} All genes: #{all_genes.length}, cancer genes: #{mut_genes.length}, unmodified genes: #{unmod_genes.length} ##"
    f.write unmod_genes.sort_by{|e| e[1] }.map{|e| e.join("\t")}.join("\n")
  }

  sleep(3)
end

puts "Done."

