require 'yaml'
require 'find'


def get_filehandle(filename, cols = ["Cancer", "Chr", "StartPosition", "EndPosition", "VarClass", "VarType", "NormalAllele1", "NormalAllele2", "TumorAllele1", "TumorAllele2"])
  f = File.open(filename, 'w')
  f.write(cols.join("\t") + "\n")
  return f
end


mutation_dir = "/Users/sarah.killcoyne/Data/TCGA/mutations"
maf_files = []
Find.find(mutation_dir) do |f|
  maf_files << f if f.match(/\.maf$/)
end

cancers = Dir["#{mutation_dir}/*"].map { |d| d.sub("#{mutation_dir}/", "") }
cancers.reject! { |c| c.match(/\.\w/) }

maf_by_cancer = {}
maf_files.each do |maf|
  cancer = "Unknown"
  cancers.each do |c|
    if maf.match(/#{c}/)
      (maf_by_cancer[c] ||= []) << maf
      break
    end
  end
end



fout = get_filehandle("#{Dir.home}/Data/TCGA/all_variants.txt")

maf_by_cancer.each_pair do |cancer, files|
  puts "Reading #{cancer} files..."
  cout = get_filehandle("#{Dir.home}/Data/TCGA/#{cancer}_variants.txt")

  files.each do |file|
    puts "\tReading #{file}..."
    File.open(file, 'r').each_with_index do |line, i|
      line.chomp!
      next if line.start_with? "#" or line.start_with? "Hugo"
      cols = line.split("\t")

      (build, chr, startpos, endpos, strand, variant_class, variant_type, ref_allele, tumor_a1, tumor_a2) = cols[3..11]
      (match_norm_a1, match_norm_a2) = cols[17..18]

      output = [cancer, chr, startpos, endpos, variant_class, variant_type, match_norm_a1, match_norm_a2, tumor_a1, tumor_a2]
      fout.write(output.join("\t") + "\n")
      cout.write(output.join("\t") + "\n")
    end
  end
  cout.close
end
fout.close