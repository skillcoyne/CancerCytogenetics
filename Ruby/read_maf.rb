require 'yaml'
require 'find'


def get_filehandle(filename,
    cols = ["Center", "Patient", "Cancer", "Chr", "Start", "End", "VarClass", "VarType"])
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

fh = Hash[Array(1..22).map { |e| [e.to_s, get_filehandle("#{Dir.home}/Data/TCGA/chr#{e}_variants.txt")] }]
fh['X'] = get_filehandle("#{Dir.home}/Data/TCGA/chrX_variants.txt")
fh['Y'] = get_filehandle("#{Dir.home}/Data/TCGA/chrY_variants.txt")

maf_by_cancer.each_pair do |cancer, files|
  puts "Reading #{cancer} files..."
  files.each do |file|
    # WU made calls in regions of the genome that they could not have. There are previously known issues with their sequencing quality so just filtering them out
    next if file.match(/genome.wustl.edu/)
    puts "\tReading #{file}..."

    File.open(file, 'r').each_with_index do |line, i|
      line.chomp!
      next if line.start_with? "#" or line.start_with? "Hugo" or line.eql? ""
      #cols = line.split("\t")
      #cols.each_with_index do |c, i|
      #  puts "#{i} #{c}"
      #end
      cols = line.split("\t")

      cols.each { |e| raise "ERROR: #{cols.join(' ')}" if e.match(/\t/) }

      (center, build, chr, startpos, endpos, strand, variant_class, variant_type) = cols[2..9]
      patient = cols[15].split("-")[2]

      output = [center, patient, cancer, chr, startpos, endpos, variant_class, variant_type]

      if fh.has_key? chr
        fh[chr].write(output.join("\t") + "\n")
      end
    end
  end
end
fh.values.each { |e| e.close }