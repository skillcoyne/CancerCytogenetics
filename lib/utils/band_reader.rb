module BandReader


  def self.bands(chr, file)
    bands = self.read_file(file)
    return bands[chr]
  end

  def self.read_file(file)
    band_by_chr = {}
    file.each_line do |line|
      line.chomp!
      line.match(/^(\d+|X|Y)([p|q].*)/)
      band_by_chr[$1] = [] unless band_by_chr.has_key? $1
      band_by_chr[$1].push($2)
    end
    return band_by_chr
  end

end