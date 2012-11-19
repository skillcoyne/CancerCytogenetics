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
      c = $1; b = $2
      band_by_chr[c] = [] unless band_by_chr.has_key? c
      band_by_chr[c].push(b)
      band_by_chr[c] = $3 if b.match(/([p|q]\d+)\.\d+/)
    end
    return band_by_chr
  end

end