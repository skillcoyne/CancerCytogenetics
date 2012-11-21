module BandReader
  @bands_by_chr = {}

  def bands(chr, file)
    file = File.open(file, 'r') unless file.is_a?File
    bands = read_file(file)
    bds = bands[chr]
    bds.uniq!
    return bds
  end

  def read_file(file)
    band_by_chr = {}
    file.each_line do |line|
      line.chomp!
      line.match(/^(\d+|X|Y)([p|q].*)/)
      c = $1; b = $2
      band_by_chr[c] = Array.new unless band_by_chr.has_key? c
      band_by_chr[c] << "#{c}#{b}"
      band_by_chr[c] << "#{c}#{$1}" if b.match(/([p|q]\d+)\.\d+/)
    end
    return band_by_chr
  end

end