require_relative '../karyotype_error'
require_relative '../aberration'
require_relative '../../lib/logging'

module KaryotypeReader
  include Logging

  def self.cleanup(abr)
    new_abr = []

    # +t(13;X)(q13;p12) doesn't need a +
    abr.sub!(/^[\+|-]/, "") unless abr.match(/^[\+|-][\d|X|Y]+$/)

    # not going to bother with aberrations that are unclear/unknown '?' or with '**'
    if (abr.match(/\?|\*\*/))
      log.warn("Removing aberration with unknown/unclear information: #{abr}")
      return new_abr
    end


    # 13x2 is normal, 13x3 is a duplicate and should read +13
    if abr.match(/^([\d+|X|Y]+)x(\d+)/)
      chr = $1; dups = $2.to_i
      if dups.eql? 0 # managed to lose both chromosomes in a diploidy karyotype
        (Array(1..dups)).map { new_abr.push("-#{chr}") }
      elsif dups > 2 # sometimes you have 13x3, really just means 1 additional chr 13 since normal ploidy is 2
        dups -= 2
        (Array(1..dups)).map { new_abr.push("+#{chr}") }
      elsif dups.eql?(1)
        new_abr.push("-#{chr}")
      end
      # add(9)(p21)x2 or add(7)x2 should indicate that this "additional material of unk origin" happened twice
    elsif abr.match(/(.*)x(\d+)$/)
      a = $1; dups = $2.to_i
      (Array(1..dups)).map { new_abr.push(a) }
      # del(7) should be -7  but not del(7)(q12)
    else # everything else
      new_abr.push(abr)
    end

    return new_abr
  end

  def self.determine_sex(str)

    raise KaryotypeError, "Definition of gender incorrect (#{str})" unless str.match(/X|Y/)

    sex_chr = {}
    # ploidy number makes no difference since this string will tell us how many or at least what the gender should be
    ['X', 'Y'].each { |c| sex_chr[c] = 0 }

    chrs = str.match(/([X|Y]+)/).to_s.split(//)
    chrs.each { |c| sex_chr[c] +=1 }


    # assume this was an XY karyotype that may have lost the Y, have only seen this in
    # severely affected karyotypes
    sex_chr['Y'] += 1 if (chrs.length.eql?(1) and chrs[0].eql?('X'))

    return sex_chr
  end

  def self.calculate_ploidy(str, haploid)
    diploid = haploid*2
    triploid = haploid*3
    quadraploid = haploid*4

    # typically see di- tri- quad- if more than that it should be noted
    ploidy = nil
    min = diploid
    max = diploid
    if str.match(/<\+(\d)n>/) # sometimes see this odd configuration: 46<+3n>
      ploidy = $1
    elsif str.match(/(\d+)-(\d+)/) # num and range or just range: 46-53
      min = $1.to_i; max = $2.to_i
    elsif str.match(/^(\d+)/) # single num:  72
      min = $1.to_i; max = $1.to_i
    end

    if ploidy.nil?
      case
        when (min.eql? diploid and max.eql? diploid)
          log.info("Normal ploidy: #{str}")
          ploidy = 2
        when ((min >= haploid and max <= diploid) or (min <= diploid and max < triploid))
          log.info("Relatively normal ploidy #{str}")
          ploidy = 2
        when (min >= haploid and max < quadraploid)
          log.info("Triploid #{str}")
          ploidy = 3
        when (min >= diploid and max >= quadraploid)
          log.info("Quadraploid #{str}")
          ploidy = 4
        else
          raise KaryotypeError, "Failed to determine ploidy for #{str}"
      end
    end
    return ploidy
  end

  # find chromosome in aberration strings ex. der(19)t(19:2)(q10;q32)
  def self.find_chr(str)
    chr_s = str.index(/\(/, 0)
    chr_e = str.index(/\)/, chr_s)
    chr = str[chr_s+1..chr_e-1]
    raise KaryotypeError, "No chromosome parsed from #{str}." unless chr.match(/\d+|X|Y/)
    return {:start_index => chr_s, :end_index => chr_e, :chr => chr.split(/;|:/)}
  end

  # find bands in aberration strings ex. der(19)t(19:2)(q10;q32)
  def self.find_bands(str, index)
    if str.length.eql?(index+1)
      log.warn("No bands defined in #{str}")
      return {:start_index => nil, :end_index => nil, :bands => []}
    end

    ei = str.index(/\(/, index)
    if str.match(/(q|p)(\d+|\?)/) and str[ei-1..ei].eql?(")(") # has bands and is not a translocation
      band_s = str.index(/\(/, index)
      band_e = str.index(/\)/, band_s)
      band_e = str.length-1 if band_e.nil?

      bands = str[band_s+1..band_e-1]

      return {:start_index => band_s, :end_index => band_e, :bands => bands.split(/:|;/)}
    end
  end

  # sometimes bands are defined for a single chr as p13q22
  def self.find_fragments(str)
    return str.scan(/[p|q][\d+][\.\d]/)
  end

end