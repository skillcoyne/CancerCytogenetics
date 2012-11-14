require_relative 'utils/karyotype_reader'
require_relative 'aberration'
require_relative 'chromosome'
require_relative 'chromosome_aberrations'
require_relative 'breakpoint'

require 'yaml'

class Karyotype
  include KaryotypeReader

  @@haploid = 23

  attr_reader :aberrations, :karyotype, :breakpoints, :ploidy, :sex

  class<<self
    attr_accessor :normal_chr, :abnormal_chr, :aberration_objs
  end

  def initialize(str)
    @karyotype = str.gsub(/\s/, "")
    @normal_chr = {}; @abnormal_chr = {}; @aberrations = {}; @breakpoints = {}
    setup_abberation_objs()
    prep_karyotype()
    find_breaks()
  end


  def analyze
    Aberration.breakpoint_regex.each do |abr_type|
      #puts abr_type
      next unless @aberrations.has_key? abr_type
      @aberrations[abr_type].each do |abr|
        chr_i = find_chr(abr)
        band_i = find_bands(abr, chr_i[:end_index])

        # just a check
        raise KaryotypeError, "Bands and chromosomes don't match up: #{abr}" unless chr_i[:chr].length.eql? band_i[:bands].length

        ## TODO deal with the case of 2 chromosomes defined in the aberration
        chr = Chromosome.new(chr_i[:chr][0])
        chr.aberration( @aberration_obj[abr_type].new(abr) )

        @abnormal_chr[chr.name] = [] unless @abnormal_chr.has_key?chr.name
        @abnormal_chr[chr.name].push(chr)
      end

    end

  end


  :private

  def setup_abberation_objs
    @aberration_obj = {}
    ChromosomeAberrations.constants.each do |ca|
      abr_obj = ChromosomeAberrations.const_get(ca)
      @aberration_obj[abr_obj.type.to_sym] = abr_obj
    end
  end


  # determine ploidy & gender, clean up each aberration and drop any "unknown"
  def prep_karyotype
    clones = @karyotype.scan(/(\[\d+\])/).collect { |a| a[0] }
    raise KaryotypeError, "Karyotype is an collection of clones, needs to be curated first. #{@karyotype} " if clones.length > 3

    @karyotype.gsub!(/\[\d+\]/, "") # don't care about numbers of cells: [5]

    (pl, sc) = @karyotype.split(",")[0..1]
    @ploidy = KaryotypeReader.calculate_ploidy(pl, @@haploid)
    sex_chr = KaryotypeReader.determine_sex(sc)
    @sex = sex_chr.keys.join("")

    (Array(1..23)).each { |c| @normal_chr[c.to_s] = @ploidy }
    sex_chr.each_pair { |c, p| @normal_chr[c] = p }

    # deal with the most common karyotype string inconsistencies
    cleaned_karyotype = []
    @karyotype.split(",")[2..-1].each do |abr|
      cleaned_karyotype = [cleaned_karyotype, KaryotypeReader.cleanup(abr)].flatten
    end
    @karyotype = cleaned_karyotype

    # classify each type of aberration in the karyotype
    @karyotype.each do |k|
      abrclass = Aberration.classify_aberration(k)
      @aberrations[abrclass] = [] unless @aberrations.has_key? abrclass
      @aberrations[abrclass].push(k)
    end
  end

  # pull out breakpoints from an aberration
  def find_breaks
    Aberration.breakpoint_regex.each do |abr_type|
      next unless @aberrations.has_key? abr_type
      @aberrations[abr_type].each do |abr|
        chrs = find_chr(abr)
        bands = find_bands(abr, chrs[:end_index])
        # just a check
        raise KaryotypeError, "Bands and chromosomes don't match up: #{abr}" unless chrs[:chr].length.eql? bands[:bands].length
        chrs[:chr].each_with_index do |c, i|
          bp = Breakpoint.new(c, bands[:bands][i], abr_type.to_s)
          @breakpoints[bp] = 0 unless @breakpoints.has_key? bp
          @breakpoints[bp] += 1
        end
      end
    end
  end

  # find chromosome in aberration strings ex. der(19)t(19:2)(q10;q32)
  def find_chr(str)
    chr_s = str.index(/\(/, 0)
    chr_e = str.index(/\)/, chr_s)
    chr = str[chr_s+1..chr_e-1]
    raise KaryotypeError, "No chromosome parsed from #{str}." unless chr.match(/\d+|X|Y/)
    return {:start_index => chr_s, :end_index => chr_e, :chr => chr.split(/;|:/)}
  end

  # find bands in aberration strings ex. der(19)t(19:2)(q10;q32)
  def find_bands(str, index)
    raise KaryotypeError, "No bands defined in #{str}" if str.length.eql?(index+1)
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
  def find_fragments(str)
    return str.scan(/[p|q][\d+][\.\d]/)
  end

end