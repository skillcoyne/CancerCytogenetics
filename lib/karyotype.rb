require_relative 'utils/karyotype_reader'
require_relative 'aberration'
require_relative 'chromosome'
require_relative 'chromosome_aberrations'
require_relative 'breakpoint'
require_relative 'logging'

require 'yaml'

class Karyotype
  include KaryotypeReader
  include Logging
#  Logging.configure(:class => 'Karyotype')  # doesn't work right now

  @@haploid = 23

  attr_reader :aberrations, :karyotype, :ploidy, :sex

  class<<self
    attr_accessor :normal_chr, :abnormal_chr, :aberration_objs, :unclear_aberrations
  end

  def initialize(str)
    raise ArgumentError, "#{str} is not a karyotype." unless str.is_a? String

    @karyotype = str.gsub(/\s/, "")
    @normal_chr = {}; @abnormal_chr = {}; @aberrations = {}; @unclear_aberrations = [];
    setup_abberation_objs()
    prep_karyotype()
    handle_ploidy_diff()
    analyze()
  end


  def analyze
    Aberration.aberration_type.each do |abr_type|
      next unless @aberrations.has_key? abr_type
      @aberrations[abr_type].each do |abr|
        chr_i = KaryotypeReader.find_chr(abr)
        #band_i = KaryotypeReader.find_bands(abr, chr_i[:end_index])

        # just a check
        #unless chr_i[:chr].length.eql? band_i[:bands].length
        #  log.error("Bands and chromosomes don't match up: #{abr}. Skipping.")
        #  next
        #end

        ## TODO deal with the case of 2 chromosomes defined in the aberration
        chr = Chromosome.new(chr_i[:chr][0])
        chr.aberration(@aberration_obj[abr_type].new(abr))

        @abnormal_chr[chr.name] = [] unless @abnormal_chr.has_key? chr.name
        @abnormal_chr[chr.name] << chr
      end
    end
  end


  # get breakpoints for the karyotype
  def report_breakpoints
    bps = []
    @abnormal_chr.each_pair do |c, chr_list|
      chr_list.each do |chr|
        bps << chr.breakpoints
      end
    end
    return bps.flatten!
  end

  def report_ploidy_change
    pd = []
    pd << @aberrations[:loss].map {|e| "-#{e}" } if @aberrations[:loss]
    pd << @aberrations[:gain].map {|e| "+#{e}" } if @aberrations[:gain]
    return pd.flatten!
  end



  def summarize
    puts "NORMAL CHROMOSOMES:"
    @normal_chr.each_pair do |chr, count|
      puts "#{chr}: #{count}"
    end

    puts "ABNORMAL:"
    @abnormal_chr.each_pair do |chr, list|
      puts "#{chr}"
      list.each do |c|
        puts c.aberrations
        puts c.breakpoints
      end
    end
  end


  :private

  def setup_abberation_objs
    @aberration_obj = Aberration.aberration_objs
  end


  def handle_ploidy_diff
    @aberrations[:loss].each { |c| @normal_chr[c] -= 1 } if @aberrations[:loss]
    @aberrations[:gain].each { |c| @normal_chr[c] += 1 } if @aberrations[:gain]
  end

# determine ploidy & gender, clean up each aberration and drop any "unknown"
  def prep_karyotype
    clones = @karyotype.scan(/(\[\d+\])/).collect { |a| a[0] }
    log.warn("Karyotype is a collection of clones, analysis may be inaccurate. #{@karyotype}") if clones.length > 3

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
      cleaned_karyotype |= [cleaned_karyotype, KaryotypeReader.cleanup(abr)].flatten
    end
    @karyotype = cleaned_karyotype

    # classify each type of aberration in the karyotype
    @karyotype.each do |k|
      abrclass = Aberration.classify_aberration(k)
      @aberrations[abrclass] = [] unless @aberrations.has_key? abrclass
      @aberrations[abrclass] << k.sub(/^(\+|-)?/, "")
    end


    @aberrations.each_pair do |abrclass, abrlist|
      next if (abrclass.eql? ChromosomeAberrations::ChromosomeGain.type or abrclass.eql? ChromosomeAberrations::ChromosomeLoss.type)
      # aberrations other than chromosome gains/losses should be uniquely represented

      counts = abrlist.inject(Hash.new(0)) { |h, i| h[i] += 1; h }
      counts.each_pair { |k, v| log.warn("#{k} was seen multiple times. Analyzed only once.") if v > 1 }

      @aberrations[abrclass] = abrlist.uniq
    end

  end


end