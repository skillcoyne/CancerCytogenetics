require 'yaml'

require_relative 'logging'
require_relative 'genome_structure'
require_relative 'utils/karyotype_reader'

class Karyotype
  include GenomeStructure
  include KaryotypeReader
  include Logging

  Logging.configure(:class => 'Karyotype')  # doesn't work right now

  @@haploid = 23

  attr_reader :aberrations, :karyotype, :ploidy, :sex, :abnormal_chr, :normal_chr

  class<<self
    attr_accessor :aberration_objs, :unclear_aberrations
  end

  def initialize(str)
    raise ArgumentError, "#{str} is not a karyotype." unless (str.is_a? String and str.length > 1)
    log.info("Reading karyotype #{str}")

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
      regex = @aberration_obj[abr_type].regex

      @aberrations[abr_type].each do |abr|
#        if abr_type
        abr.match(regex)
        log.warn("Aberration has two chromosomes #{abr} but only the first one is handled.") unless ($2.nil? or $1.eql?$2 )

        ## TODO deal with the case of 2 chromosomes defined in the aberration
        chr = Chromosome.new($1, true)
        chr.aberration(@aberration_obj[abr_type].new(abr))

        @abnormal_chr[chr.name] = [] unless @abnormal_chr.has_key? chr.name
        @abnormal_chr[chr.name] << chr
      end
    end
  end


  # get breakpoints for the karyotype
  def report_breakpoints
    bps = Array.new
    @abnormal_chr.each_pair do |c, chr_list|
      chr_list.each do |chr|
        bps << chr.breakpoints
      end
    end
    bps.delete_if { |c| c.empty? }
    bps.flatten!
    return bps
  end

  #def report_ploidy_change
  #  ploidy_report = {}
  #  @normal_chr.each_pair {|chr, count| ploidy_report[chr] = count-2 if count > 2 }
  #  return ploidy_report
  #end


  def report_fragments
    frags = []
    @abnormal_chr.each_pair do |c, chr_list|
      chr_list.each do |chr|
        frags << chr.fragments
      end
    end
    frags.delete_if {|c| c.empty?}
    frags.flatten!
    return frags
  end

  def report_ploidy_change
    pd = []
    pd << @aberrations[:loss].map {|e| "-#{e}" } if @aberrations[:loss]
    pd << @aberrations[:gain].map {|e| "+#{e}" } if @aberrations[:gain]
    pd.flatten!
    return pd
  end



  def summarize
    summary = "NORMAL CHROMOSOMES\n"
    @normal_chr.each_pair do |chr, count|
      summary = "#{summary} #{chr}: #{count}\n"
    end

    summary = "#{summary}\nABNORMAL:"
    @abnormal_chr.each_pair do |chr, list|
      summary = "#{summary}\n#{chr}"
      list.each do |c|
        summary = "#{summary}\n#{c.aberrations}\n"
        summary = "#{summary}\n#{c.breakpoints}\n"
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
    @karyotype.gsub!(/\s/, "")
    clones = @karyotype.scan(/(\[\d+\])/).collect { |a| a[0] }
    log.warn("Karyotype is a collection of clones, analysis may be inaccurate.") if clones.length > 3

    @karyotype.gsub!(/\[\d+\]/, "") # don't care about numbers of cells: [5]

    (pl, sc) = @karyotype.split(",")[0..1]
    @ploidy = KaryotypeReader.calculate_ploidy(pl, @@haploid)
    sex_chr = KaryotypeReader.determine_sex(sc)

    st = sex_chr.values.inject {|sum,v| sum+v }
    @sex = nil
    karyotype_index = 1 # sometimes the sex is not indicated and there's no case information to figure it out
    if st > 0
      @sex = sex_chr.keys.join("")
      karyotype_index = 2
    end

    (Array(1..23)).each { |c| @normal_chr[c.to_s] = @ploidy.to_i }

    sex_chr.each_pair { |c, p| @normal_chr[c] = p.to_i }

    # deal with the most common karyotype string inconsistencies
    cleaned_karyotype = []

    @karyotype.split(",")[karyotype_index..-1].each do |abr|
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