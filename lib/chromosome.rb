require_relative 'chromosome_aberrations'
require_relative 'logging'
require_relative 'utils/band_reader'

class Chromosome
  include ChromosomeAberrations
  include BandReader
  include Logging

  class<<self
    attr_accessor :normal_bands
  end

  attr_reader :name, :aberrations


  def initialize(*args)
    chr = args[0]
    raise ArgumentError, "#{chr} is not a valid chromosome identifier." unless (chr.is_a? String and chr.match(/^\d+|X|Y$/))
    @name = chr
    @aberrations = []
    ## TODO get bands
    load_bands() if (args.length > 1 and args[1].eql?true)
  end

  def aberration(obj)
    raise ArgumentError, "Not an Aberration object" unless obj.is_a? Aberration
    @aberrations << obj
  end

  def breakpoints
    bps = []
    @aberrations.each do |a|
      a.breakpoints.each do |bp|
        bps << bp.to_s
      end
    end
    return bps
  end


  class Fragment
    attr_reader :start, :end, :gene

    def initialize(*args)
      if args.size == 2
        #raise ArgumentError, "Arguments should be 'Band'." unless (args[0].kind_of?(Band) and args[1].kind_of? Band)
        @start = args[0]
        @end = args[1]
      #elsif args.size == 3
      #  @start = Band.new(args[0], args[1])
      #  @end = Band.new(args[0], args[2])
      #else
      #  raise ArgumentError, "Incorrect number of arguments, expected 2 (from, to) or 3 (chromosome, from, to)"
      end
    end

    def add_gene(gene)
      @gene = gene
    end

    def as_string
      return "#{@start.to_s} --> #{@end.to_s}"
    end
  end


  :private
  def load_bands # TODO quit hardcoding
    @normal_bands = []
    File.open("../resources/HsBands.txt", 'r').each_line do |line|
      line.chomp!
      if line.match(/^#{self.name}/)
        @normal_bands << line
      end
    end
  end

end