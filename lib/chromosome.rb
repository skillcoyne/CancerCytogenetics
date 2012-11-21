require_relative 'utils/band_reader'

class Chromosome
  include BandReader


  class<<self
    attr_accessor :normal_bands
  end

  attr_reader :name, :aberrations


  def initialize(*args)
    Logging.configure(:class => self.class.name)

    chr = args[0]
    raise ArgumentError, "#{chr} is not a valid chromosome identifier." unless (chr.is_a? String and chr.match(/^\d+|X|Y$/))
    @name = chr
    @aberrations = []
    @normal_bands = bands(@name, "/Users/sarah.killcoyne/workspace/CancerCytogenetics/resources/HsBands.txt") if (args.length > 1 and args[1].eql? true) ## TODO quit hardcoding
  end

  def to_s
    "#{@name}"
  end

  def aberration(obj)
    raise ArgumentError, "Not an Aberration object" unless obj.is_a? Aberration

    #obj.breakpoints.each do |bp|
    #  log.warn("Band #{bp.to_s} doesn't exist. Removing.") if @normal_bands.index(bp.to_s).nil?
    #end

    ## TODO Deal with bands, HOWEVER because the chromosome has aberration objects breakpoints can include
    ## bands for which no chromosome object is created

    #obj.breakpoints.reject {|bp|
    #  @normal_bands.index(bp.to_s).nil?
    #}

    @aberrations << obj
  end

  def breakpoints
    bps = []
    @aberrations.each { |a| bps << a.breakpoints }
    return bps
  end

  def fragments
    frags = []
    @aberrations.each do |a|
      frags << a.fragments
    end
    frags
  end

end