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
    ## TODO do something with the bands
    load_bands() if (args.length > 1 and args[1].eql? true)
  end

  def to_s
    "#{@name}"
  end

  def aberration(obj)
    raise ArgumentError, "Not an Aberration object" unless obj.is_a? Aberration

    obj.breakpoints.delete_if {|bp|
      @normal_bands.index(bp.to_s).nil?
      log.warn("Band #{bp.to_s} doesn't exist. Removing.")
    }
    @aberrations << obj
  end

  def breakpoints
    bps = []
    @aberrations.each {|a| bps << a.breakpoints }
    return bps
  end

  def fragments
    frags = []
    @aberrations.each do |a|
      frags << a.fragments
    end
    frags
  end


  :private

  def load_bands # TODO quit hardcoding



    @normal_bands = []
    File.open("/Users/sarah.killcoyne/workspace/CancerCytogenetics/resources/HsBands.txt", 'r').each_line do |line|
      line.chomp!
      if line.match(/^#{self.name}/)
        @normal_bands << line
      end
    end
  end


end