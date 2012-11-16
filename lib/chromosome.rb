#require_relative 'chromosome_aberrations'
#require_relative 'logging'
#require_relative 'breakpoint'
require_relative 'utils/band_reader'

class Chromosome
  #include ChromosomeAberrations
  #include BandReader
  #include Logging




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

    def aberration(obj)
      raise ArgumentError, "Not an Aberration object" unless obj.is_a? Aberration
      @aberrations << obj
    end

    def breakpoints
      bps = []
      @aberrations.each do |a|
        bps << a.breakpoints
        #a.breakpoints.each do |bp|
        #  bps << bp
        #  bps << bp.to_s
        #end
      end
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
      File.open("../resources/HsBands.txt", 'r').each_line do |line|
        line.chomp!
        if line.match(/^#{self.name}/)
          @normal_bands << line
        end
      end
    end


end