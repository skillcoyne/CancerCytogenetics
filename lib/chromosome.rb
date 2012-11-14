require_relative 'chromosome_aberrations'
require_relative 'utils/band_reader'

class Chromosome
  include ChromosomeAberrations
  include BandReader

  class<<self
    attr_accessor :normal_bands, :aberrations
  end

  attr_reader :name


  def initialize(chr)
    raise ArgumentError, "#{chr} is not a valid chromosome identifier." unless (chr.is_a?String and chr.match(/^\d+|X|Y$/))
    @name = chr
    @aberrations = []
    ## TODO get bands
  end

  def aberration(obj)
    raise ArgumentError, "Not an Aberration object" unless obj.is_a?Aberration
    @aberrations.push(obj)
  end


end