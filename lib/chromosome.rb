require_relative 'chromosome_aberrations'
require_relative 'utils/band_reader'

class Chromosome
  include ChromosomeAberrations
  include BandReader

  class<<self
    attr_accessor :normal_bands
  end

  attr_reader :name


  def initialize(*args)
    raise ArgumentError, "#{args}" unless args[:chr]
  end


end