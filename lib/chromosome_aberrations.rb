require_relative 'aberration'
require_relative 'karyotype_error'

module ChromosomeAberrations

  ## INVERSION
  class Inversion < Aberration
    @kt = 'inv'
    @rx = /^inv\(/
  end

  ## DUPLICATION
  class Duplication < Aberration
    @kt = 'dup'
    @rx = /^dup\(/
  end

  ## INSERTION
  class Insertion < Aberration
    @kt = 'ins'
    @rx = /^ins\(/
  end

  ## DELETION
  class Deletion < Aberration
    @kt = 'del'
    @rx = /^del\(/
  end

  ## ADD (addition of unknown material)
  class Addition < Aberration
    @kt = 'add'
    @rx = /^add\(/
  end

  ## ISOCHROMOSOME
  class Isochromosome < Aberration
    @kt = 'iso'
    @rx = /^i\(/
  end

  ## DICENTRIC
  class DicentricChromosome < Aberration
    @kt = 'dic'
    @rx = /^dic\(/
  end

  ## RING
  class RingChromosome < Aberration
    @kt = 'ring'
    @rx = /^r\(/
  end

  ## ROBERTSONIAN
  class Robertsonian < Aberration
    @kt = 'rob'
    @rx = /^rob\(/
  end

  ## DERIVATIVE
  class Derivative < Aberration
    @kt = 'der'
    @rx = /^der\(/


    def get_breakpoints(str)
      @aberrations = []

      ab_objs = Aberration.aberration_objs

      chr_i = find_chr(str)
      derivative_abr = str[chr_i[:end_index]+1..str.length]

      # separate different abnormalities within the derivative chromosome and clean it up to make it parseable
      abnormalities = derivative_abr.scan(/([^\(\)]+\(([^\(\)]|\)\()*\))/).collect { |a| a[0] }

      abnormalities.each do |abn|
        abrclass = Aberration.classify_aberration(abn)
        ab_obj = ab_objs[abrclass].new(abn)
        @aberrations << ab_obj
        @breakpoints << ab_obj.breakpoints
      end
      #@breakpoints.flatten!
    end

  end

  ## TRANSLOCATION ... this is typically a subset of Derivative chromosomes, but have seen it on it's own
  class Translocation < Aberration
    @kt = 'trans'
    @rx = /^t\(/

    def get_breakpoints(str)
      chr_i = find_chr(str)
      band_i = find_bands(str, chr_i[:end_index])

      chr_i[:chr].each_with_index do |c,i|
        @breakpoints << Breakpoint.new(c, band_i[:bands][i], 'trans')
      end
      #@breakpoints.flatten!
    end

    :private
    # :last_chr, :last_band, :abr
    def recursive_trans_read(*args)

    end


  end

  ## FRAGMENT
  class Fragment < Aberration
    @kt = 'frag'
    @rx = /^frag\(/
  end

  ## CHROMOSOME GAIN
  class ChromosomeGain < Aberration
    @kt = 'gain'
    @rx = /^\+\d+|X|Y$/

    def initialize
      @abr = str.sub("+", "")
      @breakpoints = []
    end
  end

  ## CHROMOSOME LOSS
  class ChromosomeLoss < Aberration
    @kt = 'loss'
    @rx = /^-\d+|X|Y$/

    def initialize
      @abr = str.sub("-", "")
      @breakpoints = []
    end

  end


end