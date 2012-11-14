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
  end

  ## TRANSLOCATION ... this is typically a subset of Derivative chromosomes, but have seen it on it's own
  class Translocation < Aberration
    @kt = 'trans'
    @rx = /^t\(/
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
  end

  ## CHROMOSOME LOSS
  class ChromosomeLoss < Aberration
    @kt = 'loss'
    @rx = /^-\d+|X|Y$/
  end


end