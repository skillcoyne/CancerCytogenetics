require_relative 'chromosome'
require_relative 'chromosome_aberrations'
require_relative 'breakpoint'
require_relative 'fragment'
require_relative 'aberration'
require_relative 'genome_structure_error'

require_relative 'logging'

module GenomeStructure
  include Logging

  Chromosome = ::Chromosome
  Aberration = ::Aberration
  ChromosomeAberrations = ::ChromosomeAberrations
  Breakpoint = ::Breakpoint
  Fragment = ::Fragment

end