require_relative 'aberration'

module ChromosomeAberrations

  ## INVERSION
  class Inversion < Aberration
    @kt = 'inv'
    @rx = /^inv\((\d+|X|Y)\)/

  end

  ## DUPLICATION
  class Duplication < Aberration
    @kt = 'dup'
    @rx = /^dup\((\d+|X|Y)\)/
  end

  ## INSERTION
  class Insertion < Aberration
    @kt = 'ins'
    @rx = /^ins\((\d+|X|Y)\)/
  end

  ## DELETION
  class Deletion < Aberration
    @kt = 'del'
    @rx = /^del\((\d+|X|Y)\)/
  end

  ## ADD (addition of unknown material)
  class Addition < Aberration
    @kt = 'add'
    @rx = /^add\((\d+|X|Y)\)/
  end

  ## ISOCHROMOSOME
  class Isochromosome < Aberration
    @kt = 'iso'
    @rx = /^i\((\d+|X|Y)\)/
  end

  ## DICENTRIC
  class DicentricChromosome < Aberration
    @kt = 'dic'
    @rx = /^dic\((\d+|X|Y)[;|:](\d+|X|Y)\)/

    #def get_breakpoints
    #  chr_i = find_chr(@abr)
    #  band_i = find_bands(@abr, chr_i[:end_index])
    #  chr_i[:chr].each_with_index do |c, i|
    #    @breakpoints << Breakpoint.new(c, band_i[:bands][i], 'dic')
    #  end
    #  # TODO am not sure how the dic rearrangment works, see this in CyDas dic(13;13)(q14;q32)
    #  #@fragments << Fragment.new( Breakpoint.new(@breakpoints[0].chr, "pter"), @breakpoints[0])
    #  #@fragments << Fragment.new( @breakpoints[1], Breakpoint.new(@breakpoints[1].chr, "#{@breakpoints[1].arm}ter"))
    #end

  end

  ## RING  ## TODO figure out the right regex for this
  #class RingChromosome < Aberration
  #  @kt = 'ring'
  #  @rx = /^r\(/
  #end

  ## ROBERTSONIAN
  #class Robertsonian < Aberration
  #  @kt = 'rob'
  #  @rx = /^rob\(/
  #end

  ## DERIVATIVE
  class Derivative < Aberration
    @kt = 'der'
    @rx = /^der\((\d+|X|Y)\)/

    def get_breakpoints
      @aberrations = []

      ab_objs = Aberration.aberration_objs

      chr_i = find_chr(@abr)
      derivative_abr = @abr[chr_i[:end_index]+1..@abr.length]

      # separate different abnormalities within the derivative chromosome and clean it up to make it parseable
      abnormalities = derivative_abr.scan(/([^\(\)]+\(([^\(\)]|\)\()*\))/).collect { |a| a[0] }

      trans_bps = []
      abnormalities.each do |abn|
        abrclass = Aberration.classify_aberration(abn)

        if abrclass.to_s.eql? 'unk' # not dealing with unknowns
          log.warn("Cannot handle #{abn}, incorrect format.")
          next
        end

        # special handling because translocations are written as a sliding window
        # translocations should also only every have 2 breakpoints...
        if abrclass.to_s.eql? ChromosomeAberrations::Translocation.type
          trans = ChromosomeAberrations::Translocation.new(abn)
          trans_bps << trans.breakpoints
          @breakpoints << trans.breakpoints
        else
          ab_obj = ab_objs[abrclass].new(abn)
          if ab_obj.breakpoints.length > 0
            @aberrations << ab_obj
            @breakpoints << ab_obj.breakpoints
          end
        end
      end
      trans_bps.delete_if {|c| c.empty? }
      add_fragments(trans_bps.flatten!) if trans_bps.length > 0
    end

    :private
    # have to reorder the array and then turn Breakpoints into fragments
    def add_fragments(tbp_list)
      sorted = []
      tbp_list.each_with_index do |e, i|
        if i <= 1
          sorted << Breakpoint.new(e.chr, "#{e.arm}ter") if i.eql? 0
          sorted << e
        elsif i%2 == 0
          sorted << tbp_list[i+1]
          sorted << tbp_list[i]
        end
      end
      sorted << Breakpoint.new(sorted[-1].chr, "#{sorted[-1].arm}ter")
      sorted.each_slice(2).to_a.each do |pair|
        @fragments << Fragment.new(pair[0], pair[1])
      end
    end

  end

  ## TRANSLOCATION ... this is typically a subset of Derivative chromosomes, but have seen it on it's own
  class Translocation < Aberration
    @kt = 'trans'
    @rx = /^t\((\d+|X|Y)[;|:](\d+|X|Y)\)/


    ## TWo ways of defining translocations:
    ## 1)  t(1;3)(p31;p13)
    def get_breakpoints
      chr_i = find_chr(@abr)
      band_i = find_bands(@abr, chr_i[:end_index])
      unless band_i
        log.warn("No bands defined in #{@abr}")
      else
        chr_i[:chr].each_with_index do |c, i|
          @breakpoints << Breakpoint.new(c, band_i[:bands][i], 'trans')
        end
      end
    end

  end

  ## FRAGMENT
  class ChromosomeFragment < Aberration
    @kt = 'frag'
    @rx = /^frag\((\d+|X|Y)\)/
  end

  ## CHROMOSOME GAIN
  class ChromosomeGain < Aberration
    @kt = 'gain'
    @rx = /^\+(\d+|X|Y)$/

    def initialize
      Logging.configure(:class => self.class.name)
      @abr = str.sub("+", "")
      @breakpoints = []
    end
  end

  ## CHROMOSOME LOSS
  class ChromosomeLoss < Aberration
    @kt = 'loss'
    @rx = /^-(\d+|X|Y)$/

    def initialize
      Logging.configure(:class => self.class.name)
      @abr = str.sub("-", "")
      @breakpoints = []
    end

  end


end