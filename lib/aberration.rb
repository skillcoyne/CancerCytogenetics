require 'yaml'

class Aberration


  attr_reader :abr, :ab_objs


  class<<self
    @kt = 'unk'
    @rx = /\?/

    def instantiate_aberrations
      aberration_obj = {}
      ChromosomeAberrations.constants.each do |ca|
        abr_obj = ChromosomeAberrations.const_get(ca)
        aberration_obj[abr_obj.type.to_sym] = abr_obj
      end
      return aberration_obj
    end

  end

  def self.type
    return @kt
  end

  def self.regex
    return @rx
  end

  def self.all_regex
    rx = {}
    ChromosomeAberrations.constants.each do |ca|
      ca_obj = ChromosomeAberrations.const_get(ca)
      rx[ca_obj.type.to_sym] = ca_obj.regex
    end
    return rx
  end

  # instantiate these
  def self.aberration_objs
    @ab_objs ||= self.instantiate_aberrations
  end

  def self.aberration_type
    abr_breaks = Aberration.all_regex.keys
    abr_breaks.delete_if { |c| c.to_s.match(/gain|loss/) }
    return abr_breaks
  end

  def self.classify_aberration(abr)
    Aberration.all_regex.each_pair do |k, regex|
      return k if abr.match(regex)
    end
    return "unknown".to_sym
  end

  def initialize(str)
    @abr = str
    @breakpoints = []

    #regex = Aberration.regex[@type.to_sym]
    # make sure it really is an inversion first
    #raise KaryotypeError, "#{str} does not appear to be a #{self.class}" unless str.match(self.regex)

    breakpoints()
  end

  def breakpoints
    @breakpoints += get_breakpoints(@abr)
  end


  :private


  def get_breakpoints(str)
    bps = []
    chr_i = find_chr(str)
    #raise KaryotypeError, "Too many chromosomes #{str}" if chr_i[:chr].length > 1

    band_i = find_bands(str, chr_i[:end_index])

    chr_i[:chr].each_with_index do |c, i|
      if band_i
        b = band_i[:bands][i]
        fragments = find_fragments(b)
        fragments.each { |f| bps.push(Breakpoint.new(c, f, @type)) }
      else
        bps.push(Breakpoint.new(c, "", @type))
      end
    end
    return bps
  end

  # Parsing aberration strings to pull out the chromosome and band definitions
  # These will result in breakpoint information
  def find_chr(str)
    chr_s = str.index(/\(/, 0)
    chr_e = str.index(/\)/, chr_s)
    chr = str[chr_s+1..chr_e-1]
    raise KaryotypeError, "No chromosome parsed from #{str}." unless chr.match(/\d+|X|Y/)
    return {:start_index => chr_s, :end_index => chr_e, :chr => chr.split(/;|:/)}
  end

  def find_bands(str, index)
    #raise KaryotypeError, "No bands defined in #{str}" if str.length.eql?(index+1)
    if str.length.eql?(index+1)
      warn("No bands defined in #{str}")
      return
    end

    ei = str.index(/\(/, index)
    if str.match(/(q|p)(\d+|\?)/) and str[ei-1..ei].eql?(")(") # has bands and is not a translocation
      band_s = str.index(/\(/, index)
      band_e = str.index(/\)/, band_s)
      band_e = str.length-1 if band_e.nil?

      bands = str[band_s+1..band_e-1]

      return {:start_index => band_s, :end_index => band_e, :bands => bands.split(/:|;/)}
    end
  end

# sometimes bands are defined for a single chr as p13q22
  def find_fragments(str)
    return str.scan(/([p|q]\d+)/).collect { |a| a[0] }
  end

end