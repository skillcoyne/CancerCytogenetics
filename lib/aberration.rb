require 'yaml'

class Aberration

  attr_reader :abr, :ab_objs, :breakpoints, :fragments

  class<<self

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
    return "unk".to_sym
  end

  def initialize(str)
    Logging.configure(:class => self.class.name)

    @abr = str
    @breakpoints = []; @fragments = []

    #regex = Aberration.regex[@type.to_sym]
    # make sure it really is an inversion first
    #raise KaryotypeError, "#{str} does not appear to be a #{self.class}" unless str.match(self.regex)
    get_breakpoints() #(@abr)
    @breakpoints.flatten!
  end

  def to_s
    "#{@abr}: #{@breakpoints.join(',')}"
  end

  :private

  def get_breakpoints
    chr_i = find_chr(@abr)
    return if chr_i.nil?

    band_i = find_bands(@abr, chr_i[:end_index])

    unless band_i.nil? # breakpoints aren't added if there is no band information
      chr_i[:chr].each_with_index do |c, i|
        fragments = find_fragments(band_i[:bands][i])
        fragments.each { |f| @breakpoints << Breakpoint.new(c, f, @type) }
      end
    else
      ## No band --> TODO add this as information somewhere but not as a breakpoint
      #@breakpoints << Breakpoint.new(c, "", @type)
    end
  end

  # Parsing aberration strings to pull out the chromosome and band definitions
  # These will result in breakpoint information
  def find_chr(str)
    chr_s = str.index(/\(/, 0)
    chr_e = str.index(/\)/, chr_s)
    chrs = str[chr_s+1..chr_e-1].split(/;|:/)
    chrs.each do |chr|
      unless chr.match(/^\d+|X|Y$/)
        log.warn("No chromosome defined from #{str}, skipped.")
        return
      end
    end
    return {:start_index => chr_s, :end_index => chr_e, :chr => chrs}
  end

  def find_bands(str, index)
    band_info = nil
    #raise KaryotypeError, "No bands defined in #{str}" if str.length.eql?(index+1)
    if str.length.eql?(index+1)
      log.warn("No bands defined in #{str}, skipped.")
      return
    end

    ei = str.index(/\(/, index)
    if str.match(/(q|p)(\d+|\?)/) and str[ei-1..ei].eql?(")(") # has bands and is not a translocation
      band_s = str.index(/\(/, index)
      band_e = str.index(/\)/, band_s)
      band_e = str.length-1 if band_e.nil?
      bands = str[band_s+1..band_e-1].split(/;|:/)

      if str[band_s+1..band_e-1].match(/::/)
        log.warn("Aberration defined using different language, not currently parsed skipping: #{@abr}")
        return band_info
      else
        bands.map! {|b| b.sub(/-[q|p]\d+$/, "")} # sometimes bands are given a range, for our purposes we'll take the first one (CyDas appears to do this as well)
        bands.each do |b|
          unless b.match(/^[p|q]\d+(\.\d)?$/)
            log.warn("Bands incorrectly defined in #{str}")
            return band_info
          end
        end
        band_info = {:start_index => band_s, :end_index => band_e, :bands => bands}
      end
    end
    return band_info
  end

# sometimes bands are defined for a single chr as p13q22
  def find_fragments(str)
    return str.scan(/([p|q]\d+)/).collect { |a| a[0] }
  end

end