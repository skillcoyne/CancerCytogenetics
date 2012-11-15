require_relative 'logging'

class Breakpoint
  include Logging

  attr_accessor :chr, :band, :type

  def initialize(*args)
    c = args[0]; b = args[1]
    @type = args[2] if args.length > 2

    unless ((c.is_a? String and c.match(/\d+|X|Y/)) and (b.is_a? String))
      log.error("#{c}#{b} is not a valid breakpoint")
      raise ArgumentError, "#{c}#{b} is not a valid breakpoint"
    end
    @chr = c; @band = b
  end

  def to_s
    return "#{@chr}#{@band}"
  end


end