
class Breakpoint
  attr_accessor :chr, :band, :type

  def initialize(*args)
    c = args[0]; b = args[1]
    @type = args[2] if args.length > 2

    unless ((c.is_a? String and c.match(/\d+|X|Y/)) and (b.is_a? String and b.length > 0))
      log.error("#{c}#{b} is not a valid breakpoint")
      raise ArgumentError, "#{c}#{b} is not a valid breakpoint"
    end
    @chr = c; @band = b
  end

  def arm
    @band.match(/(q|p)\d+/)
    return $1
  end

  def to_s
    return "#{@chr}#{@band}"
  end

end