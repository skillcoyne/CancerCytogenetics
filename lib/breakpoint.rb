class Breakpoint

  attr_accessor :chr, :band, :type

  def initialize(*args)
    c = args[0]; b = args[1]
    @type = args[2] if args.length > 2
    raise ArgumentError, "#{c}#{b} is not a valid breakpoint" unless ( (c.is_a?String and c.match(/\d+|X|Y/)) and (b.is_a?String) )
    @chr = c; @band = b
  end

end