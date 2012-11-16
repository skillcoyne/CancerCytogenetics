class Fragment
  attr_reader :chr, :start, :end, :genes

  def initialize(*args)
    unless (args.length.eql?2 and (args[0].is_a?Breakpoint and args[1].is_a?Breakpoint))
      raise ArgumentError, "Expected arguments are missing or are not Breakpoints: #{args}"
    end

    #@genes = []
    @start = args[0]
    @end = args[1]
    @chr = @start.chr

    unless @start.chr.eql?@end.chr
      raise GenomeStructureError, "Fragments must be within the same chromosome: #{args}"
    end
  end

  def add_gene(gene)
    @genes << gene
  end

  def to_s
    return "#{@start.to_s} --> #{@end.to_s}"
  end

  ## TODO this will require length in basepairs of each band
  #def length
  #
  #end

end

