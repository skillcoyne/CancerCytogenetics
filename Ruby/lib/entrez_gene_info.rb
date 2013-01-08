require 'yaml'


class EntrezGeneInfo

  attr_reader :tax_id, :gene_id, :official_symbol, :tag, :synonyms, :alt_ids, :chr, :band, :other_bands, :desc, :type, :auth_sym, :auth_full_name, :date

  def initialize(line)
    return nil if line.length <= 0
    parse(line)
  end

  def self.parse_line(line)
    return nil if line.length <= 0
    self.new(line)
  end

  # Format: (tab is used as a separator, pound sign - start of a comment)
  #tax_id
  #GeneID
  #Symbol
  #LocusTag
  #Synonyms
  #dbXrefs
  #chromosome
  #map_location
  #description
  #type_of_gene
  #Symbol_from_nomenclature_authority
  #Full_name_from_nomenclature_authority
  #Nomenclature_status
  #Other_designations
  #Modification_date
  :private

  def parse(line)
    line.chomp!
    (@tax_id, @gene_id, @official_symbol, @tag, @synonyms, alt_ids, @chr, bands, @desc, @type, @auth_sym, @auth_full_name, name_status, other, @date) = line.split("\t")

    @synonyms = @synonyms.split("|")
    alt_ids = alt_ids.split("|")
    @alt_ids = {}
    alt_ids.each do |e|
      (db, id) = e.split(":")
      @alt_ids[db] = id
    end


    bands = bands.split(/\||-/)
    @band = bands[0]
    if @band
      @band.sub!(@chr, "")
      @band = nil if (@band.eql? '-' or @band.eql?'p' or @band.eql?'q')
      @other_bands = bands[1, bands.length]
    end
  end


end