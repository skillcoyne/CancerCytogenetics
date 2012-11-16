require 'yaml'
require 'date'

require_relative 'lib/karyotype'

require_relative 'lib/logging'
include Logging

dir = "/Users/sarah.killcoyne/Data/sky-cgh"


Logging.configure(:out => "#{dir}/karyotype.txt")
log.datetime_format = "%M"


# Start with ESI
esidir = "#{dir}/ESI/karyotype"


Dir.foreach(esidir) do |entry|
  file = "#{esidir}/#{entry}"
  next if entry.start_with?(".")
  next if File.directory?(file)

  puts "Reading #{entry}..."

  File.open("#{esidir}/#{entry}", 'r').each_with_index do |line, i|
    next if i.eql?0
    line.chomp!

    (kcase, diag, stage, karyotype) = line.split("\t")
    next if kcase.match(/mouse/)

    begin
      kt = Karyotype.new(karyotype)
        #puts kt.report_breakpoints
        #puts kt.report_breakpoints

    rescue KaryotypeError => error
      log.error("Failed to read karyotype: #{error.message}")
    end


  end



end


