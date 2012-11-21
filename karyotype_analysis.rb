require 'yaml'
require 'date'

require_relative 'lib/karyotype'
require_relative 'lib/logging'
include Logging

def write_fragments(file, fragments)
  fragments.each do |f|
    file.write("#{f.chr}\t#{f.start}\t#{f.end}\n")
  end
end

def write_breakpoints(file, breakpoints)
  breakpoints.each do |bp|
    file.write("#{bp.type}\t#{bp.to_s}\t#{bp.chr}\n")
  end
end

## NCBI SKY-FISH
def ncbi_skyfish(args)
  bpf = args[:bpf]; fragf = args[:fragf]; dir = args[:dir]

  esidir = "#{dir}/ESI/karyotype"
  Dir.foreach(esidir) do |entry|
    file = "#{esidir}/#{entry}"
    next if entry.start_with?(".")
    next if File.directory?(file)

    puts "Reading #{entry}..."

    File.open("#{esidir}/#{entry}", 'r').each_with_index do |line, i|
      next if i.eql? 0
      line.chomp!

      (kcase, diag, stage, karyotypes) = line.split("\t")
      next if kcase.match(/mouse/)
      log.info("Reading #{file} karyotype #{i}")

      karyotypes.split(/\//).each do |karyotype|
        begin
          kt = Karyotype.new(karyotype)
          write_breakpoints(bpf, kt.report_breakpoints)
          write_fragments(fragf, kt.report_fragments)
        rescue GenomeStructureError => gse
          log.info("#{gse.message}: #{entry}")
          #rescue => error
          #  log.error("Failed to parse karyotype from #{entry}:#{karyotype}: #{error.message}")
          #  log.error(error.backtrace)
          #  puts error.backtrace
        end
      end
    end
  end
end

## Mitelman karyotypes
def mitelman(args)
  bpf = args[:bpf]; fragf = args[:fragf]; dir = args[:dir]

  File.open("#{dir}/mm-karyotypes.txt", 'r').each_with_index do |line, i|
    line.chomp!
    next if line.start_with? "#"
    puts "Reading  Mitelman karyotype # #{i}"
    log.info("Reading  Mitelman karyotype # #{i}: #{dir}/mm-karyotypes.txt")
    begin
      kt = Karyotype.new(line)
      write_breakpoints(bpf, kt.report_breakpoints)
      write_fragments(fragf, kt.report_fragments)
    rescue GenomeStructureError => gse
      log.info("#{gse.message}: Mitelman line #{i}")
      #rescue => error
      #  log.error("Failed to parse karyotype from Mitelman line #{i}: #{error.message}")
      #  log.error(error.backtrace)
      #  puts error.backtrace
    end
  end
end

def cam_tissues(args)
  bpf = args[:bpf]; fragf = args[:fragf]; dir = args[:dir]

  ## Cambridge
  camdir = "#{dir}/path.cam.ac.uk"

  Dir.foreach(camdir) do |tissuedir|
#  file = "#{esidir}/#{entry}"
    next if tissuedir.start_with?(".")
    next unless File.directory? "#{camdir}/#{tissuedir}"

    Dir.foreach("#{camdir}/#{tissuedir}") do |entry|
      next if entry.start_with?(".")
      next if entry.eql? "url.txt"
      file = "#{camdir}/#{tissuedir}/#{entry}"

      puts "Reading #{file}..."
      log.info("Reading #{file}")
      File.open(file, 'r').each_line do |karyotype|
        karyotype.chomp!
        next if karyotype.length <= 1

        begin
          kt = Karyotype.new(karyotype)
          write_breakpoints(bpf, kt.report_breakpoints)
          write_fragments(fragf, kt.report_fragments)
        rescue GenomeStructureError => gse
          log.info("#{gse.message}: #{file}")
          #rescue => error
          #  log.error("Failed to parse karyotype from #{file}: #{error.message}")
          #  log.error(error.backtrace)
          #  puts error.backtrace
        end
      end
    end
  end
end


##### ------ MAIN ------ ####
dir = "/Users/sarah.killcoyne/Data/sky-cgh"

time = Time.new
date = "#{time.day}#{time.month}#{time.year}"

Logging.configure(:out => "#{dir}/karyotype-parse-errors.#{date}.txt")
log.datetime_format = "%M"

bpf = File.open("#{dir}/breakpoints.#{date}.txt", '2')
bpf.write("Event\tBreakpoint\tChr\n")

fragf = File.open("#{dir}/fragments.#{date}.txt", '2')
fragf.write("Chr\tFrom\tTo\n")


args = {:bpf => bpf, :fragf => fragf, :dir => dir}
ncbi_skyfish(args)
mitelman(args)
cam_tissues(args)

