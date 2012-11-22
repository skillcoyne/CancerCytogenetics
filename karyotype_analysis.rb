require 'yaml'
require 'date'

require_relative 'lib/karyotype'
require_relative 'lib/logging'
include Logging

def write_fragments(file, fragments, cancer)
  fragments.each do |f|
    file.write("#{f.chr}\t#{f.start}\t#{f.end}\t#{cancer}\n")
  end
end

def write_breakpoints(file, breakpoints, cancer)
  breakpoints.each do |bp|
    file.write("#{bp.type}\t#{bp.to_s}\t#{bp.chr}\t#{cancer}\n")
  end
end

def write_ploidy(file, ploidy, cancer)
  ploidy.each {|e| file.write("#{e}\t#{cancer}\n") }
end

## NCBI SKY-FISH
def ncbi_skyfish(dir, args)
  bpf = args[:bpf]; fragf = args[:fragf]; pf = args[:pf]

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
          write_breakpoints(bpf, kt.report_breakpoints, diag)
          write_fragments(fragf, kt.report_fragments, diag)
          write_ploidy(pf, kt.report_ploidy_change, diag)
        rescue GenomeStructureError => gse
          log.info("#{gse.message}: #{entry}")
        end
      end
    end
  end
end

## Mitelman karyotypes
def mitelman(dir, args)
  bpf = args[:bpf]; fragf = args[:fragf]; pf = args[:pf]

  File.open("#{dir}/mm-kary-cancer.txt", 'r').each_with_index do |line, i|
    line.chomp!
    next if line.start_with? "#"
    puts "Reading  Mitelman karyotype # #{i}"
    log.info("Reading  Mitelman karyotype # #{i}: #{dir}/mm-karyotypes.txt")
    begin
      (karyotype, morph, shortmorph) = line.split(/\t/)
      kt = Karyotype.new(karyotype)
      write_breakpoints(bpf, kt.report_breakpoints, morph)
      write_fragments(fragf, kt.report_fragments, morph)
      write_ploidy(pf, kt.report_ploidy_change, morph)
    rescue GenomeStructureError => gse
      log.info("#{gse.message}: Mitelman line #{i}")
    end
  end
end

def cam_tissues(dir, args)
  bpf = args[:bpf]; fragf = args[:fragf]; pf = args[:pf]

  ## Cambridge
  camdir = "#{dir}/path.cam.ac.uk"

  Dir.foreach(camdir) do |tissuedir|
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
          write_breakpoints(bpf, kt.report_breakpoints, tissuedir)
          write_fragments(fragf, kt.report_fragments, tissuedir)
          write_ploidy(pf, kt.report_ploidy_change, tissuedir)
        rescue GenomeStructureError => gse
          log.info("#{gse.message}: #{file}")
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

comment = "## Includes mitelman/ncbi/cam karyotypes\n"

bpf = File.open("#{dir}/breakpoints.#{date}.txt", 'w')
bpf.write(comment)
bpf.write("Event\tBreakpoint\tChr\tCancer\n")

fragf = File.open("#{dir}/fragments.#{date}.txt", 'w')
fragf.write(comment)
fragf.write("Chr\tStart\tEnd\tCancer\n")

pf = File.open("#{dir}/ploidy.#{date}.txt", 'w')
pf.write(comment)
pf.write("Change\tCancer\n")

files = {:bpf => bpf, :fragf => fragf, :pf => pf}

## Data readers
ncbi_skyfish(dir, files)
mitelman(dir, files)
cam_tissues(dir, files)

