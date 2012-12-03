require 'yaml'
require 'date'
require 'fileutils'

require 'cytogenetics'
require 'logger'


def write_aberrations(file, aber, cancer, source)
  aber.each do |abr|
    file.write("#{abr}\t#{cancer}\t#{source}\n")
  end
end

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
      $LOG.info("Reading #{file} karyotype #{i}")

      karyotypes.split(/\//).each do |karyotype|
        begin
          kt = Cytogenetics.karyotype(karyotype)
          write_breakpoints(args[:bpf], kt.report_breakpoints, diag)
          write_fragments(args[:fragf], kt.report_fragments, diag)
          write_ploidy(args[:pf], kt.report_ploidy_change, diag)
          write_aberrations(args[:abr], kt.karyotype, diag, "ncbi")
        rescue Cytogenetics::StructureError => gse
          $LOG.info("#{gse.message}: #{entry}")
        end
      end
    end
  end
end

## Mitelman karyotypes
def mitelman(dir, args)
  File.open("#{dir}/mitelman/mm-kary_cleaned.txt", 'r').each_with_index do |line, i|
    line.chomp!
    next if line.start_with? "#"

    puts "Reading  Mitelman karyotype # #{i}"
    $LOG.info("Reading  Mitelman karyotype # #{i}: #{dir}/mm-karyotypes.txt")
    (karyotype, morph, shortmorph, refno, caseno) = line.split(/\t/)
    begin
      kt = Cytogenetics.karyotype(karyotype)
      write_breakpoints(args[:bpf], kt.report_breakpoints, morph)
      write_fragments(args[:fragf], kt.report_fragments, morph)
      write_ploidy(args[:pf], kt.report_ploidy_change, morph)
      write_aberrations(args[:abr], kt.karyotype, morph, "mitelman")
    rescue Cytogenetics::StructureError => gse
      $LOG.info("#{gse.message}: Mitelman line #{i}")
    end
  end
end

def cam_tissues(dir, args)
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
      $LOG.info("Reading #{file}")
      File.open(file, 'r').each_line do |karyotype|
        karyotype.chomp!
        next if karyotype.length <= 1

        begin
          kt = Cytogenetics.karyotype(karyotype)
          write_breakpoints(args[:bpf], kt.report_breakpoints, tissuedir)
          write_fragments(args[:fragf], kt.report_fragments, tissuedir)
          write_ploidy(args[:pf], kt.report_ploidy_change, tissuedir)
          write_aberrations(args[:abr], kt.karyotype, tissuedir, "cambridge")
        rescue Cytogenetics::StructureError => gse
          $LOG.info("#{gse.message}: #{file}")
        end
      end
    end
  end
end


##### ------ MAIN ------ ####
dir = "/Users/sarah.killcoyne/Data/sky-cgh"

time = Time.new
date = "#{time.day}#{time.month}#{time.year}"

outdir = "#{dir}/output/#{date}"
logdir = "#{dir}/logs/#{date}"

FileUtils.mkpath(outdir)
FileUtils.mkpath(logdir)

$LOG = Logger.new("#{logdir}/karyotype-parse-errors.txt")
$LOG.datetime_format = "%M"
$LOG.level = Logger::INFO
Cytogenetics.logger = $LOG

comment = "## Includes mitelman/ncbi/cam karyotypes\n"
bpf = File.open("#{outdir}/breakpoints.txt", 'w')
bpf.write(comment)
bpf.write("Event\tBreakpoint\tChr\tCancer\n")

fragf = File.open("#{outdir}/fragments.txt", 'w')
fragf.write(comment)
fragf.write("Chr\tStart\tEnd\tCancer\n")

pf = File.open("#{outdir}/ploidy.txt", 'w')
pf.write(comment)
pf.write("Change\tCancer\n")

abr = File.open("#{outdir}/aberrations.txt", 'w')
abr.write(comment)
abr.write("Aberration\tCancer\tSource\n")

files = {:bpf => bpf, :fragf => fragf, :pf => pf, :abr => abr}

## Data readers
ncbi_skyfish(dir, files)
mitelman(dir, files)
cam_tissues(dir, files)

