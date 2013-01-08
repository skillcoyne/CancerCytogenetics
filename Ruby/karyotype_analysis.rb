require 'yaml'
require 'date'
require 'fileutils'

require 'cytogenetics'
require 'logger'


def translate_cancer(cancer)
  if $CANCER_TRANSLATION.has_key? cancer
    return $CANCER_TRANSLATION[cancer]
  else
    puts "#{cancer} not found"
    return cancer
  end
end


def cancer_name(cancer)
  cancer.chomp!
  cancer = cancer.split(",")[0]
  cancer = 'unknown' if (cancer.nil? or cancer.match(/N\/a|NA|N\/A/) or cancer.eql? "")
  cancer.downcase!
  cancer.sub!(/nos .*/, "") if cancer.match(/ nos /)
  cancer.lstrip!
  cancer.rstrip!
  cancer = translate_cancer(cancer)
  return cancer.capitalize
end

def quote_str(terms)
  terms.map! { |e| "\"#{e}\"" }
  return terms.join("\t")
end


def write_aberrations(file, aber, cancer, source)
  aber.each do |abr|
    str = quote_str([abr, Cytogenetics.classify_aberration(abr), cancer, source])
    file.write("#{str}\n")
  end
end

def write_fragments(file, fragments, cancer)
  fragments.each do |f|
    str = quote_str([f.chr, f.start, f.end, cancer])
    file.write("#{str}\n")
  end
end

def write_breakpoints(file, breakpoints, cancer)
  breakpoints.each do |bp|
    str = quote_str([bp.type, bp.to_s, bp.chr, cancer])
    file.write("#{str}\n")
  end
end

def write_ploidy(file, ploidy, cancer)
  ploidy.each { |e|
    file.write("#{quote_str([e, cancer])}\n")
  }
end

## NCBI SKY-FISH
def ncbi_skyfish(dir, args)
  count = 0
  esidir = "#{dir}/ESI/karyotype"
  Dir.foreach(esidir) do |entry|
    file = "#{esidir}/#{entry}"
    next if entry.start_with?(".")
    next if File.directory?(file)

#    puts "Reading #{entry}..."

    File.open("#{esidir}/#{entry}", 'r').each_with_index do |line, i|
      next if i.eql? 0
      line.chomp!

      (kcase, diag, stage, karyotypes) = line.split("\t")
      diag = cancer_name(diag)

      next if kcase.match(/mouse/)
      $LOG.info("Reading #{file} karyotype #{i}")

      karyotypes.split(/\//).each do |karyotype|
        begin
          kt = Cytogenetics.karyotype(karyotype)
          write_breakpoints(args[:bpf], kt.report_breakpoints, diag)
          write_fragments(args[:fragf], kt.report_fragments, diag)
          write_ploidy(args[:pf], kt.report_ploidy_change, diag)
          write_aberrations(args[:abr], kt.karyotype, diag, "ncbi")
          count += 1
        rescue Cytogenetics::StructureError => gse
          $LOG.info("#{gse.message}: #{entry}")
        end
      end
    end
  end
  count
end

## Mitelman karyotypes
def mitelman(dir, args)
  count = 0
  File.open("#{dir}/mitelman/mm-kary_cleaned.txt", 'r').each_with_index do |line, i|
    line.chomp!
    next if line.start_with? "#"

#    puts "Reading  Mitelman karyotype # #{i}"
    $LOG.info("Reading  Mitelman karyotype # #{i}: #{dir}/mm-karyotypes.txt")
    (karyotype, morph, shortmorph, refno, caseno) = line.split(/\t/)
    morph = cancer_name(morph)

    begin
      kt = Cytogenetics.karyotype(karyotype)
      write_breakpoints(args[:bpf], kt.report_breakpoints, morph)
      write_fragments(args[:fragf], kt.report_fragments, morph)
      write_ploidy(args[:pf], kt.report_ploidy_change, morph)
      write_aberrations(args[:abr], kt.karyotype, morph, "mitelman")
      count += 1
    rescue Cytogenetics::StructureError => gse
      $LOG.info("#{gse.message}: Mitelman line #{i}")
    end
  end
  count
end

def cam_tissues(dir, args)
  ## Cambridge
  camdir = "#{dir}/path.cam.ac.uk"
  count = 0
  Dir.foreach(camdir) do |tissuedir|
    next if tissuedir.start_with?(".")
    next unless File.directory? "#{camdir}/#{tissuedir}"

    Dir.foreach("#{camdir}/#{tissuedir}") do |entry|
      next if entry.start_with?(".")
      next if entry.eql? "url.txt"
      file = "#{camdir}/#{tissuedir}/#{entry}"

#      puts "Reading #{file}..."
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
          count += 1
        rescue Cytogenetics::StructureError => gse
          $LOG.info("#{gse.message}: #{file}")
        end
      end
    end
  end
  count
end

def create_file(name, comment="#\n", columns)
  file = File.open(name, 'w')
  file.write(comment)
  file.write(columns.join("\t"))
  file.write("\n")
  return file
end

##### ------ MAIN ------ ####


dir = "#{Dir.home}/Data/sky-cgh"

time = Time.new
date = time.strftime("%d%m%Y")

outdir = "#{dir}/output/#{date}"
logdir = "#{dir}/logs/#{date}"


FileUtils.mkpath(outdir)
FileUtils.mkpath(logdir)

$LOG = Logger.new("#{logdir}/karyotype-parse-errors.txt")
$LOG.datetime_format = "%M"
$LOG.level = Logger::INFO
Cytogenetics.logger = $LOG


$CANCER_TRANSLATION = {}
File.open("resources/cancers.csv", 'r').each_line do |line|
  line.chomp!
  (cancer, translation) = line.split(",")
  translation = cancer if translation.nil?
  $CANCER_TRANSLATION[cancer.downcase] = translation.downcase
end

comment = "## Includes mitelman/ncbi/cam karyotypes\n"
files = {
    :bpf => create_file("#{outdir}/breakpoints.txt", comment, ['Event', 'Breakpoint', 'Chr', 'Cancer']),
    :fragf => create_file("#{outdir}/fragments.txt", comment, ['Chr', 'Start', 'End', 'Cancer']),
    :pf => create_file("#{outdir}/ploidy.txt", comment, ['Ploidy', 'Cancer']),
    :abr => create_file("#{outdir}/aberrations.txt", comment, ['Aberration', 'Class', 'Cancer', 'Source'])
}

## Data readers
total = ncbi_skyfish(dir, files)
total += mitelman(dir, files)

comment = "## cam karyotypes only\n"
files = {
    :bpf => create_file("#{outdir}/cam_breakpoints.txt", comment, ['Event', 'Breakpoint', 'Chr', 'Cancer']),
    :fragf => create_file("#{outdir}/cam_fragments.txt", comment, ['Chr', 'Start', 'End', 'Cancer']),
    :pf => create_file("#{outdir}/cam_ploidy.txt", comment, ['Ploidy', 'Cancer']),
    :abr => create_file("#{outdir}/cam_aberrations.txt", comment, ['Aberration', 'Cancer', 'Source'])
}
total += cam_tissues(dir, files)

print "#{total} karyotypes."
