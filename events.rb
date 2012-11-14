#require 'mysql2'
require 'yaml'
require_relative 'lib/event_matrix'


def parse_events(karyotype, event_matrix)
  aberrations = karyotype.split(",")
  events = []

  aberrations.each_with_index do |a, i|
    next if i <= 1
    a.gsub!(/\s/, "")
    a.gsub!(/\[\d+\]/, "")
    next if (a.length <= 1 or a.match(/\?/))

    a.sub!(/^[\+|-]/, "") unless a.match(/^[\+|-][\d|X|Y]+$/) # +t(13;X)(q13;p12) doesn't need a +

    # 13x2 is normal, 13x3 is a duplicate and should read +13
    if a.match(/^([\d+|X|Y]+)x(\d+)/)
      chr = $1; dups = $2.to_i
      if dups.eql? 0
        a = "-#{chr}"
        event_matrix.add_event(a, 2)
        events.push(a)
      elsif dups > 2
        dups -= 2
        a = "+#{chr}"
        event_matrix.add_event(a, dups)
        events.push(a)
      end
      # add(7)x2 should read +7, presume 'add' does not indicate normal diploid even when x2
    elsif a.match(/^add\(([\d|X|Y]+)\)x(\d+)/)
      chr = $1; dups = $2.to_i
      a = "+#{chr}"
      event_matrix.add_event(a, dups)
      events.push(a)
      # add(9)(p21)x2 should indicate that this happened twice
    elsif a.match(/(.*)x(\d+)$/)
      a = $1; dups = $2.to_i
      event_matrix.add_event(a, dups)
      events.push(a)
      # del(7) should be -7  but not del(7)(q12)
    elsif a.match(/^del\(([\d|X|Y]+)\)$/)
      chr = $1
      a = "-#{chr}"
      event_matrix.add_event(a)
      events.push(a)
    else # everything else
      event_matrix.add_event(a)
      events.push(a)
    end
  end

  event_matrix.link_events(events)

  return event_matrix
end


if ARGV.length <= 0
  print "Root directory for ESI/cam karyotypes required"
  exit
end
dir = ARGV[0]


em = EventMatrix.new()

#em.add_event('a')
#em.add_event('b', 5)
#em.add_event('c')
#
#em.add_event('b')
#
##
##em.link_events(['a', 'c'])
##em.link_events(['a', 'b'])
#events = em.events
#puts "\t" + events.join("\t")
#em.matrix.each_with_index do |column, i|
#  puts "#{events[i]}\t" + column.join("\t")
#end
#
#em.output(File.new("testmatrix.txt", 'w'))
#
#exit
#

ktsql = File.open("#{dir}/sql/karyotypes.txt", 'w')
src = 'mitelman'
uid_kt = 1
## mitelman karyotypes
File.open("#{dir}/mm-karyotypes.txt").each_line do |karyotype|
  karyotype.chomp!
  next if karyotype.start_with? "#"
  em = parse_events(karyotype, em)
  #ktsql.write("#{uid_kt}\t#{src}\t#{karyotype}\n")
  uid_kt += 1
end

## NCBI sky-fish karyotypes
esidir = "#{dir}/ESI/karyotype"
src = 'ncbi'
Dir.foreach(esidir) do |entry|
  file = "#{esidir}/#{entry}"
  next if entry.start_with?(".")
  next if File.directory?(file)

  unless (File.basename(entry).match(/\.karyotype/) or File.basename(entry).match(/\.kt/))
    puts "#{entry} is not a karyotype file"
    next
  end

  puts "#{file}"

  kts = 0
  File.open(file, 'r').each_line do |line|
    line.chomp
    next if line.length <= 0
    next if line.match(/mouse/)
    next if line.match(/Case/) # column names
    karyotype = line.split(/\t/)[-1].gsub!(/\s/, "")

    em = parse_events(karyotype, em)

    #ktsql.write("#{uid_kt}\t#{src}\t#{karyotype}\n")
    kts += 1
  end
  puts "#{entry}: #{kts}"

end

## Cambridge karyotypes
camdir = "#{dir}/path.cam.ac.uk"
src = 'cam'
Dir.foreach(camdir) do |cd|
  ktdir = "#{camdir}/#{cd}"
  next if cd.start_with?(".")
  next unless File.directory? ktdir

  Dir.foreach(ktdir) do |entry|
    next if entry.start_with?(".")
    file = "#{ktdir}/#{entry}"

    unless (File.basename(entry).match(/\.karyotype/) or File.basename(entry).match(/\.kt/))
      puts "#{entry} is not a karyotype file"
      next
    end

    File.open(file, 'r').each_line do |karyotype|
      karyotype.chomp!
      em = parse_events(karyotype, em)
      #ktsql.write("#{uid_kt}\t#{src}\t#{karyotype}\n")
      uid_kt += 1
    end
  end
end

#ktsql.close


puts em.events.length

em.output(File.new("#{dir}/events.txt", 'w'))

File.new("#{dir}/codes.txt", 'w') { |f|
  f.write(em.events.join("\n"))
}

File.new("#{dir}/events-coded.txt", 'w') { |f|
  em.matrix.each_with_index do |column, i|
    f.write(column.join("\t") + "\n")
  end
}

#File.open("#{dir}/events.txt", 'w') { |f|
#  f.write("event\tfrequency\n")
#  events.each_pair { |k, v| f.write("#{k}\t#{v}\n") }
#}



