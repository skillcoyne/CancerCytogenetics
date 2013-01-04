

fout = File.open("bands_by_chr.txt", 'w') 
fout.write("Chr\tBand\n")
File.open("bands.txt", 'r').each_line do |line|
  line.chomp!

	line.match(/(\d+|X|Y)[q|p].*/)
  chr = $1

  fout.write("#{chr}\t#{line}\n")
end
