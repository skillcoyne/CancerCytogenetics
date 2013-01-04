
require 'yaml'

cancers = {}
File.open("../resources/cancers.csv", 'r').each_line do |line|
  line.chomp!
  (cancer, translation) = line.split(",")
  translation = cancer if translation.nil?
  cancers[cancer] = translation
end

#puts YAML::dump cancers


str = "Acute monoblastic leukemia with differentiation (fab type m5b)"


puts cancers.has_key?str



#cancers.each_key do |c|
#
#  if str.match(/^#{c}/)
#    puts "#{str}: #{c}"
#  end
#
#end


