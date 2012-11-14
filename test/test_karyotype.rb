require 'yaml'
require_relative '../lib/karyotype'
require_relative '../lib/chromosome_aberrations'


k = [
    "43-45,XY,add(2;3)(q1?3),-3,-5,der(5)t(5;17)(q1?3;q2?1),-7,i(8)(q10),-11,-17,ider(19)(q10)add(19)(q13)",
     "46(45-49), XY, 1x1, der(1)t(1;4), 2x2, 3x1, der(3)t(5;3;7), der(3)t(3;5)(p12;p12), der(4)t(4;1;4;7), der(4)del(4)t(4;7), der(5)t(4;5), i(5p), der(5?)t(5;18), 6x2, der(6)t(6;8), 7x1, 8x2, 9x2, 10x2, 11x2, 12x2, 13x2, der(13)t(8;13), 14x2, 15x2, dup(15), 16x2, 17x2, 18x1, der(18)t(14;18;8), 19x2, 20x1, 21x2, 22x2",
     "66~185<5n>,XXXX,+X[3],-X[3],del(X)(q22)[7],+del(X)(q22)[5],+del(X)(q22)[3],-1[5],+2[2],-2[4],del(2)(p21)[2],+3[2],-3[4],+4[3],-4[4],+5[6],+5[4],+5[2],+5[2],+der(5;21)(p10;q10)[4],-6[6],-6[4],+7[6],+7[4],+7[2],+del(7)(p14)[3],+8[5],+8[3],+8[2],-9[6],-9[6],-10[6],-10[5],del(10)(q21)[3],-11[5],del(11)(p14)[2],del(11)(q11)[4],del(11)(q11)[2],der(11)del(11)(p14)del(11)q(23)[2],-12[4],-12[3],-13[4],-13[3],+del(13)(q22)[3],-14[7],-14[5],-14[3],-15[4],-15[3],del(15)(q22)x2[5],+16[3],-16[4],-16[3],+17[6],del(17)(p11)x2[8],+del(17)(p11)x2[6],-18[5],-18[4],der(18)t(14;18)(q1?3;p11.3)[2],-19[6],-19[3],der(19;21)(?p10;q10)[2],+20[2],-20[5],+21[3],-21[3],-21[3],+22[2],-22[5],3~8min[5]"
     ]


#puts k[0]
kt = Karyotype.new(k[0])
puts kt.sex
puts kt.ploidy
puts kt.karyotype.join(",")
puts  kt.aberrations
kt.analyze

#puts kt.breakpoints

#puts ChromosomeAberrations::Inversion.type
#puts ChromosomeAberrations::Inversion.regex
#exit

exit



ChromosomeAberrations.constants.each do |ca|
  puts "----#{ca}-----"

  #puts "*****INV*****" if ChromosomeAberrations.const_get(ca).type.eql? "inv"

  obj = ChromosomeAberrations.const_get(ca)
  puts obj.type
  puts obj.regex
  obj = obj.new("der(3)t(3;5)(p12;p12)")

  #puts ChromosomeAberrations.const_get(ca).name
  #puts ChromosomeAberrations.const_get(ca).type
  #
  #puts ChromosomeAberrations.const_get(ca).regex
end




#iv = ChromosomeAberrations::Inversion.new("inv(4)(q13p4)")
#puts YAML::dump iv.breakpoints


