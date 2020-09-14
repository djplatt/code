# compute Dusart epsilon for 20 < log (x) < 21

def small_table(k, step):
  b = 20
  s = 0.86
  m = 5
  d = 1.595e-5
  T = 1132492
  b_list = []
  e_list = []
  while b < 21:
    b_list.append(b)
    e_list.append(eps(b, m, s, d, T).n())
    b = b + step
    print "{}% finished.".format(round(100*(b - 20), 2))
  o = open('small_table_{}_{}.txt'.format(k, step), 'a')
  # o = open('../data/small_table_{}_{}.txt'.format(k, step), 'a')
  o.write("Values of epsilon, 20 < log(x) < 21, quick and dirty version. \r\n\r\n")
  for i in range(len(b_list)):
    o.write("For b = {}, epsilon = {}. \r\n".format(b_list[i], e_list[i]))
  o.write("\r\n\r\nValues of eta_{}.\r\n\r\n".format(k))
  for i in range(len(b_list) - 1):
    o.write("For x in (e^{},e^{}), \eta_{} >= {}. \r\n".format(b_list[i], b_list[i+1], k, b_list[i+1]^k*e_list[i]))
  o.close
  print "Done."
