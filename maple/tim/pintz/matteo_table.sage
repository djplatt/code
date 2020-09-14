def matteo_table():
    o = open("epsilon_table.txt", "a")
    o.write("Values of epsilon, 1000 < log(x) < 3500, quick and dirty version. \r\n\r\n")
    s = 0.97
    d = 1.574e-11
    m = 9
    T = 1012519281279
    for b in range (1000,1500,5):
        o.write("For b = {}, epsilon = {}. \r\n".format(b, eps(b, m, s, d, T).n()))
        print "{}% finished.".format(((b - 1000)/25).n())
    #quick end
    o.close
    return "done!"
    s = 0.98
    d = 8.852e-12
    T = 1019509030546
    m = 5
    for b in range (1500,2000,5):
        o.write("For b = {}, epsilon = {}. \r\n".format(b, eps(b, m, s, d, T).n()))
        print "{}% finished.".format(((b - 1000)/25).n())
    m = 2
    d = 3.381e-12
    T = 1364832983117
    for b in range (2000,2500,5):
        o.write("For b = {}, epsilon = {}. \r\n".format(b, eps(b, m, s, d, T).n()))
        print "{}% finished.".format(((b - 1000)/25).n())
    d = 1.193e-12
    T = 2445999556029
    for b in range (2500,3000,5):
        o.write("For b = {}, epsilon = {}. \r\n".format(b, eps(b, m, s, d, T).n()))
        print "{}% finished.".format(((b - 1000)/25).n())
    d = 4.209e-13
    T = 2445999556030
    for b in range (3000,3501,5):
        o.write("For b = {}, epsilon = {}. \r\n".format(b, eps(b, m, s, d, T).n()))
        print "{}% finished.".format(((b - 1000)/25).n())
    o.close
