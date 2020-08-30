#define hashmodulus ((1L<<61)-1)
#define hashlength 464
#define hashmin 4096
#define hashmax 8192
static struct {
	unsigned long p,c;
} hashcoeff[hashlength]={
{4099,326490430436040986},
{4111,559705121321738418},
{4127,1027143540648291608},
{4129,1614463795034667624},
{4133,455689193399227776},
{4139,812966537786397194},
{4153,2073755909783565705},
{4157,1309198521558998535},
{4159,486216762465058766},
{4177,1847926951704044964},
{4201,2093254198748441889},
{4211,566490051061970630},
{4217,150232564538834691},
{4219,1356728749119613735},
{4229,987635478950895264},
{4231,1799657824218406718},
{4241,1921416341351658148},
{4243,1827423174086145070},
{4253,1750068052304413179},
{4259,1335382038239683171},
{4261,1126485543032891129},
{4271,2189612557070775619},
{4273,1588425437950725921},
{4283,1906114970148501816},
{4289,1748768066189739575},
{4297,1105553649419536683},
{4327,41823880976068680},
{4337,2246893936121305098},
{4339,680675478232219153},
{4349,1096492737570930588},
{4357,1064058600983463886},
{4363,2124681778445686677},
{4373,1153253523204927664},
{4391,1624949962307564146},
{4397,884760591894578064},
{4409,722684359584774534},
{4421,469294503455959899},
{4423,1078657853083538250},
{4441,497558833647780143},
{4447,430880240537243608},
{4451,1008306263808672160},
{4457,871262219757043849},
{4463,1895004365215350114},
{4481,553114791335573273},
{4483,928282405904509326},
{4493,1298199971472090520},
{4507,1361509731647030069},
{4513,426420832006230709},
{4517,750020119738494463},
{4519,892950654076632414},
{4523,1225464410814600612},
{4547,1911848480297925904},
{4549,842847261377168671},
{4561,836690411740782547},
{4567,36595684701041066},
{4583,57074465600036538},
{4591,35391454785304773},
{4597,1027606372000412697},
{4603,858149375821293895},
{4621,1216392374703443454},
{4637,59308853655414224},
{4639,1030486962058546718},
{4643,382910609159726501},
{4649,768789926722341438},
{4651,762735381378628682},
{4657,1005758989771948074},
{4663,1657009501638647508},
{4673,1783661361016004740},
{4679,796798233969021059},
{4691,1658520847567437422},
{4703,502975179223457818},
{4721,2063998821801160708},
{4723,2126598223478603304},
{4729,817551008849328795},
{4733,1793074162393397615},
{4751,1287596263315727892},
{4759,1629305847068896729},
{4783,2282065591485919335},
{4787,1280388906297308209},
{4789,173159035165825250},
{4793,1203194438340773505},
{4799,2146825320332825798},
{4801,847076010454532974},
{4813,2132606604399767971},
{4817,865350797130078274},
{4831,421223214903942949},
{4861,2202859852828864983},
{4871,1627572340776304831},
{4877,1301036100621122535},
{4889,2151172683210770621},
{4903,555918022010940381},
{4909,1195820575245311406},
{4919,2060813966122583132},
{4931,824196499832939747},
{4933,1252314214858834971},
{4937,380498114208176064},
{4943,621869463771460120},
{4951,1487674193901485781},
{4957,1569074147090699661},
{4967,1723498454689514459},
{4969,1489838779667276265},
{4973,607626788325788389},
{4987,93543108859195056},
{4993,1874271115035734974},
{4999,1456016012787031897},
{5003,619764822731213939},
{5009,1812449603590865741},
{5011,808484663842074461},
{5021,2009697952400734786},
{5023,1525933978789885248},
{5039,343887624789001682},
{5051,1182376379945660137},
{5059,1982314473921546769},
{5077,1109549848371395693},
{5081,1037594154159590924},
{5087,1071053104849367160},
{5099,1322181949714913233},
{5101,1516660949039528341},
{5107,960526604699918173},
{5113,1729904691101240134},
{5119,261117919934717464},
{5147,2271784899875479358},
{5153,756802274277310875},
{5167,1289220444092802136},
{5171,474369139841197116},
{5179,1716815258254385285},
{5189,103716246685267192},
{5197,543779117105835462},
{5209,1645057139707767457},
{5227,895800586311529398},
{5231,1255427590538696616},
{5233,152478208398822237},
{5237,59235267842928844},
{5261,1502771737122401274},
{5273,1149578551939377903},
{5279,1470772656511184950},
{5281,1546086255370076952},
{5297,1723497785943073942},
{5303,778240149963762596},
{5309,240870114509877266},
{5323,394305328258085500},
{5333,2102620516550230799},
{5347,1039820873553197464},
{5351,979798654654721830},
{5381,880027557663442629},
{5387,1676981816531131145},
{5393,1802107305139241263},
{5399,1972433293052973713},
{5407,2107405063590590043},
{5413,1798917982073452520},
{5417,1369268024301602286},
{5419,867033797562981667},
{5431,1038357135783187942},
{5437,758476292223849603},
{5441,1948092882600628075},
{5443,2207529277533454374},
{5449,1821419918118374849},
{5471,1231889908299259230},
{5477,566310110224392380},
{5479,1609356725483962542},
{5483,280378617804444931},
{5501,1072662998681271815},
{5503,116308709127642766},
{5507,1193169610307430309},
{5519,866966243748392804},
{5521,166237193327216135},
{5527,1077013023941018041},
{5531,404884253921467160},
{5557,786088301434511589},
{5563,1383535122407493085},
{5569,2280658829488325172},
{5573,101154688442168806},
{5581,186007322364504054},
{5591,132651484623670765},
{5623,2214024743056683473},
{5639,2082072212962344576},
{5641,1527055902872993253},
{5647,914904768868572390},
{5651,828477094595207304},
{5653,1020679050708770534},
{5657,482636846586846145},
{5659,1930865547754160712},
{5669,1593671129282272719},
{5683,1493198467868909485},
{5689,729902645271416500},
{5693,275540268357558312},
{5701,164114802119030362},
{5711,788447619988896953},
{5717,1762740703670330645},
{5737,660855577878083177},
{5741,1651988416921493024},
{5743,740652833177384429},
{5749,1112201596451006206},
{5779,415698847934834932},
{5783,1211582319647132127},
{5791,1146510220721650373},
{5801,1849436445614060470},
{5807,2087092872652683432},
{5813,2118502348483502728},
{5821,1356524772912098481},
{5827,1199384942357517449},
{5839,172551026757446140},
{5843,578031956729941707},
{5849,523340081847222890},
{5851,1076777027268874489},
{5857,504399020033657060},
{5861,1278551106709414382},
{5867,2159465951497451565},
{5869,1178157191616537256},
{5879,204263226455195995},
{5881,1056341819781968292},
{5897,183521353142147658},
{5903,2188450004032853736},
{5923,815413180157425263},
{5927,1872285744226329343},
{5939,959184959959358956},
{5953,473007083155872003},
{5981,655761716995053547},
{5987,1131460430873190185},
{6007,2139124645518872072},
{6011,511733859594496686},
{6029,15198510254334311},
{6037,1224323599606986326},
{6043,717867206610437778},
{6047,2091512354759023324},
{6053,372342232752868676},
{6067,1361511712413436237},
{6073,1389190973283340505},
{6079,394349220142131124},
{6089,2079377585202309849},
{6091,353365880305796299},
{6101,2032166139485738617},
{6113,1890917131797951728},
{6121,242865361432353437},
{6131,1418792507475867019},
{6133,2119099350463010017},
{6143,1014188227490285243},
{6151,479492624224518275},
{6163,1303029569429482669},
{6173,517247294593876834},
{6197,1554557044656123283},
{6199,750281115903727536},
{6203,2167122262389919937},
{6211,760554688782332821},
{6217,2023636030598854916},
{6221,1790146557619247357},
{6229,386163722563943194},
{6247,1515274606763521578},
{6257,2204179368292080266},
{6263,964158696771121369},
{6269,303439105863023359},
{6271,8182230548124380},
{6277,1750434984088519049},
{6287,1725590414598766182},
{6299,1265114980378421064},
{6301,1015227773830014864},
{6311,229929992560423398},
{6317,764214183816115409},
{6323,538352539450824188},
{6329,1941773060895353999},
{6337,1068434172733967371},
{6343,1355790773646160387},
{6353,459324502245141234},
{6359,609129328626229402},
{6361,1241119177010491262},
{6367,1783576433920437207},
{6373,1523680846139002895},
{6379,882824005398680507},
{6389,413096479776864968},
{6397,522865969927243974},
{6421,1858351603281690756},
{6427,1968585526421383793},
{6449,2178118415854385403},
{6451,2071714574011626742},
{6469,2075065799199309684},
{6473,2276241901353008033},
{6481,303400925906664587},
{6491,1426227202230524239},
{6521,1930606302598963877},
{6529,249953308414640146},
{6547,611228839507773914},
{6551,1672745442514341102},
{6553,467604306000306674},
{6563,1474554813214382459},
{6569,1601661712875312382},
{6571,614840167992924737},
{6577,1228071177654928913},
{6581,527816710270111610},
{6599,2217787042387174521},
{6607,639805394326521740},
{6619,222549283662427192},
{6637,1360905187147121266},
{6653,2218130040969485290},
{6659,1295851844663939225},
{6661,563784543912533038},
{6673,1995338666855533310},
{6679,1570565903061390324},
{6689,1421390998286027062},
{6691,1394318358097125191},
{6701,1259069656723159936},
{6703,782274544912671248},
{6709,727119931274242152},
{6719,461373271832281770},
{6733,431218333850664873},
{6737,1192819027123234430},
{6761,2078764559709872649},
{6763,185598300798682005},
{6779,753027393642717163},
{6781,39457098005678485},
{6791,1334017905593361063},
{6793,2208208003949042369},
{6803,995759906937041788},
{6823,1045940157364976040},
{6827,194824647782216037},
{6829,550631184874398695},
{6833,1360200364068800381},
{6841,1357865448826768161},
{6857,1831861326200370539},
{6863,942093021910086667},
{6869,1640270910790040055},
{6871,186615109286328085},
{6883,1330440696074470319},
{6899,499018273810238035},
{6907,502274974614414055},
{6911,1207335215870481547},
{6917,2013999866627689866},
{6947,1419916425046140717},
{6949,191559056573160841},
{6959,1328802988676857752},
{6961,1405960078185023606},
{6967,227507798797399340},
{6971,1637526486952132401},
{6977,1076968863810265335},
{6983,944510191997220613},
{6991,1301386330701215932},
{6997,285779824044017183},
{7001,1429750858521890899},
{7013,1618865668058420542},
{7019,841779507635076338},
{7027,2271885690336656780},
{7039,1950830875641497149},
{7043,2020789551919109899},
{7057,975546679421148460},
{7069,1197104163269028491},
{7079,1270315990156796392},
{7103,748604252817308486},
{7109,816129261753596233},
{7121,384118410847738091},
{7127,2113266006559319391},
{7129,1338854039375748358},
{7151,1361143499198430117},
{7159,633423014922436774},
{7177,1290791779633361737},
{7187,81273616335831288},
{7193,734007502359373790},
{7207,1803343551649794557},
{7211,178160046107106100},
{7213,1669700173018758407},
{7219,1829836142710185153},
{7229,1253431308749847288},
{7237,70019619094993502},
{7243,939065521645602191},
{7247,571602252457140250},
{7253,26887212485499413},
{7283,984459396949257361},
{7297,852773633209386873},
{7307,2289526104158020696},
{7309,756333221468239349},
{7321,478223842701987702},
{7331,2004947225028724200},
{7333,526770890233938212},
{7349,1661268713623636486},
{7351,1595033927594418161},
{7369,1532663022957125952},
{7393,364955822302609296},
{7411,603258635519191127},
{7417,371859597962583054},
{7433,94282227629658712},
{7451,2160611809915747887},
{7457,27000232625437140},
{7459,22687777651226933},
{7477,734430233276692626},
{7481,1127699960534747774},
{7487,346857527391478939},
{7489,399588948728484631},
{7499,1369575405845760568},
{7507,2217803687757289581},
{7517,2206814713288610370},
{7523,130141496340768529},
{7529,861110681132541840},
{7537,230850531138610791},
{7541,1780590839341524422},
{7547,1923534983071749673},
{7549,1055631719355441015},
{7559,1222514615506258219},
{7561,937915311769343786},
{7573,852868812961957254},
{7577,718656592041199719},
{7583,2250542267067365785},
{7589,2169537354592688250},
{7591,1568074419444165342},
{7603,853778925104674827},
{7607,105031681250262217},
{7621,1204393036537417656},
{7639,592755100672600484},
{7643,1509207668054427766},
{7649,1409630039748867615},
{7669,433329873945170157},
{7673,168130078420452894},
{7681,701434349299435396},
{7687,1736119639406860361},
{7691,1801042332324036889},
{7699,82826810621030003},
{7703,581092394588713697},
{7717,1513323039712657034},
{7723,2086339870071049553},
{7727,512802587457892537},
{7741,1294754943443095033},
{7753,1486581673100914879},
{7757,930909063370627411},
{7759,2280060915913643774},
{7789,219424962331973086},
{7793,118156503193461485},
{7817,743557822870685069},
{7823,1997655344719642813},
{7829,393161419260218815},
{7841,1086985962035808335},
{7853,2119375969747368461},
{7867,1650489163325525904},
{7873,1967094695603069467},
{7877,916149623124228391},
{7879,1122737829960120900},
{7883,144869810337397940},
{7901,2261458899736343261},
{7907,1226838560319847571},
{7919,897743852062682980},
{7927,45750188043851908},
{7933,1858576614171946931},
{7937,1568041120172266851},
{7949,289541060457747472},
{7951,1539585379217609033},
{7963,866887122666596526},
{7993,6060188892447452},
{8009,1707684831658632807},
{8011,1062812350167057257},
{8017,887626467969935688},
{8039,1968363468003847982},
{8053,2169897168216456361},
{8059,217716763626832970},
{8069,413611451367326769},
{8081,336255814660537144},
{8087,1464084696245397914},
{8089,1902174501258288151},
{8093,1440415903059217441},
{8101,302153101507069755},
{8111,1558366710940453537},
{8117,717776382684618355},
{8123,1206381076465295510},
{8147,1308718247292688437},
{8161,555835170043458452},
{8167,1029518900794643490},
{8171,1034197980057311552},
{8179,131258234416495689},
{8191,260799345029896943}};