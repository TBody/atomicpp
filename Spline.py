# from scipy.interpolate import RectBivariateSpline #this is the Anaconda version
from interpolate import RectBivariateSpline #this is the local version
import numpy as np

x = np.array([1,2,3,4,5])
y = np.array([1,3,4,7,10])

z = np.zeros((len(x), len(y)))

for i in range(len(x)):
	for j in range(len(y)):
		z[i,j] = x[i]**2 + y[j]**2 + x[i]*y[j]


# x = np.array([ -6.98770000e-01,-5.22680000e-01,-3.00830000e-01,-1.54700000e-01,2.00000000e-04,1.76290000e-01,3.01230000e-01,4.77320000e-01,6.99170000e-01,8.45300000e-01,1.00020000e+00,1.17629000e+00,1.30123000e+00,1.47732000e+00,1.69917000e+00,1.84530000e+00,2.00020000e+00,2.17629000e+00,2.30123000e+00,2.47732000e+00,2.69917000e+00,2.84530000e+00,3.00020000e+00,3.17629000e+00,3.30123000e+00,3.47732000e+00,3.69917000e+00,3.84530000e+00,4.00020000e+00,4.17629000e+00])
# y = np.array([13.69897,14.0,14.30103,14.69897,15.0,15.30103,15.69897,16.0,16.30103,16.69897,17.0,17.30103,17.69897,18.0,18.30103,18.69897,19.0,19.30103,19.69897,20.0,20.30103,20.69897,21.0,21.30103])
# z = np.array([[-37.96838,-37.93698,-37.91227,-37.9028,-37.88125,-37.80117,-37.65822
# 	,-37.52789,-37.32822,-36.98097,-36.67732,-36.34809,-35.89351,-35.55016
# 	,-35.21716,-34.79272,-34.48058,-34.17345,-33.772,-33.46974,-33.16805
# 	,-32.77547,-32.48308,-32.19131]
# 	,[-29.84573,-29.80925,-29.76818,-29.69848,-29.64733,-29.615,-29.56207
# 	,-29.50993,-29.46273,-29.39264,-29.31933,-29.21915,-29.03446,-28.85587
# 	,-28.6504,-28.34568,-28.09438,-27.82795,-27.4594,-27.17236,-26.88111
# 	,-26.49646,-26.20735,-25.91849]
# 	,[-23.37955,-23.35104,-23.3223,-23.28389,-23.25113,-23.21166,-23.15281
# 	,-23.10371,-23.04844,-22.96667,-22.90184,-22.83939,-22.76673,-22.71579
# 	,-22.65796,-22.54525,-22.41653,-22.24392,-21.96146,-21.71788,-21.45691
# 	,-21.10333,-20.83593,-20.56842]
# 	,[-20.62358,-20.59714,-20.56919,-20.53087,-20.49953,-20.46425,-20.41104
# 	,-20.36668,-20.3196,-20.25276,-20.20054,-20.15049,-20.09281,-20.0533
# 	,-20.00974,-19.92596,-19.82803,-19.69163,-19.45797,-19.24826,-19.01876
# 	,-18.71544,-18.49349,-18.27216]
# 	,[-18.56734,-18.54403,-18.51945,-18.48526,-18.45585,-18.42096,-18.36887
# 	,-18.32701,-18.28342,-18.22215,-18.17553,-18.13324,-18.08734,-18.05687
# 	,-18.02245,-17.95336,-17.87101,-17.75595,-17.5589,-17.38245,-17.19175
# 	,-16.95616,-16.79552,-16.63657]
# 	,[-16.91803,-16.90007,-16.88054,-16.85062,-16.82371,-16.79305,-16.74828
# 	,-16.7131,-16.67737,-16.62839,-16.59175,-16.55876,-16.52279,-16.49816
# 	,-16.4692,-16.41007,-16.34008,-16.24339,-16.08059,-15.93906,-15.79213
# 	,-15.62276,-15.51482,-15.40897]
# 	,[-16.04895,-16.03551,-16.02052,-15.99601,-15.97316,-15.94733,-15.90982
# 	,-15.88076,-15.85185,-15.81225,-15.78258,-15.75619,-15.72693,-15.70603
# 	,-15.68087,-15.62816,-15.56601,-15.48156,-15.34053,-15.21978,-15.09772
# 	,-14.96203,-14.87861,-14.79738]
# 	,[-15.16732,-15.15826,-15.14753,-15.12886,-15.11089,-15.09074,-15.06288
# 	,-15.04182,-15.02051,-14.99098,-14.96894,-14.94953,-14.92757,-14.91054
# 	,-14.88833,-14.84167,-14.78786,-14.7158,-14.59595,-14.49531,-14.39757
# 	,-14.29158,-14.22709,-14.16448]
# 	,[-14.42328,-14.41807,-14.41122,-14.39847,-14.38633,-14.37332,-14.35649
# 	,-14.34369,-14.32964,-14.30949,-14.29441,-14.28075,-14.26371,-14.24943
# 	,-14.23091,-14.19103,-14.14455,-14.08263,-13.98037,-13.89488,-13.81211
# 	,-13.72377,-13.67135,-13.62068]
# 	,[-14.07152,-14.06729,-14.06231,-14.05358,-14.0453,-14.03631,-14.0244
# 	,-14.01521,-14.00518,-13.98969,-13.97758,-13.967,-13.95347,-13.94138
# 	,-13.92527,-13.89042,-13.84945,-13.79415,-13.70237,-13.62447,-13.54708
# 	,-13.46552,-13.41911,-13.37461]
# 	,[-13.78514,-13.78225,-13.77819,-13.77129,-13.76575,-13.76068,-13.75333
# 	,-13.74656,-13.73871,-13.72689,-13.71783,-13.70974,-13.69937,-13.68963
# 	,-13.67569,-13.64566,-13.61039,-13.56196,-13.47875,-13.40664,-13.33533
# 	,-13.25926,-13.21475,-13.17182]
# 	,[-13.53813,-13.53584,-13.53312,-13.52908,-13.52583,-13.52244,-13.51781
# 	,-13.51374,-13.50851,-13.49991,-13.49312,-13.48712,-13.47903,-13.47091
# 	,-13.45911,-13.43434,-13.40518,-13.36402,-13.2908,-13.22512,-13.15828
# 	,-13.08647,-13.04481,-13.0047,]
# 	,[-13.4025,-13.4005,-13.39888,-13.39635,-13.39386,-13.39121,-13.38759
# 	,-13.38451,-13.38079,-13.37436,-13.36863,-13.36287,-13.35595,-13.34954
# 	,-13.33928,-13.31739,-13.29176,-13.25524,-13.18843,-13.12708,-13.06376
# 	,-12.99505,-12.95492,-12.91626]
# 	,[-13.25488,-13.25416,-13.25267,-13.25039,-13.249,-13.24789,-13.24563
# 	,-13.2431,-13.24022,-13.23586,-13.2319,-13.22731,-13.22155,-13.21639
# 	,-13.20816,-13.19005,-13.16842,-13.13729,-13.07845,-13.02266,-12.96363
# 	,-12.89871,-12.8607,-12.82412]
# 	,[-13.12536,-13.12477,-13.1242,-13.12307,-13.12206,-13.12121,-13.11983
# 	,-13.11832,-13.11639,-13.11341,-13.11088,-13.10799,-13.10313,-13.09838
# 	,-13.0922,-13.07868,-13.0617,-13.03669,-12.98856,-12.94079,-12.88701
# 	,-12.82591,-12.78973,-12.75487]
# 	,[-13.06797,-13.0675,-13.06716,-13.06636,-13.06556,-13.0649,-13.06374
# 	,-13.06252,-13.06116,-13.059,-13.05683,-13.05408,-13.05018,-13.0466
# 	,-13.04137,-13.02981,-13.01559,-12.99457,-12.95263,-12.90955,-12.85995
# 	,-12.80248,-12.76781,-12.73429]
# 	,[-13.02256,-13.02236,-13.02224,-13.02177,-13.02116,-13.02051,-13.0197
# 	,-13.01898,-13.01791,-13.01597,-13.01406,-13.01172,-13.00854,-13.0056
# 	,-13.00111,-12.9911,-12.97886,-12.96079,-12.92395,-12.88516,-12.83952
# 	,-12.78436,-12.74936,-12.71524]
# 	,[-12.97917,-12.97901,-12.97865,-12.97823,-12.97796,-12.97756,-12.97699
# 	,-12.97646,-12.97561,-12.97413,-12.97279,-12.97123,-12.96865,-12.96598
# 	,-12.96223,-12.95412,-12.94396,-12.92867,-12.89801,-12.8646,-12.82251
# 	,-12.76951,-12.73521,-12.70167]
# 	,[-12.95164,-12.9513,-12.95108,-12.95077,-12.9505,-12.95031,-12.94988
# 	,-12.94938,-12.9488,-12.94768,-12.94654,-12.94525,-12.94312,-12.94082
# 	,-12.93746,-12.9304,-12.92154,-12.90788,-12.88056,-12.85007,-12.81028
# 	,-12.75988,-12.72755,-12.69598]
# 	,[-12.91679,-12.91665,-12.91674,-12.91668,-12.9164,-12.91613,-12.91564
# 	,-12.91524,-12.91493,-12.91391,-12.91278,-12.91179,-12.91006,-12.9081
# 	,-12.90547,-12.89949,-12.89177,-12.88021,-12.8577,-12.83271,-12.79935
# 	,-12.75271,-12.71917,-12.68585]
# 	,[-12.87687,-12.87703,-12.87751,-12.87776,-12.87747,-12.87704,-12.87643
# 	,-12.87618,-12.87618,-12.87525,-12.8741,-12.8735,-12.87222,-12.87068
# 	,-12.86901,-12.86431,-12.85792,-12.84902,-12.83265,-12.81475,-12.78984
# 	,-12.74793,-12.71244,-12.67639]
# 	,[-12.85212,-12.85249,-12.85323,-12.85368,-12.85338,-12.85285,-12.85216
# 	,-12.852,-12.85221,-12.85134,-12.85017,-12.84982,-12.84884,-12.84758
# 	,-12.84655,-12.84269,-12.83716,-12.83002,-12.81772,-12.80448,-12.78509
# 	,-12.74629,-12.70956,-12.67177]
# 	,[-12.82679,-12.82737,-12.82838,-12.82905,-12.82873,-12.82809,-12.82732
# 	,-12.82726,-12.82769,-12.82687,-12.8257,-12.82562,-12.82495,-12.82399
# 	,-12.82364,-12.82066,-12.81603,-12.81078,-12.80281,-12.79448,-12.78091
# 	,-12.74539,-12.70739,-12.66777]
# 	,[-12.79876,-12.79958,-12.80091,-12.80182,-12.80148,-12.80071,-12.79986
# 	,-12.79991,-12.80059,-12.79984,-12.79866,-12.79888,-12.79857,-12.79794
# 	,-12.79837,-12.79639,-12.79276,-12.78968,-12.78663,-12.78388,-12.77687
# 	,-12.74507,-12.70566,-12.66401]
# 	,[-12.77923,-12.78022,-12.78177,-12.78285,-12.7825,-12.78164,-12.78072
# 	,-12.78085,-12.7817,-12.781,-12.77982,-12.78026,-12.78021,-12.7798
# 	,-12.78079,-12.77951,-12.7766,-12.77505,-12.7755,-12.77669,-12.7743
# 	,-12.74513,-12.70475,-12.66169]
# 	,[-12.75203,-12.75327,-12.75514,-12.75646,-12.75608,-12.75509,-12.75408
# 	,-12.75433,-12.75543,-12.7548,-12.75361,-12.75435,-12.75466,-12.75458
# 	,-12.75635,-12.75606,-12.75414,-12.75477,-12.76015,-12.76687,-12.77095
# 	,-12.74546,-12.70376,-12.65874]
# 	,[-12.71814,-12.71969,-12.72195,-12.72358,-12.72318,-12.72202,-12.7209
# 	,-12.72128,-12.7227,-12.72216,-12.72096,-12.72209,-12.72285,-12.72318
# 	,-12.72595,-12.7269,-12.72619,-12.72958,-12.74118,-12.75482,-12.76696
# 	,-12.74608,-12.70279,-12.65536]
# 	,[-12.69596,-12.69771,-12.70024,-12.70207,-12.70165,-12.70039,-12.69919
# 	,-12.69967,-12.70129,-12.70081,-12.6996,-12.70099,-12.70205,-12.70264
# 	,-12.70607,-12.70783,-12.70792,-12.71313,-12.72883,-12.747,-12.7644
# 	,-12.74653,-12.70223,-12.65325]
# 	,[-12.67253,-12.6745,-12.67731,-12.67936,-12.67892,-12.67754,-12.67626
# 	,-12.67683,-12.67868,-12.67826,-12.67704,-12.67871,-12.68008,-12.68096
# 	,-12.68508,-12.6877,-12.68862,-12.69577,-12.71581,-12.73877,-12.76169
# 	,-12.747,-12.70166,-12.65108]
# 	,[-12.64597,-12.64818,-12.65132,-12.65361,-12.65314,-12.65164,-12.65026
# 	,-12.65095,-12.65305,-12.65269,-12.65147,-12.65345,-12.65518,-12.65638
# 	,-12.6613,-12.66489,-12.66673,-12.6761,-12.70108,-12.72944,-12.75857
# 	,-12.74747,-12.70101,-12.64865]])

interp = RectBivariateSpline(x, y, z)

print(interp(1.5,5))
