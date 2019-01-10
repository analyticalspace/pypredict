from pprint import pprint
import predict

tle = '0 OBJECT NY\n1 43550U 98067NY  19009.55938219 +.00013482 +00000-0 +17279-3 0  9995\n2 43550 051.6378 066.3469 0004106 279.6394 080.4135 15.59492665028120'

x = predict.observe(tle, [0,0,0])

pprint(x)
