function mpc = SimBench_Urban6
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 0.63;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [  %% (Pd and Qd are specified in kW & kVAr here and then converted to MW)
1 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
2 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
3 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
4 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
5 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
6 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
7 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
8 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
9 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
11 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
12 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
13 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
14 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
15 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
16 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
17 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
18 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
19 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
20 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
21 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
22 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
23 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
24 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
25 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
26 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
27 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
28 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
29 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
30 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
31 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
32 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
33 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
34 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
35 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
36 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
37 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
38 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
39 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
40 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
41 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
42 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
43 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
44 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
45 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
46 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
47 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
48 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
49 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
50 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
51 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
52 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
53 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
54 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
55 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
56 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
57 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
58 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
59 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
60 3 0 0 0 0 1 1 0 20 1 1 0.9;


];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
1 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
2 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
3 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
4 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
5 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
6 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
7 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
8 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
9 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
11 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
12 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
13 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
14 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
15 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
16 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
17 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
18 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
19 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
20 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
21 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
22 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
23 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
24 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
25 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
26 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
27 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
28 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
29 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
30 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
31 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
32 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
33 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
34 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
35 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
36 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
37 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
38 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
39 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
40 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
41 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
42 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
43 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
44 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
45 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
46 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
47 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
48 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
49 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
50 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
51 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
52 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
53 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
54 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
55 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
56 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
57 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
58 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
59 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
60 0.63 0 0 0 1 100 1 0.63 0 0 0 0 0 0 0 0 0 0 0 0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
60 9 0.001 0.001 0 0.63 0.63 0.63 1 0 1 -360 360;
9 17 0.00256761351 0.004174652 0 0.247 0.247 0.247 0 0 1 -360 360;
12 18 0.00186924311 0.00303918 0 0.247 0.247 0.247 0 0 1 -360 360;
5 34 0.0000274049566 0.0000445574 0 0.247 0.247 0.247 0 0 1 -360 360;
3 24 0.00486095953 0.007903375 0 0.247 0.247 0.247 0 0 1 -360 360;
20 43 0.000028810313 0.0000468423 0 0.247 0.247 0.247 0 0 1 -360 360;
20 37 0.00462595637 0.007521287 0 0.247 0.247 0.247 0 0 1 -360 360;
13 33 0.000384871522 0.000625758 0 0.247 0.247 0.247 0 0 1 -360 360;
15 2 0.00313991741 0.005105154 0 0.247 0.247 0.247 0 0 1 -360 360;
26 45 0.00228277924 0.003711543 0 0.247 0.247 0.247 0 0 1 -360 360;
30 13 0.0018361364 0.002985352 0 0.247 0.247 0.247 0 0 1 -360 360;
6 39 0.000310815372 0.000505351 0 0.247 0.247 0.247 0 0 1 -360 360;
9 29 0.00240129442 0.003904236 0 0.247 0.247 0.247 0 0 1 -360 360;
38 27 0.000542007396 0.000881243 0 0.247 0.247 0.247 0 0 1 -360 360;
2 19 0.00305932354 0.004974117 0 0.247 0.247 0.247 0 0 1 -360 360;
18 11 0.000471088338 0.000765937 0 0.247 0.247 0.247 0 0 1 -360 360;
28 6 0.00237677797 0.003864375 0 0.247 0.247 0.247 0 0 1 -360 360;
24 4 0.00131572882 0.002139228 0 0.247 0.247 0.247 0 0 1 -360 360;
40 46 0.00128871638 0.002095308 0 0.247 0.247 0.247 0 0 1 -360 360;
46 1 0.00358868881 0.005834806 0 0.247 0.247 0.247 0 0 1 -360 360;
31 30 0.000469840343 0.000763908 0 0.247 0.247 0.247 0 0 1 -360 360;
34 44 0.0000273099316 0.0000444029 0 0.247 0.247 0.247 0 0 1 -360 360;
35 41 0.00235400998 0.003827356 0 0.247 0.247 0.247 0 0 1 -360 360;
33 38 0.00192918488 0.003136638 0 0.247 0.247 0.247 0 0 1 -360 360;
29 25 0.0000237420596 0.0000386019 0 0.247 0.247 0.247 0 0 1 -360 360;
7 12 0.00166232934 0.002702761 0 0.247 0.247 0.247 0 0 1 -360 360;
36 6 0.00226132893 0.003676667 0 0.247 0.247 0.247 0 0 1 -360 360;
39 7 0.000949297216 0.001543451 0 0.247 0.247 0.247 0 0 1 -360 360;
3 35 0.00257406254 0.004185137 0 0.247 0.247 0.247 0 0 1 -360 360;
22 8 0.00290722019 0.004726814 0 0.247 0.247 0.247 0 0 1 -360 360;
27 4 0.00272606453 0.004432275 0 0.247 0.247 0.247 0 0 1 -360 360;
45 31 0.00183366575 0.002981335 0 0.247 0.247 0.247 0 0 1 -360 360;
44 22 0.00175891275 0.002859795 0 0.247 0.247 0.247 0 0 1 -360 360;
17 28 0.00508695432 0.008270818 0 0.247 0.247 0.247 0 0 1 -360 360;
43 3 0.00219246748 0.003564706 0 0.247 0.247 0.247 0 0 1 -360 360;
42 1 0.00985040553 0.016015655 0 0.247 0.247 0.247 0 0 1 -360 360;
11 23 0.00180884522 0.00294098 0 0.247 0.247 0.247 0 0 1 -360 360;
25 2 0.00172094076 0.002798057 0 0.247 0.247 0.247 0 0 1 -360 360;
4 32 0.00654591749 0.010642928 0 0.247 0.247 0.247 0 0 1 -360 360;
23 40 0.000966174923 0.001570892 0 0.247 0.247 0.247 0 0 1 -360 360;
21 26 0.000428662843 0.000696958 0 0.247 0.247 0.247 0 0 1 -360 360;
16 14 0.00465874633 0.007574599 0 0.247 0.247 0.247 0 0 1 -360 360;
16 15 0.00659356936 0.010720405 0 0.247 0.247 0.247 0 0 1 -360 360;
19 47 0.0032942 0.005356 0 0.247 0.247 0.247 0 0 1 -360 360;
48 49 0.0021539 0.003502 0 0.247 0.247 0.247 0 0 1 -360 360;
50 51 0.0020272 0.003296 0 0.247 0.247 0.247 0 0 1 -360 360;
51 52 0.0024073 0.003914 0 0.247 0.247 0.247 0 0 1 -360 360;
49 50 0.0027874 0.004532 0 0.247 0.247 0.247 0 0 1 -360 360;
52 53 0.0015204 0.002472 0 0.247 0.247 0.247 0 0 1 -360 360;
53 54 0.002534 0.00412 0 0.247 0.247 0.247 0 0 1 -360 360;
9 48 0.0027874 0.004532 0 0.247 0.247 0.247 0 0 1 -360 360;
9 55 0.0057015 0.00927 0 0.247 0.247 0.247 0 0 1 -360 360;
56 9 0.0017738 0.002884 0 0.247 0.247 0.247 0 0 1 -360 360;
56 57 0.0020272 0.003296 0 0.247 0.247 0.247 0 0 1 -360 360;
57 58 0.0019005 0.00309 0 0.247 0.247 0.247 0 0 1 -360 360;
58 59 0.001267 0.00206 0 0.247 0.247 0.247 0 0 1 -360 360;
37 9 0.00374762129 0.006093212 0 0.247 0.247 0.247 0 0 1 -360 360;
9 5 0.00230927221 0.003754618 0 0.247 0.247 0.247 0 0 1 -360 360;
];


%% cost data
%ID rateA cost (€)
mpc.cos = [    
1 0.63 61730;
2 0.247 2143.24786708861;
3 0.247 1560.30153797468;
4 0.247 22.8755669620253;
5 0.247 4057.55815822785;
6 0.247 24.0486512658228;
7 0.247 3861.39544936709;
8 0.247 321.261383544304;
9 0.247 2620.9635;
10 0.247 1905.48994936709;
11 0.247 1532.66658227848;
12 0.247 259.444959493671;
13 0.247 2004.41737974684;
14 0.247 452.42642278481;
15 0.247 2553.68988607595;
16 0.247 393.228603797468;
17 0.247 1983.95291772152;
18 0.247 1098.27003797468;
19 0.247 1075.72211392405;
20 0.247 2995.56362658228;
21 0.247 392.18687278481;
22 0.247 22.7962473417722;
23 0.247 1964.94793670886;
24 0.247 1610.33635443038;
25 0.247 19.81806;
26 0.247 1387.58570886076;
27 0.247 1887.58486708861;
28 0.247 792.400891139241;
29 0.247 2148.63102531646;
30 0.247 2426.72561392405;
31 0.247 2275.51068987342;
32 0.247 1530.6042721519;
33 0.247 1468.20617088608;
34 0.247 4246.20136708861;
35 0.247 1830.10458227848;
36 0.247 8222.36701898734;
37 0.247 1509.88598734177;
38 0.247 1436.51005063291;
39 0.247 5464.03251265823;
40 0.247 806.489113291139;
41 0.247 357.815037341772;
42 0.247 3888.76600632911;
43 0.247 5503.80865822785;
44 0.247 2749.74683544304;
45 0.247 1797.91139240506;
46 0.247 1692.15189873418;
47 0.247 2009.43037974684;
48 0.247 2326.70886075949;
49 0.247 1269.11392405063;
50 0.247 2115.18987341772;
51 0.247 2326.70886075949;
52 0.247 4759.17721518987;
53 0.247 1480.63291139241;
54 0.247 1692.15189873418;
55 0.247 1586.39240506329;
56 0.247 1057.59493670886;
57 0.247 3128.22833544304;
58 0.247 1927.60425949367;
];

%% convert branch impedances from Ohms to p.u.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
Vbase = mpc.bus(3, BASE_KV) * 1e3;      %% in Volts
Sbase = mpc.baseMVA * 1e6;              %% in VA
mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);

