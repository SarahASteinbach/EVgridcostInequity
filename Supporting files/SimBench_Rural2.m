function mpc = SimBench_Rural2
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 0.25;

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
10 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
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
49 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
50 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
51 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
52 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
53 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
55 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
56 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
57 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
58 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
59 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
61 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
62 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
63 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
65 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
66 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
67 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
68 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
69 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
70 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
71 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
72 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
73 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
74 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
75 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
77 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
78 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
79 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
80 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
81 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
82 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
83 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
84 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
85 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
86 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
87 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
88 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
89 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
90 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
91 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
92 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
93 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
94 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
95 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
96 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
97 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
99 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
100 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
101 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
102 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
103 1 0.001 0 0 0 1 1 0 0.4 1 1 0.9;
104 3 0 0 0 0 1 1 0 20 1 1 0.9;
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
10 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
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
49 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
50 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
51 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
52 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
53 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
55 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
56 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
57 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
58 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
59 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
61 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
62 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
63 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
65 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
66 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
67 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
68 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
69 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
70 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
71 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
72 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
73 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
74 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
75 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
77 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
78 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
79 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
80 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
81 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
82 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
83 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
84 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
85 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
86 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
87 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
88 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
89 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
90 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
91 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
92 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
93 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
94 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
95 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
96 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
97 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
99 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
100 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
101 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
102 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
103 0 0 0 0 1 100 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
104 0.25 0 0.25 0 1 100 1 0.25 0 0 0 0 0 0 0 0 0 0 0 0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
104 19 0.001 0.001 0 0.25 0.25 0.25 1 0 1 -360 360;
15 66 0.0010839617 0.000420956 0 0.187 0.187 0.187 0 0 1 -360 360;
72 80 0.003689357 0.00143276 0 0.187 0.187 0.187 0 0 1 -360 360;
78 57 0.0061701738 0.002396184 0 0.187 0.187 0.187 0 0 1 -360 360;
21 94 0.0043748838 0.001698984 0 0.187 0.187 0.187 0 0 1 -360 360;
99 15 0.0032373106 0.001257208 0 0.187 0.187 0.187 0 0 1 -360 360;
52 20 0.0010418038 0.000404584 0 0.187 0.187 0.187 0 0 1 -360 360;
7 4 0.00061794232 0.0002399776 0 0.187 0.187 0.187 0 0 1 -360 360;
94 29 0.0031676414 0.001230152 0 0.187 0.187 0.187 0 0 1 -360 360;
75 53 0.0053682158 0.002084744 0 0.187 0.187 0.187 0 0 1 -360 360;
2 21 0.003819034 0.00148312 0 0.187 0.187 0.187 0 0 1 -360 360;
102 39 0.003468422 0.00134696 0 0.187 0.187 0.187 0 0 1 -360 360;
61 95 0.0036704668 0.001425424 0 0.187 0.187 0.187 0 0 1 -360 360;
103 83 0.0004875505 0.00018934 0 0.187 0.187 0.187 0 0 1 -360 360;
30 46 0.0040073386 0.001556248 0 0.187 0.187 0.187 0 0 1 -360 360;
96 92 0.004507486 0.00175048 0 0.187 0.187 0.187 0 0 1 -360 360;
68 19 0.00134091168 0.0005207424 0 0.187 0.187 0.187 0 0 1 -360 360;
10 26 0.0067597046 0.002625128 0 0.187 0.187 0.187 0 0 1 -360 360;
83 34 0.0044457684 0.001726512 0 0.187 0.187 0.187 0 0 1 -360 360;
22 89 0.005010435 0.0019458 0 0.187 0.187 0.187 0 0 1 -360 360;
13 93 0.00108167098 0.0004200664 0 0.187 0.187 0.187 0 0 1 -360 360;
12 72 0.0036520504 0.001418272 0 0.187 0.187 0.187 0 0 1 -360 360;
3 51 0.00062581358 0.0002430344 0 0.187 0.187 0.187 0 0 1 -360 360;
57 2 0.00139069364 0.0005400752 0 0.187 0.187 0.187 0 0 1 -360 360;
5 24 0.0034212274 0.001328632 0 0.187 0.187 0.187 0 0 1 -360 360;
24 45 0.003316909 0.00128812 0 0.187 0.187 0.187 0 0 1 -360 360;
92 65 0.0053659292 0.002083856 0 0.187 0.187 0.187 0 0 1 -360 360;
8 9 0.0085634818 0.003325624 0 0.187 0.187 0.187 0 0 1 -360 360;
77 33 0.001064711 0.00041348 0 0.187 0.187 0.187 0 0 1 -360 360;
4 70 0.00037883606 0.0001471208 0 0.187 0.187 0.187 0 0 1 -360 360;
41 87 0.00127915906 0.0004967608 0 0.187 0.187 0.187 0 0 1 -360 360;
70 59 0.00123818978 0.0004808504 0 0.187 0.187 0.187 0 0 1 -360 360;
80 91 0.0039992428 0.001553104 0 0.187 0.187 0.187 0 0 1 -360 360;
23 88 0.0003654749 0.000141932 0 0.187 0.187 0.187 0 0 1 -360 360;
95 75 0.0048125514 0.001868952 0 0.187 0.187 0.187 0 0 1 -360 360;
44 71 0.0050296754 0.001953272 0 0.187 0.187 0.187 0 0 1 -360 360;
37 41 0.0057455048 0.002231264 0 0.187 0.187 0.187 0 0 1 -360 360;
91 1 0.0060148498 0.002335864 0 0.187 0.187 0.187 0 0 1 -360 360;
34 73 0.00032097066 0.0001246488 0 0.187 0.187 0.187 0 0 1 -360 360;
101 37 0.0037434114 0.001453752 0 0.187 0.187 0.187 0 0 1 -360 360;
6 85 0.00137795666 0.0005351288 0 0.187 0.187 0.187 0 0 1 -360 360;
20 40 0.0040643594 0.001578392 0 0.187 0.187 0.187 0 0 1 -360 360;
19 96 0.0031465676 0.001221968 0 0.187 0.187 0.187 0 0 1 -360 360;
51 82 0.0032273608 0.001253344 0 0.187 0.187 0.187 0 0 1 -360 360;
97 99 0.00123830308 0.0004808944 0 0.187 0.187 0.187 0 0 1 -360 360;
62 100 0.00065312918 0.0002536424 0 0.187 0.187 0.187 0 0 1 -360 360;
32 56 0.00125593256 0.0004877408 0 0.187 0.187 0.187 0 0 1 -360 360;
35 63 0.00140800588 0.0005467984 0 0.187 0.187 0.187 0 0 1 -360 360;
65 43 0.0014114399 0.000548132 0 0.187 0.187 0.187 0 0 1 -360 360;
36 22 0.0004580719 0.000177892 0 0.187 0.187 0.187 0 0 1 -360 360;
82 35 0.00031347638 0.0001217384 0 0.187 0.187 0.187 0 0 1 -360 360;
100 17 0.0048172276 0.001870768 0 0.187 0.187 0.187 0 0 1 -360 360;
59 81 0.0010599318 0.000411624 0 0.187 0.187 0.187 0 0 1 -360 360;
56 27 0.00049047364 0.0001904752 0 0.187 0.187 0.187 0 0 1 -360 360;
89 69 0.0032817036 0.001274448 0 0.187 0.187 0.187 0 0 1 -360 360;
31 103 0.0006050117 0.000234956 0 0.187 0.187 0.187 0 0 1 -360 360;
85 68 0.0045480886 0.001766248 0 0.187 0.187 0.187 0 0 1 -360 360;
38 36 0.0045217618 0.001756024 0 0.187 0.187 0.187 0 0 1 -360 360;
14 58 0.003898035 0.0015138 0 0.187 0.187 0.187 0 0 1 -360 360;
29 32 0.004079624 0.00158432 0 0.187 0.187 0.187 0 0 1 -360 360;
53 42 0.0056605298 0.002198264 0 0.187 0.187 0.187 0 0 1 -360 360;
28 10 0.0045186306 0.001754808 0 0.187 0.187 0.187 0 0 1 -360 360;
73 8 0.00067986798 0.0002640264 0 0.187 0.187 0.187 0 0 1 -360 360;
11 67 0.0045263144 0.001757792 0 0.187 0.187 0.187 0 0 1 -360 360;
55 86 0.00125431958 0.0004871144 0 0.187 0.187 0.187 0 0 1 -360 360;
66 101 0.0052434004 0.002036272 0 0.187 0.187 0.187 0 0 1 -360 360;
43 3 0.005256605 0.0020414 0 0.187 0.187 0.187 0 0 1 -360 360;
93 14 0.00106558856 0.0004138208 0 0.187 0.187 0.187 0 0 1 -360 360;
39 23 0.0045432476 0.001764368 0 0.187 0.187 0.187 0 0 1 -360 360;
17 79 0.0003284155 0.00012754 0 0.187 0.187 0.187 0 0 1 -360 360;
79 77 0.00051526986 0.0002001048 0 0.187 0.187 0.187 0 0 1 -360 360;
40 61 0.0057973344 0.002251392 0 0.187 0.187 0.187 0 0 1 -360 360;
22 12 0.0041334724 0.001605232 0 0.187 0.187 0.187 0 0 1 -360 360;
49 102 0.00064554632 0.0002506976 0 0.187 0.187 0.187 0 0 1 -360 360;
84 16 0.0035273998 0.001369864 0 0.187 0.187 0.187 0 0 1 -360 360;
67 28 0.0045025626 0.001748568 0 0.187 0.187 0.187 0 0 1 -360 360;
86 50 0.0012550653 0.000487404 0 0.187 0.187 0.187 0 0 1 -360 360;
25 90 0.0035300572 0.001370896 0 0.187 0.187 0.187 0 0 1 -360 360;
18 38 0.0003480061 0.000135148 0 0.187 0.187 0.187 0 0 1 -360 360;
26 78 0.005925693 0.00230124 0 0.187 0.187 0.187 0 0 1 -360 360;
74 62 0.0045374384 0.001762112 0 0.187 0.187 0.187 0 0 1 -360 360;
58 7 0.0040957332 0.001590576 0 0.187 0.187 0.187 0 0 1 -360 360;
27 49 0.004978711 0.00193348 0 0.187 0.187 0.187 0 0 1 -360 360;
1 84 0.00063453768 0.0002464224 0 0.187 0.187 0.187 0 0 1 -360 360;
45 97 0.0037039006 0.001438408 0 0.187 0.187 0.187 0 0 1 -360 360;
88 13 0.0003361611 0.000130548 0 0.187 0.187 0.187 0 0 1 -360 360;
71 6 0.0042826164 0.001663152 0 0.187 0.187 0.187 0 0 1 -360 360;
81 44 0.000697413 0.00027084 0 0.187 0.187 0.187 0 0 1 -360 360;
63 74 0.00132367154 0.0005140472 0 0.187 0.187 0.187 0 0 1 -360 360;
90 31 0.0031040698 0.001205464 0 0.187 0.187 0.187 0 0 1 -360 360;
50 5 0.00046572686 0.0001808648 0 0.187 0.187 0.187 0 0 1 -360 360;
69 25 0.0035988612 0.001397616 0 0.187 0.187 0.187 0 0 1 -360 360;
26 52 0.010186494 0.00395592 0 0.187 0.187 0.187 0 0 1 -360 360;
19 55 0.01545 0.006 0 0.187 0.187 0.187 0 0 1 -360 360;
16 19 0.0013802 0.000536 0 0.187 0.187 0.187 0 0 1 -360 360;
30 11 0.006633406 0.00257608 0 0.187 0.187 0.187 0 0 1 -360 360;

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

